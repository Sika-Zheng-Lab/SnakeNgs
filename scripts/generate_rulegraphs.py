#!/usr/bin/env python3
"""Build documentation rulegraphs from real Snakefiles and isolated dummy inputs."""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET

import yaml

ROOT = Path(__file__).resolve().parents[1]


def run(args, **kwargs):
    result = subprocess.run(args, text=True, capture_output=True, timeout=180, **kwargs)
    if result.returncode:
        raise RuntimeError(f"Command failed: {args}\n{result.stdout}\n{result.stderr}")
    return result.stdout


def fixtures(base):
    run(['bash', str(ROOT / '.github/scripts/generate_test_data.sh')],
        cwd=base, env={**os.environ, 'GITHUB_WORKSPACE': str(base)})
    data = base / 'tests/testdata'
    cases = []
    for path in sorted((data / 'configs').glob('*.yaml')):
        name = path.stem.removeprefix('config_')
        pipeline = name
        for prefix in ['preprocessing_RNAseq', 'preprocessing_ChIPseq', 'preprocessing_CLIPseq']:
            if name.startswith(prefix + '_') and name != 'preprocessing_RNAseq_long':
                pipeline = prefix
        if name == 'kb_nac':
            pipeline = 'kb-nac'
        cases.append((name, pipeline, path))
    for layout in ['paired', 'single']:
        name = f'denovo_assembly_RNAseq_{layout}'
        work = data / 'workdirs' / name
        (work / 'fastq').mkdir(parents=True)
        for suffix in (['_1', '_2'] if layout == 'paired' else ['']):
            (work / 'fastq' / f'sample1{suffix}.fastq.gz').touch()
        cfg = dict(workdir=str(work), samples=['sample1'], layout=layout,
                   swissprot_db=str(data / 'swissprot'), max_memory='2G')
        path = data / 'configs' / f'config_{name}.yaml'
        path.write_text(yaml.safe_dump(cfg))
        cases.append((name, 'denovo_assembly_RNAseq', path))
    cfg = yaml.safe_load((ROOT / 'examples/scrna_seq/config.yaml').read_text())
    cfg['workdir'] = str(data / 'workdirs/scrna_seq')
    cfg['sample_table'] = str(data / 'scrna_samples.tsv')
    cfg['reference'].update(fasta=str(data / 'dummy.fa'), gtf=str(data / 'dummy.gtf'))
    cfg['cellranger']['transcriptome'] = str(data / 'transcriptome')
    cfg['celltypist']['model_path'] = str(data / 'model.pkl')
    (data / 'model.pkl').touch()
    (data / 'scrna_samples.tsv').write_text(
        'sample\tbatch\tdonor\tcondition\tR1\tR2\n'
        f'sample1\tb1\td1\tcontrol\t{data}/fastqs/sample1_1.fastq.gz\t{data}/fastqs/sample1_2.fastq.gz\n')
    for mode in ['fastq', 'counts']:
        cfg['input_mode'] = mode
        if mode == 'counts':
            cfg['sample_table'] = str(data / 'scrna_counts.tsv')
            for f in ['filtered.h5ad', 'unfiltered.h5ad', 'metrics.csv']:
                (data / f).touch()
            (data / 'scrna_counts.tsv').write_text(
                'sample\tbatch\tdonor\tcondition\tkb_filtered\tkb_unfiltered\tcellranger_metrics\n'
                f'sample1\tb1\td1\tcontrol\t{data}/filtered.h5ad\t{data}/unfiltered.h5ad\t{data}/metrics.csv\n')
        name = f'scrna_seq_{mode}'
        path = data / 'configs' / f'config_{name}.yaml'
        path.write_text(yaml.safe_dump(cfg))
        cases.append((name, 'scrna_seq', path))
    expected = {p.stem for p in (ROOT / 'snakefile').glob('*.smk')}
    covered = {pipeline for _, pipeline, _ in cases}
    if expected != covered:
        raise ValueError(f'Rulegraph coverage mismatch: missing={expected-covered}, obsolete={covered-expected}')
    return sorted(cases)


def main():
    output = ROOT / 'docs/assets/rulegraphs'
    with tempfile.TemporaryDirectory(prefix='snakengs-rulegraphs-') as tmp:
        base = Path(tmp)
        cases = fixtures(base)
        staged = base / 'graphs'
        staged.mkdir()
        groups = {}
        for name, pipeline, cfg in cases:
            print(f'Rendering {name}', flush=True)
            work = base / 'tests/testdata/workdirs/bam2cram' if pipeline == 'bam2cram' else base
            dot = run(['snakemake', '-s', str(ROOT / 'snakefile' / f'{pipeline}.smk'),
                       '--configfile', str(cfg), '--cores', '1', '--forceall',
                       '--rulegraph', 'dot', 'all'], cwd=work)
            # Do not silently accept status messages or empty/incomplete graphs.
            if not dot.lstrip().startswith('digraph') or 'label = "all"' not in dot:
                raise ValueError(f'Invalid rulegraph for {name}: {dot[:500]}')
            dotfile = staged / f'{name}.dot'
            dotfile.write_text('\n'.join(line.rstrip() for line in dot.splitlines()) + '\n')
            for suffix, skip in [('', []), ('.slim', ['-s', 'all', 'multiqc'])]:
                svg = staged / f'{name}{suffix}.svg'
                run(['snakevision', str(dotfile), '-o', str(svg), '-y', 'label_font_size=14', 'label_font_color=#1f2937', *skip], cwd=base)
                tree = ET.parse(svg).getroot()
                if not tree.tag.endswith('svg'):
                    raise ValueError(f'Invalid SVG: {svg}')
                labels = {node.text for node in tree.iter() if node.tag.endswith('}text')}
                expected_labels = set(re.findall(r'label\s*=\s*"([^"]+)"', dot))
                if suffix:
                    expected_labels -= {'all', 'multiqc'}
                if labels != expected_labels:
                    raise ValueError(f'Incorrect SVG rule labels: {svg}')
            groups.setdefault(pipeline, []).append(name)
        output.mkdir(parents=True, exist_ok=True)
        # Publish only after all pipelines have rendered successfully.
        for p in staged.iterdir():
            (output / p.name).write_bytes(p.read_bytes())
        gallery = ['# Pipeline rulegraphs\n',
                   'Generated from the current Snakefiles using Snakemake and '
                   '[snakevision](https://github.com/OpenOmics/snakevision). '
                   'The slim view hides the aggregate `all` and `multiqc` nodes; '
                   'the full graph retains every rule. These are representative configurations, '
                   'not a union of every optional setting.\n']
        for pipeline, names in groups.items():
            gallery.append(f'## {pipeline}\n')
            for name in names:
                gallery.append(f'### {name}\n\n![{name}](assets/rulegraphs/{name}.slim.svg)\n\n'
                               f'[Full SVG](assets/rulegraphs/{name}.svg) · [DOT source](assets/rulegraphs/{name}.dot)\n')
        (ROOT / 'docs/rulegraphs.md').write_text('\n'.join(gallery))
        for pipeline, names in groups.items():
            page = ROOT / 'docs/usage' / f'{pipeline}.md'
            if not page.exists():
                continue
            content = page.read_text()
            block = '<!-- rulegraphs:start -->\n\n'
            for name in names:
                block += (f'**{name}**\n\n![{name} rulegraph](../assets/rulegraphs/{name}.slim.svg)\n\n'
                          f'[Full graph](../assets/rulegraphs/{name}.svg) · '
                          f'[Graph generation](../rulegraphs-maintenance.md)\n\n')
            block += '<!-- rulegraphs:end -->'
            if '<!-- rulegraphs:start -->' in content:
                content = re.sub(r'<!-- rulegraphs:start -->.*?<!-- rulegraphs:end -->',
                                 lambda _: block, content, flags=re.S)
            else:
                title, rest = content.split('\n', 1)
                content = title + '\n\n' + block + '\n' + rest
            page.write_text(content)

        print(f'Generated {len(cases)} configurations covering {len(groups)} pipelines.')


if __name__ == '__main__':
    main()
