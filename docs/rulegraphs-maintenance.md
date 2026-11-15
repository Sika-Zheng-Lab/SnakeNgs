# Maintaining pipeline rulegraphs

Run the following from the repository root using Python 3.12:

```bash
pip install -r docs/requirements.txt
python scripts/generate_rulegraphs.py
mkdocs build --strict
```

The generator uses each real top-level `snakefile/*.smk`, requests the explicit
`all` target with `--forceall --rulegraph dot`, and renders it with
[snakevision](https://github.com/OpenOmics/snakevision). It creates isolated
placeholder inputs and configs in a temporary directory. No sequencing jobs,
container pulls, model downloads or reference-index builds are executed.
The placeholder files are suitable only for constructing graphs.

The [gallery](rulegraphs.md) covers every top-level pipeline. Paired/single RNA,
ChIP and de novo assembly, both CLIP examples, and FASTQ/counts scRNA-seq have
separate views. Other workflows use the existing CI example configuration.
These views represent those configurations; optional branches such as alternative
annotation backends or grouped FASTQ merging may differ in a user's run.

Each configuration produces three files under `docs/assets/rulegraphs/`:

- `.slim.svg`: the default documentation view, hiding `all` and `multiqc`.
- `.svg`: the complete rulegraph, including reporting and aggregate targets.
- `.dot`: the original Snakemake rulegraph for inspection or alternate rendering.

The slim view uses snakevision's `--skip-rules` option. It does not change the
workflow. MultiQC remains available in the full graph. Single-node pipelines
may look identical in both views after removal of `all`.

The generator validates pipeline coverage and SVG syntax, and stops if any graph
cannot be produced. Add a representative config when adding a new pipeline;
missing pipeline coverage fails the documentation build. Generated blocks in
usage pages are delimited by `<!-- rulegraphs:start -->` and
`<!-- rulegraphs:end -->`; edit prose outside these blocks.

The Publish MkDocs workflow regenerates graphs before building the site on
pull requests and pushes to `develop`. Pull requests upload a site preview
artifact; pushes to `develop` also deploy the site. SVGs and links are served
from the documentation itself and do not depend on raw GitHub image URLs.

## Animated versions

Each graph also has `.animated.svg` (full) and `.slim.animated.svg` variants.
Follow **Animated graph** on a usage page or the animated links in the gallery.
They are standalone SVG files with embedded CSS; no JavaScript or video player
is required. The static images remain the default documentation view.

Nodes pulse according to their dependency level in the original DOT graph.
Parallel branches at the same level pulse together. The timing illustrates
workflow structure, not measured runtime, execution scheduling, or live progress.
Open the animated SVG directly and hover over a node or focus the SVG with Tab to pause.
The OS/browser reduced-motion preference disables animation automatically.

Animation is added by `scripts/animate_rulegraph.py` after snakevision renders
the graph. Static SVG and DOT files remain available alongside each animation.
