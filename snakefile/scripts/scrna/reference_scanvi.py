import sys
sys.path.insert(0, snakemake.scriptdir)
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scipy import sparse
from core import aggregate_class_probabilities, composition_intervals, map_major_cell_types, uncertainty, validate_probabilities
from model_utils import configure_training, save_history
from reference_utils import harmonize_counts, select_count_source, unique_batch_ids
from runtime import logged, sha256, versions, write_json


def _subsample_reference(reference, cfg):
    p = cfg['annotation_reference']
    limit = p['max_cells_per_cluster']
    if limit is None:
        return reference
    key = p['cluster_key'] or p['label_key']
    if key not in reference.obs:
        raise ValueError(f'Reference cluster_key {key!r} is absent from reference.obs')
    rng = np.random.default_rng(cfg['seed'])
    keep = []
    for _, indices in reference.obs.groupby(key, observed=True, sort=True).indices.items():
        indices = np.asarray(indices)
        keep.extend(rng.choice(indices, size=min(limit, len(indices)), replace=False))
    return reference[np.sort(np.asarray(keep))].copy()


def _reference_metadata(reference, cfg):
    p = cfg['annotation_reference']
    for key in [p['label_key'], p['batch_key']]:
        if key not in reference.obs:
            raise ValueError(f'Reference is missing required obs column {key!r}')
    labels = reference.obs[p['label_key']].astype(str)
    batches = reference.obs[p['batch_key']].astype(str)
    if labels.isin(['', 'nan', 'None']).any() or batches.isin(['', 'nan', 'None']).any():
        raise ValueError('Reference labels and batches must be nonempty for every cell')
    return labels, batches


with logged(snakemake):
    cfg = snakemake.params.cfg
    p = cfg['annotation_reference']
    query = ad.read_h5ad(snakemake.input.query)
    reference = ad.read_h5ad(snakemake.input.reference)
    reference = _subsample_reference(reference, cfg)
    reference_labels, reference_batches = _reference_metadata(reference, cfg)
    query_source = cfg['depth']['model_layer']
    reference_source = select_count_source(reference, p['count_source'])
    ref_counts, query_counts, genes, gene_stats = harmonize_counts(
        reference, query, reference_source, query_source, p['gene_symbol_column'],
        'gene_symbol', p['min_shared_genes'])

    query_batches = query.obs[cfg['integration']['batch_key']].astype(str).to_numpy()
    ref_batch_ids, query_batch_ids = unique_batch_ids(reference_batches, query_batches)
    ref_obs = pd.DataFrame(index=pd.Index(['reference::' + value for value in reference.obs_names.astype(str)]))
    ref_obs['annotation_source'] = 'reference'
    ref_obs['annotation_label'] = reference_labels.to_numpy()
    ref_obs['annotation_batch'] = ref_batch_ids
    q_obs = query.obs.copy()
    q_obs['annotation_source'] = 'query'
    q_obs['annotation_label'] = p['unlabeled_category']
    q_obs['annotation_batch'] = query_batch_ids
    combined_obs = pd.concat([ref_obs, q_obs], axis=0)
    combined = ad.AnnData(sparse.vstack([ref_counts, query_counts], format='csr'), obs=combined_obs,
                          var=pd.DataFrame(index=genes))
    combined.layers['counts'] = combined.X.copy()

    configure_training(cfg, snakemake.threads)
    if combined.n_vars <= p['n_hvg']:
        combined.var['highly_variable'] = True
    else:
        sc.pp.highly_variable_genes(combined, layer='counts', flavor='seurat_v3',
                                    n_top_genes=p['n_hvg'], batch_key='annotation_batch')
    hvg = combined[:, combined.var['highly_variable']].copy()
    if hvg.n_vars < 2:
        raise ValueError('Reference annotation requires at least two highly variable genes')
    scvi.model.SCVI.setup_anndata(hvg, layer='counts', batch_key='annotation_batch', labels_key='annotation_label')
    scvi_model = scvi.model.SCVI(hvg, n_hidden=p['n_hidden'], n_layers=p['n_layers'],
                                 n_latent=p['n_latent'], gene_likelihood='nb')
    scvi_model.train(max_epochs=p['scvi_max_epochs'], accelerator=cfg['integration']['accelerator'],
                     devices=cfg['integration']['devices'], batch_size=p['batch_size'], early_stopping=False)
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, adata=hvg, labels_key='annotation_label', unlabeled_category=p['unlabeled_category'])
    scanvi_model.train(max_epochs=p['scanvi_max_epochs'], accelerator=cfg['integration']['accelerator'],
                       devices=cfg['integration']['devices'], batch_size=p['batch_size'], early_stopping=False)

    query_rows = np.flatnonzero(combined.obs['annotation_source'].to_numpy() == 'query')
    total = squares = conditional_entropy = None
    for _ in range(p['posterior_samples']):
        prediction = scanvi_model.predict(soft=True, use_posterior_mean=False)
        if not isinstance(prediction, pd.DataFrame) or not prediction.index.equals(hvg.obs_names):
            raise ValueError('Unexpected scANVI reference prediction index/type')
        values = validate_probabilities(prediction.to_numpy())
        if total is None:
            classes = prediction.columns.to_numpy(dtype=str)
            total = np.zeros_like(values)
            squares = np.zeros_like(values)
            conditional_entropy = np.zeros(values.shape[0])
        elif not np.array_equal(classes, prediction.columns.to_numpy(dtype=str)):
            raise ValueError('scANVI class ordering changed between posterior draws')
        total += values
        squares += values * values
        conditional_entropy += uncertainty(values)['entropy']
    probabilities = total / p['posterior_samples']
    sd = np.sqrt(np.maximum(0, squares / p['posterior_samples'] - probabilities ** 2))
    query_probabilities = probabilities[query_rows]
    metrics = uncertainty(query_probabilities)
    output = query.copy()
    for name, values in metrics.items():
        output.obs['scanvi_' + name] = values
    output.obs['scanvi_latent_information'] = np.maximum(
        0, metrics['entropy'] - conditional_entropy[query_rows] / p['posterior_samples'])
    output.obsm['scanvi_probabilities'] = query_probabilities.astype(np.float32)
    output.obsm['scanvi_probability_sd'] = sd[query_rows].astype(np.float32)
    output.uns['scanvi_classes'] = classes
    output.obs['scanvi_label'] = pd.Categorical(classes[query_probabilities.argmax(axis=1)])
    confident = ((metrics['max_probability'] >= p['min_probability'])
                 & (metrics['normalized_entropy'] <= p['max_normalized_entropy']))
    output.obs['annotation_uncertain'] = ~confident
    output.obs['cell_type'] = output.obs['scanvi_label'].astype(str).where(confident, p['unlabeled_category'])
    major_key = p['major_label_key']
    if major_key and major_key in reference.obs:
        label_to_major = dict(zip(reference_labels, reference.obs[major_key].astype(str)))
        major_map = {label: label_to_major[label] for label in classes if label in label_to_major}
    else:
        major_map = cfg['annotation']['major_celltype_map']
    output.obs['scanvi_major_label'] = pd.Categorical(map_major_cell_types(
        output.obs['scanvi_label'].astype(str), major_map, cfg['annotation']['unmapped_major_celltype'], p['unlabeled_category']))
    output.obs['major_cell_type'] = pd.Categorical(map_major_cell_types(
        output.obs['cell_type'].astype(str), major_map, cfg['annotation']['unmapped_major_celltype'], p['unlabeled_category']))
    latent = scanvi_model.get_latent_representation()
    hvg.obsm['X_scANVI'] = latent
    output.obsm['X_scANVI_reference'] = latent[query_rows]
    sc.pp.neighbors(hvg, use_rep='X_scANVI', n_neighbors=min(cfg['integration']['n_neighbors'], hvg.n_obs - 1),
                    key_added='scanvi_reference', random_state=cfg['seed'])
    sc.tl.umap(hvg, neighbors_key='scanvi_reference', random_state=cfg['seed'])
    output.obsm['X_umap_scANVI_reference'] = hvg.obsm['X_umap'][query_rows].copy()
    output.uns['annotation_reference'] = {
        'path': str(Path(snakemake.input.reference).resolve()), 'sha256': sha256(snakemake.input.reference),
        'n_reference_cells': int(reference.n_obs), 'n_reference_labels': int(reference_labels.nunique()),
        'reference_batches': sorted(reference_batches.unique().tolist()),
        'query_batches': sorted(set(query_batches.tolist())), 'count_source': reference_source,
        **gene_stats, 'n_hvg': int(hvg.n_vars), 'model_parameters': {
            key: p[key] for key in ['n_hidden', 'n_layers', 'n_latent', 'scvi_max_epochs', 'scanvi_max_epochs', 'batch_size']}}
    if isinstance(reference.uns.get('snakengs_reference'), dict):
        output.uns['annotation_reference_metadata'] = reference.uns['snakengs_reference']

    composition = composition_intervals(query_probabilities, output.obs, classes, cfg['scanvi']['composition_draws'], cfg['seed'])
    major_probabilities, major_classes = aggregate_class_probabilities(
        query_probabilities, classes, major_map, cfg['annotation']['unmapped_major_celltype'])
    composition_major = composition_intervals(major_probabilities, output.obs, major_classes,
                                              cfg['scanvi']['composition_draws'], cfg['seed'])
    scanvi_model.save(snakemake.output.model, overwrite=True, save_anndata=False)
    save_history(scanvi_model, snakemake.output.history)
    for frame in [composition, composition_major]:
        for key in ['donor', 'condition', 'batch']:
            metadata = output.obs[['sample', key]].drop_duplicates().set_index('sample')[key]
            frame[key] = frame['sample'].map(metadata)
    output.uns['uncertainty_scope'] = 'Conditional on reference labels, observed query cells and fitted weights; excludes reference-model uncertainty.'
    output.obs.to_csv(snakemake.output.cells, sep='\t', index_label='cell_id')
    composition.to_csv(snakemake.output.composition, sep='\t', index=False)
    composition_major.to_csv(snakemake.output.composition_major, sep='\t', index=False)
    output.write_h5ad(snakemake.output.adata, compression='gzip')
    write_json(snakemake.output.summary, {
        'enabled': True, 'backend': 'scanvi_reference', 'uncertain_cells': int((~confident).sum()),
        'n_classes': len(classes), 'classes': classes.tolist(), 'reference': output.uns['annotation_reference'],
        'versions': versions()})
