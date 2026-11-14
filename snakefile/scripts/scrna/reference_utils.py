"""Validation and deterministic harmonization for external AnnData references."""
import numpy as np
import pandas as pd
from scipy import sparse

from core import integer_counts


def _matrix(adata, source):
    if source == 'X':
        return adata.X
    if source == 'raw':
        if adata.raw is None:
            return None
        return adata.raw.X
    return adata.layers.get(source)


def _valid_counts(matrix, min_library_size=1):
    if matrix is None:
        return False
    try:
        counts = integer_counts(matrix, 'reference raw counts')
    except (TypeError, ValueError):
        return False
    libraries = np.asarray(counts.sum(axis=1)).ravel()
    return bool(len(libraries) and np.all(np.isfinite(libraries))
                and np.all(libraries >= min_library_size)
                and np.median(libraries) >= min_library_size)


def select_count_source(adata, requested='auto'):
    """Return a validated raw-count source name without silently using log data."""
    available = ['X'] + list(adata.layers.keys())
    if adata.raw is not None:
        available.append('raw')
    if requested != 'auto':
        if requested not in available:
            raise ValueError(f"Reference count_source {requested!r} is unavailable; choose X, raw, or a layer")
        if not _valid_counts(_matrix(adata, requested)):
            raise ValueError(f"Reference count_source {requested!r} is not a finite nonnegative integer-count matrix")
        return requested
    valid = [name for name in available if _valid_counts(_matrix(adata, name))]
    if not valid:
        raise ValueError('Could not identify raw integer counts in reference .X, .raw, or any layer; normalized/scaled data are invalid')
    preferred = ['counts', 'counts_raw', 'raw', 'X']
    return next((name for name in preferred if name in valid), valid[0])


def gene_identifiers(adata, symbol_column=None):
    """Build stable identifiers, preferring configured symbols and falling back to var names."""
    var_names = pd.Index(adata.var_names.astype(str))
    if symbol_column and symbol_column in adata.var:
        symbols = adata.var[symbol_column].astype(str)
        valid = symbols.notna() & symbols.str.strip().ne('') & symbols.ne('nan')
        identifiers = pd.Index(np.where(valid, symbols, var_names))
    else:
        identifiers = var_names
    identifiers = identifiers.str.replace(r'\.[0-9]+$', '', regex=True)
    return identifiers


def _gene_owner(adata, source):
    return adata.raw if source == 'raw' and adata.raw is not None else adata


def aggregate_by_identifier(matrix, identifiers):
    """Sum duplicate identifiers in first-seen order while retaining sparse storage."""
    identifiers = pd.Index(identifiers.astype(str))
    codes, names = pd.factorize(identifiers, sort=False)
    mapping = sparse.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)),
                                shape=(len(codes), len(names)))
    return (sparse.csr_matrix(matrix) @ mapping).tocsr(), pd.Index(names)


def harmonize_counts(reference, query, reference_source, query_source, reference_symbol_column=None,
                     query_symbol_column='gene_symbol', min_shared_genes=1000):
    """Align count matrices by identifiers and return shared genes in query order."""
    ref_owner = _gene_owner(reference, reference_source)
    query_owner = _gene_owner(query, query_source)
    ref_x, ref_ids = aggregate_by_identifier(
        _matrix(reference, reference_source), gene_identifiers(ref_owner, reference_symbol_column))
    query_x, query_ids = aggregate_by_identifier(
        _matrix(query, query_source), gene_identifiers(query_owner, query_symbol_column))
    ref_positions = {gene: i for i, gene in enumerate(ref_ids)}
    query_positions = {gene: i for i, gene in enumerate(query_ids)}
    shared = pd.Index([gene for gene in query_ids if gene in ref_positions])
    if len(shared) < min_shared_genes:
        raise ValueError(f'Only {len(shared)} shared genes are available; need at least {min_shared_genes}')
    ref_cols = [ref_positions[gene] for gene in shared]
    query_cols = [query_positions[gene] for gene in shared]
    stats = {'reference_genes': int(len(ref_ids)), 'query_genes': int(len(query_ids)),
             'shared_genes': int(len(shared)), 'shared_fraction_query': float(len(shared) / len(query_ids)),
             'reference_duplicates_removed': int(ref_owner.n_vars - len(ref_ids)),
             'query_duplicates_removed': int(query_owner.n_vars - len(query_ids)),
             'genes_used_by_model': int(len(shared))}
    return ref_x[:, ref_cols].tocsr(), query_x[:, query_cols].tocsr(), shared, stats


def unique_batch_ids(reference_batches, query_batches):
    """Namespace reference and query batches so equal labels cannot collapse."""
    return (np.array([f'reference::{value}' for value in reference_batches], dtype=str),
            np.array([f'query::{value}' for value in query_batches], dtype=str))
