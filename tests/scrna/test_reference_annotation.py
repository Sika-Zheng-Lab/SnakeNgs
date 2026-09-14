import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from reference_utils import harmonize_counts, select_count_source, unique_batch_ids


def make_reference(matrix, layers=None):
    obs = pd.DataFrame({'cell_type': ['A', 'B'], 'reference_batch': ['b1', 'b2']}, index=['r1', 'r2'])
    var = pd.DataFrame({'gene_symbol': ['G2', 'G1', 'G1']}, index=['id2', 'id1', 'id1b'])
    return ad.AnnData(sparse.csr_matrix(matrix), obs=obs, var=var, layers=layers or {})


def test_raw_counts_in_x_are_detected():
    assert select_count_source(make_reference([[1, 0, 2], [0, 3, 1]])) == 'X'


def test_raw_counts_in_layer_are_detected():
    reference = make_reference(np.log1p([[1, 0, 2], [0, 3, 1]]), {'counts': sparse.csr_matrix([[1, 0, 2], [0, 3, 1]])})
    assert select_count_source(reference) == 'counts'


def test_normalized_only_reference_fails():
    with pytest.raises(ValueError, match='raw integer counts'):
        select_count_source(make_reference(np.array([[.2, 0, .8], [0, .4, .6]])))


def test_genes_align_by_identity_and_sum_duplicate_symbols():
    reference = make_reference([[1, 2, 3], [4, 5, 6]])
    query = ad.AnnData(sparse.csr_matrix([[10, 20, 30], [40, 50, 60]]),
                       obs=pd.DataFrame(index=['q1', 'q2']),
                       var=pd.DataFrame({'gene_symbol': ['G1', 'G2', 'G3']}, index=['x1', 'x2', 'x3']))
    ref, qry, genes, stats = harmonize_counts(reference, query, 'X', 'X', 'gene_symbol', 'gene_symbol', 2)
    assert list(genes) == ['G1', 'G2']
    np.testing.assert_array_equal(ref.toarray(), [[5, 1], [11, 4]])
    np.testing.assert_array_equal(qry.toarray(), [[10, 20], [40, 50]])
    assert stats['reference_duplicates_removed'] == 1


def test_reference_and_query_batch_ids_cannot_collide():
    reference, query = unique_batch_ids(['shared', 'r2'], ['shared', 'q2'])
    assert set(reference).isdisjoint(query)
    assert reference[0] == 'reference::shared'
    assert query[0] == 'query::shared'
