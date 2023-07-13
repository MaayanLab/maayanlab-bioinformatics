import os
import contextlib
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

class _DevNull:
    def write(self, *args, **kwargs): pass
    def flush(self, *args, **kwargs): pass

def deseq2_differential_expression(
  controls_mat: pd.DataFrame,
  cases_mat: pd.DataFrame,
  n_cpus=os.cpu_count() or 1,
  remove_na=True,
  sorted=True,
  stdout=_DevNull(),
):
    ''' Use pydeseq2 for differential expression.

    Note that this function expects the original raw gene counts.

    :param controls_mat: (pd.DataFrame) the control samples (samples as columns and genes as rows)
    :param cases_mat: (pd.DataFrame) the case samples (samples as columns and genes as rows)
    :param n_cpus: (int) number of CPUs used (default: number of cpus)
    :param remove_na: (bool) remove genes with NAN values (default: True)
    :param sorted: (bool) sort genes from most significant to least significant (default: True)
    :param stdout: (writeable stream) direct deseq's output, e.g. sys.stdout (default: suppress)
    :return: A data frame with the results
    '''
    # Check if controls_mat and cases_mat have the same number of rows
    if controls_mat.shape[0] != cases_mat.shape[0]:
        raise ValueError("controls_mat and cases_mat must have the same number of rows.")
    if (controls_mat.shape[1] < 2) | (cases_mat.shape[1] < 2):
        raise ValueError("controls_mat and cases_mat must have at least two samples.")
    with contextlib.redirect_stdout(stdout):
        exp = pd.concat([controls_mat, cases_mat], axis=1)
        condition_labels = ['C'] * controls_mat.shape[1] + ['RS'] * cases_mat.shape[1]
        sample_names = controls_mat.columns.tolist() + cases_mat.columns.tolist()
        metadata = pd.DataFrame({'Sample': sample_names, 'Condition': condition_labels}).set_index("Sample")
        dds = DeseqDataSet(counts=exp.T, clinical=metadata, design_factors="Condition")
        dds.deseq2()
        stat_res = DeseqStats(dds, n_cpus=n_cpus, contrast = ('Condition', 'RS', 'C'))
        stat_res.summary()
        if sorted:
            stat_res.results_df = stat_res.results_df.sort_values("pvalue")
    if remove_na:
        return stat_res.results_df.dropna()
    else:
        return stat_res.results_df
