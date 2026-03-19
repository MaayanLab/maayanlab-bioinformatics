import numpy as np
import pandas as pd
import typing as t

def quantile_normalize(data: pd.DataFrame, axis: int = 1, target: t.Union[None, np.ndarray] = None, ncpus: int = 1):
  from qnorm import quantile_normalize as qnorm
  return qnorm(data, axis=axis, target=target, ncpus=ncpus)
