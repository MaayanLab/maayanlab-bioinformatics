import numpy as np
import pandas as pd

def bridge_plot(select: pd.Series, weights: pd.Series = None):
  ''' Use the filter to construct a bridge plot.

  ```python
  import numpy as np
  from matplotlib import pyplot as plt
  from maayanlab_bioinformatics.plotting import bridge_plot

  x, y = bridge_plot(select)
  plt.plot(x, y)
  plt.vlines(np.argwhere(select.values)[:, 0], ymin=-1, ymax=0)
  plt.show()
  ```

  :param select: (pd.Series) selection of hits (i.e. `df['gene'] == 'my_target'`) in ranked order
  :param weights: (pd.Series) optional weights for each hit in the same order
  :return: (Tuple[np.array, np.array]) x and y arrays ready to be plotted.
  '''
  if weights is None:
    weights = pd.Series(np.ones(select.shape[0]), index=select.index)
  max_es = weights[select].abs().sum() # maximum enrichment score if we were to hit everything (positively)
  up = select * weights / max_es # go up/dn by normalized weight on each hit
  dn = - (1 - select) * up.sum() / (~select).sum()
  x = np.arange(select.shape[0]+1)
  y = np.concatenate([
    np.zeros(1),
    np.cumsum(up + dn),
  ])
  return x, y
