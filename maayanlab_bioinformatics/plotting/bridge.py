import numpy as np
import pandas as pd

def bridge_plot(select: pd.Series):
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

  :param select: (pd.Series) selection of hits (i.e. `df['gene'] == 'my_target'`)
  :return: (Tuple[np.array, np.array]) x and y arrays ready to be plotted.
  '''
  up = (~select).sum() / select.shape[0]
  dn = -(1-up)
  x = np.arange(select.shape[0]+1)
  y = np.concatenate([
    np.zeros(1),
    np.cumsum(select.apply(lambda s: up if s else dn)),
  ])
  return x, y
