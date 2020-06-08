import itertools
import pandas as pd
from typing import Dict, Set, Hashable

def upset_from_dict_of_sets(inputs: Dict[Hashable, Set[Hashable]]):
  ''' Given a dictionary of sets, produce input ready for `upsetplot` python package

  We produce this input by computing set intersections of all relevant combinations
   of sets interacting with one another.

  Example:
  ```python
  import upsetplot
  from maayanlab_bioinformatics.plotting import upset_from_dict_of_sets
  upsetplot.plot(upset_from_dict_of_sets({
    'A': {'a', 'b', 'c'},
    'B': {'b', 'c', 'd'},
    'C': {'d', 'e', 'f'},
  }))
  ```
  :param inputs: (Dict[Hashable, Set[Hashable]]) Several named sets
  :return: (pd.DataFrame) in a form ready for `upsetplot.plot`
  '''
  sets = []
  for n in range(1, len(inputs)+1):
    if n == 1:
      it = [[k] for k in inputs.keys()]
    else:
      it = map(list, itertools.combinations(inputs.keys(), n))
    for V in it:
      size = len(inputs[V[0]] if n == 1 else set.intersection(*[inputs[v] for v in V]))
      if size > 0:
        sets.append(dict({vv: vv in V for vv in inputs.keys()}, size=size))
  return pd.DataFrame(sets).groupby(list(inputs.keys()))['size'].sum()
