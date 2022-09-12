import uuid
from collections import Counter

class IDMapper:
  ''' Stores id mappings and makes it easy to use many of them in tandem.

  ```python
  mapper = IDMapper()

  mapper.update({ 'a': {'A', 'C'} }, namespace='source_1')
  mapper.update({ 'b': {'A', 'B'} }, namespace='source_2')
  mapper.get('C', namespace='source_2') == 'b'

  Because of the overlap in synonyms it is inferred that source_1's 'a' and source_2's 'b' correspond to the same
    id, we can get using any of the synyonyms to retreive the id in a given namespace.
  Since this can be problematic when synonyms are malformed, mapper.conflicts_summary() and mapper.conflicts_counts()
    provide ways of debugging excess synonym applications.
  ```
  '''
  def __init__(self):
    # { uuid1: {id1: 1, id2: 1, ...} }
    self._forward = {}
    # { id1: uuid1, id2, uuid1, ... }
    self._reverse = {}
    # { uuid1: { ns1: id1 }, ... }
    self._namespaces = {}
    # { ns1: { shared_synonym: { conflictid1: origid1 }, ... } } }
    self._conflicts = {}

  def summary(self):
    ''' Return counts of overlapping namespaces (like a venn diagram)
    '''
    return Counter(
      frozenset(ns_ids.keys())
      for ns_ids in self._namespaces.values()
    )
    
  def conflicts_summary(self):
    ''' Return counts of conflicts in each namespace
    '''
    return Counter({
      ns: len(conflicts)
      for ns, conflicts in self._conflicts.items()
    })

  def top_conflicts(self):
    ''' Return conflicting synonym counts
    '''
    return Counter({
      (ns, conflict): len(cases)
      for ns, cc in self._conflicts.items()
      for conflict, cases in cc.items()
    })

  def get_id(self, id, namespace=None):
    if id is None: return None
    if namespace is None:
      return dict(
        id=id,
        refs=self._namespaces[id],
        synonyms=self._forward[id],
      )
    else:
      return self._namespaces[id].get(namespace)

  def get(self, term, namespace=None):
    id = self._reverse.get(term)
    return self.get_id(id, namespace=namespace)

  def find(self, term):
    potential_ids = {
      id
      for k, id in self._reverse.items()
      if str(term).lower().strip() in str(k).lower().strip() or str(k).lower().strip() in str(term).lower().strip()
    }
    return {
      id: self.get_id(id)
      for id in potential_ids
    }

  def update(self, mappings, namespace=None):
    ''' Add mappings of the form:
    { identifier: { synonyms } }
    '''
    for key, synonyms in (mappings.items() if type(mappings) == dict else mappings):
      id = uuid.uuid4()
      self._forward[id] = Counter()
      self._namespaces[id] = {namespace: key}
      for synonym in {key, *synonyms}:
        if synonym not in self._reverse:
          self._forward[id].update([synonym])
          self._reverse[synonym] = id
        else:
          orig_id = self._reverse[synonym]
          if orig_id == id:
            self._forward[id].update([synonym])
          else:
            for ns, k in self._namespaces.pop(id, {}).items():
              if orig_id not in self._namespaces:
                self._namespaces[orig_id] = {}
              orig_k = self._namespaces[orig_id].get(ns)
              if orig_k is not None:
                if orig_k != k:
                  if ns not in self._conflicts:
                    self._conflicts[ns] = {}
                  if synonym not in self._conflicts[ns]:
                    self._conflicts[ns][synonym] = {}
                  self._conflicts[ns][synonym][k] = orig_k
              else:
                self._namespaces[orig_id][ns] = k
            new_cnt = self._forward.pop(id)
            self._forward[orig_id] += new_cnt
            self._reverse.update({s: orig_id for s in new_cnt.keys()})
            id = orig_id
