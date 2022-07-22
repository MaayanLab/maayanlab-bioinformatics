import uuid

class IDMapper:
  ''' Stores id mappings and makes it easy to use many of them in tandem.

  ```python
  mapper = IDMapper()

  mapper.update({ 'a': {'A', 'C'} }, namespace='source_1')
  mapper.update({ 'b': {'A', 'B'} }, namespace='source_2')
  mapper.get('C', namespace='source_2') == 'b'

  # Because of the overlap in synonyms it is inferred that source_1's 'a' and source_2's 'b' correspond to the same
  #  id, we can get using any of the synyonyms to retreive the id in a given namespace.
  ```
  '''
  def __init__(self):
    self._forward = {}
    self._reverse = {}
    self._namespaces = {}
    self._counts = {}

  def summary(self):
    ''' Return summaries about all namespaces, and specifically id overlap counts between them.
    '''
    summary = {}
    for synonyms in self._forward.values():
      ns = frozenset(
        ns
        for synonym in synonyms
        for ns in self._namespaces[synonym]
      )
      if ns not in summary:
        summary[ns] = 0
      summary[ns] += 1
    return summary
    
  def get(self, term, namespace=None):
    ''' Given a term corresponding to some synonym, find the corresponding identifier in the provided namespace
    (or in all namespaces if not provided).
    '''
    id = self._reverse.get(term)
    if id is None: return None
    refs = {
      namespace: keys
      for synonym in self._forward[id]
      for namespace, keys in self._namespaces[synonym].items()
    }
    return dict(
      id=id,
      refs=refs,
      synonyms=self._forward[id],
    ) if namespace is None else refs.get(namespace)

  def update(self, mappings, namespace=None):
    ''' Add mappings of the form for a specific namespace:
    { identifier: { synonyms, ... }, ... }
    '''
    for key, synonyms in mappings.items():
      id = uuid.uuid4()
      self._forward[id] = set()
      for synonym in (key, *synonyms):
        if synonym not in self._reverse:
          self._forward[id].add(synonym)
          self._reverse[synonym] = id
          if synonym not in self._namespaces:
            self._namespaces[synonym] = {}
          if namespace not in self._namespaces[synonym]:
            self._namespaces[synonym][namespace] = set()
          self._namespaces[synonym][namespace].add(key)
        else:
          orig_id = self._reverse[synonym]
          if orig_id != id:
            for s in self._forward[id]:
              self._forward[orig_id].add(s)
              self._reverse[s] = orig_id
            if synonym not in self._namespaces:
              self._namespaces[synonym] = {}
            if namespace not in self._namespaces[synonym]:
              self._namespaces[synonym][namespace] = set()
            self._namespaces[synonym][namespace].add(key)
            del self._forward[id]
            id = orig_id
