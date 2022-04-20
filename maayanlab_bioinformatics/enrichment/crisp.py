# import fisher
import scipy.stats
from typing import Union, Dict, Set, Iterable, Tuple, Hashable, Any, TypeVar, Optional
from dataclasses import dataclass

@dataclass(frozen=True)
class FisherOverlap:
  pvalue: float
  odds_ratio: float
  n_overlap: int
  overlap: Optional[Set[Hashable]]

T = TypeVar('T')
DictOrIterableTuple = Union[Dict[Hashable, T], Iterable[Tuple[Hashable, T]]]
CompatibleSignature = Union[DictOrIterableTuple[Any], Set[Hashable]]
CompatibleSignatures = DictOrIterableTuple[CompatibleSignature]
EnrichmentResult = Iterable[Tuple[Hashable, FisherOverlap]]

def _dict_or_iterable_tuple(it: DictOrIterableTuple[T]) -> Iterable[Tuple[Hashable, T]]:
  if callable(getattr(it, 'items', None)):
    return it.items()
  else:
    return it

def safe_odds_ratio(a, b, c, d):
  ''' Compute the odds ratio returning helpful answers in the case of division by zero issues..
  '''
  # numerator
  if a == 0 and c == 0:
    ac = float('nan')
  elif c == 0: # a != 0
    ac = float('inf')
  else:
    ac = float(a / c)
  # denominator
  if b == 0 and d == 0:
    bd = float('nan')
  elif d == 0: # b != 0
    bd = float('inf')
  else:
    bd = float(b / d)
  # odds ratio (numerator / denominator)
  if ac == float('nan') or bd == float('nan'):
    # not going to bother.. this would only happen if you had empty signatures
    return float('nan')
  elif ac == float('inf') and bd == float('inf'):
    # this would mean *everything* is in the input set..
    # inf probably makes sense given that the occurrence
    # of the event would be *guaranteed* in this case
    return float('inf')
  elif ac == float('inf'): # bd != float('inf')
    # inf / number = inf
    return float('inf')
  elif bd == float('inf'): # ac != float('inf')
    # number / inf = 0
    return 0.0
  elif bd == 0:
    return float('inf')
  else:
    return ac / bd

def fisher_overlap(
  input_signature: Set[Hashable],
  background_signature: Set[Hashable],
  n_background_entities: int,
  preserve_overlap: bool = False,
) -> Optional[FisherOverlap]:
  ''' Given input and background set, compute the overlap, fisher significance, and odds ratio.
  In the case of no overlap, will return None.
  '''
  overlap = input_signature & background_signature
  n_overlap = len(overlap)
  n_input_signature = len(input_signature)
  n_background_signature = len(background_signature)
  if n_overlap == 0:
    return None
  #
  a = n_overlap
  b = n_input_signature - n_overlap
  c = n_background_signature - n_overlap
  d = n_background_entities - n_background_signature - n_input_signature + n_overlap
  if d < 0:
    raise Exception('The total population cannot be smaller than the current overlap..')
  #
  # pvalue = fisher.pvalue(a, b, c, d).right_tail
  pvalue = scipy.stats.fisher_exact([[a, b], [c, d]], 'greater')[1]
  odds_ratio = safe_odds_ratio(a, b, c, d)
  #
  return FisherOverlap(
    pvalue=pvalue,
    odds_ratio=odds_ratio,
    n_overlap=n_overlap,
    overlap=overlap if preserve_overlap else None,
  )

def enrich_crisp(
  input_signature: CompatibleSignature,
  background_signatures: CompatibleSignatures,
  n_background_entities: int,
  preserve_overlap: bool = False,
) -> Iterable[Tuple[Hashable, FisherOverlap]]:
  ''' Perform crisp set enrichment analysis using fisher overlap.
  Eriches the signature in input_signature against signatures in background_signatures.

  :param n_background_entities: should correspond to the approximate number of entities exist, in the case of Human Genes for instance this might be 21000.
  '''
  input_signature = set(input_signature)
  for background_signature_term, background_signature in _dict_or_iterable_tuple(background_signatures):
    background_signature = set(background_signature)
    result = fisher_overlap(
      input_signature,
      background_signature,
      n_background_entities=n_background_entities,
      preserve_overlap=preserve_overlap,
    )
    if result is not None:
      yield background_signature_term, result
