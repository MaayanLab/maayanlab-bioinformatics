import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def silhouette_analysis(mat: pd.DataFrame, min_clusters=2, max_clusters=25, metric='cosine', random_state=None, **kwargs):
  ''' Compute KMeans repeatedly on the matrix with different cluster
  values between min_clusters and max_clusters, compute the silhouette_score,
  and return the best kmeans model/predictions.
  '''
  silhouette_scores = {}
  best = None
  for n in range(min_clusters, max_clusters+1):
    km = KMeans(n_clusters=n, random_state=random_state)
    y_pred = km.fit_predict(mat.values)
    score = silhouette_score(mat.values, y_pred, metric='cosine')
    silhouette_scores[n] = score
    if best is None or score > best[0]:
      best = (score, km, y_pred)
  #
  silhouette_scores = pd.DataFrame([
    {'N Clusters': k, 'Silhouette Score': v}
    for k, v in silhouette_scores.items()
  ])
  #
  score, km, y_pred = best
  y_pred = pd.DataFrame({
    'Cluster': [
      'Cluster {c}'.format(c=c)
      for c in km.fit_predict(mat.values)
    ]
  }, index=mat.index)
  return type('SilhouetteAnalysis', tuple(), dict(
    silhouette_scores=silhouette_scores,
    best_score=score,
    best_km=km,
    best_preds=y_pred,
  ))
