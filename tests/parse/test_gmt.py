import os

def test_gmt_read_dict():
  from maayanlab_bioinformatics.parse.gmt import gmt_read_dict
  gmt = gmt_read_dict(open(os.path.join(os.path.dirname(__file__), '..', 'test_gmt.gmt'), 'r'))
  assert 'ChIP-seq\ttest desc' in gmt
  assert gmt['ChIP-seq\ttest desc']['LRRC37A3'] == 2.0
  assert 'LRRC37A3' in gmt['ChIP-seq\ttest desc']
  assert 'TLE3' in gmt['ChIP-seq\ttest desc']
  assert gmt['ChIP-seq\ttest desc']['TLE3'] == 5.0
  assert 'Data aggregation' in gmt
  assert 'CD24L4' in gmt['Data aggregation']
  assert gmt['Data aggregation']['CD24L4'] == 1

def test_gmt_read_pd():
  from maayanlab_bioinformatics.parse.gmt import gmt_read_pd
  gmt = gmt_read_pd(os.path.join(os.path.dirname(__file__), '..', 'test_gmt.gmt'))
  assert 'ChIP-seq\ttest desc' in gmt.columns
  assert 'LRRC37A3' in gmt.index
  assert 'TLE3' in gmt.index
  assert gmt.loc['LRRC37A3', 'ChIP-seq\ttest desc'] == 2.0
  assert gmt.loc['TLE3', 'ChIP-seq\ttest desc'] == 5.0
  assert 'Data aggregation' in gmt.columns
  assert 'CD24L4' in gmt.index
  assert gmt.loc['CD24L4', 'Data aggregation'] == 1
