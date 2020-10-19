def display_clustergrammer(net):
  ''' This function displays clustergrammer in a jupyter notebook without dependencies
  on ipywidgets or any locally installed jupyter extensions. This is convenient for
  static exports, colab, and appyters.

  Example:
  ```python
  from maayanlab_bioinformatics.plotting import display_clustergrammer
  from clustergrammer import Network
  net = Network()
  net.load_df(df)
  net.cluster()
  display_clustergrammer(net)
  ```
  '''
  from IPython.display import HTML
  import uuid, json
  id = '_' + str(uuid.uuid4())
  return HTML(f"""
    <div id='{id}' style="width: 950px; height: 800px"></div>
    <script src="/static/components/requirejs/require.js"></script>
    <script>
    requirejs.config({json.dumps({
      'paths': {
        'jquery': 'https://cdnjs.cloudflare.com/ajax/libs/jquery/1.11.2/jquery.min',
        'd3': 'https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.15/d3.min',
        f"clustergrammer": 'https://raw.githack.com/MaayanLab/clustergrammer/0024d8cd245dc597113a860db9f1dc989a8876c2/clustergrammer',
      },
      'shim': {
        f"clustergrammer": {
          'exports': 'Clustergrammer',
          'deps': ['d3', 'jquery'],
        }
      }
    })})
    require(['clustergrammer'], function (clustergrammer) {{
      clustergrammer({json.dumps({
        'root': f"#{id}",
        'network_data': json.loads(net.export_net_json()),
      })})
    }})
    </script>
  """)
