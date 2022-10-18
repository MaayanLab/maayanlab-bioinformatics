from setuptools import setup, find_packages

setup(
  name='maayanlab_bioinformatics',
  python_requires='>=3.8.0',
  version='0.5.2',
  url='https://github.com/maayanLab/maayanlab-bioinformatics/',
  author='Daniel J. B. Clarke',
  author_email='danieljbclarkemssm@gmail.com',
  long_description=open('README.md', 'r').read(),
  license='Apache-2.0',
  install_requires=list(map(str.strip, open('requirements.txt', 'r').readlines())),
  extras_require={
    'h5py': ["h5py"],
    'limma_voom': ["rpy2==3.3.6"],
    'progress': ["tqdm"],
    'enrichr_user_list': ['bs4', 'lxml'],
  },
  packages=find_packages(),
  include_package_data=True,
)
