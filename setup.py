from setuptools import setup, find_packages
from os import path


this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), 'r', encoding='utf-8') as f:
    long_description = f.read()

DESCRIPTION = 'Tools to filter sam o bam files by percent identity or percent of matched sequence',
LONG_DESCRIPTION = long_description,
LONG_DESCRIPTION_CONTENT_TYPE = 'text/markdown'
NAME = 'filtersam'
AUTHOR = "Semidán Robaina Estévez, 2021-2022"
AUTHOR_EMAIL = "srobaina@ull.edu.es"
MAINTAINER = "Semidán Robaina Estévez"
MAINTAINER_EMAIL = "srobaina@gmail.com"
DOWNLOAD_URL = 'http://github.com/robaina/filterSAM'
LICENSE = 'Creative Commons Attribution 4.0 International'
VERSION = '0.0.9'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=DOWNLOAD_URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=find_packages(),
      install_requires=['numpy', 'pysam', 'parallelbam'],
      entry_points ={
            'console_scripts': [
                'filtersam = filtersam.cli:main'
            ]
        }
      )
