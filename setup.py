from setuptools import setup
from distutils.extension import Extension
import numpy as np

import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


if __name__ == '__main__':
    setup(
        name='scptseq',
        packages=['scptseq'],
        version=get_version("scptseq/__init__.py"),
        entry_points={
          'console_scripts': [
              'scptseq = scptseq.main:main'
          ]
        },
        include_package_data=True,
        include_dirs=[np.get_include()],
        install_requires=[
            "numpy>=1.20.0",
            "docopt>=0.6.2",
            "biopython>=1.79",
            "matplotlib>=3.5.2",
            ],
        zip_safe=False,
        author='John Hawkins',
        author_email='hawkjo@gmail.com',
        description='scPT-seq preprocessing package',
        url='https://github.com/hawkjo/scptseq',
        download_url='',
        keywords=['DNA', 'RNA', 'CRISPR', 'NGS', 'bioinformatics'],
        python_requires='>=3.0',
        classifiers=['Development Status :: 3 - Alpha',
                     'Natural Language :: English',
                     'Intended Audience :: Science/Research',
                     'Operating System :: POSIX :: Linux',
                     'Programming Language :: Python :: 3',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     ]
    )
