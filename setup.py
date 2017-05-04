#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import setuptools

with open('README.rst') as readme_file:
    README = readme_file.read()

with open('HISTORY.rst') as history_file:
    HISTORY = history_file.read()

INSTALL_REQUIRES = [
    'future', 'pandas>=0.19.0', 'numexpr', 'pysam>=0.9.1', 'toolz', 'pyfaidx',
    'scipy', 'intervaltree', 'pathlib2', 'htseq>=0.7.2', 'matplotlib',
    'seaborn', 'typing; python_version < "3.5"'
]

EXTRAS_REQUIRE = {
    'de_single': ['rpy2'],
    'dev': [
        'sphinx', 'sphinx-autobuild', 'sphinx-rtd-theme', 'bumpversion',
        'pytest>=2.7', 'pytest-mock', 'pytest-helpers-namespace', 'pytest-cov',
        'python-coveralls', 'seaborn'
    ]
}

setuptools.setup(
    name='imfusion',
    version='0.3.0',
    description=('Tool for identifying transposon insertions in '
                 'Insertional Mutagenesis screens from gene-transposon '
                 'fusions using single- and paired-end RNA-sequencing data.'),
    long_description=README + '\n\n' + HISTORY,
    url='https://github.com/jrderuiter/im-fusion',
    author='Julian de Ruiter',
    author_email='julianderuiter@gmail.com',
    license='MIT license',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'imfusion = imfusion.main.dispatch:main',
            'imfusion-build = imfusion.main.build:main',
            'imfusion-insertions = imfusion.main.insertions:main',
            'imfusion-ctg = imfusion.main.ctg:main',
            'imfusion-expression = imfusion.main.expression:main',
            'imfusion-merge = imfusion.main.merge:main'
        ]
    },
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    zip_safe=False,
    classifiers=[])
