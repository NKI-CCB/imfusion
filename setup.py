#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import setuptools

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

# General requirements.
requirements = ['future', 'pandas', 'numexpr', 'pysam', 'toolz', 'seaborn',
                'pyfaidx', 'scipy', 'intervaltree', 'frozendict']

extra_requirements = {
    'de_single': ['rpy2'],
    'dev': ['sphinx', 'pytest', 'pytest-mock', 'pytest-datafiles',
            'pytest-cov', 'pytest-helpers-namespace']
}

# Check setuptools version, as recommended by:
# https://hynek.me/articles/conditional-python-dependencies/.
if int(setuptools.__version__.split('.', 1)[0]) < 18:
    assert 'bdist_wheel' not in sys.argv

    # Add pathlib for Pythons before 3.4.
    if sys.version_info[0:2] < (3, 4):
        requirements.append('pathlib2')
else:
    extra_requirements[":python_version<'3.4'"] = ['pathlib2']

setuptools.setup(
    name='imfusion',
    version='0.2.0.dev0',
    description=('Tool for identifying transposon insertions in '
                 'Insertional Mutagenesis screens from gene-transposon '
                 'fusions using single- and paired-end RNA-sequencing data.'),
    long_description=readme + '\n\n' + history,
    url='https://github.com/jrderuiter/im-fusion',
    author='Julian de Ruiter',
    author_email='julianderuiter@gmail.com',
    license='MIT license',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    entry_points={'console_scripts': [
        'imfusion-build = imfusion.main.build:main',
        'imfusion-insertions = imfusion.main.insertions:main',
        'imfusion-ctg = imfusion.main.ctg:main',
        'imfusion-expression = imfusion.main.expression:main',
        'imfusion-merge = imfusion.main.merge:main'
    ]},
    install_requires=requirements,
    extras_require=extra_requirements,
    zip_safe=False,
    classifiers=[])
