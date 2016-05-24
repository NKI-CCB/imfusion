IM-Fusion
=========

IM-Fusion is a tool for identifying transposon insertion sites in insertional mutagenesis screens using single- and paired-end RNA-sequencing data. It essentially identifies insertion sites from gene-transposon fusions in the RNA-sequencing data, which represent splicing events between the transposon and endogeneous genes.

IM-Fusion also identifies candidate genes for a given screen using a
statistical test (based on the Poisson distribution) that selects Commonly
Targeted Genes (CTGs) – genes that are more frequently affected by insertions
than would be expected by chance. To futher narrow down a list of CTGs, which
may contain hundreds of genes, IM-Fusion also tests if insertions in a CTG have
a significant effect on the expression of the gene, which is a strong indicator
of them having an actual biological effect.

IM-Fusion has the following key features:

- It identifies transposon insertion sites from both single- and paired-end
  RNA-sequencing data, without having any special sequencing requirements.
- Its gene-centric approach – both in the identification of insertions and in
  testing for differential expression of identified candidate genes – greatly reduces the number of false positive candidate genes.
- Its differential expression tests provide detailed insight into the effects
  of insertions on the expression of their target gene(s). By providing both a
  groupwise and a single-sample version of the test, IM-Fusion can identify
  effects for a single insertion in a specific sample, or determine the
  general effect of insertions on a given gene within the tumor cohort.

For more details on the approach and a comparison with existing DNA-sequencing
approaches, please see our preprint "Identifying transposon insertions and
their effects from RNA-sequencing data" on bioRxiv.


Documentation
=============

IM-Fusion's documentation is available on ReadTheDocs.


Installation
============

IM-Fusion
---------

Using Pip
~~~~~~~~~

IM-Fusion can be installed using pip:

.. code:: bash

    pip install git+https://github.com/jrderuiter/im-fusion.git

The single-sample differential expression test requires rpy2, which is not
installed by default as it can be complicated to install. To use the test,
ensure a working version of rpy2 is installed or install the required
rpy2 dependencies using:

.. code:: bash

    pip install git+https://github.com/jrderuiter/im-fusion.git[de_single]

To install the development dependencies of IM-Fusion, use:

.. code:: bash

    pip install git+https://github.com/jrderuiter/im-fusion.git[dev]

Using Conda for dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although we don't yet provide a conda package for IM-Fusion, some of
IM-Fusion's dependencies may be easier to install using conda, which
provides pre-compiled binaries of many Python packages. Using conda,
IM-Fusion can be installed in an environment named 'imfusion' as follows:

.. code:: bash

    conda create -n imfusion pip future toolz pandas numexpr scipy seaborn
    pip install git+https://github.com/jrderuiter/im-fusion.git

Alternatively, you can use the provided *environment.yml* file to
create the environment:

.. code:: bash

   conda create -f ./conda/environment.yml

To include the development dependencies use:

.. code:: bash

   conda create -f ./conda/environment-dev.yml

Python dependencies
-------------------

IM-Fusion is compatible with both Python 2.7 and Python 3.3+. However,
identifying insertions unforuntately only works in Python 2.7
as Tophat2 currently does not run under Python 3.

IM-Fusion was tested using the following versions for its dependencies:

- future -- 0.15.2 or higher
- toolz -- 0.7.4 or higher
- pandas -- 0.18 or higher
- numexpr -- 2.5 or higher
- scipy -- 0.17 or higher
- seaborn -- 0.7 or higher
- pysam -- 0.8.4
- intervaltree -- 2.1.0

The following additional dependencies are required for running tests
and building the docs:

- pytest
- pytest-mock
- pytest-datafiles
- pytest-helpers-namespace
- sphinx

External dependencies
---------------------

IM-Fusion has two external dependencies:

- Tophat2 -- for the alignment and detection of gene fusions
- featureCounts -- for generating expression counts

Bothof these tools must be available in the current PATH. For help with
their installation, see below.

TopHat2
~~~~~~~

Tophat2 is a RNA-seq read aligner that also includes Tophat-Fusion,
an aligner that can identify gene fusions from single- and paired-end
RNA-sequencing data. A binary version of Tophat2 can be downloaded from the
Tophat `website <https://ccb.jhu.edu/software/tophat/index.shtml>`_.
After unpacking, make sure that the binary *Tophat2* is available in PATH.

Alternatively, Tophat2 can also be installed using `homebrew <http://brew.sh>`_
or `linuxbrew <http://linuxbrew.sh>`_. To do so, first install homebrew or
linuxbrew. Then install TopHat2 using:

.. code:: bash

    brew update
    brew install tophat

featureCounts
~~~~~~~~~~~~~

featureCounts is an efficient tool for generating expression counts from
RNA-seq alignments. It is part of the subread package, which can be downloaded
from `sourceforge <http://subread.sourceforge.net>`_. After downloading and
unpacking the correct binary version for your machine, make sure the binary
featureCounts is available in PATH.

Subread can also be installed using `homebrew <http://brew.sh>`_
or `linuxbrew <http://linuxbrew.sh>`_ using a custom tap:

.. code:: bash

    brew tap jrderuiter/tap
    brew update
    brew install subread

License
=======

This software is released under the MIT license.
