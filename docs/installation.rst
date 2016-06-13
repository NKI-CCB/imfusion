Installation
============

Using conda (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~

The recommended way to install IM-Fusion is using conda, as this allows us to
install IM-Fusion together with its external dependencies (tophat2, bowtie and
subread) into an isolated environment using a single command:

.. code:: bash

  conda create -n imfusion -c jrderuiter -c bioconda -c r im-fusion

Note that this uses the bioconda channel to install the external
dependencies (tophat2, bowtie and subread) and the r channel to install R
and rpy2.

Alternatively, IM-Fusion can be installed in an existing environent using:

.. code:: bash

  conda install -c jrderuiter -c bioconda -c r im-fusion

Using pip
~~~~~~~~~

IM-Fusion can also be installed using pip:

.. code:: bash

    pip install git+https://github.com/jrderuiter/im-fusion.git

Note that in this case, rpy2 is not installed by default as it can be
difficult to install and is only required for the single-sample differential
expression test. To use the test, ensure a working version of rpy2 is installed or install the required rpy2 dependencies using:

.. code:: bash

    pip install git+https://github.com/jrderuiter/im-fusion.git[de_single]

To also install the development dependencies of IM-Fusion use:

.. code:: bash

    pip install git+https://github.com/jrderuiter/im-fusion.git[dev]

Dependencies
~~~~~~~~~~~~

Python dependencies
-------------------

IM-Fusion is compatible with both Python 2.7 and Python 3.3+. However,
identifying insertions only works in Python 3.3+ if a Python
3 compatible version of Tophat2 is installed. Such a version can be obtained
from bioconda or by applying the 2to3 conversion to Tophats scripts.

IM-Fusion was tested using the following versions for its dependencies:

- future -- 0.15.2
- toolz -- 0.7.4
- pandas -- 0.18
- numexpr -- 2.5
- scipy -- 0.17
- seaborn -- 0.7
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

IM-Fusion has three external dependencies:

- Tophat2 -- for the alignment and detection of gene fusions
- Bowtie1 -- used by Tophat2 for the alignment
- featureCounts -- for generating expression counts

These tools must be available in the current PATH for the insertion
detection to work properly. For help with their installation, see below.

Tophat2
*******

Tophat2 is a RNA-seq read aligner that also includes Tophat-Fusion,
an aligner that can identify gene fusions from single- and paired-end
RNA-sequencing data. A binary version of Tophat2 can be downloaded from the
Tophat `website <https://ccb.jhu.edu/software/tophat/index.shtml>`_.
After unpacking, make sure that the binary *Tophat2* is available in PATH.

Alternatively, Tophat2 can also be installed using bioconda, or using
`homebrew <http://brew.sh>`_ or `linuxbrew <http://linuxbrew.sh>`_. For
homebrew/linuxbrew, first install homebrew or linuxbrew and then install
Tophat2 using:

.. code:: bash

    brew update
    brew install tophat

Bowtie
******

Bowtie can be installed in the same fashion as Tophat, by downloading a binary
version from the `bowtie website
<http://bowtie-bio.sourceforge.net/index.shtml>`_ and ensuring that bowtie
is available in PATH. Bowtie can also be installed using bioconda or using
homebrew/linuxbrew:

.. code:: bash

    brew update
    brew install bowtie

featureCounts
*************

featureCounts is an efficient tool for generating expression counts from
RNA-seq alignments. It is part of the subread package, which can be downloaded
from `sourceforge <http://subread.sourceforge.net>`_. After downloading and
unpacking the correct binary version for your machine, make sure the binary
featureCounts is available in PATH.

Subread can also be installed using bioconda, or by using
`homebrew <http://brew.sh>`_ or `linuxbrew <http://linuxbrew.sh>`_
with a custom tap:

.. code:: bash

    brew tap jrderuiter/tap
    brew update
    brew install subread
