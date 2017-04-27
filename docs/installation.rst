Installation
============

Using conda (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~

The recommended way to install IM-Fusion is using conda, as this allows us to
install IM-Fusion together with its external dependencies (Tophat2/Bowtie,
STAR and Subread) into an isolated environment using a single command:

.. code:: bash

  conda create -n imfusion -c bioconda -c r -c jrderuiter imfusion

Note that this uses the bioconda channel to install the external
dependencies (Tophat2/STAR, Bowtie and Subread) and the r channel to install R
and rpy2.

Alternatively, IM-Fusion can be installed in an existing environment using:

.. code:: bash

  conda install -c bioconda -c r -c jrderuiter imfusion

Using pip
~~~~~~~~~

IM-Fusion can also be installed using pip:

.. code:: bash

    pip install git+https://github.com/jrderuiter/imfusion.git

Note that in this case, rpy2 is not installed by default as it can be
difficult to install and is only required for the single-sample differential
expression test. To use the test, ensure a working version of rpy2 is
installed or install the required rpy2 dependencies using:

.. code:: bash

    pip install git+https://github.com/jrderuiter/imfusion.git[de_single]

To also install the development dependencies of IM-Fusion use:

.. code:: bash

    pip install git+https://github.com/jrderuiter/imfusion.git[dev]

Dependencies
~~~~~~~~~~~~

Python dependencies
-------------------

IM-Fusion is compatible with both Python 2.7 and Python 3.4+. However,
identifying insertions using Tophat-Fusion only works in Python 2.7, as
Tophat2 is not compatible with Python 3.

External dependencies
---------------------

IM-Fusion has three external dependencies:

- STAR - for the alignment and detection of gene fusions.
- Tophat2/Bowtie1 -- for the alignment and detection of gene fusions.
- Stringtie - For performing reference-guided transcript assembly.
- featureCounts -- for generating expression counts.

These tools must be available in the current PATH for the insertion
detection to work properly. Note that either STAR or Tophat2 needs to be
installed, depending on your preference for the aligner. For help with
their installation, see below.

STAR
****

STAR is an RNA-seq aligner that can also perform chimeric read alignment,
enabling the detection of gene-transposon fusions. STAR can be downloaded
from the `STAR GitHub repository`_. After unpacking, make sure that the
STAR binary is available in ``$PATH``.

STAR can also be installed using bioconda:

.. code:: bash

    conda install -c bioconda star

Or using homebrew/linuxbrew:

.. code:: bash

    brew tap homebrew/science
    brew update
    brew install star

.. _`STAR GitHub repository`: https://github.com/alexdobin/STAR


StringTie
*********

StringTie is a fast and highly efficient assembler of RNA-Seq alignments into
potential transcripts. A binary version of StringTie can be downloaded from the
StringTie `website <https://ccb.jhu.edu/software/stringtie>`_. After unpacking,
make sure that the stringtie binary is available in ``$PATH``.

STAR can also be installed using bioconda:

.. code:: bash

    conda install -c bioconda stringtie

Or using homebrew/linuxbrew:

.. code:: bash

    brew tap homebrew/science
    brew update
    brew install stringtie

Tophat2/Bowtie1
***************

Tophat2 is an RNA-seq read aligner that also includes Tophat-Fusion,
an aligner that can identify gene fusions from single- and paired-end
RNA-sequencing data. A binary version of Tophat2 can be downloaded from the
Tophat `website <https://ccb.jhu.edu/software/tophat/index.shtml>`_.
After unpacking, make sure that the Tophat2 binary is available in ``$PATH``.

Alternatively, Tophat2 can also be installed using bioconda:

.. code:: bash

    conda install -c bioconda tophat

Tophat2 can also be installed using `homebrew <http://brew.sh>`_ or
`linuxbrew <http://linuxbrew.sh>`_. For this option, first install
homebrew or linuxbrew and then install Tophat2 using:

.. code:: bash

    brew tap homebrew/science
    brew update
    brew install tophat

Bowtie can be installed in the same fashion as Tophat, by downloading a binary
version from the `bowtie website
<http://bowtie-bio.sourceforge.net/index.shtml>`_ and ensuring that bowtie
is available in ``$PATH``.

Bowtie can also be installed using bioconda:

.. code:: bash

    conda install -c bioconda bowtie

Or using homebrew/linuxbrew:

.. code:: bash

    brew tap homebrew/science
    brew update
    brew install bowtie

featureCounts
*************

featureCounts is an efficient tool for generating expression counts from
RNA-seq alignments. It is part of the subread package, which can be downloaded
from `sourceforge <http://subread.sourceforge.net>`_. After downloading and
unpacking the correct binary version for your machine, make sure the binary
featureCounts is available in PATH.

Subread can also be installed using bioconda:

.. code:: bash

    conda install -c bioconda subread

Or by using `homebrew <http://brew.sh>`_ or `linuxbrew <http://linuxbrew.sh>`_
with a custom tap:

.. code:: bash

    brew tap jrderuiter/tap
    brew update
    brew install subread
