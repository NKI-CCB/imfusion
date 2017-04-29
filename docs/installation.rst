Installation
============

Using bioconda (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The recommended way to install IM-Fusion is using bioconda, as this installs
IM-Fusion together with its external dependencies into an isolated environment
using a single command:

.. code:: bash

  conda create -n imfusion imfusion

Note that this assumes that condas channels have been configured as
recommended_ by bioconda.

.. _recommended: https://bioconda.github.io/#set-up-channels

Alternatively, IM-Fusion can be installed in an existing environment using:

.. code:: bash

  conda install imfusion

Note that the bioconda package does not install rpy2, as R/rpy2 can be
difficult to install. If you wish to use the single-sample differential
expression test (which requires rpy2), you can install R and rpy2 can be
installed using:

.. code:: bash

    conda install rpy2

Rpy2 can also be installed using pip if you want to use an existing R
installation:

.. code:: bash

    pip install rpy2

Using pip
~~~~~~~~~

IM-Fusion can also be installed using pip:

.. code:: bash

    pip install git+https://github.com/jrderuiter/imfusion.git

Note that this does not install any of the required external dependencies,
which must therefore be installed manually (see below for more details).
This does also not install rpy2, which must also be installed using pip
if you plan to use single-sample differential expression test.

Hardware dependencies
~~~~~~~~~~~~~~~~~~~~~

To build references and align reads using the STAR aligner, you will need to
have ~30GB RAM. For analyses with Tophat2, you will typically need ~5-15GB of
RAM, depending on the number of reads in your sample.

Software dependencies
~~~~~~~~~~~~~~~~~~~~~

IM-Fusion has been tested on macOS and linux, and is compatible with both
Python 2.7 and Python 3.4+. However, identifying insertions using Tophat-Fusion
only works in Python 2.7, as Tophat2 is unfortunately not compatible with
Python 3.

IM-Fusion has the following external dependencies that must be installed
for full functionality:

- STAR_ - for insertion detection using STAR.
- Tophat2_/Bowtie_ (Python 2.7 only) -- for insertion detection using
  Tophat-Fusion.
- featureCounts_ -- for generating expression counts.
- StringTie_ (optional) - for performing reference-guided transcript assembly.
- STAR-Fusion_ (optional) -- for identifying endogenous gene fusions using STAR.
- Sambamba_ (optional) -- for faster sorting of BAM files.

Note that either STAR or Tophat2 needs to be installed, depending on your
preference for the RNA-seq aligner. If you are unsure, we recommend using STAR.
StringTie is only required if you plan to perform reference-guided assembly
(used to detect novel transcripts resulting from insertions). Similarly,
STAR-Fusion is only needed if you plan to do endogenous fusion detection.
Sambamba is used for faster sorting of BAM files in the STAR pipeline, but is
also not a strict requirement.

.. _Bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _Tophat2: https://ccb.jhu.edu/software/tophat/index.shtml
.. _StringTie: https://ccb.jhu.edu/software/stringtie
.. _featureCounts: http://subread.sourceforge.net
.. _STAR: https://github.com/alexdobin/STAR
.. _STAR-Fusion: https://github.com/STAR-Fusion/STAR-Fusion
.. _Sambamba: http://lomereiter.github.io/sambamba

Each of these packages can be installed manually by downloading binaries from
the respective websites and making sure the unpacked binaries are available in
the current ``$PATH``. Alternatively, the dependencies can easily be installed
using (bio)conda:

.. code:: bash

    conda install -c bioconda star bowtie tophat stringtie \
        subread star-fusion sambamba

Note that Tophat will not function properly in Python 3. Bioconda also does not
provide any macOS builds of STAR-Fusion, meaning that STAR-Fusion will have to
be installed manually on macOS.

If you don't use conda, a number of these tools can also be installed using
`homebrew <http://brew.sh>`_ or `linuxbrew <http://linuxbrew.sh>`_:

.. code:: bash

    brew tap homebrew/science
    brew update
    brew install bowtie star stringtie tophat sambamba

Unfortunately, STAR-Fusion and subread (featureCounts) are currently not
available in homebrew and will therefore need to be installed manually.
