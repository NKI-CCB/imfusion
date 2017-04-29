======
Extras
======

Using STAR-fusion to identify endogenous fusions
------------------------------------------------

Motivation
~~~~~~~~~~

An important advantage of using RNA-sequencing to identify insertions, rather
than targeted DNA-sequencing, is that RNA-seq data provides more information
about other events that may be contributing to tumorigenesis. For example,
besides the gene-transposon insertions, we may expect somatic gene-fusions
(between endogenous genes) to occur in these tumors, which may also contribute
to tumorigenesis.

In the IM-Fusion manuscript, we demonstrate this idea for both the SB and the
B-ALL datasets by identifying endogenous fusions from STAR RNA-seq alignments
using `STAR-Fusion`_. Using this approach, we identified the engineered
*Etv6-Runx1* fusion in the B-ALL tumors and identified several
known oncogenic fusions in the SB tumors.

Here we provide a brief description of how to apply STAR-Fusion to identify
endogenous gene fusions. In the first approach, we describe how to run
STAR-Fusion as part of IM-Fusion (when STAR is being used to identify
insertions). In the second, we describe how to run STAR-Fusion from
scratch, which requires extra computation time and disk space but has the
advantage of allowing you to use other parameter settings (specifically geared
towards identifying endogenous fusions) for STAR during the alignment.

.. _STAR-Fusion: http://star-fusion.github.io

Installing STAR-Fusion
~~~~~~~~~~~~~~~~~~~~~~

STAR-Fusion can be downloaded from the `STAR-Fusion GitHub
repository`_. After unpacking, make sure that the ``STAR-Fusion`` binary
is available in ``$PATH``. Note that STAR-Fusion requires STAR and a few
non-standard `perl modules`_ to be installed.

Alternatively, STAR-Fusion can also be installed using bioconda:

.. code:: bash

    conda install -c bioconda star-fusion

This should also install any dependencies of STAR-Fusion (such as STAR). Note
that the bioconda installation only works on linux at the moment, as bioconda
currently does not provide any STAR-Fusion builds for macOS.

.. _STAR-Fusion GitHub repository: https://github.com/STAR-Fusion/STAR-Fusion/releases
.. _perl modules: https://github.com/STAR-Fusion/STAR-Fusion/wiki#tools-required

Preparing the reference
~~~~~~~~~~~~~~~~~~~~~~~

Building a STAR-Fusion reference is quite involved as it
requires a number of additional tools and resources to be installed. We
therefore recommend using a pre-built reference if possible. Pre-built
references are available for human_ reference genomes and for the mouse mm10_
reference genome (using M9/Gencode).

.. _human: https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB
.. _mm10: ftp://ftp.broadinstitute.org/pub/users/bhaas/tmp/MouseGencodeM9.tar.gz

Custom references for other organisms or other reference genome versions can be
built using FusionFilter. For more details on how to do so, see the
`FusionFilter wiki`_. For convenience, we also provide a `Python script`_
that executes the required commands for building a STAR-Fusion reference. Note
that this script does expect the required dependencies to be installed
(FusionFilter, Repeatmasker and Blast).

.. _FusionFilter wiki: https://github.com/FusionFilter/FusionFilter/wiki/Building-a-Custom-FusionFilter-Dataset
.. _Python script: https://github.com/jrderuiter/imfusion/blob/develop/scripts/starfusion_build_reference.py

Identifying fusions
~~~~~~~~~~~~~~~~~~~

Using IM-Fusion
===============

Endogenous gene fusions can be identified using IM-Fusion by supplying a
STAR-Fusion reference when identifying insertions using STAR:

.. code:: bash

    imfusion-insertions star
        --fastq sample_s1.R1.fastq.gz \
        --reference references/GRCm38.76.t2onc.star \
        --star_fusion_reference /path/to/star_fusion/reference \
        --output_dir output/sample_s1

This will produce a file called ``gene_fusions.txt`` in the output
directory, which describes the identified endogenous fusions. See the
`STAR-Fusion documentation`_ for a detailed description of the structure of
this file. Note that the STAR-Fusion reference must use the same reference
genome as the  IM-Fusion reference (also using the same chromosome naming
convention), otherwise no fusions will be found.

.. _STAR-Fusion documentation: https://github.com/STAR-Fusion/STAR-Fusion/wiki#output-from-star-fusion

From scratch
============

STAR-Fusion can also be run from scratch (starting with fastq files), which
has the advantage of performing the STAR alignment with settings
specifically tuned for the identification of endogenous gene fusions.
Depending on the dataset, this may provide better results than the settings
used by IM-Fusion, which are geared towards identifying novel fusion events
(with the transposon), rather than focussing on detecting fusions between
known reference genes. As a result, analyses based on the IM-Fusion alignments
may be more prone to false positives than analyses performed using the STAR
settings used by STAR-Fusion.

This type of analysis can be performed using STAR-Fusion with the
following command:

.. code:: bash

    STAR-Fusion --genome_lib_dir /path/to/reference \
                --left_fq sample_s1.R1.fastq.gz \
                --right_fq sample_s1.R2.fastq.gz \
                --output_dir ./out/sample_s1

Note that the ``--right_fq`` argument is optional and can be omitted for
single-end sequencing data.

Alternatively, STAR can first be called separately to perform the alignment,
so that you have more control over the supplied options or the location
of the generated alignments (which can be useful for further analyses). In this
case, the authors of STAR-Fusion recommend using the following parameters for
the alignment:

.. code:: bash

     STAR --genomeDir ${star_index_dir} \
          --readFilesIn ${left_fq_filename} ${right_fq_filename} \
          --twopassMode Basic \
          --outReadsUnmapped None \
          --chimSegmentMin 12 \
          --chimJunctionOverhangMin 12 \
          --alignSJDBoverhangMin 10 \
          --alignMatesGapMax 200000 \
          --alignIntronMax 200000 \
          --chimSegmentReadGapMax parameter 3 \
          --alignSJstitchMismatchNmax 5 -1 5 5 \
          --runThreadN ${THREAD_COUNT} \
          --limitBAMsortRAM 31532137230 \
          --outSAMtype BAM SortedByCoordinate

After the alignment, the produced ``Chimeric.out.junction`` file(s) can be
analyzed to identify fusions with the following command:

.. code:: bash

    STAR-Fusion --genome_lib_dir /path/to/your/CTAT_resource_lib \
                -J Chimeric.out.junction \
                --output_dir star_fusion_outdir
