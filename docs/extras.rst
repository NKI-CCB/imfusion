======
Extras
======

Using STAR-fusion to identify endogenous fusions
------------------------------------------------

An important advantage of using RNA-sequencing to identify insertions, rather
than targeted DNA-sequencing, is that RNA-seq data provides more information
about other events that may be contributing to tumorigenesis. For example,
besides the gene-transposon insertions, we may expect somatic gene-fusions
(between endogenous genes) to occur in these tumors, which may also contribute
to tumorigenesis.

In our manuscript, we show a proof of concept of this idea for both the SB
and the B-ALL datasets by identifying endogenous fusions from STAR RNA-seq
alignments using `STAR-Fusion <http://star-fusion.github.io>`_.
Using this approach, we identified the engineered *Etv6-Runx1* fusion in the
B-ALL tumors and identified several known oncogenic fusions in the SB tumors.

Here we provide a brief description of how to apply STAR-Fusion to identify
endogenous gene fusions. In the first approach, we describe how to apply
STAR-Fusion to output from IM-Fusion, in which we assume the STAR aligner was
to identify fusions. In the second, we describe how to run STAR-Fusion from
scratch.

The main advantage of the first approach is that we only have to
perform a single alignment (saving computational power and disk space).
However, IM-Fusion uses different settings for the alignment than STAR-Fusion,
as IM-Fusion is aimed at identifying novel fusion events (between the
transposon and potentially novel gene transcripts), whilst STAR-Fusion is
mainly designed to identify fusions between known reference genes. As a result,
the settings used by STAR-Fusion are likely more conservative and running
STAR-Fusion from scratch may therefore be less likely to identify
false positives.

Installation
~~~~~~~~~~~~

STAR-Fusion can be downloaded from its GitHub
`repository <https://github.com/STAR-Fusion/STAR-Fusion/releases>`_. After
unpacking, make sure that the ``STAR-Fusion`` binary is available in ``$PATH``.
Note that STAR-Fusion also requires STAR to be installed if run from scratch.

Alternatively, STAR-Fusion can also be installed using bioconda:

.. code:: bash

    conda install -c bioconda star-fusion

This should also install any dependencies of STAR-Fusion (such as STAR). Note
that the bioconda installation only works on linux at the moment, as bioconda
currently does not provide any STAR-Fusion builds for macOS.

Preparing a reference
~~~~~~~~~~~~~~~~~~~~~

Unfortunately, building a STAR-Fusion reference is quite involved as it
requires a number of additional tools and resources to be installed. We
therefore recommend using a pre-built reference if possible. Pre-built
references (for human) are available
`here <https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB>`_. A
reference for the mouse mm10 genome (using M9/Gencode) is available
`here <ftp://ftp.broadinstitute.org/pub/users/bhaas/tmp/MouseGencodeM9.tar.gz>`_.

Custom references for other organisms or other reference genome versions can be
built using FusionFilter. For more details on how to do so, see the FusionFilter
`wiki <https://github.com/FusionFilter/FusionFilter/wiki/Building-a-Custom-FusionFilter-Dataset>`_.
For convenience, we also provide a Python
`script <https://github.com/jrderuiter/imfusion/blob/develop/scripts/starfusion_build_reference.py>`_
that executes the required commands for building a STAR-Fusion reference. Note
that this script does expect the required dependencies to be installed
(FusionFilter, Repeatmasker and Blast).

Identifying fusions
~~~~~~~~~~~~~~~~~~~

Using IM-Fusion
===============

The alignments produced by IM-Fusion (assuming was used STAR as aligner) can be
used to directly identify fusions without performing another alignment (with
the above mentioned caveats). This analysis can be performed by supplying the
``Chimeric.out.junction`` files from the alignments to STAR-Fusion as follows:

.. code:: bash

    STAR-Fusion --genome_lib_dir /path/to/reference \
                -J out/sample_s1/_star/Chimeric.out.junction \
                --output_dir fusions/sample_s1

This produces a tab-delimited file named
``star-fusion.fusion_candidates.final.abridged`` in the specified output
directory, which contains the gene fusions identified by STAR-Fusion.

From scratch
============

STAR-Fusion can also be run from scratch (starting with fastq files), which
has the advantage of performing the RNA-seq alignment with settings
specifically tuned for the identification of endogenous gene fusions. This type
of analysis can be performed using the following command:

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
analyzed using STAR-Fusion in the same fashion as for the IM-Fusion alignments.
