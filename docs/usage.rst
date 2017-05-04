=====
Usage
=====

Building a reference
--------------------

Before detecting insertion sites, IM-Fusion first needs to build an
augmented version of the host reference genome that contains the
sequence of the transposon as additional sequence. This augmented reference
is created using the ``imfusion-build`` command, which essentially concatenates
the host reference genome and the transposon sequence (both in Fasta
format) into a new Fasta file and builds the indices needed for alignment.
Separate sub-commands are provided for each supported aligner (currently STAR
and Tophat-Fusion).

The basic command for building a (STAR-based) reference is as follows:

.. code:: bash

    imfusion-build star \
        --reference_seq Mus_musculus.GRCm38.dna.primary_assembly.fa \
        --reference_gtf Mus_musculus.GRCm38.76.gtf \
        --transposon_seq transposon.fa \
        --transposon_features transposon.features.txt \
        --output_dir references/GRCm38.76.t2onc.star

Here, the ``--reference_seq`` argument points to the fasta file containing
the sequence of the reference genome, ``--reference_gtf`` points to the gtf
file containing the reference gene features, ``--transposon_seq`` points to
the transposon sequence and ``--transposon_features`` points to a file
describing the transposon features. The reference files are written to the
path specified by ``--output_dir``.

Optionally, specific genomic regions can be blacklisted using the
``--blacklist_genes`` and ``--blacklist_regions`` arguments. The former can
be used to blacklist entire genes (specified by their Ensembl ID), whilst the
latter can be used to blacklist specific regions (specified as chr:start-end).
Sequences of blacklisted regions are replaced by 'N' nucleotides in the
generated reference.

The command for building a Tophat-Fusion-based reference is nearly identical:

.. code:: bash

    imfusion-build tophat \
        --reference_seq Mus_musculus.GRCm38.dna.primary_assembly.fa \
        --reference_gtf Mus_musculus.GRCm38.76.gtf \
        --transposon_seq transposon.fa \
        --transposon_features transposon.features.txt \
        --output_dir references/GRCm38.76.t2onc.tophat

However, both aligners do have some aligner-specific arguments for building
their references. See the help of the respective sub-commands for more details.
For STAR, special attention should be paid to memory usage, as STAR requires
approximately 30GB of memory for building the reference genome.

Detecting insertions (per sample)
---------------------------------

After building the augmented reference, we can detect transposon insertions
in each sample using the ``imfusion-insertions`` command. This command
essentially runs a gene-fusion aware aligner (either STAR or Tophat-Fusion,
with seperate sub-commands for each aligner) to align the RNA-seq reads and
identify gene fusions, extracts gene-transposon fusions from the results and
derives the corresponding insertion sites. The insertions are written as the
tab-separated file ``insertions.txt`` in the output directory.

The basic command for ``imfusion-insertions`` using the STAR aligner is:

.. code:: bash

    imfusion-insertions star
        --fastq sample_s1.R1.fastq.gz \
        --reference references/GRCm38.76.t2onc.star \
        --output_dir output/sample_s1

In this command, the ``--fastq`` argument specifies a path to the fastq
file containing RNA-seq reads for the given sample. For paired-end samples, the
second pair should be provided using the optional ``--fastq2`` argument.
The ``--reference`` argument should point to the previously built augmented
reference, whilst the ``--output_dir`` argument specifies where the
sample output should be written.

An optional ``--assemble`` argument indicates whether IM-Fusion should perform
a reference-guided transcript assembly. If given, IM-Fusion runs Stringtie
after the RNA-seq alignment to detect novel gene transcripts based on the
RNA-seq alignment. The results of this assembly are subsequently used in the
insertion detection step to annotate insertions that involve novel transcripts.

The command for using Tophat-Fusion is nearly identical:

.. code:: bash

    imfusion-insertions tophat
        --fastq sample_s1.R1.fastq.gz \
        --reference references/GRCm38.76.t2onc.tophat \
        --output_dir output/sample_s1

However, both aligners do have some aligner-specific arguments concerning the
alignment. See the help of the respective sub-commands for more details. For
STAR, special attention should be paid to memory usage, as STAR requires
approximately 30GB of memory (per process) for loading the reference genome.

Quantifying expression (per sample)
-----------------------------------

After detecting insertions, we use the generated RNA-seq alignment to quantify
exon expression counts for the given sample. These counts are later used to
test for differential expression when identifying candidate genes from a group
of samples.

The expression counts are generated using the ``imfusion-expression`` command:

.. code:: bash

    imfusion-expression \
        --sample_dir ./output/s1 \
        --exon_gtf exons.gtf

Here, the ``--sample_dir`` argument should point to a sample directory (which
was previously generated by ``imfusion-insertions``). The generated counts are
written to this sample directory as the TSV file ``exon_counts.txt``.
The ``--exon_gtf`` argument refers to a gtf file containing the flattened
exon representation of previously used reference gtf.

The optional arguments ``--paired`` and ``--stranded`` are used to indicate
whether the alignment contains paired-end sequencing data and to indicate
the strandedness of this data (0 = unstranded, 1 = stranded
and 2 = reversely stranded).

Merging sample results
----------------------

To detect genes that are recurrently mutated across samples, we first merge
the individual sample results into a combined dataset using ``imfusion-merge``.
This command effectively concatenates the individual results into combined
``insertions.txt`` and ``exon_counts.txt`` files.

The basic command is as follows:

.. code:: bash

    im-fusion merge --sample_dirs ./output/sample_s1 \
                                  ./output/sample_s2 \
                    --output ./output/merged.insertions.txt \
                    --output_expression ./output/merged.exon_counts.txt

In this command, the ``--sample_dirs`` argument points to the sample
directories that should be merged and ``output`` indicates that path to
which the merged insertion file should be written. The ``--output_expression``
argument indicates where merged expression counts should be written. This
argument may be omitted if no expression counts were generated.

Selecting (DE) CTGs
-------------------

To identify genes that are commonly targeted by insertions (commonly targeted
genes, or CTGs), IM-Fusion uses the Poisson distribution to test whether a
given gene has more insertions than may be expected by chance. This test is
performed on the merged dataset using the ``imfusion-ctg`` command:

.. code:: bash

    im-fusion ctg --insertions ./output/merged.insertions.txt  \
                  --expression ./output/merged.exon_counts.txt \
                  --reference references/GRCm38.76.t2onc.star \
                  --output ./output/merged.ctgs.txt

Here, ``--insertions`` and ``--expression`` should point to the merged
insertions and expression files generated by ``imfusion-merge``. The
``--reference`` argument refers to the same reference as used for the
alignment, whilst ``--output`` specifies the path where the CTG output
should be written.

The parameters for the CTG test can be changed using the ``--window``,
``--pattern`` and ``--chromosomes`` arguments. The ``--window`` parameter
specifies the size of the window around genes within which insertions should
be included. The ``--pattern`` argument can be used to account for integration
biases of the transposon, if the transposon is known to integrate at specific
nucleotide sequences. The ``--chromosome`` argument specifies which chromosomes
should be included in the test and is mainly used to omit chromosomes
containing the transposon donor loci.

Finally, the significance thresholds for the CTG and DE tests can be specified
using the ``--threshold`` and ``--de_threshold`` arguments:

.. code:: bash

    im-fusion ctg --insertions ./output/merged.insertions.txt  \
                  --expression ./output/merged.exon_counts.txt \
                  --reference references/GRCm38.76.t2onc.star \
                  --output ./output/merged.ctgs.txt \
                  --threshold 0.05 \
                  --de_threshold 0.05

Optionally, the differential expression test can be skipped by not providing
the expression data. In this case, only the CTG test is performed:

.. code:: bash

    im-fusion ctg --insertions ./merged.insertions.txt  \
                  --reference references/GRCm38.76.t2onc.star \
                  --output ./output/merged.ctgs.txt
                  --threshold 0.05 \
