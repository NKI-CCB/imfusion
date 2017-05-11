Getting started
===============

Overview
--------

An IM-Fusion analysis consists of the following steps:

    1. **Building a reference** - An augmented reference genome is created,
       which contains both the host reference sequence and the transposon
       sequence. This reference genome is also indexed for use in alignment.
    2. **Identifying insertions** - Insertions are identified (per sample)
       using one of the fusion-aware aligners (currently STAR or
       Tophat-Fusion), which are used to align reads to the reference genome
       and identify gene-transposon fusions from chimeric read alignments.
       These fusions are analyzed to identify the corresponding transposon
       insertions.
    3. **Quantifying exon expression** - Exon level expression counts are
       generated (per sample) using the RNA-seq alignment from the previous
       step. These counts can be used to (interactively) visualize insertions
       together with the expression of affected genes, or to test for the
       differential expression of candidate genes identified in Step 5.
    4. **Merging samples** - Insertions and expression counts from multiple
       samples are merged into a single dataset, which is more convenient for
       downstream analyses.
    5. **Selecting (DE) CTGs** - The merged dataset is analyzed to identify
       Commonly Targeted Genes (CTGs) – genes that are more frequently affected
       by insertions (across samples) than would be expected by chance. The
       exon expression counts are used to test whether the expression of CTGs
       is significantly affected by their insertions. This provides an
       additional 'biological' filter, allowing us to specifically select for
       these differentially expressed CTGs (DE CTGs), which are more likely to
       have a true biological effect.

Each of these steps is implemented using a dedicated command. For more details
on the individual steps and their commands, see :doc:`usage`. Additionally,
supporting code is provided for interactive analyses (such as plotting
differential expression) and manually running certain steps of the pipeline.
For more details, see :doc:`api`.

For the STAR aligner, we also support identifying endogenous gene fusions
using STAR-Fusion. See :doc:`extras` for more details about this type of
analysis.

Required files
--------------

IM-Fusion does not require any special preparation for sequencing and should
be applicable to standard RNA-sequencing data, as long as these data are not
somehow depleted for the transposon sequences. Both single- and paired-end
sequencing data is supported, although paired-end data is likely to provide
more sensitivity in the results.

For building a reference, the following files are required:

    - The reference genome sequence (fasta file)
    - The reference gene features (GTF file, adhering to Ensembl's GTF format)
    - The transposon sequence (fasta file)
    - The transposon features file (tab-separated file)

The transposon feature file describes the locations of the splice-donor and
splice-acceptor sites that are present in the transposon. The file should
contain the following columns: ``name``, ``start``, ``end``, ``strand`` and
``type``; with each row describing a single feature. The first four of these
columns describe the name, location and orientation of each feature,
respectively. The last column specifies the type of the feature, which should
be either ‘SD’ or ‘SA’ for splice-donor or splice-acceptor sites respectively.
The field may also be left empty for other types of features, however these
features will not be used by IM-Fusion.

For generating the exon-level expression counts, IM-Fusion needs a
flattened exon representation of the reference gene features (in GTF format).
Example reference files
-----------------------

Reference files for various genomes can be downloaded here:

    - mm10 (Ensembl, current): `Fasta <ftp://ftp.ensembl.org//pub/release-88/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz>`_, `GTF <ftp://ftp.ensembl.org//pub/release-88/gtf/mus_musculus/Mus_musculus.GRCm38.88.gtf.gz>`_
    - mm9 (Ensembl, archive): `Fasta <ftp://ftp.ensembl.org//pub/release-67/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz>`_, `GTF <ftp://ftp.ensembl.org//pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz>`_

Example sequence and feature files for the Sleeping Beauty
*T2/Onc* and Piggybac transposons are available in our GitHub `respository
<https://github.com/jrderuiter/imfusion/tree/develop/data>`_. For other
transposon,s these files need to be created manually using knowledge of the
transposon sequence and the locations of its splice acceptor/donor sites.
