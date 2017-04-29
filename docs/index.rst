IM-Fusion
=========

IM-Fusion is a tool for identifying transposon insertion sites in
insertional mutagenesis screens using single- and paired-end RNA-sequencing
data. It essentially identifies insertion sites from gene-transposon fusions
in the RNA-sequencing data, which represent splicing events between the
transposon and endogeneous genes.

IM-Fusion also identifies candidate genes for a given screen using a
statistical test (based on the Poisson distribution) that identifies Commonly
Targeted Genes (CTGs) -- genes that are more frequently affected by insertions
than would be expected by chance. To further narrow down a list of CTGs, which
may contain hundreds of genes, IM-Fusion also tests if insertions in a CTG have
a significant effect on the expression of the gene, which is a strong indicator
of them having an actual biological effect.

IM-Fusion has the following key features:

- It identifies transposon insertion sites from both single- and paired-end
  RNA-sequencing data, without having any special sequencing requirements.
- It uses a gene-centric approach -- both for the identification of insertions
  and for testing of differential expression for identified candidate genes --
  which greatly reduces the number of false positive candidate genes.
- It implements several exon-level and gene-level differential expression
  tests, which provide detailed insight into the effects of insertions on
  the expression of their target  gene(s). By providing both a group-wise and
  a single-sample version of the test, IM-Fusion can identify effects for a
  single insertion in a specific sample, or determine the general
  effect of insertions on a given gene within the tumor cohort.

Additionally, by integrating STAR-Fusion into its STAR-based insertion
identification pipeline, IM-Fusion enables simultaneous identification of
transposon insertions and endogenous gene fusions. As shown in our manuscript,
this approach can identify endogenous gene-fusions driving tumorigenesis, which
are impossible to identify using targeted DNA-based transposon-sequencing
approaches.

.. _STAR-Fusion: https://github.com/STAR-Fusion/STAR-Fusion/wiki

For more details on the approach and a comparison with existing DNA-sequencing
approaches, please see our paper **"Identifying transposon insertions and
their effects from RNA-sequencing data"** (*Currently under revision*).

.. toctree::
   :maxdepth: 2
   :hidden:

   self
   installation
   getting_started
   usage
   extras
   api
   contributing
   authors
   history
