Home
====

IM-Fusion
---------

IM-Fusion is a tool for identifying transposon insertion sites in
insertional mutagenesis screens using single- and paired-end RNA-sequencing
data. It essentially identifies insertion sites from gene-transposon fusions
in the RNA-sequencing data, which represent splicing events between the
transposon and endogeneous genes.

IM-Fusion also identifies candidate genes for a given screen using a
statistical test (based on the Poisson distribution) that identifies Commonly
Targeted Genes (CTGs) -- genes that are more frequently affected by insertions
than would be expected by chance. To futher narrow down a list of CTGs, which
may contain hundreds of genes, IM-Fusion also tests if insertions in a CTG have
a significant effect on the expression of the gene, which is a strong indicator
of them having an actual biological effect.

IM-Fusion has the following key features:

- It identifies transposon insertion sites from both single- and paired-end
  RNA-sequencing data, without having any special sequencing requirements.
- Its gene-centric approach -- both in the identification of insertions and
  in testing for differential expression of identified candidate genes --
  greatly reduces the number of false positive candidate genes.
- Its differential expression tests provide detailed insight into the effects
  of insertions on the expression of their target gene(s). By providing both
  a groupwise and a single-sample version of the test, IM-Fusion can identify
  effects for a single insertion in a specific sample, or determine the general
  effect of insertions on a given gene within the tumor cohort.

For more details on the approach and a comparison with existing
DNA-sequencing approaches, please see our preprint *"Identifying transposon
insertions and their effects from RNA-sequencing data"* on bioRxiv.

Where to start?
---------------

For details on installing IM-Fusion and its dependencies, see
:doc:`installation`.

The main interface of IM-Fusion is its command-line interface, which provides
various commands for performing each step of its analysis. For an overview
of the analysis and for a guide to preparing the required input files,
see :doc:`getting_started`. For more details on the different options available
for each command, see :doc:`command_line_usage`.

An API interface is also provided for programmatic use of IM-Fusions various
steps, or for creating custom pipelines based on our code. For details on the
provided functions, see the :doc:`api` documentation.

Sources
-------

The source code of IM-Fusion is available on Github in the repository
`jrderuiter/im-fusion <http://github.com/jrderuiter/im-fusion>`_. Any
questions or issues should also be addressed there.

License
-------

IM-Fusion is released under the MIT license.
