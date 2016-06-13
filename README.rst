.. image:: https://img.shields.io/travis/jrderuiter/im-fusion/develop.svg?style=flat-square
    :target: https://travis-ci.org/jrderuiter/im-fusion

.. image:: https://img.shields.io/coveralls/jrderuiter/imfusion/develop.svg?style=flat-square
    :target: https://coveralls.io/github/jrderuiter/imfusion

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
approaches, please see our paper "Identifying transposon insertions and
their effects from RNA-sequencing data" (*submitted*).

Documentation
=============

IM-Fusion's documentation is available on ReadTheDocs.

Installation
============

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

For details on how to install IM-Fusion without using conda, please see
the installation section of the documentation.

License
=======

This software is released under the MIT license.
