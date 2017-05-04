=======
History
=======

0.3.0 (2017-05-04)
------------------

* Refactored external tools into the ``imfusion.external`` module.
* Use docker/tox for testing against multiple Python versions locally.
* Added additional checks for inputs and improved error messages.
* Added support for DataFrame insertion inputs to DE testing functions.
* Added building of exon gtf as part of imfusion-build.
* Added identification of endogenous fusions using STAR-Fusion as part
  of imfusion-insertions (using STAR). Also adds script for building
  (murine) STAR-Fusion references.
* Made matplotlib/seaborn lazy imports that are only required when actually
  using the plotting functions. This makes IM-Fusion easier to use on
  headless servers/HPCs.

0.2.0 (2017-03-09)
------------------

* Added support for the STAR aligner.
* Added detection of novel transcripts using StringTie.
* Changed reference building to generate a self-contained reference.
* Refactored differential expression tests + added gene-level test.

0.1.0 (2016-03-26)
------------------

* First release on GitHub.
