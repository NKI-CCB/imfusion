=======
History
=======

0.3.2 (2017-05-11)
------------------

* Properly added star-fusion support to star aligner (was previously not
  fully merged).
* Changed documentation URLs to new repository.

0.3.1 (2017-05-09)
------------------

* Several small fixes for Python 2.7 compatibility.
* Fixed issue in CTG test that occurs when no insertions are within the
  gene windows.
* Replaced usage of deprecated .ix indexer for pandas DataFrames.
* Updated documentation.

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
