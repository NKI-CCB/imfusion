API
===

imfusion.aligners
-----------------

The **aligners** module contains the different aligners that can be used to
identify insertion sites from gene-transposon fusions. Currently only
Tophat2 is supported, though other aligners may be added in the future. The
aligner modules only implement functionality that is specific to the aligner.
Generic code for building reference genomes and identifying insertions are
contained in the util modules **util.fusions**, **util.insertions** and
**util.reference**.

imfusion.aligners.tophat2
~~~~~~~~~~~~~~~~~~~~~~~~~

The **tophat2** submodule contains the functionality that is used to
identify insertions from gene-transposon fusions using Tophat2. The module
contains two main functions: `build_reference` and `identify_insertions`.
The function `build_reference` is used to build an augmented reference genome
and the corresponding indices for Tophat2. This augmented reference is can be
supplied to `identify_insertions` to perform the actual
insertion identification.

.. autofunction:: imfusion.aligners.tophat2.build_reference
.. autofunction:: imfusion.aligners.tophat2.identify_insertions
.. autofunction:: imfusion.aligners.tophat2.get_version
.. autofunction:: imfusion.aligners.tophat2.get_bowtie_version

imfusion.ctg
------------

The **ctg** module contains the functions for identifying commonly targeted
genes (CTGs) from a collection of insertions from multiple samples. The module
contains two main functions: `test_ctgs` and `test_de`. The `test_ctgs`
function performs the actual enrichment test and returns CTGs and their
corresponding (corrected) p-values. The function `test_de` takes the CTG
frame and tests each of the genes for differential expression, filtering
out CTGs that are not significantly differentially expressed.

.. autofunction:: imfusion.ctg.test_ctgs
.. autofunction:: imfusion.ctg.test_de

imfusion.expression
-------------------

The **expression** module contains functionality for generating expression
counts and testing for differential expression. The **counts** submodule
handles the count generation using the featureCounts tool. The **de_test**
submodule implements the various differential expression tests.

imfusion.expression.counts
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: imfusion.expression.counts.exon_counts
.. autofunction:: imfusion.expression.counts.feature_counts

imfusion.expression.de_test
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: imfusion.expression.de_test.de_exon
.. autoclass:: imfusion.expression.de_test.DeExonResult
    :members:

.. autofunction:: imfusion.expression.de_test.de_exon_single
.. autoclass:: imfusion.expression.de_test.DeExonSingleResult
    :members:

imfusion.merge
--------------

The **merge** module contains functions for merging the results of
the individual sample analyses (the insertions and expression counts) into
a single combined dataset. This combined dataset is used as input
for the CTG and differential expression analysis.

.. autofunction:: imfusion.merge.merge_samples

imfusion.model
--------------

Two model classes, `Fusion` and `Insertion`, are used to represent fusions and
insertions respectively. These classes are mainly used to track which
attributes fusions and insertions have and to convert between model instances
and DataFrame representations.

.. autoclass:: imfusion.model.Fusion
.. autoclass:: imfusion.model.Insertion

imfusion.util
-------------

The **util** module contains various helper modules shared between
different parts of im-fusion. The most important submodules are **fusions**,
**insertions** and  **reference**, which contain functions that are used by the
aligners (currently only Tophat2) to generate augmented reference genomes and
identify insertion sites.

imfusion.util.check
~~~~~~~~~~~~~~~~~~~

Utility functions for checking the validity of inputs.

.. autofunction:: imfusion.util.check.check_features

imfusion.util.fusions
~~~~~~~~~~~~~~~~~~~~~

Utility functions for annotating fusions and converting fusions into
insertions by determining the approximate position of the corresponding
insertion in the genome ('placing' the fusion).

.. autofunction:: imfusion.util.fusions.annotate_fusions
.. autofunction:: imfusion.util.fusions.place_fusions

imfusion.util.insertions
~~~~~~~~~~~~~~~~~~~~~~~~

Utility functions for filtering invalid/unwanted insertions.

.. autofunction:: imfusion.util.insertions.filter_invalid_insertions
.. autofunction:: imfusion.util.insertions.filter_wrong_orientation
.. autofunction:: imfusion.util.insertions.filter_unexpected_sites
.. autofunction:: imfusion.util.insertions.filter_blacklist

imfusion.util.reference
~~~~~~~~~~~~~~~~~~~~~~~

Utility functions used for generating augmented reference genomes.

.. autofunction:: imfusion.util.reference.concatenate_fastas
.. autofunction:: imfusion.util.reference.mask_regions
.. autofunction:: imfusion.util.reference.blacklist_for_regions
.. autofunction:: imfusion.util.reference.blacklist_for_genes

imfusion.util.tabix
~~~~~~~~~~~~~~~~~~~

Utility classes used for fast access to Gtf and Bed files.

.. autoclass:: imfusion.util.tabix.GtfFile
    :members:
    :inherited-members:
