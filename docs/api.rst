API
===

Building references
-------------------

Indexers
~~~~~~~~

IM-Fusion uses ``Indexer`` classes to generate references for the
different RNA-seq aligners, with one class provided for each aligner. The main
methods provided by the Indexers are the ``check_dependencies`` and ``build``
methods. The former is used to check whether required external dependencies
are present in ``$PATH``, whilst the latter should be called with the required
reference files to build the actual reference.

.. autoclass:: imfusion.build.indexers.Indexer
    :members:

.. autoclass:: imfusion.build.indexers.StarIndexer
    :members:

.. autoclass:: imfusion.build.indexers.TophatIndexer
    :members:

References
~~~~~~~~~~

``Indexer`` classes return ``Reference`` objects after building
references, which define paths to specific files within the reference.
This includes, for example, paths to the augmented fasta file, to the file
describing transposon features and to the reference index. The ``Reference``
base class defines paths that should be available for all references.
Subclasses may define additional paths that are specific to the
corresponding aligner.

.. autoclass:: imfusion.build.indexers.Reference
    :members:

.. autoclass:: imfusion.build.indexers.StarReference
    :members:

.. autoclass:: imfusion.build.indexers.TophatReference
    :members:


Identifying insertions
----------------------

IM-Fusion uses ``Aligner`` classes to identify insertions using the different
RNA-seq aligners, with one class provided for each supported aligner. The main
methods provided by the Indexers are the ``check_dependencies`` and
``identify_insertions`` methods. The latter should be called with a reference
instance and the sequence read files to identify insertions.

.. autoclass:: imfusion.insertions.aligners.Aligner
    :members:

.. autoclass:: imfusion.insertions.aligners.StarAligner
    :members:

.. autoclass:: imfusion.insertions.aligners.TophatAligner
    :members:


Differential expression
-----------------------

Quantification
~~~~~~~~~~~~~~

The ``imfusion.expression`` module provides several functions for generating,
reading and normalizing exon expression counts. These functions are used by
IM-Fusion to prepare counts for the differential expression tests.

.. autofunction:: imfusion.expression.generate_exon_counts
.. autofunction:: imfusion.expression.read_exon_counts
.. autofunction:: imfusion.expression.normalize_counts
.. autofunction:: imfusion.expression.estimate_size_factors

DE tests
~~~~~~~~

IM-Fusion provides three differential expression tests, which are implemented
in their own functions (``test_de_exon``, ``test_de_exon_single`` and
``test_de_gene``). Each of these functions returns a corresponding test result
instance (see below for more details). The ``test_de`` function is a
convenient wrapper that applies the group-wise exon-level test for multiple
genes, and provides an optional fallback to the gene-level test for genes that
cannot be tested at the exon-level.

.. autofunction:: imfusion.expression.test_de
.. autofunction:: imfusion.expression.test_de_exon
.. autofunction:: imfusion.expression.test_de_exon_single
.. autofunction:: imfusion.expression.test_de_gene

DE results
~~~~~~~~~~

Each of the three differential expression tests returns a result instance that
contains the results of the test and provides several functions for
visualizing the test result.

.. autoclass:: imfusion.expression.test.DeResult
    :members:

.. autoclass:: imfusion.expression.test.DeSingleResult
    :members:

.. autoclass:: imfusion.expression.test.DeGeneResult
    :members:
