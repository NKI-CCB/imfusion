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

``Aligner`` classes are used to identify insertions using the different RNA-seq
aligners, with one class provided for each supported aligner. The main
methods provided by the Indexers are the ``check_dependencies`` and
``identify_insertions`` methods. The latter should be called with a reference
instance and the sequence read files to identify insertions.

.. autoclass:: imfusion.insertions.aligners.Aligner
    :members:

.. autoclass:: imfusion.insertions.aligners.StarAligner
    :members:

.. autoclass:: imfusion.insertions.aligners.TophatAligner
    :members:


Testing for differential expression
-----------------------------------

TODO
