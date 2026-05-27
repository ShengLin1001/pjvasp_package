Structure Building
------------------

.. automodule:: mymetal.build
   :members:

Bulk structures
---------------

``mymetal.build.bulk.create`` contains bulk and slab construction helpers.
Some functions import the runtime ``myvasp`` helper package, so this module is
described here but not imported during the online documentation build.

Important entry points include ``create_fcc_111`` and ``create_hcp_basal``.

Film construction
-----------------

.. automodule:: mymetal.build.film.stretch

Principal entry point: ``mymetal.build.film.stretch.generate_film``.

.. automodule:: mymetal.build.film.hydroxyl

.. automodule:: mymetal.build.film.findprim
   :members:
   :show-inheritance:

.. automodule:: mymetal.build.film.findhetero

Principal entry points include ``find_hetero`` and ``build_supercells``.
