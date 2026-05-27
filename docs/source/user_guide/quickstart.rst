Quick Start
===========

Generate a film
---------------

.. code-block:: python

   from mymetal.build.film.stretch import generate_film

   slab = generate_film(
       symbols="Au",
       structure="fcc",
       num_layers=12,
       my_vacuum=20.0,
       slice_plane=(1, 1, 1),
   )

Read and write VASP files
-------------------------

.. code-block:: python

   from mymetal.io.vasp import my_read_vasp, my_write_vasp

   atoms, scale = my_read_vasp("CONTCAR")
   my_write_vasp("POSCAR", atoms, lattice_scale_factor=scale)

Calculate surface energy
------------------------

.. code-block:: python

   from mymetal.calculate.calenergy.surfenergy import cal_surface_energy

   gamma = cal_surface_energy(
       bulk_energy=-100.0,
       bulk_atoms_number=4,
       relaxed_surface_energy=-190.0,
       surface_atoms_number=8,
       area=25.0,
       energy_unit="eV",
   )

Export n2p2 data
----------------

.. code-block:: python

   from mymetal.ml.n2p2.dataset import nnpdata

   data = nnpdata()
   data.load_from_outcar("./OUTCAR", index=":", tag="train")
   data.write("./input.data")
