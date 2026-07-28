Calculations
============

.. automodule:: mymetal.calculate

Surface energy
--------------

.. automodule:: mymetal.calculate.calenergy.surfenergy

.. autofunction:: mymetal.calculate.calenergy.surfenergy.cal_surface_energy

Contract
~~~~~~~~

* ``bulk_energy`` 与 ``relaxed_surface_energy``：eV；
* ``bulk_atoms_number`` 与 ``surface_atoms_number``：正整数；
* ``area``：一个表面的面积，Å²，必须大于 0；
* ``factor``：excess energy 对应的等价表面数，正整数；
* ``energy_unit="eV"`` 返回 eV/Å²；``"J"`` 返回 J/m²。输入能量始终按
  eV 解释。

函数实现的公式为：

.. math::

   \frac{E_\mathrm{slab}
   - E_\mathrm{bulk}N_\mathrm{slab}/N_\mathrm{bulk}}
   {fA}.

返回 ``float``/NumPy scalar。非法 atom count、area、factor 或
``energy_unit`` 抛 ``ValueError``；无法相除的非数值输入会抛 Python 原生
``TypeError``。

最小示例
~~~~~~~~

.. code-block:: python

   from mymetal.calculate.calenergy.surfenergy import cal_surface_energy

   gamma = cal_surface_energy(
       bulk_energy=-47.04775262,
       bulk_atoms_number=12,
       relaxed_surface_energy=-46.24580704,
       surface_atoms_number=12,
       area=7.026888773593904,
       factor=2,
   )
   assert abs(gamma - 0.05706263510343279) < 1e-12

Related tutorials
~~~~~~~~~~~~~~~~~

See :doc:`../tutorials/surface_energy`.

Mechanics and mismatch
----------------------

.. automodule:: mymetal.calculate.calmechanics.deformation
   :members:
   :show-inheritance:

.. automodule:: mymetal.calculate.calmechanics.strain
   :members:
   :show-inheritance:

.. automodule:: mymetal.calculate.calmechanics.stress
   :members:
   :show-inheritance:

.. automodule:: mymetal.calculate.calmismatch.calhetero

Principal entry point: ``mymetal.calculate.calmismatch.calhetero.cal_mismatch``.

Electronic and materials analysis
---------------------------------

.. automodule:: mymetal.calculate.calqm.kpoints
   :members:
   :show-inheritance:

.. automodule:: mymetal.calculate.material_science.schmid
