Migration Guide for cclib 1.9
==============================

This guide helps you migrate your code from cclib 1.8.x to cclib 1.9.

Overview
--------

cclib 1.9 introduces several API improvements that make the library more user-friendly
and scientifically accurate. The main change is that **all energy attributes now use
atomic units (hartree) instead of eV** for consistency.

Breaking Changes
----------------

Energy Units: eV → Hartree
~~~~~~~~~~~~~~~~~~~~~~~~~~

**What changed:**

All energy attributes that were previously in eV are now in atomic units (hartree):

- ``scfenergies`` (HF/DFT energies)
- ``ccenergies`` (Coupled-Cluster energies)
- ``mpenergies`` (Møller-Plesset energies)
- ``dispersionenergies`` (dispersion corrections)
- ``moenergies`` (molecular orbital energies)
- ``scanenergies`` (potential energy scan)

Attributes that were already in hartree remain unchanged:

- ``enthalpy``
- ``freeenergy``
- ``zpve``

**Why this was changed:**

1. **Scientific consistency:** Atomic units are the natural units for quantum chemistry
2. **Eliminates confusion:** Thermodynamic properties (enthalpy, free energy) were already in hartree
3. **Easier calculations:** No need to convert between eV and hartree when computing correlation energies

**How to migrate:**

.. code-block:: python

    # OLD CODE (cclib 1.8.x - energies in eV):
    import cclib
    data = cclib.io.ccread("water.out")
    scf_energy_ev = data.scfenergies[0]

    # NEW CODE (cclib 1.9 - energies in hartree):
    import cclib
    data = cclib.parse("water.out")  # or cclib.io.ccread()
    scf_energy_hartree = data.scfenergies[0]

    # Convert to eV if needed:
    scf_energy_ev = data.convert('scfenergies', 'eV')[0]

    # Or use the conversion factor (1 hartree = 27.21138505 eV):
    scf_energy_ev = data.scfenergies[0] * 27.21138505

Example Migration Scenarios
----------------------------

Scenario 1: Simple energy extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # OLD (1.8.x):
    data = cclib.io.ccread("calculation.out")
    energy_ev = data.scfenergies[-1]  # Last SCF energy in eV

    # NEW (1.9):
    data = cclib.parse("calculation.out")
    energy_hartree = data.scfenergies[-1]  # Last SCF energy in hartree
    energy_ev = data.convert('scfenergies', 'eV')[-1]  # Convert to eV

Scenario 2: Energy differences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # OLD (1.8.x):
    data = cclib.io.ccread("optimization.out")
    delta_e_ev = data.scfenergies[-1] - data.scfenergies[0]

    # NEW (1.9):
    data = cclib.parse("optimization.out")
    delta_e_hartree = data.scfenergies[-1] - data.scfenergies[0]
    delta_e_ev = (data.scfenergies[-1] - data.scfenergies[0]) * 27.21138505

    # Or convert after calculation:
    delta_e_kcal = delta_e_hartree * 627.5094741  # hartree to kcal/mol

Scenario 3: Molecular orbital energies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # OLD (1.8.x):
    data = cclib.io.ccread("molecule.out")
    homo_energy_ev = data.moenergies[0][data.homos[0]]

    # NEW (1.9):
    data = cclib.parse("molecule.out")
    homo_energy_hartree = data.moenergies[0][data.homos[0]]
    homo_energy_ev = data.convert('moenergies', 'eV')[0][data.homos[0]]

Scenario 4: Correlation energies (MP2, CCSD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # OLD (1.8.x):
    data = cclib.io.ccread("mp2_calculation.out")
    hf_energy = data.scfenergies[0]  # in eV
    mp2_energy = data.mpenergies[0]  # in eV
    correlation_ev = mp2_energy - hf_energy

    # NEW (1.9) - Now consistent in hartree:
    data = cclib.parse("mp2_calculation.out")
    hf_energy = data.referenceenergies[0]  # Clearer name, in hartree
    mp2_energy = data.electronicenergies[0]  # Auto-selects MP2, in hartree
    correlation_hartree = mp2_energy - hf_energy

    # Convert if needed:
    correlation_ev = correlation_hartree * 27.21138505

New Features
------------

1. Simple parse() function
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A new convenience function combines file opening and parsing:

.. code-block:: python

    # Instead of:
    parser = cclib.io.ccopen("file.out")
    data = parser.parse()

    # You can now use:
    data = cclib.parse("file.out")

2. has() method
~~~~~~~~~~~~~~~

Check attribute availability before accessing:

.. code-block:: python

    if data.has('vibfreqs'):
        print(f"Found {len(data.vibfreqs)} frequencies")
    else:
        print("No vibrational frequencies")

3. convert() method
~~~~~~~~~~~~~~~~~~~

Easy unit conversion for any attribute:

.. code-block:: python

    # Convert energies
    data.convert('scfenergies', 'eV')
    data.convert('scfenergies', 'kcal/mol')
    data.convert('scfenergies', 'kJ/mol')
    data.convert('moenergies', 'wavenumber')

    # Also works for other attributes
    data.convert('atomcoords', 'bohr')
    data.convert('time', 'time_au')

Supported units:

- Energy: ``hartree``, ``eV``, ``kcal/mol``, ``kJ/mol``, ``wavenumber``
- Length: ``angstrom``, ``bohr``
- Time: ``fs``, ``time_au``
- Dipole: ``Debye``, ``ebohr``

4. New energy properties
~~~~~~~~~~~~~~~~~~~~~~~~~

Clearer terminology for energy attributes:

.. code-block:: python

    # referenceenergies - clearer than scfenergies for post-HF methods
    data.referenceenergies  # Alias for scfenergies

    # electronicenergies - automatically selects best available energy
    data.electronicenergies  # CCSD > MP > SCF (auto-selected)

5. Better error messages
~~~~~~~~~~~~~~~~~~~~~~~~

More helpful messages when accessing unavailable attributes:

.. code-block:: python

    >>> data.vibfreqs  # If not available
    AttributeError: 'ccData' object does not have attribute 'vibfreqs'.
    This attribute was not parsed from the output file. This usually means
    the calculation did not produce this property, or it is not supported
    by the parser for this QM package. Use has('vibfreqs') to check if an
    attribute is available before accessing it.

Quick Reference: Unit Conversion
---------------------------------

Common conversion factors:

.. code-block:: python

    # Energy conversions FROM hartree:
    eV = hartree * 27.21138505
    kcal_mol = hartree * 627.5094741
    kJ_mol = hartree * 2625.4996398
    wavenumber = hartree * 219474.6313708

    # Energy conversions TO hartree:
    hartree = eV / 27.21138505
    hartree = kcal_mol / 627.5094741
    hartree = kJ_mol / 2625.4996398
    hartree = wavenumber / 219474.6313708

Or simply use the ``convert()`` method:

.. code-block:: python

    # No need to remember conversion factors!
    data.convert('scfenergies', 'eV')
    data.convert('scfenergies', 'kcal/mol')

Testing Your Migration
-----------------------

To verify your code works correctly with cclib 1.9:

1. **Check energy values:** Make sure you're not accidentally using hartree values where eV was expected
2. **Look for hardcoded conversions:** Search for ``27.21`` or ``27.2114`` in your code - these may need updating
3. **Test calculations:** Verify that energy differences, reaction energies, etc. are calculated correctly
4. **Update documentation:** Update any documentation or comments that mention "eV" for energy units

Common Pitfalls
---------------

1. **Mixing units:** Don't mix hartree and eV values in calculations:

   .. code-block:: python

       # WRONG:
       barrier = scf_ts_eV - scf_reactant_hartree  # Mixed units!

       # CORRECT:
       barrier = scf_ts_hartree - scf_reactant_hartree
       # or
       barrier = scf_ts_eV - (scf_reactant_hartree * 27.21138505)

2. **Plotting without conversion:** If you're plotting energies, you may want to convert to familiar units:

   .. code-block:: python

       import matplotlib.pyplot as plt

       # Convert for better readability on plots
       energies_kcal = [e * 627.51 for e in data.scfenergies]
       # or
       energies_kcal = data.convert('scfenergies', 'kcal/mol')

       plt.plot(energies_kcal)
       plt.ylabel('Energy (kcal/mol)')

3. **Output formatting:** Update any print statements or file output:

   .. code-block:: python

       # OLD:
       print(f"Energy: {data.scfenergies[0]:.2f} eV")

       # NEW:
       print(f"Energy: {data.scfenergies[0]:.6f} hartree")
       # or
       print(f"Energy: {data.convert('scfenergies', 'eV')[0]:.2f} eV")

Getting Help
------------

If you encounter issues during migration:

1. Check the updated `documentation <https://cclib.github.io>`_
2. Review the `API reference <data.html>`_ for attribute details
3. File an issue on `GitHub <https://github.com/cclib/cclib/issues>`_

Summary
-------

**Key takeaways:**

✅ All energy attributes now in hartree (was eV)

✅ Use ``data.convert(attr, 'eV')`` to get eV values

✅ New ``cclib.parse()`` function for simpler parsing

✅ New ``has()`` method to check attribute availability

✅ New ``referenceenergies`` and ``electronicenergies`` properties

✅ Better error messages for missing attributes

The changes in cclib 1.9 make the library more consistent, scientifically accurate,
and easier to use. While there is a breaking change in energy units, the migration
is straightforward using the new ``convert()`` method.
