# ORCA and XTB Parser Test Suite

## Overview

This document describes the comprehensive sanity test suite for ORCA and XTB parsers (`test_orca_xtb_sanity.py`).

## Test Coverage

### Total: 38 Tests

#### ORCA Parser (20 tests)

**Basic Parsing Tests:**
- `test_dvb_sp_dft_basic_attributes`: Validates natom, charge, mult, nbasis
- `test_dvb_sp_dft_scf_energy`: Tests exact SCF energy (-382.055108 Eh)
- `test_dvb_sp_dft_metadata`: Package version, success status, methods
- `test_dvb_sp_dft_mulliken_charges`: Charge values and conservation

**Post-HF Methods:**
- `test_water_ccsd_basic`: CCSD calculation parsing
- `test_water_ccsd_t_basic`: CCSD(T) calculation parsing
- `test_water_mp2_basic`: MP2 perturbation theory
- `test_water_mp3_basic`: MP3 perturbation theory
- `test_water_ccsd_vs_ccsd_t`: Energy ordering (CCSD(T) < CCSD)

**Geometry and Dynamics:**
- `test_dvb_gopt_basic`: Geometry optimization attributes
- `test_dvb_gopt_convergence`: Optimization convergence status
- `test_dvb_scan_relaxed`: Relaxed surface scan

**Spectroscopy:**
- `test_dvb_ir_frequencies`: 54 vibrational modes for 20-atom system
- `test_dvb_ir_thermochemistry`: Temperature, pressure, enthalpy, free energy
- `test_dvb_raman_intensities`: Raman spectrum data
- `test_dvb_td_excited_states`: TDDFT excited states and oscillator strengths
- `test_dvb_nmr_shieldings`: NMR calculation results

**Advanced Features:**
- `test_dvb_dispersion_energy`: Dispersion correction parsing
- `test_Trp_polar_polarizability`: Polarizability calculations
- `test_dvb_sp_hf_vs_dft`: Energy differences between methods

#### XTB Parser (14 tests)

**Basic Parsing:**
- `test_dvb_sp_basic_attributes`: Core attributes validation
- `test_dvb_sp_charges`: Mulliken atomic charges
- `test_dvb_sp_metadata`: Package version 6.6.1, success status
- `test_dvb_sp_energy_value`: Exact energy (-26.425939 Eh)
- `test_dvb_sp_energy_range`: Energy reasonableness checks

**Geometry Optimization:**
- `test_dvb_opt_basic`: Optimization attributes and SCF energies
- `test_dvb_opt_convergence`: Successful convergence verification
- `test_dvb_opt_energy_decreases`: Energy minimization validation
- `test_dvb_opt_final_energy_lower`: Final energy at minimum

**Vibrational Analysis:**
- `test_dvb_ir_frequencies`: All 60 modes (3N for 20 atoms)
- `test_dvb_ir_intensities`: IR intensity data
- `test_dvb_ir_zero_point_energy`: ZPE calculation
- `test_dvb_ir_has_negative_frequencies`: Frequency sign validation

**Data Quality:**
- `test_charge_conservation`: Atomic charges sum to total charge

#### Data Consistency (4 tests)

**Cross-Parser Validation:**
- `test_orca_atom_count_consistency`: natom matches atomnos/atomcoords length
- `test_xtb_atom_count_consistency`: Same for XTB parser
- `test_orca_coordinates_3d`: All coordinates are 3D vectors
- `test_xtb_coordinates_3d`: Same for XTB parser

## Running the Tests

### Run all tests:
```bash
python -m pytest test_orca_xtb_sanity.py -v
```

### Run specific test class:
```bash
python -m pytest test_orca_xtb_sanity.py::TestORCABasicORCA60 -v
python -m pytest test_orca_xtb_sanity.py::TestXTBBasicXTB661 -v
python -m pytest test_orca_xtb_sanity.py::TestDataConsistency -v
```

### Run specific test:
```bash
python -m pytest test_orca_xtb_sanity.py::TestORCABasicORCA60::test_dvb_sp_dft_scf_energy -v
```

## Test Data

Tests use actual output files from:
- `data/ORCA/basicORCA6.0/`: ORCA 6.0 calculations
- `data/XTB/basicXTB6.6.1/`: xTB 6.6.1 calculations

## Validation Strategy

### Exact Value Tests
Some tests verify exact values extracted from raw outputs:
- ORCA SP energy: -382.05510839256243 Eh
- XTB SP energy: -26.425939358406 Eh
- ORCA Mulliken charge (atom 0): -0.004225

### Range Tests
Other tests verify values are within reasonable ranges:
- Energies are negative
- Charges sum to total charge (±0.01 e)
- Temperatures match standard conditions (298.15 K)

### Consistency Tests
- Array lengths match declared sizes
- Coordinates are 3D
- Energy ordering follows quantum chemistry principles

### Physical Correctness
- CCSD(T) energy < CCSD energy
- Optimization decreases energy
- Real vibrational modes are positive
- HF and DFT give different results

## Dependencies

- `pytest`: Test framework
- `numpy`: Numerical comparisons
- `cclib`: Parser library being tested

## Future Enhancements

Potential areas for expansion:
1. Add tests for ORCA 5.0 outputs
2. Test error handling with malformed files
3. Add performance benchmarks
4. Test additional calculation types (TD-DFT, CASSCF, etc.)
5. Cross-version compatibility tests
6. Basis set parsing validation
7. Molecular orbital data checks

## Maintenance

When adding new test data:
1. Extract expected values from raw output
2. Add test with descriptive docstring
3. Use `np.testing.assert_allclose` for floating point
4. Document tolerance choices in comments
5. Run full suite to ensure no regressions

## Success Criteria

All 38 tests should pass. Current status: **38/38 PASSED** ✅

Last verified: 2025-11-16
