# Negative Tests Summary

## Purpose
Negative tests verify that our positive tests are actually validating the correct behavior by intentionally checking for wrong values, missing attributes, or incorrect data types. **All negative tests should FAIL** if the parser implementation is correct.

## Test Results
✅ **All 44 negative tests are FAILING** (as expected)

This confirms that our positive tests are properly validating:
- Exact reference values (not accepting wrong values)
- Attribute existence (not passing when attributes are missing)
- Data types and shapes
- Physical constraints (e.g., dispersion energy must be negative, entropy must be positive)

## Test Categories

### XTB Single Point Tests (17 tests)
Tests verify parsing of:
- SCF energies (exact values, not approximations)
- Dispersion energies (wrong values, wrong sign, wrong array length)
- Basis functions (nbasis, nmo)
- Dipole moments
- Polarizabilities (values and matrix shape)
- Molecular charge
- HOMO indices
- Timing data (wall time, CPU time, and their relationship)

### XTB Geometry Optimization Tests (7 tests)
Tests verify:
- Initial and final energies
- Number of optimization steps
- Convergence status (optdone attribute)
- Energy monotonic decrease
- Molecular charge

### XTB IR/Frequency Tests (12 tests)
Tests verify:
- Raman intensities (all zeros for XTB)
- Atomic masses (XTB-specific values)
- Imaginary frequency count
- Dispersion energy
- Thermodynamic properties (entropy, temperature)
- Entropy sign (must be positive)

### ORCA IR Tests (3 tests)
Tests verify:
- Wall time and CPU time with exact reference values
- Entropy values

### ORCA Raman Tests (3 tests)
Tests verify:
- Raman intensity existence
- Maximum Raman intensity value
- Specific Raman intensity values

### ORCA Single Point Tests (2 tests)
Tests verify:
- Wall time and CPU time with exact reference values

## Key Validation Patterns

### 1. Wrong Exact Values
```python
# Should FAIL: checking for wrong energy
wrong_energy = -26.999999  # Actual: -26.425939358406
assert pytest.approx(data.scfenergies[0], abs=1e-6) == wrong_energy
```

### 2. Missing Attributes
```python
# Should FAIL: nbasis should exist
assert not hasattr(data, "nbasis")
```

### 3. Wrong Data Types/Shapes
```python
# Should FAIL: polarizability should be 3x3, not 2x2
assert data.polarizabilities[0].shape == (2, 2)
```

### 4. Physical Constraints
```python
# Should FAIL: dispersion energy should be negative (attractive)
assert data.dispersionenergies[0] > 0

# Should FAIL: entropy should be positive
assert data.entropy < 0
```

### 5. Relationship Constraints
```python
# Should FAIL: CPU time should be >= wall time for parallel execution
assert cpu_seconds < wall_seconds
```

## Files
- `test_negative_xtb.py`: Contains all 44 negative tests
- Run with: `pytest test_negative_xtb.py -v`

## Expected Outcome
When the parser is working correctly, running these tests should produce:
```
============================== 44 failed in ~2s ===============================
```

If any negative test PASSES, it indicates a problem with either:
1. The parser implementation (not validating correctly)
2. The positive test (too lenient or not checking the right thing)
3. The negative test itself (checking for something that legitimately shouldn't exist)
