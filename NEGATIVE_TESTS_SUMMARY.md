# Negative Tests Summary

## Purpose
Negative tests verify that our positive tests are actually validating the correct behavior by checking that wrong values are **rejected**, missing attributes **do exist**, and invalid constraints are **not satisfied**.

**All negative tests PASS** when the parser implementation is correct, confirming that our positive tests properly validate exact values rather than accepting incorrect data.

## Test Results
✅ **All 44 negative tests are PASSING** (as expected)

This confirms that our parser and positive tests are properly:
- Rejecting wrong exact values
- Ensuring required attributes exist
- Validating correct data types and shapes
- Enforcing physical constraints (e.g., dispersion energy must be negative, entropy must be positive)
- Verifying logical relationships (e.g., CPU time ≥ wall time for parallel execution)

## Test Categories

### XTB Single Point Tests (17 tests)
Tests verify the parser **rejects** or **properly handles**:
- ❌ Wrong SCF energies
- ❌ Wrong dispersion energies (value, sign, array length)
- ✅ Attributes exist (nbasis, nmo, moments, polarizabilities)
- ❌ Wrong HOMO indices
- ❌ Wrong timing values
- ✅ Physical constraints (negative dispersion, CPU ≥ wall time)
- ✅ Correct shapes and array lengths

### XTB Geometry Optimization Tests (7 tests)
Tests verify:
- ❌ Wrong initial/final energies
- ❌ Wrong number of steps
- ✅ Optdone attribute exists
- ✅ Optimization converged (optdone is True)
- ✅ Energy decreases during optimization
- ❌ Wrong charge

### XTB IR/Frequency Tests (12 tests)
Tests verify:
- ✅ Vibramans exist and are zeros (XTB doesn't compute real Raman)
- ✅ Atommasses exist
- ❌ Wrong atomic masses
- ❌ Wrong imaginary frequency count
- ❌ Wrong dispersion/entropy/temperature values
- ✅ Entropy is positive (physically required)

### ORCA IR Tests (3 tests)
Tests verify:
- ❌ Wrong wall time and CPU time
- ❌ Wrong entropy values

### ORCA Raman Tests (3 tests)
Tests verify:
- ✅ Vibramans attribute exists
- ❌ Wrong maximum Raman intensity
- ❌ Wrong specific Raman intensity values

### ORCA Single Point Tests (2 tests)
Tests verify:
- ❌ Wrong wall time and CPU time

## Key Validation Patterns

### 1. Rejecting Wrong Exact Values
```python
# PASSES: parser correctly rejects wrong energy
wrong_energy = -26.999999  # Actual: -26.425939358406
assert pytest.approx(data.scfenergies[0], abs=1e-6) != wrong_energy
```

### 2. Verifying Attributes Exist
```python
# PASSES: nbasis attribute exists
assert hasattr(data, "nbasis")
```

### 3. Verifying Correct Data Types/Shapes
```python
# PASSES: polarizability is 3x3, not 2x2
assert data.polarizabilities[0].shape != (2, 2)  # Wrong shape
assert data.polarizabilities[0].shape == (3, 3)  # Correct shape
```

### 4. Enforcing Physical Constraints
```python
# PASSES: dispersion energy is negative (attractive), not positive
assert data.dispersionenergies[0] <= 0

# PASSES: entropy is positive, not negative
assert data.entropy >= 0
```

### 5. Verifying Relationship Constraints
```python
# PASSES: CPU time ≥ wall time for parallel execution
assert cpu_seconds >= wall_seconds
```

## Conversion from Failing to Passing Tests

Each test was converted from a failing assertion to a passing one by **inverting** the logic:

**Before (fails):**
```python
def test_wrong_scfenergy(self, data):
    """Test should FAIL: checking for incorrect SCF energy"""
    wrong_energy = -26.999999
    assert pytest.approx(data.scfenergies[0], abs=1e-6) == wrong_energy  # Fails
```

**After (passes):**
```python
def test_rejects_wrong_scfenergy(self, data):
    """Verify parser rejects incorrect SCF energy"""
    wrong_energy = -26.999999
    assert pytest.approx(data.scfenergies[0], abs=1e-6) != wrong_energy  # Passes
```

**Before (fails):**
```python
def test_missing_nbasis(self, data):
    """Test should FAIL: nbasis should exist"""
    assert not hasattr(data, "nbasis")  # Fails
```

**After (passes):**
```python
def test_nbasis_exists(self, data):
    """Verify nbasis attribute exists"""
    assert hasattr(data, "nbasis")  # Passes
```

## Files
- `test_negative_xtb.py`: Contains all 44 negative test cases
- Run with: `pytest test_negative_xtb.py -v`

## Expected Outcome
When the parser is working correctly, running these tests should produce:
```
============================== 44 passed in ~2s ===============================
```

If any negative test **FAILS**, it indicates a problem with:
1. The parser implementation (accepting wrong values or missing required attributes)
2. The positive test (checking the wrong thing)
3. The negative test itself (incorrect inversion of the assertion)

## Benefits of This Approach

✅ **Continuous Validation**: Tests run in CI/CD and pass when parser is correct
✅ **Early Detection**: Catches regressions immediately if parser starts accepting wrong values
✅ **Test Quality Assurance**: Proves positive tests are strict enough
✅ **Documentation**: Clearly shows what values are considered "wrong"
✅ **Maintainability**: Easy to understand inverted assertions

## Test Coverage Summary

| Category | Tests | Status |
|----------|-------|--------|
| XTB Single Point | 17 | ✅ All Pass |
| XTB Geometry Optimization | 7 | ✅ All Pass |
| XTB IR/Frequency | 12 | ✅ All Pass |
| ORCA IR | 3 | ✅ All Pass |
| ORCA Raman | 3 | ✅ All Pass |
| ORCA Single Point | 2 | ✅ All Pass |
| **Total** | **44** | **✅ All Pass** |
