"""
Negative tests for XTB parser - these tests SHOULD FAIL if the parser is working correctly.

These tests verify that our positive tests are actually checking the right things
by intentionally checking for wrong values, missing attributes, or incorrect data.
"""

import numpy
import pytest
from cclib.io import ccread


class TestXTBNegativeSP:
    """Negative tests for XTB single point calculations - these should all FAIL"""

    @pytest.fixture
    def data(self):
        return ccread("data/XTB/basicXTB6.6.1/dvb_sp/dvb_sp.out")

    def test_wrong_scfenergy(self, data):
        """Test should FAIL: checking for incorrect SCF energy"""
        wrong_energy = -26.999999  # Wrong value
        assert pytest.approx(data.scfenergies[0], abs=1e-6) == wrong_energy

    def test_wrong_dispersion_energy(self, data):
        """Test should FAIL: checking for incorrect dispersion energy"""
        wrong_dispersion = -0.999999  # Wrong value
        assert pytest.approx(data.dispersionenergies[0], abs=1e-8) == wrong_dispersion

    def test_missing_nbasis(self, data):
        """Test should FAIL: nbasis should exist"""
        assert not hasattr(data, "nbasis")

    def test_wrong_nbasis_value(self, data):
        """Test should FAIL: checking for incorrect nbasis"""
        assert data.nbasis == 999  # Wrong value, should be 50

    def test_wrong_nmo_value(self, data):
        """Test should FAIL: checking for incorrect nmo"""
        assert data.nmo == 100  # Wrong value, should be 50

    def test_missing_moments(self, data):
        """Test should FAIL: moments should exist"""
        assert not hasattr(data, "moments")

    def test_wrong_dipole_nonzero(self, data):
        """Test should FAIL: dipole should be zero for symmetric molecule"""
        # DVB is symmetric, dipole should be ~zero
        assert abs(data.moments[0][0]) > 1.0  # Wrong: should be near zero

    def test_missing_polarizabilities(self, data):
        """Test should FAIL: polarizabilities should exist"""
        assert not hasattr(data, "polarizabilities")

    def test_wrong_polarizability_value(self, data):
        """Test should FAIL: checking for incorrect polarizability"""
        wrong_alpha = 999.999
        assert pytest.approx(data.polarizabilities[0][0, 0], abs=1e-5) == wrong_alpha

    def test_wrong_charge(self, data):
        """Test should FAIL: charge should be 0, not 1"""
        assert data.charge == 1  # Wrong value

    def test_wrong_homo_index(self, data):
        """Test should FAIL: checking for incorrect HOMO index"""
        assert 99 in data.homos  # Wrong value, should be 24

    def test_positive_dispersion_energy(self, data):
        """Test should FAIL: dispersion energy should be negative (attractive)"""
        assert data.dispersionenergies[0] > 0  # Wrong: should be negative

    def test_wrong_wall_time(self, data):
        """Test should FAIL: checking for incorrect wall time"""
        wall_seconds = data.metadata["wall_time"][0].total_seconds()
        assert pytest.approx(wall_seconds, abs=0.001) == 999.999  # Wrong value

    def test_wrong_cpu_time(self, data):
        """Test should FAIL: checking for incorrect CPU time"""
        cpu_seconds = data.metadata["cpu_time"][0].total_seconds()
        assert pytest.approx(cpu_seconds, abs=0.001) == 999.999  # Wrong value

    def test_cpu_less_than_wall(self, data):
        """Test should FAIL: CPU time should be >= wall time for parallel execution"""
        wall_seconds = data.metadata["wall_time"][0].total_seconds()
        cpu_seconds = data.metadata["cpu_time"][0].total_seconds()
        assert cpu_seconds < wall_seconds  # Wrong: CPU should be >= wall

    def test_wrong_polarizability_shape(self, data):
        """Test should FAIL: polarizability should be 3x3"""
        assert data.polarizabilities[0].shape == (2, 2)  # Wrong shape

    def test_wrong_dispersion_array_length(self, data):
        """Test should FAIL: should have exactly 1 dispersion energy"""
        assert len(data.dispersionenergies) == 5  # Wrong length


class TestXTBNegativeGeoOpt:
    """Negative tests for XTB geometry optimization - these should all FAIL"""

    @pytest.fixture
    def data(self):
        return ccread("data/XTB/basicXTB6.6.1/dvb_opt/dvb_opt.out")

    def test_wrong_final_energy(self, data):
        """Test should FAIL: checking for incorrect final energy"""
        wrong_energy = -99.999999
        assert pytest.approx(data.scfenergies[-1], abs=1e-6) == wrong_energy

    def test_wrong_initial_energy(self, data):
        """Test should FAIL: checking for incorrect initial energy"""
        wrong_energy = -99.999999
        assert pytest.approx(data.scfenergies[0], abs=1e-6) == wrong_energy

    def test_wrong_num_steps(self, data):
        """Test should FAIL: should have 6 optimization steps"""
        assert len(data.scfenergies) == 99  # Wrong number

    def test_missing_optdone(self, data):
        """Test should FAIL: optdone should exist"""
        assert not hasattr(data, "optdone")

    def test_wrong_convergence_step(self, data):
        """Test should FAIL: optimization should have converged"""
        # optdone should be True (it converged), not False
        assert data.optdone is False  # Wrong value

    def test_energy_increases(self, data):
        """Test should FAIL: energy should decrease during optimization"""
        # Check if energy increases (it shouldn't)
        energy_diffs = numpy.diff(data.scfenergies)
        assert numpy.any(energy_diffs > 0.1)  # Wrong: should decrease

    def test_wrong_charge(self, data):
        """Test should FAIL: charge should be 0, not -1"""
        assert data.charge == -1  # Wrong value


class TestXTBNegativeIR:
    """Negative tests for XTB IR calculations - these should all FAIL"""

    @pytest.fixture
    def data(self):
        return ccread("data/XTB/basicXTB6.6.1/dvb_ir/dvb_ir.out")

    def test_missing_vibramans(self, data):
        """Test should FAIL: vibramans should exist"""
        assert not hasattr(data, "vibramans")

    def test_wrong_vibramans_all_nonzero(self, data):
        """Test should FAIL: vibramans should be all zeros for XTB"""
        # XTB doesn't compute real Raman intensities
        assert numpy.any(data.vibramans > 1.0)  # Wrong: should be all zeros

    def test_missing_atommasses(self, data):
        """Test should FAIL: atommasses should exist"""
        assert not hasattr(data, "atommasses")

    def test_wrong_atommass_value(self, data):
        """Test should FAIL: checking for incorrect atomic mass"""
        # Carbon mass should be ~12, not 999
        carbon_indices = numpy.where(data.atomnos == 6)[0]
        assert pytest.approx(data.atommasses[carbon_indices[0]], abs=0.01) == 999.0

    def test_wrong_imaginary_freqs(self, data):
        """Test should FAIL: should have 0 imaginary frequencies"""
        assert data.metadata["imaginary_freqs"] == 99  # Wrong value

    def test_missing_dispersion(self, data):
        """Test should FAIL: dispersionenergies should exist"""
        assert not hasattr(data, "dispersionenergies")

    def test_wrong_dispersion_value(self, data):
        """Test should FAIL: checking for incorrect dispersion energy"""
        wrong_value = -0.999999
        assert pytest.approx(data.dispersionenergies[0], abs=1e-8) == wrong_value

    def test_missing_entropy(self, data):
        """Test should FAIL: entropy should exist"""
        assert not hasattr(data, "entropy")

    def test_wrong_entropy_value(self, data):
        """Test should FAIL: checking for incorrect entropy"""
        wrong_entropy = 0.999999
        assert pytest.approx(data.entropy, abs=1e-10) == wrong_entropy

    def test_negative_entropy(self, data):
        """Test should FAIL: entropy should be positive"""
        assert data.entropy < 0  # Wrong: entropy is always positive

    def test_missing_temperature(self, data):
        """Test should FAIL: temperature should exist"""
        assert not hasattr(data, "temperature")

    def test_wrong_temperature(self, data):
        """Test should FAIL: temperature should be 298.15 K"""
        assert pytest.approx(data.temperature, abs=0.01) == 999.99  # Wrong value


class TestORCANegativeIR:
    """Negative tests for ORCA IR calculations - these should all FAIL"""

    @pytest.fixture
    def data(self):
        return ccread("data/ORCA/basicORCA5.0/dvb_ir.out")

    def test_wrong_wall_time(self, data):
        """Test should FAIL: checking for incorrect wall time"""
        wall_seconds = data.metadata["wall_time"][0].total_seconds()
        assert pytest.approx(wall_seconds, abs=0.001) == 999.999  # Wrong value

    def test_wrong_cpu_time(self, data):
        """Test should FAIL: checking for incorrect CPU time"""
        cpu_seconds = data.metadata["cpu_time"][0].total_seconds()
        assert pytest.approx(cpu_seconds, abs=0.001) == 999.999  # Wrong value

    def test_wrong_entropy(self, data):
        """Test should FAIL: checking for incorrect entropy"""
        wrong_entropy = 0.999999
        assert pytest.approx(data.entropy, abs=1e-10) == wrong_entropy


class TestORCANegativeRaman:
    """Negative tests for ORCA Raman calculations - these should all FAIL"""

    @pytest.fixture
    def data(self):
        return ccread("data/ORCA/basicORCA5.0/dvb_raman.out")

    def test_missing_vibramans(self, data):
        """Test should FAIL: vibramans should exist"""
        assert not hasattr(data, "vibramans")

    def test_wrong_max_vibramans(self, data):
        """Test should FAIL: checking for incorrect max Raman intensity"""
        wrong_max = 9999.999
        assert pytest.approx(data.vibramans.max(), abs=0.001) == wrong_max

    def test_wrong_vibramans_values(self, data):
        """Test should FAIL: checking for incorrect Raman intensity values"""
        expected_first_five = [99.0, 99.0, 99.0, 99.0, 99.0]  # Wrong values
        for i in range(5):
            assert pytest.approx(data.vibramans[i], abs=0.001) == expected_first_five[i]


class TestORCANegativeSP:
    """Negative tests for ORCA single point - these should all FAIL"""

    @pytest.fixture
    def data(self):
        return ccread("data/ORCA/basicORCA5.0/dvb_sp.out")

    def test_wrong_wall_time(self, data):
        """Test should FAIL: checking for incorrect wall time for ORCA 5.0"""
        wall_seconds = data.metadata["wall_time"][0].total_seconds()
        assert pytest.approx(wall_seconds, abs=0.001) == 999.999  # Wrong value

    def test_wrong_cpu_time(self, data):
        """Test should FAIL: checking for incorrect CPU time for ORCA 5.0"""
        cpu_seconds = data.metadata["cpu_time"][0].total_seconds()
        assert pytest.approx(cpu_seconds, abs=0.001) == 999.999  # Wrong value


if __name__ == "__main__":
    print("Running negative tests - these should all FAIL if the parser works correctly")
    print("=" * 80)
    pytest.main([__file__, "-v", "--tb=short"])
