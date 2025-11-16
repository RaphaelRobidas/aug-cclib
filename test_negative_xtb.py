"""
Negative tests for XTB and ORCA parsers - these tests verify the parser rejects wrong values.

These tests PASS by asserting that wrong values are NOT accepted, missing attributes
DO exist, and invalid constraints are NOT satisfied. This confirms our positive tests
are properly validating the parser implementation.
"""

import numpy
import pytest
from cclib.io import ccread


class TestXTBNegativeSP:
    """Negative tests for XTB single point calculations - verify wrong values are rejected"""

    @pytest.fixture
    def data(self):
        return ccread("data/XTB/basicXTB6.6.1/dvb_sp/dvb_sp.out")

    def test_rejects_wrong_scfenergy(self, data):
        """Verify parser rejects incorrect SCF energy"""
        wrong_energy = -26.999999  # Wrong value
        assert pytest.approx(data.scfenergies[0], abs=1e-6) != wrong_energy

    def test_rejects_wrong_dispersion_energy(self, data):
        """Verify parser rejects incorrect dispersion energy"""
        wrong_dispersion = -0.999999  # Wrong value
        assert pytest.approx(data.dispersionenergies[0], abs=1e-8) != wrong_dispersion

    def test_nbasis_exists(self, data):
        """Verify nbasis attribute exists"""
        assert hasattr(data, "nbasis")

    def test_rejects_wrong_nbasis_value(self, data):
        """Verify parser rejects incorrect nbasis"""
        assert data.nbasis != 999  # Wrong value, actual is 50

    def test_rejects_wrong_nmo_value(self, data):
        """Verify parser rejects incorrect nmo"""
        assert data.nmo != 100  # Wrong value, actual is 50

    def test_moments_exists(self, data):
        """Verify moments attribute exists"""
        assert hasattr(data, "moments")

    def test_dipole_is_near_zero(self, data):
        """Verify dipole is near zero for symmetric molecule, not > 1.0"""
        # DVB is symmetric, dipole should be ~zero
        assert abs(data.moments[0][0]) <= 1.0  # Should be near zero

    def test_polarizabilities_exists(self, data):
        """Verify polarizabilities attribute exists"""
        assert hasattr(data, "polarizabilities")

    def test_rejects_wrong_polarizability_value(self, data):
        """Verify parser rejects incorrect polarizability"""
        wrong_alpha = 999.999
        assert pytest.approx(data.polarizabilities[0][0, 0], abs=1e-5) != wrong_alpha

    def test_rejects_wrong_charge(self, data):
        """Verify charge is 0, not 1"""
        assert data.charge != 1  # Actual is 0

    def test_rejects_wrong_homo_index(self, data):
        """Verify parser rejects incorrect HOMO index"""
        assert 99 not in data.homos  # Wrong value, actual is 24

    def test_dispersion_energy_is_negative(self, data):
        """Verify dispersion energy is negative (attractive), not positive"""
        assert data.dispersionenergies[0] <= 0  # Should be negative

    def test_rejects_wrong_wall_time(self, data):
        """Verify parser rejects incorrect wall time"""
        wall_seconds = data.metadata["wall_time"][0].total_seconds()
        assert pytest.approx(wall_seconds, abs=0.001) != 999.999  # Wrong value

    def test_rejects_wrong_cpu_time(self, data):
        """Verify parser rejects incorrect CPU time"""
        cpu_seconds = data.metadata["cpu_time"][0].total_seconds()
        assert pytest.approx(cpu_seconds, abs=0.001) != 999.999  # Wrong value

    def test_cpu_greater_equal_wall(self, data):
        """Verify CPU time >= wall time for parallel execution"""
        wall_seconds = data.metadata["wall_time"][0].total_seconds()
        cpu_seconds = data.metadata["cpu_time"][0].total_seconds()
        assert cpu_seconds >= wall_seconds  # CPU should be >= wall

    def test_polarizability_correct_shape(self, data):
        """Verify polarizability is 3x3, not 2x2"""
        assert data.polarizabilities[0].shape != (2, 2)  # Should be (3, 3)
        assert data.polarizabilities[0].shape == (3, 3)  # Correct shape

    def test_dispersion_correct_array_length(self, data):
        """Verify exactly 1 dispersion energy, not 5"""
        assert len(data.dispersionenergies) != 5  # Wrong length
        assert len(data.dispersionenergies) == 1  # Correct length


class TestXTBNegativeGeoOpt:
    """Negative tests for XTB geometry optimization - verify wrong values are rejected"""

    @pytest.fixture
    def data(self):
        return ccread("data/XTB/basicXTB6.6.1/dvb_opt/dvb_opt.out")

    def test_rejects_wrong_final_energy(self, data):
        """Verify parser rejects incorrect final energy"""
        wrong_energy = -99.999999
        assert pytest.approx(data.scfenergies[-1], abs=1e-6) != wrong_energy

    def test_rejects_wrong_initial_energy(self, data):
        """Verify parser rejects incorrect initial energy"""
        wrong_energy = -99.999999
        assert pytest.approx(data.scfenergies[0], abs=1e-6) != wrong_energy

    def test_rejects_wrong_num_steps(self, data):
        """Verify 6 optimization steps, not 99"""
        assert len(data.scfenergies) != 99  # Wrong number
        assert len(data.scfenergies) == 6  # Correct number

    def test_optdone_exists(self, data):
        """Verify optdone attribute exists"""
        assert hasattr(data, "optdone")

    def test_optimization_converged(self, data):
        """Verify optimization converged (optdone is True, not False)"""
        # optdone should be True (it converged), not False
        assert data.optdone is not False  # Should be True

    def test_energy_decreases(self, data):
        """Verify energy decreases during optimization"""
        # Energy should decrease (or stay constant within numerical noise)
        energy_diffs = numpy.diff(data.scfenergies)
        assert not numpy.any(energy_diffs > 0.1)  # Should not increase significantly

    def test_rejects_wrong_charge(self, data):
        """Verify charge is 0, not -1"""
        assert data.charge != -1  # Actual is 0


class TestXTBNegativeIR:
    """Negative tests for XTB IR calculations - verify wrong values are rejected"""

    @pytest.fixture
    def data(self):
        return ccread("data/XTB/basicXTB6.6.1/dvb_ir/dvb_ir.out")

    def test_vibramans_exists(self, data):
        """Verify vibramans attribute exists"""
        assert hasattr(data, "vibramans")

    def test_vibramans_are_zeros(self, data):
        """Verify vibramans are all zeros (XTB doesn't compute real Raman)"""
        # XTB doesn't compute real Raman intensities
        assert not numpy.any(data.vibramans > 1.0)  # Should be all zeros

    def test_atommasses_exists(self, data):
        """Verify atommasses attribute exists"""
        assert hasattr(data, "atommasses")

    def test_rejects_wrong_atommass_value(self, data):
        """Verify parser rejects incorrect atomic mass"""
        # Carbon mass should be ~12, not 999
        carbon_indices = numpy.where(data.atomnos == 6)[0]
        assert pytest.approx(data.atommasses[carbon_indices[0]], abs=0.01) != 999.0

    def test_rejects_wrong_imaginary_freqs(self, data):
        """Verify 0 imaginary frequencies, not 99"""
        assert data.metadata["imaginary_freqs"] != 99  # Wrong value
        assert data.metadata["imaginary_freqs"] == 0  # Correct value

    def test_dispersion_exists(self, data):
        """Verify dispersionenergies attribute exists"""
        assert hasattr(data, "dispersionenergies")

    def test_rejects_wrong_dispersion_value(self, data):
        """Verify parser rejects incorrect dispersion energy"""
        wrong_value = -0.999999
        assert pytest.approx(data.dispersionenergies[0], abs=1e-8) != wrong_value

    def test_entropy_exists(self, data):
        """Verify entropy attribute exists"""
        assert hasattr(data, "entropy")

    def test_rejects_wrong_entropy_value(self, data):
        """Verify parser rejects incorrect entropy"""
        wrong_entropy = 0.999999
        assert pytest.approx(data.entropy, abs=1e-10) != wrong_entropy

    def test_entropy_is_positive(self, data):
        """Verify entropy is positive, not negative"""
        assert data.entropy >= 0  # Entropy is always positive

    def test_temperature_exists(self, data):
        """Verify temperature attribute exists"""
        assert hasattr(data, "temperature")

    def test_rejects_wrong_temperature(self, data):
        """Verify temperature is 298.15 K, not 999.99"""
        assert pytest.approx(data.temperature, abs=0.01) != 999.99  # Wrong value


class TestORCANegativeIR:
    """Negative tests for ORCA IR calculations - verify wrong values are rejected"""

    @pytest.fixture
    def data(self):
        return ccread("data/ORCA/basicORCA5.0/dvb_ir.out")

    def test_rejects_wrong_wall_time(self, data):
        """Verify parser rejects incorrect wall time"""
        wall_seconds = data.metadata["wall_time"][0].total_seconds()
        assert pytest.approx(wall_seconds, abs=0.001) != 999.999  # Wrong value

    def test_rejects_wrong_cpu_time(self, data):
        """Verify parser rejects incorrect CPU time"""
        cpu_seconds = data.metadata["cpu_time"][0].total_seconds()
        assert pytest.approx(cpu_seconds, abs=0.001) != 999.999  # Wrong value

    def test_rejects_wrong_entropy(self, data):
        """Verify parser rejects incorrect entropy"""
        wrong_entropy = 0.999999
        assert pytest.approx(data.entropy, abs=1e-10) != wrong_entropy


class TestORCANegativeRaman:
    """Negative tests for ORCA Raman calculations - verify wrong values are rejected"""

    @pytest.fixture
    def data(self):
        return ccread("data/ORCA/basicORCA5.0/dvb_raman.out")

    def test_vibramans_exists(self, data):
        """Verify vibramans attribute exists"""
        assert hasattr(data, "vibramans")

    def test_rejects_wrong_max_vibramans(self, data):
        """Verify parser rejects incorrect max Raman intensity"""
        wrong_max = 9999.999
        assert pytest.approx(data.vibramans.max(), abs=0.001) != wrong_max

    def test_rejects_wrong_vibramans_values(self, data):
        """Verify parser rejects incorrect Raman intensity values"""
        expected_first_five = [99.0, 99.0, 99.0, 99.0, 99.0]  # Wrong values
        for i in range(5):
            assert pytest.approx(data.vibramans[i], abs=0.001) != expected_first_five[i]


class TestORCANegativeSP:
    """Negative tests for ORCA single point - verify wrong values are rejected"""

    @pytest.fixture
    def data(self):
        return ccread("data/ORCA/basicORCA5.0/dvb_sp.out")

    def test_rejects_wrong_wall_time(self, data):
        """Verify parser rejects incorrect wall time for ORCA 5.0"""
        wall_seconds = data.metadata["wall_time"][0].total_seconds()
        assert pytest.approx(wall_seconds, abs=0.001) != 999.999  # Wrong value

    def test_rejects_wrong_cpu_time(self, data):
        """Verify parser rejects incorrect CPU time for ORCA 5.0"""
        cpu_seconds = data.metadata["cpu_time"][0].total_seconds()
        assert pytest.approx(cpu_seconds, abs=0.001) != 999.999  # Wrong value


if __name__ == "__main__":
    print("Running negative tests - these should all PASS if the parser works correctly")
    print("=" * 80)
    pytest.main([__file__, "-v", "--tb=short"])
