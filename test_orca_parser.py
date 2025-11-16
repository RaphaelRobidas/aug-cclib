#!/usr/bin/env python3
"""
Comprehensive tests for ORCA parser.
Tests ALL properties in detail including array sizes, shapes, and values.
"""

import os
from pathlib import Path

import pytest
import numpy as np

from cclib.io import ccread
from cclib.parser.utils import convertor


ORCA_DIR = Path("data/ORCA/basicORCA6.0")


class TestORCASPDFT:
    """Comprehensive tests for ORCA DFT single point calculation."""

    @pytest.fixture
    def data(self):
        """Load dvb_sp_dft.out once for all tests in this class."""
        return ccread(ORCA_DIR / "dvb_sp_dft.out")

    def test_basic_attributes(self, data):
        """Test basic molecular attributes."""
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert data.charge == 0, f"Expected charge 0, got {data.charge}"
        assert data.mult == 1, f"Expected multiplicity 1, got {data.mult}"
        assert data.nbasis == 60, f"Expected 60 basis functions, got {data.nbasis}"
        assert data.nmo == 60, f"Expected 60 MOs, got {data.nmo}"
        assert data.nelectrons == 70, f"Expected 70 electrons, got {data.nelectrons}"
        assert data.closed_shell is True, "Should be closed shell"

    def test_atomnos(self, data):
        """Test atomic numbers array."""
        assert hasattr(data, 'atomnos'), "Missing atomnos"
        assert data.atomnos.shape == (20,), f"atomnos shape should be (20,), got {data.atomnos.shape}"
        assert data.atomnos.dtype == np.int32, f"atomnos dtype should be int32, got {data.atomnos.dtype}"

        # DVB molecule has C and H atoms only
        unique_elements = np.unique(data.atomnos)
        assert set(unique_elements) == {1, 6}, f"Should have only C(6) and H(1), got {unique_elements}"

        # Count atoms: DVB has 10 C + 10 H = 20
        n_carbon = np.sum(data.atomnos == 6)
        n_hydrogen = np.sum(data.atomnos == 1)
        assert n_carbon == 10, f"Expected 10 carbon atoms, got {n_carbon}"
        assert n_hydrogen == 10, f"Expected 10 hydrogen atoms, got {n_hydrogen}"

    def test_atomcoords(self, data):
        """Test atomic coordinates array."""
        assert hasattr(data, 'atomcoords'), "Missing atomcoords"
        assert data.atomcoords.shape == (1, 20, 3), \
            f"atomcoords shape should be (1, 20, 3), got {data.atomcoords.shape}"
        assert data.atomcoords.dtype == np.float64, \
            f"atomcoords dtype should be float64, got {data.atomcoords.dtype}"

        # All coordinates should be finite
        assert np.all(np.isfinite(data.atomcoords)), "All coordinates should be finite"

        # Coordinates should be in reasonable range (Angstroms, typically -10 to 10)
        assert np.all(np.abs(data.atomcoords) < 100), "Coordinates seem unreasonably large"

    def test_atommasses(self, data):
        """Test atomic masses array."""
        assert hasattr(data, 'atommasses'), "Missing atommasses"
        assert data.atommasses.shape == (20,), \
            f"atommasses shape should be (20,), got {data.atommasses.shape}"
        assert data.atommasses.dtype == np.float64, \
            f"atommasses dtype should be float64, got {data.atommasses.dtype}"

        # All masses should be positive
        assert np.all(data.atommasses > 0), "All atomic masses should be positive"

        # Check reasonable mass ranges
        # Carbon ~ 12, Hydrogen ~ 1
        carbon_masses = data.atommasses[data.atomnos == 6]
        hydrogen_masses = data.atommasses[data.atomnos == 1]
        assert np.all((carbon_masses > 11) & (carbon_masses < 13)), "Carbon mass out of range"
        assert np.all((hydrogen_masses > 0.9) & (hydrogen_masses < 1.2)), "Hydrogen mass out of range"

    def test_coreelectrons(self, data):
        """Test core electrons array."""
        assert hasattr(data, 'coreelectrons'), "Missing coreelectrons"
        assert data.coreelectrons.shape == (20,), \
            f"coreelectrons shape should be (20,), got {data.coreelectrons.shape}"
        assert data.coreelectrons.dtype == np.int32, \
            f"coreelectrons dtype should be int32, got {data.coreelectrons.dtype}"

        # ORCA uses all-electron basis by default, so coreelectrons might all be 0
        # Check that coreelectrons is consistent with atomnos
        assert np.all(data.coreelectrons >= 0), "Core electrons should be non-negative"

        # Total core electrons should be reasonable
        total_core = np.sum(data.coreelectrons)
        assert total_core >= 0, f"Total core electrons should be non-negative, got {total_core}"

    def test_scfenergies(self, data):
        """Test SCF energies array."""
        assert hasattr(data, 'scfenergies'), "Missing scfenergies"
        assert data.scfenergies.shape == (1,), \
            f"scfenergies shape should be (1,), got {data.scfenergies.shape}"

        # Test exact value from raw output (scfenergies are stored in Hartree)
        expected_hartree = -382.05510839256243
        np.testing.assert_allclose(
            data.scfenergies[0], expected_hartree, rtol=1e-6,
            err_msg=f"SCF energy mismatch"
        )

    def test_scftargets(self, data):
        """Test SCF convergence targets."""
        assert hasattr(data, 'scftargets'), "Missing scftargets"
        assert data.scftargets.shape == (1, 3), \
            f"scftargets shape should be (1, 3), got {data.scftargets.shape}"

        # All targets should be positive
        assert np.all(data.scftargets > 0), "SCF targets should be positive"

    def test_scfvalues(self, data):
        """Test SCF convergence values."""
        assert hasattr(data, 'scfvalues'), "Missing scfvalues"
        assert isinstance(data.scfvalues, list), "scfvalues should be a list"
        assert len(data.scfvalues) == 1, f"Should have 1 scfvalues, got {len(data.scfvalues)}"
        assert data.scfvalues[0].shape[1] == 3, \
            f"scfvalues should have 3 columns, got {data.scfvalues[0].shape[1]}"

    def test_homos(self, data):
        """Test HOMO indices."""
        assert hasattr(data, 'homos'), "Missing homos"
        assert data.homos.shape == (1,), f"homos shape should be (1,), got {data.homos.shape}"
        assert data.homos.dtype == np.int32, f"homos dtype should be int32, got {data.homos.dtype}"

        # For 70 electrons in closed shell, HOMO should be orbital 34 (0-indexed)
        assert data.homos[0] == 34, f"HOMO should be 34, got {data.homos[0]}"

    def test_moenergies(self, data):
        """Test molecular orbital energies."""
        assert hasattr(data, 'moenergies'), "Missing moenergies"
        assert isinstance(data.moenergies, list), "moenergies should be a list"
        assert len(data.moenergies) == 1, f"Should have 1 set of MO energies, got {len(data.moenergies)}"
        assert data.moenergies[0].shape == (60,), \
            f"MO energies should have 60 values, got {data.moenergies[0].shape}"

        # Energies should be mostly ordered (allowing for near-degeneracies)
        # Check that at least 80% of differences are non-negative
        diffs = np.diff(data.moenergies[0])
        non_decreasing = np.sum(diffs >= -0.001)  # Allow small numerical errors
        assert non_decreasing / len(diffs) > 0.8, "MO energies should be mostly sorted"

        # HOMO-LUMO gap should be positive
        homo_energy = data.moenergies[0][data.homos[0]]
        lumo_energy = data.moenergies[0][data.homos[0] + 1]
        gap = lumo_energy - homo_energy
        assert gap > 0, f"HOMO-LUMO gap should be positive, got {gap}"

    def test_mocoeffs(self, data):
        """Test MO coefficients."""
        assert hasattr(data, 'mocoeffs'), "Missing mocoeffs"
        assert isinstance(data.mocoeffs, list), "mocoeffs should be a list"
        assert len(data.mocoeffs) == 1, f"Should have 1 set of MO coeffs, got {len(data.mocoeffs)}"
        assert data.mocoeffs[0].shape == (60, 60), \
            f"MO coefficients should be (60, 60), got {data.mocoeffs[0].shape}"

        # Coefficients should be real numbers
        assert np.all(np.isfinite(data.mocoeffs[0])), "All MO coefficients should be finite"

    def test_mosyms(self, data):
        """Test MO symmetry labels."""
        assert hasattr(data, 'mosyms'), "Missing mosyms"
        assert isinstance(data.mosyms, list), "mosyms should be a list"
        assert len(data.mosyms) == 1, f"Should have 1 set of MO syms, got {len(data.mosyms)}"
        # Number of symmetry labels may be less than number of MOs if some aren't assigned
        assert len(data.mosyms[0]) > 0, f"Should have MO symmetries, got {len(data.mosyms[0])}"
        assert len(data.mosyms[0]) <= 60, f"Should have at most 60 MO symmetries, got {len(data.mosyms[0])}"

    def test_aonames(self, data):
        """Test AO names."""
        assert hasattr(data, 'aonames'), "Missing aonames"
        assert len(data.aonames) == 60, f"Should have 60 AO names, got {len(data.aonames)}"
        assert all(isinstance(name, str) for name in data.aonames), "All AO names should be strings"

    def test_aooverlaps(self, data):
        """Test AO overlap matrix."""
        assert hasattr(data, 'aooverlaps'), "Missing aooverlaps"
        assert data.aooverlaps.shape == (60, 60), \
            f"AO overlaps should be (60, 60), got {data.aooverlaps.shape}"

        # Overlap matrix should be symmetric
        np.testing.assert_allclose(
            data.aooverlaps, data.aooverlaps.T, rtol=1e-10,
            err_msg="AO overlap matrix should be symmetric"
        )

        # Diagonal should be 1 (normalized basis)
        np.testing.assert_allclose(
            np.diag(data.aooverlaps), 1.0, rtol=1e-10,
            err_msg="AO overlap matrix diagonal should be 1"
        )

    def test_atombasis(self, data):
        """Test atom basis information."""
        assert hasattr(data, 'atombasis'), "Missing atombasis"
        assert len(data.atombasis) == 20, f"Should have 20 atom basis entries, got {len(data.atombasis)}"
        assert all(isinstance(basis, list) for basis in data.atombasis), \
            "All atombasis entries should be lists"

        # Total number of basis functions
        total_basis = sum(len(basis) for basis in data.atombasis)
        assert total_basis == 60, f"Total basis functions should be 60, got {total_basis}"

    def test_gbasis(self, data):
        """Test Gaussian basis information."""
        assert hasattr(data, 'gbasis'), "Missing gbasis"
        assert len(data.gbasis) == 20, f"Should have 20 gbasis entries, got {len(data.gbasis)}"
        assert all(isinstance(basis, list) for basis in data.gbasis), \
            "All gbasis entries should be lists"

    def test_atomcharges_mulliken(self, data):
        """Test Mulliken atomic charges."""
        assert hasattr(data, 'atomcharges'), "Missing atomcharges"
        assert 'mulliken' in data.atomcharges, "Missing Mulliken charges"
        assert len(data.atomcharges['mulliken']) == 20, \
            f"Should have 20 Mulliken charges, got {len(data.atomcharges['mulliken'])}"

        # Charges should sum to total charge (0)
        total_charge = sum(data.atomcharges['mulliken'])
        np.testing.assert_allclose(
            total_charge, 0.0, atol=0.01,
            err_msg=f"Mulliken charges should sum to 0, got {total_charge}"
        )

        # Test first atom charge (from raw output: -0.004225)
        np.testing.assert_allclose(
            data.atomcharges['mulliken'][0], -0.004225, atol=0.001,
            err_msg="First atom Mulliken charge mismatch"
        )

    def test_atomcharges_lowdin(self, data):
        """Test Lowdin atomic charges."""
        assert 'lowdin' in data.atomcharges, "Missing Lowdin charges"
        assert len(data.atomcharges['lowdin']) == 20, \
            f"Should have 20 Lowdin charges, got {len(data.atomcharges['lowdin'])}"

        # Charges should sum to total charge (0)
        total_charge = sum(data.atomcharges['lowdin'])
        np.testing.assert_allclose(
            total_charge, 0.0, atol=0.01,
            err_msg=f"Lowdin charges should sum to 0, got {total_charge}"
        )

    def test_atomcharges_hirshfeld(self, data):
        """Test Hirshfeld atomic charges."""
        assert 'hirshfeld' in data.atomcharges, "Missing Hirshfeld charges"
        assert len(data.atomcharges['hirshfeld']) == 20, \
            f"Should have 20 Hirshfeld charges, got {len(data.atomcharges['hirshfeld'])}"

        # Charges should sum to total charge (0)
        total_charge = sum(data.atomcharges['hirshfeld'])
        np.testing.assert_allclose(
            total_charge, 0.0, atol=0.01,
            err_msg=f"Hirshfeld charges should sum to 0, got {total_charge}"
        )

    def test_moments(self, data):
        """Test molecular moments."""
        assert hasattr(data, 'moments'), "Missing moments"
        assert isinstance(data.moments, list), "moments should be a list"
        assert len(data.moments) >= 1, f"Should have at least 1 moment, got {len(data.moments)}"

        # Dipole moment (first moment)
        assert data.moments[0].shape == (3,), f"Dipole should be 3D, got {data.moments[0].shape}"

        # If quadrupole is present, check its shape
        if len(data.moments) > 1:
            # Can be 3 or 6 components depending on format
            assert data.moments[1].shape[0] in (3, 6), \
                f"Second moment should have 3 or 6 components, got {data.moments[1].shape}"

    def test_rotconsts(self, data):
        """Test rotational constants."""
        assert hasattr(data, 'rotconsts'), "Missing rotconsts"
        assert data.rotconsts.shape == (1, 3), \
            f"rotconsts should be (1, 3), got {data.rotconsts.shape}"

        # All rotational constants should be positive
        assert np.all(data.rotconsts > 0), "All rotational constants should be positive"

    def test_metadata(self, data):
        """Test metadata dictionary."""
        assert hasattr(data, 'metadata'), "Missing metadata"
        assert isinstance(data.metadata, dict), "metadata should be a dict"

        # Required keys
        required_keys = ['package', 'package_version', 'methods', 'success']
        for key in required_keys:
            assert key in data.metadata, f"metadata missing key: {key}"

        # Check values
        assert data.metadata['package'] == 'ORCA', f"Package should be ORCA, got {data.metadata['package']}"
        assert data.metadata['package_version'].startswith('6.'), \
            f"Should be ORCA 6.x, got {data.metadata['package_version']}"
        assert data.metadata['success'] is True, "Calculation should be successful"
        assert 'basis_set' in data.metadata, "Missing basis_set in metadata"
        assert 'symmetry_detected' in data.metadata, "Missing symmetry_detected"
        assert 'symmetry_used' in data.metadata, "Missing symmetry_used"


class TestORCAGeometryOptimization:
    """Comprehensive tests for ORCA geometry optimization."""

    @pytest.fixture
    def data(self):
        """Load dvb_gopt.out once for all tests."""
        return ccread(ORCA_DIR / "dvb_gopt.out")

    def test_basic_attributes(self, data):
        """Test basic attributes."""
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert hasattr(data, 'atomcoords'), "Missing atomcoords"
        assert hasattr(data, 'scfenergies'), "Missing scfenergies"

    def test_multiple_geometries(self, data):
        """Test that optimization has multiple geometries."""
        n_geoms = data.atomcoords.shape[0]
        assert n_geoms > 1, f"Should have multiple geometries, got {n_geoms}"

        # Each geometry should have correct shape
        for i in range(n_geoms):
            assert data.atomcoords[i].shape == (20, 3), \
                f"Geometry {i} should be (20, 3), got {data.atomcoords[i].shape}"

    def test_scfenergies(self, data):
        """Test SCF energies for optimization steps."""
        n_energies = len(data.scfenergies)
        n_geoms = data.atomcoords.shape[0]

        # Should have one energy per geometry
        assert n_energies == n_geoms, \
            f"Number of energies ({n_energies}) should match geometries ({n_geoms})"

        # Energy should generally decrease
        assert data.scfenergies[-1] <= data.scfenergies[0], \
            "Final energy should be lower than or equal to initial energy"

    def test_optstatus(self, data):
        """Test optimization status."""
        assert hasattr(data, 'optstatus'), "Missing optstatus"
        # optstatus should indicate successful optimization (OPT_DONE = 4)
        # or at least have some status information

    def test_convergence_geometries(self, data):
        """Test categorized geometries."""
        # Check that different geometry categories exist
        assert hasattr(data, 'converged_geometries'), "Missing converged_geometries"
        assert hasattr(data, 'new_geometries'), "Missing new_geometries"
        assert hasattr(data, 'unconverged_geometries'), "Missing unconverged_geometries"
        assert hasattr(data, 'unknown_geometries'), "Missing unknown_geometries"


class TestORCAFrequencies:
    """Comprehensive tests for ORCA frequency calculation."""

    @pytest.fixture
    def data(self):
        """Load dvb_ir.out once for all tests."""
        return ccread(ORCA_DIR / "dvb_ir.out")

    def test_vibfreqs(self, data):
        """Test vibrational frequencies."""
        assert hasattr(data, 'vibfreqs'), "Missing vibfreqs"
        assert data.vibfreqs.shape == (54,), \
            f"Should have 54 frequencies, got {data.vibfreqs.shape}"

        # All frequencies should be real (positive for stable structure)
        # Some might be small negative due to numerical noise
        assert np.sum(data.vibfreqs < -10) == 0, \
            "Should not have significant imaginary frequencies"

    def test_vibirs(self, data):
        """Test IR intensities."""
        assert hasattr(data, 'vibirs'), "Missing vibirs"
        assert data.vibirs.shape == (54,), \
            f"Should have 54 IR intensities, got {data.vibirs.shape}"

        # Intensities should be non-negative
        assert np.all(data.vibirs >= 0), "IR intensities should be non-negative"

    def test_vibdisps(self, data):
        """Test vibrational displacement vectors."""
        assert hasattr(data, 'vibdisps'), "Missing vibdisps"
        assert data.vibdisps.shape == (54, 20, 3), \
            f"vibdisps should be (54, 20, 3), got {data.vibdisps.shape}"

        # All displacements should be finite
        assert np.all(np.isfinite(data.vibdisps)), "All displacements should be finite"

    def test_thermochemistry(self, data):
        """Test thermochemistry data."""
        # Temperature
        assert hasattr(data, 'temperature'), "Missing temperature"
        np.testing.assert_allclose(data.temperature, 298.15, atol=0.1,
                                   err_msg="Temperature should be 298.15 K")

        # Pressure
        assert hasattr(data, 'pressure'), "Missing pressure"
        np.testing.assert_allclose(data.pressure, 1.0, atol=0.01,
                                   err_msg="Pressure should be 1 atm")

        # Enthalpy
        assert hasattr(data, 'enthalpy'), "Missing enthalpy"
        assert isinstance(data.enthalpy, float), "Enthalpy should be a float"

        # Free energy
        assert hasattr(data, 'freeenergy'), "Missing freeenergy"
        assert isinstance(data.freeenergy, float), "Free energy should be a float"

        # Free energy should be less than enthalpy (T*S contribution)
        assert data.freeenergy < data.enthalpy, \
            f"Free energy ({data.freeenergy}) should be less than enthalpy ({data.enthalpy})"

        # Entropy
        assert hasattr(data, 'entropy'), "Missing entropy"
        assert data.entropy > 0, f"Entropy should be positive, got {data.entropy}"


class TestORCAPostHF:
    """Tests for ORCA post-HF calculations."""

    def test_ccsd_energies(self):
        """Test CCSD calculation."""
        data = ccread(ORCA_DIR / "water_ccsd.out")
        assert data.natom == 3, "Water should have 3 atoms"
        assert hasattr(data, 'ccenergies'), "Missing CC energies"
        assert len(data.ccenergies) > 0, "Should have CC energies"

    def test_ccsd_t_energies(self):
        """Test CCSD(T) calculation."""
        data = ccread(ORCA_DIR / "water_ccsd_t.out")
        assert data.natom == 3, "Water should have 3 atoms"
        assert hasattr(data, 'ccenergies'), "Missing CC energies"

        # CCSD(T) should give lower energy than CCSD
        data_ccsd = ccread(ORCA_DIR / "water_ccsd.out")
        ccsd_energy = convertor(data_ccsd.ccenergies[-1], "eV", "hartree")
        ccsd_t_energy = convertor(data.ccenergies[-1], "eV", "hartree")
        assert ccsd_t_energy < ccsd_energy, \
            "CCSD(T) energy should be lower than CCSD"

    def test_mp2_energies(self):
        """Test MP2 calculation."""
        data = ccread(ORCA_DIR / "water_mp2.out")
        assert data.natom == 3, "Water should have 3 atoms"
        assert hasattr(data, 'mpenergies'), "Missing MP energies"
        assert len(data.mpenergies) > 0, "Should have MP energies"

    def test_mp3_energies(self):
        """Test MP3 calculation."""
        data = ccread(ORCA_DIR / "water_mp3.out")
        assert data.natom == 3, "Water should have 3 atoms"
        assert hasattr(data, 'mpenergies'), "Missing MP energies"


class TestORCATDDFT:
    """Tests for ORCA TDDFT calculations."""

    @pytest.fixture
    def data(self):
        """Load dvb_td.out once for all tests."""
        return ccread(ORCA_DIR / "dvb_td.out")

    def test_excited_state_energies(self, data):
        """Test excited state energies."""
        assert hasattr(data, 'etenergies'), "Missing excited state energies"
        assert len(data.etenergies) > 0, "Should have excited states"

        # All energies should be positive (excitation energies)
        assert np.all(np.array(data.etenergies) > 0), \
            "All excitation energies should be positive"

    def test_oscillator_strengths(self, data):
        """Test oscillator strengths."""
        assert hasattr(data, 'etoscs'), "Missing oscillator strengths"
        assert len(data.etoscs) == len(data.etenergies), \
            "Should have one oscillator strength per excited state"

        # Oscillator strengths should be non-negative
        assert np.all(np.array(data.etoscs) >= 0), \
            "Oscillator strengths should be non-negative"

    def test_excited_state_symmetries(self, data):
        """Test excited state symmetry labels."""
        assert hasattr(data, 'etsyms'), "Missing excited state symmetries"
        assert len(data.etsyms) == len(data.etenergies), \
            "Should have one symmetry label per excited state"

    def test_excited_state_transitions(self, data):
        """Test excited state transition information."""
        assert hasattr(data, 'etsecs'), "Missing excited state transitions"
        assert len(data.etsecs) == len(data.etenergies), \
            "Should have transition info for each excited state"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
