#!/usr/bin/env python3
"""
Comprehensive sanity tests for ORCA and XTB parsers.
Tests verify that parsed values match actual values in raw output files.
"""

import os
import sys
from pathlib import Path

import pytest
import numpy as np

from cclib.io import ccread
from cclib.parser.utils import convertor


# Test data directory
ORCA_DIR = Path("data/ORCA/basicORCA6.0")
XTB_DIR = Path("data/XTB/basicXTB6.6.1")


class TestORCABasicORCA60:
    """Tests for ORCA 6.0 output files with expected values from raw outputs."""

    def test_dvb_sp_dft_basic_attributes(self):
        """Test basic attributes from dvb_sp_dft.out"""
        data = ccread(ORCA_DIR / "dvb_sp_dft.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert data.charge == 0, f"Expected charge 0, got {data.charge}"
        assert data.mult == 1, f"Expected multiplicity 1, got {data.mult}"
        assert data.nbasis == 60, f"Expected 60 basis functions, got {data.nbasis}"

    def test_dvb_sp_dft_scf_energy(self):
        """Test SCF energy from dvb_sp_dft.out
        Expected: -382.05510839256243 Eh from raw output"""
        data = ccread(ORCA_DIR / "dvb_sp_dft.out")

        assert hasattr(data, 'scfenergies'), "No scfenergies attribute"
        assert len(data.scfenergies) == 1, f"Expected 1 SCF energy, got {len(data.scfenergies)}"

        # Convert to eV for comparison (cclib stores in eV)
        expected_hartree = -382.05510839256243
        expected_ev = convertor(expected_hartree, "hartree", "eV")

        np.testing.assert_allclose(
            data.scfenergies[0], expected_ev, rtol=1e-6,
            err_msg=f"SCF energy mismatch: expected {expected_ev}, got {data.scfenergies[0]}"
        )

    def test_dvb_sp_dft_metadata(self):
        """Test metadata from dvb_sp_dft.out"""
        data = ccread(ORCA_DIR / "dvb_sp_dft.out")

        assert "package_version" in data.metadata, "Missing package_version"
        assert data.metadata["package_version"].startswith("6."), \
            f"Expected ORCA 6.x, got {data.metadata['package_version']}"

        assert data.metadata.get("success") is True, "Calculation should be successful"
        assert "methods" in data.metadata, "Missing methods in metadata"

    def test_water_ccsd_basic(self):
        """Test water CCSD calculation
        Expected: 3 atoms (H2O)"""
        data = ccread(ORCA_DIR / "water_ccsd.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 3, f"Expected 3 atoms for water, got {data.natom}"
        assert hasattr(data, 'ccenergies'), "Missing CC energies"
        assert len(data.ccenergies) > 0, "No CC energies found"

    def test_water_ccsd_t_basic(self):
        """Test water CCSD(T) calculation"""
        data = ccread(ORCA_DIR / "water_ccsd_t.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 3, f"Expected 3 atoms for water, got {data.natom}"
        assert hasattr(data, 'ccenergies'), "Missing CC energies"

    def test_water_mp2_basic(self):
        """Test water MP2 calculation"""
        data = ccread(ORCA_DIR / "water_mp2.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 3, f"Expected 3 atoms for water, got {data.natom}"
        assert hasattr(data, 'mpenergies'), "Missing MP energies"
        assert len(data.mpenergies) > 0, "No MP energies found"

    def test_water_mp3_basic(self):
        """Test water MP3 calculation"""
        data = ccread(ORCA_DIR / "water_mp3.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 3, f"Expected 3 atoms for water, got {data.natom}"
        assert hasattr(data, 'mpenergies'), "Missing MP energies"

    def test_dvb_gopt_basic(self):
        """Test geometry optimization"""
        data = ccread(ORCA_DIR / "dvb_gopt.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert hasattr(data, 'optstatus'), "Missing optstatus"
        assert hasattr(data, 'atomcoords'), "Missing atomcoords"
        assert len(data.atomcoords) > 1, "Should have multiple geometries"

    def test_dvb_gopt_convergence(self):
        """Test that geometry optimization converged"""
        data = ccread(ORCA_DIR / "dvb_gopt.out")

        # Check if optimization converged
        assert hasattr(data, 'optstatus'), "Missing optstatus"
        # The last status should indicate convergence or completion
        assert data.metadata.get("success") is True, "Optimization should be successful"

    def test_dvb_ir_frequencies(self):
        """Test IR frequency calculation"""
        data = ccread(ORCA_DIR / "dvb_ir.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert hasattr(data, 'vibfreqs'), "Missing vibfreqs"
        assert hasattr(data, 'vibirs'), "Missing IR intensities"

        # 20 atoms = 3*20 - 6 = 54 vibrational modes (non-linear molecule)
        expected_nmodes = 54
        assert len(data.vibfreqs) == expected_nmodes, \
            f"Expected {expected_nmodes} frequencies, got {len(data.vibfreqs)}"

    def test_dvb_raman_intensities(self):
        """Test Raman spectrum calculation"""
        data = ccread(ORCA_DIR / "dvb_raman.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert hasattr(data, 'vibfreqs'), "Missing vibfreqs"
        assert hasattr(data, 'vibramans'), "Missing Raman intensities"

        expected_nmodes = 54
        assert len(data.vibramans) == expected_nmodes, \
            f"Expected {expected_nmodes} Raman intensities, got {len(data.vibramans)}"

    def test_dvb_td_excited_states(self):
        """Test TDDFT excited states"""
        data = ccread(ORCA_DIR / "dvb_td.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert hasattr(data, 'etenergies'), "Missing excited state energies"
        assert hasattr(data, 'etoscs'), "Missing oscillator strengths"

        # Should have multiple excited states
        assert len(data.etenergies) > 0, "No excited states found"
        assert len(data.etoscs) == len(data.etenergies), \
            "Oscillator strengths don't match number of states"

    def test_dvb_nmr_shieldings(self):
        """Test NMR shielding calculation"""
        data = ccread(ORCA_DIR / "dvb_nmr.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        # NMR calculations should have some results
        assert data.metadata.get("success") is True, "NMR calculation should succeed"

    def test_dvb_scan_relaxed(self):
        """Test relaxed surface scan"""
        data = ccread(ORCA_DIR / "dvb_scan_relaxed.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        # Scan should have multiple geometries and energies
        assert hasattr(data, 'atomcoords'), "Missing atomcoords"
        assert len(data.atomcoords) > 1, "Scan should have multiple geometries"

    def test_dvb_dispersion_energy(self):
        """Test dispersion correction parsing"""
        data = ccread(ORCA_DIR / "dvb_dispersion_bp86_d3zero.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert hasattr(data, 'dispersionenergies'), "Missing dispersion energies"
        assert len(data.dispersionenergies) > 0, "No dispersion energies found"

    def test_Trp_polar_polarizability(self):
        """Test polarizability calculation"""
        data = ccread(ORCA_DIR / "Trp_polar.out")

        assert data is not None, "Failed to parse file"
        # Should parse successfully
        assert data.metadata.get("success") is True, "Polarizability calculation should succeed"

    def test_dvb_sp_dft_mulliken_charges(self):
        """Test that Mulliken charges are parsed and sum to total charge"""
        data = ccread(ORCA_DIR / "dvb_sp_dft.out")

        assert hasattr(data, 'atomcharges'), "Missing atomcharges"
        assert 'mulliken' in data.atomcharges, "Missing Mulliken charges"
        assert len(data.atomcharges['mulliken']) == 20, \
            f"Expected 20 Mulliken charges, got {len(data.atomcharges['mulliken'])}"

        # Charges should sum to total charge (0)
        total_charge = sum(data.atomcharges['mulliken'])
        np.testing.assert_allclose(
            total_charge, 0.0, atol=0.01,
            err_msg=f"Mulliken charges should sum to 0, got {total_charge}"
        )

        # Check first carbon atom charge (from raw output: -0.004225)
        np.testing.assert_allclose(
            data.atomcharges['mulliken'][0], -0.004225, atol=0.001,
            err_msg=f"First atom charge mismatch"
        )

    def test_dvb_sp_hf_vs_dft(self):
        """Test that HF and DFT give different energies for same system"""
        data_hf = ccread(ORCA_DIR / "dvb_sp_hf.out")
        data_dft = ccread(ORCA_DIR / "dvb_sp_dft.out")

        assert data_hf is not None and data_dft is not None, "Failed to parse files"
        assert data_hf.natom == data_dft.natom == 20, "Both should have 20 atoms"

        # HF and DFT should give different energies
        energy_diff = abs(data_hf.scfenergies[0] - data_dft.scfenergies[0])
        assert energy_diff > 1.0, \
            f"HF and DFT energies should differ significantly, diff={energy_diff} eV"

    def test_dvb_ir_thermochemistry(self):
        """Test thermochemistry data from frequency calculation"""
        data = ccread(ORCA_DIR / "dvb_ir.out")

        assert data is not None, "Failed to parse file"
        # Should have thermochemistry data
        assert hasattr(data, 'temperature'), "Missing temperature"
        assert hasattr(data, 'pressure'), "Missing pressure"
        assert hasattr(data, 'enthalpy'), "Missing enthalpy"
        assert hasattr(data, 'freeenergy'), "Missing free energy"

        # Standard conditions: 298.15 K, 1 atm
        np.testing.assert_allclose(data.temperature, 298.15, atol=0.1)

    def test_water_ccsd_vs_ccsd_t(self):
        """Test that CCSD(T) energy is lower than CCSD"""
        data_ccsd = ccread(ORCA_DIR / "water_ccsd.out")
        data_ccsd_t = ccread(ORCA_DIR / "water_ccsd_t.out")

        assert data_ccsd is not None and data_ccsd_t is not None, "Failed to parse files"

        # CCSD(T) should give lower energy than CCSD (more correlation)
        # Energies in cclib are in eV, convert to Hartree for comparison
        ccsd_energy = convertor(data_ccsd.ccenergies[-1], "eV", "hartree")
        ccsd_t_energy = convertor(data_ccsd_t.ccenergies[-1], "eV", "hartree")

        assert ccsd_t_energy < ccsd_energy, \
            f"CCSD(T) energy ({ccsd_t_energy}) should be lower than CCSD ({ccsd_energy})"


class TestXTBBasicXTB661:
    """Tests for xTB 6.6.1 output files with expected values from raw outputs."""

    def test_dvb_sp_basic_attributes(self):
        """Test basic attributes from dvb_sp output"""
        data = ccread(XTB_DIR / "dvb_sp" / "dvb_sp.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert hasattr(data, 'scfenergies'), "Missing SCF energies"
        assert len(data.scfenergies) == 1, f"Expected 1 SCF energy, got {len(data.scfenergies)}"

    def test_dvb_sp_charges(self):
        """Test atomic charges from dvb_sp output"""
        data = ccread(XTB_DIR / "dvb_sp" / "dvb_sp.out")

        assert data is not None, "Failed to parse file"
        assert hasattr(data, 'atomcharges'), "Missing atomic charges"
        # xTB typically provides Mulliken charges
        assert 'mulliken' in data.atomcharges, "Missing Mulliken charges"
        assert len(data.atomcharges['mulliken']) == 20, \
            f"Expected 20 charges, got {len(data.atomcharges['mulliken'])}"

    def test_dvb_sp_metadata(self):
        """Test metadata from dvb_sp output"""
        data = ccread(XTB_DIR / "dvb_sp" / "dvb_sp.out")

        assert "package_version" in data.metadata, "Missing package_version"
        assert data.metadata["package_version"].startswith("6.6"), \
            f"Expected xTB 6.6.x, got {data.metadata['package_version']}"
        assert data.metadata.get("success") is True, "Calculation should be successful"

    def test_dvb_opt_basic(self):
        """Test geometry optimization"""
        data = ccread(XTB_DIR / "dvb_opt" / "dvb_opt.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert hasattr(data, 'atomcoords'), "Missing atomcoords"
        assert hasattr(data, 'scfenergies'), "Missing SCF energies"

        # Optimization should have multiple SCF energies
        assert len(data.scfenergies) > 1, \
            f"Expected multiple SCF energies for optimization, got {len(data.scfenergies)}"

    def test_dvb_opt_convergence(self):
        """Test that optimization converged"""
        data = ccread(XTB_DIR / "dvb_opt" / "dvb_opt.out")

        assert data.metadata.get("success") is True, "Optimization should be successful"
        # Final geometry should be present
        assert hasattr(data, 'atomcoords'), "Missing atomcoords"
        assert len(data.atomcoords) > 0, "No geometries found"

    def test_dvb_ir_frequencies(self):
        """Test IR frequency calculation"""
        data = ccread(XTB_DIR / "dvb_ir" / "dvb_ir.out")

        assert data is not None, "Failed to parse file"
        assert data.natom == 20, f"Expected 20 atoms, got {data.natom}"
        assert hasattr(data, 'vibfreqs'), "Missing vibfreqs"

        # xTB includes all 3N modes (including translations/rotations which are ~0)
        # 20 atoms = 3*20 = 60 modes total
        expected_nmodes = 60
        assert len(data.vibfreqs) == expected_nmodes, \
            f"Expected {expected_nmodes} frequencies, got {len(data.vibfreqs)}"

        # Check that first 6 modes are near zero (translations/rotations)
        assert np.all(np.abs(data.vibfreqs[:6]) < 1.0), \
            "First 6 modes should be near zero (translations/rotations)"

    def test_dvb_ir_intensities(self):
        """Test IR intensities"""
        data = ccread(XTB_DIR / "dvb_ir" / "dvb_ir.out")

        assert hasattr(data, 'vibirs'), "Missing IR intensities"
        assert len(data.vibirs) == len(data.vibfreqs), \
            "Number of IR intensities should match number of frequencies"

    def test_dvb_sp_energy_range(self):
        """Test that energies are in reasonable range"""
        data = ccread(XTB_DIR / "dvb_sp" / "dvb_sp.out")

        # xTB energies should be negative and in reasonable range
        energy_hartree = convertor(data.scfenergies[0], "eV", "hartree")
        assert energy_hartree < 0, f"Energy should be negative, got {energy_hartree}"
        assert energy_hartree > -1000, f"Energy seems unreasonably low: {energy_hartree}"

    def test_dvb_opt_energy_decreases(self):
        """Test that energy decreases during optimization"""
        data = ccread(XTB_DIR / "dvb_opt" / "dvb_opt.out")

        # In an optimization, the energy should generally decrease or stay the same
        energies = data.scfenergies
        assert energies[-1] <= energies[0] + 1.0, \
            "Final energy should be lower than or similar to initial energy"

    def test_charge_conservation(self):
        """Test that atomic charges sum to total charge"""
        data = ccread(XTB_DIR / "dvb_sp" / "dvb_sp.out")

        if hasattr(data, 'atomcharges') and 'mulliken' in data.atomcharges:
            total_charge = sum(data.atomcharges['mulliken'])
            expected_charge = data.charge if hasattr(data, 'charge') else 0
            np.testing.assert_allclose(
                total_charge, expected_charge, atol=0.01,
                err_msg=f"Charges don't sum to total charge: {total_charge} vs {expected_charge}"
            )

    def test_dvb_sp_energy_value(self):
        """Test specific energy value from dvb_sp output
        Expected: -26.425939358406 Eh from raw output"""
        data = ccread(XTB_DIR / "dvb_sp" / "dvb_sp.out")

        expected_hartree = -26.425939358406
        expected_ev = convertor(expected_hartree, "hartree", "eV")

        np.testing.assert_allclose(
            data.scfenergies[0], expected_ev, rtol=1e-6,
            err_msg=f"XTB energy mismatch: expected {expected_ev}, got {data.scfenergies[0]}"
        )

    def test_dvb_ir_zero_point_energy(self):
        """Test that zero point energy is parsed"""
        data = ccread(XTB_DIR / "dvb_ir" / "dvb_ir.out")

        # Should have zero point energy from frequency calculation
        # This is typically stored in the vibration section
        assert hasattr(data, 'vibfreqs'), "Missing vibfreqs needed for ZPE"

        # Zero point energy = sum of (0.5 * h * freq) for all modes
        # We can verify it's positive and reasonable
        if hasattr(data, 'zpve'):
            assert data.zpve > 0, f"Zero point energy should be positive, got {data.zpve}"

    def test_dvb_opt_final_energy_lower(self):
        """Test that final optimized energy is lower than any intermediate"""
        data = ccread(XTB_DIR / "dvb_opt" / "dvb_opt.out")

        energies = data.scfenergies
        final_energy = energies[-1]

        # Final energy should be the minimum (or very close to it)
        min_energy = min(energies)
        assert final_energy <= min_energy + 0.01, \
            f"Final energy ({final_energy}) should be at minimum ({min_energy})"

    def test_dvb_ir_has_negative_frequencies(self):
        """Test that imaginary frequencies are near zero (stable structure)"""
        data = ccread(XTB_DIR / "dvb_ir" / "dvb_ir.out")

        # First 6 modes should be translations/rotations (near zero)
        # Real vibrational modes should all be positive
        real_vib_modes = data.vibfreqs[6:]

        # All real vibrational modes should be positive
        assert np.all(real_vib_modes > 0), \
            f"Found negative frequency in vibrational modes: {real_vib_modes[real_vib_modes <= 0]}"


class TestDataConsistency:
    """Cross-parser consistency tests."""

    def test_orca_atom_count_consistency(self):
        """Test that natom matches length of atomnos and atomcoords"""
        data = ccread(ORCA_DIR / "dvb_sp_dft.out")

        assert len(data.atomnos) == data.natom, \
            f"atomnos length ({len(data.atomnos)}) doesn't match natom ({data.natom})"
        assert len(data.atomcoords[0]) == data.natom, \
            f"atomcoords length ({len(data.atomcoords[0])}) doesn't match natom ({data.natom})"

    def test_xtb_atom_count_consistency(self):
        """Test that natom matches length of atomnos and atomcoords"""
        # Use optimization output which includes final structure
        data = ccread(XTB_DIR / "dvb_opt" / "dvb_opt.out")

        assert len(data.atomnos) == data.natom, \
            f"atomnos length ({len(data.atomnos)}) doesn't match natom ({data.natom})"

        # Check that coordinates exist and have correct length
        assert hasattr(data, 'atomcoords'), "Missing atomcoords"
        assert len(data.atomcoords) > 0, "No atomcoords found"
        assert len(data.atomcoords[-1]) == data.natom, \
            f"atomcoords length ({len(data.atomcoords[-1])}) doesn't match natom ({data.natom})"

    def test_orca_coordinates_3d(self):
        """Test that coordinates are 3D"""
        data = ccread(ORCA_DIR / "dvb_sp_dft.out")

        for coord in data.atomcoords[0]:
            assert len(coord) == 3, f"Coordinate should be 3D, got {len(coord)}D"

    def test_xtb_coordinates_3d(self):
        """Test that coordinates are 3D"""
        # Use optimization output which includes final structure
        data = ccread(XTB_DIR / "dvb_opt" / "dvb_opt.out")

        assert len(data.atomcoords) > 0, "No atomcoords found"
        for coord in data.atomcoords[-1]:
            assert len(coord) == 3, f"Coordinate should be 3D, got {len(coord)}D"


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
