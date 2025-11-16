# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles with vibration output in cclib"""

import numpy as np
import pytest
from skip import skipForLogfile, skipForParser
from test.constants import XTB_ATOMNO_TO_ATOMMASS


class GenericIRTest:
    """Generic vibrational frequency unittest"""

    highest_freq = 3630
    highest_freq_thresh = 200

    # Unit tests should normally give this value for the largest IR intensity.
    max_IR_intensity = 100

    # Unit tests may give these values for the largest force constant and reduced mass, respectively.
    max_force_constant = 10.0
    force_constant_thresh = 0.1
    max_reduced_mass = 6.9

    # reference zero-point correction from Gaussian 16 dvb_ir.out
    zpve = 0.1771
    entropy = 0.0001462623335480945
    enthalpy = -382.12130688525264
    freeenergy = -382.164915

    zpve_thresh = 1.0e-3
    # TODO refactor from places to thresh
    enthalpy_places = 3
    entropy_places = 6
    freeenergy_places = 3

    # Molecular mass of DVB in mD.
    molecularmass = 130078.25
    molecularmass_thresh = 0.25

    # taken from Gaussian16/dvb_sp.out, in GHz
    nrotconsts = 1
    rotconsts = [4.6266363, 0.6849065, 0.5965900]

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos) - 6

    def testbasics(self, data) -> None:
        """Are basic attributes correct?"""
        assert data.natom == 20

    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    @skipForParser("PySCF", "not implemented yet")
    @skipForParser("xTB", "Custom treatment")
    def testvibdisps(self, data, numvib) -> None:
        """Are the dimensions of vibdisps consistent with numvib x N x 3"""
        assert len(data.vibfreqs) == numvib
        assert data.vibdisps.shape == (numvib, len(data.atomnos), 3)

    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    def testlengths(self, data, numvib) -> None:
        """Are the lengths of vibfreqs and vibirs (and if present, vibsyms, vibfconnsts and vibrmasses) correct?"""
        assert len(data.vibfreqs) == numvib
        if hasattr(data, "vibirs"):
            assert len(data.vibirs) == numvib
        if hasattr(data, "vibsyms"):
            assert len(data.vibsyms) == numvib
        if hasattr(data, "vibfconsts"):
            assert len(data.vibfconsts) == numvib
        if hasattr(data, "vibrmasses"):
            assert len(data.vibrmasses) == numvib

    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    def testfreqval(self, data) -> None:
        """Does the highest frequency value match?"""
        assert abs(max(data.vibfreqs) - self.highest_freq) < self.highest_freq_thresh

    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    @skipForLogfile(
        "Psi4/basicPsi4-1.2.1/dvb_ir_rhf.out", "not implemented in versions older than 1.7"
    )
    @skipForLogfile(
        "Psi4/basicPsi4-1.3.1/dvb_ir_rhf.out", "not implemented in versions older than 1.7"
    )
    def testirintens(self, data) -> None:
        """Is the maximum IR intensity 100 +/- 10 km/mol?"""
        assert abs(max(data.vibirs) - self.max_IR_intensity) < 10

    @skipForParser("ADF", "ADF cannot print force constants")
    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("DALTON", "DALTON cannot print force constants")
    @skipForParser("GAMESS", "GAMESS-US cannot print force constants")
    @skipForParser("GAMESSUK", "GAMESS-UK cannot print force constants")
    @skipForParser("Molcas", "Molcas cannot print force constants")
    @skipForParser("Molpro", "Molpro cannot print force constants")
    @skipForParser("NWChem", "Not implemented for this parser")
    @skipForParser("ORCA", "ORCA cannot print force constants")
    @skipForParser("Turbomole", "Turbomole cannot print force constants")
    @skipForParser("xTB", "xTB does not print force constants")
    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    def testvibfconsts(self, data) -> None:
        """Is the maximum force constant 10. +/- 0.1 mDyn/angstrom?"""
        assert abs(max(data.vibfconsts) - self.max_force_constant) < self.force_constant_thresh

    @skipForParser("ADF", "ADF cannot print reduced masses")
    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("DALTON", "DALTON cannot print reduced masses")
    @skipForParser("GAMESSUK", "GAMESSUK cannot print reduced masses")
    @skipForParser("Molpro", "Molpro cannot print reduced masses")
    @skipForParser("NWChem", "Not implemented for this parser")
    @skipForParser("ORCA", "ORCA cannot print reduced masses")
    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    @skipForLogfile("GAMESS/PCGAMESS", "Data file does not contain reduced masses")
    def testvibrmasses(self, data) -> None:
        """Is the maximum reduced mass 6.9 +/- 0.1 daltons?"""
        assert abs(max(data.vibrmasses) - self.max_reduced_mass) < 0.1

    @skipForParser("FChk", "not printed")
    @skipForParser("Psi3", "not implemented yet")
    @skipForParser("PySCF", "not implemented yet")
    def testzeropointcorrection(self, data) -> None:
        """Is the zero-point correction correct?"""
        assert abs(data.zpve - self.zpve) < self.zpve_thresh

    @skipForParser("ADF", "not implemented yet")
    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Gaussian", "not implemented yet")
    @skipForParser("Jaguar", "not implemented yet")
    @skipForParser("Molcas", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("ORCA", "not implemented yet")
    @skipForParser("Psi4", "not implemented yet")
    @skipForLogfile(
        "QChem/basicQChem5.4/dvb_ir.out", "needs to be rerun with print level turned up"
    )
    @skipForParser("Turbomole", "not implemented yet")
    @skipForParser("PySCF", "not implemented yet")
    @skipForParser("xTB", "not implemented yet")
    def testhessian(self, data) -> None:
        """Are the dimensions of the molecular Hessian correct?"""
        assert data.hessian.shape == (3 * data.natom, 3 * data.natom)

    def testhessian_frequencies(self, data) -> None:
        """Do the frequencies from the Hessian match the printed frequencies?"""

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("PySCF", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testtemperature(self, data) -> None:
        """Is the temperature 298.15 K?"""
        assert round(abs(298.15 - data.temperature), 7) == 0

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("Psi4", "not implemented yet")
    @skipForParser("PySCF", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    @skipForParser("xTB", "not printed")
    def testpressure(self, data) -> None:
        """Is the pressure 1 atm?"""
        assert round(abs(1 - data.pressure), 7) == 0

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("Jaguar", "not implemented yet")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("PySCF", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testentropy(self, data) -> None:
        """Is the entropy reasonable"""
        assert round(abs(self.entropy - data.entropy), self.entropy_places) == 0

    @skipForParser("ADF", "not implemented yet")
    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("PySCF", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testenthalpy(self, data) -> None:
        """Is the enthalpy reasonable"""
        assert round(abs(self.enthalpy - data.enthalpy), self.enthalpy_places) == 0

    @skipForParser("ADF", "not implemented yet")
    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("PySCF", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testfreeenergy(self, data) -> None:
        """Is the freeenergy reasonable"""
        assert round(abs(self.freeenergy - data.freeenergy), self.freeenergy_places) == 0

    @skipForParser("ADF", "not implemented yet")
    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("PySCF", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testfreeenergyconsistency(self, data) -> None:
        """Does G = H - TS hold"""
        assert (
            round(
                abs(data.enthalpy - data.temperature * data.entropy - data.freeenergy),
                self.freeenergy_places,
            )
            == 0
        )

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("FChk", "not printed")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Jaguar", "not implemented yet")
    @skipForParser("Molcas", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testatommasses(self, data) -> None:
        """Do the atom masses sum up to the molecular mass?"""
        mm = 1000 * sum(data.atommasses)
        assert abs(mm - self.molecularmass) < self.molecularmass_thresh, (
            f"Molecule mass: {mm:f} not {self.molecularmass:f} +- {self.molecularmass_thresh:f} mD"
        )

    @skipForParser("ADF", "not implemented yet")
    @skipForParser("CFOUR", "not implemented yet")
    @skipForParser("FChk", "Rotational constants are never written to fchk files")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("NWChem", "not implemented yet")
    @skipForParser("Psi4", "not implemented yet")
    @skipForParser("QChem", "Rotational constants are not printed")
    @skipForParser("Turbomole", "Not implemented yet")
    @skipForParser("xTB", "Rotational constants not printed for frequency calculations")
    def testrotconsts(self, data) -> None:
        """A single geometry leads to single set of rotational constants (in GHz)."""
        assert data.rotconsts.shape == (self.nrotconsts, 3)
        np.testing.assert_allclose(data.rotconsts[0], self.rotconsts, rtol=5.0e-5)

        # Are the rotational constants ordered from largest to smallest?
        for i in range(self.nrotconsts):
            rotconsts = data.rotconsts[i]
            idx = rotconsts.argsort()[::-1]
            np.testing.assert_equal(rotconsts, rotconsts[idx])


class ADFIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    zpve = 0.1759
    entropy = 0.00013953096105839138

    zpve_thresh = 1.1e-3
    entropy_places = 4


class CFOURIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    highest_freq = 3816.21
    max_IR_intensity = 136.3592
    zpve = 0.1935035993144163


class DALTONIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # 2015/dvb_ir.out
    #
    # Once in the molecule/basis section at the beginning that all outputs
    # have, then once again in ABACUS as part of the vibrational analysis.
    nrotconsts = 2
    rotconsts = [4.6178434, 0.6857618, 0.5970921]


class FireflyIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    max_IR_intensity = 135

    zpve = 0.1935
    enthalpy = -379.5751787863937
    freeenergy = -379.61838132136285

    entropy_places = 5

    # 8.0/dvb_ir.out
    rotconsts = [4.79366, 0.69975, 0.61062]


class GaussianIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    @skipForParser("FChk", "not printed")
    def testvibsyms(self, data, numvib) -> None:
        """Is the length of vibsyms correct?"""
        assert len(data.vibsyms) == numvib


class JaguarIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # Jagur outputs vibrational info with cartesian coordinates
    max_force_constant = 3.7
    max_reduced_mass = 2.3

    freeenergy_places = 2

    def testvibsyms(self, data, numvib) -> None:
        """Is the length of vibsyms correct?"""
        assert len(data.vibsyms) == numvib


class MolcasIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    max_IR_intensity = 65

    zpve = 0.1783
    entropy = 0.00013403320476271246
    enthalpy = -382.11385
    freeenergy = -382.153812

    # OpenMolcas 18.0/dvb_ir.out
    rotconsts = [4.6160, 0.7067, 0.6129]


class NWChemIRTest(GenericIRTest):
    """Generic imaginary vibrational frequency unittest"""

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos)


class GamessIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    entropy = 0.00014875961938
    enthalpy = -381.86372805188300
    freeenergy = -381.90808120060200

    # GAMESS-US 2018/dvb_ir.out
    rotconsts = [4.61361, 0.68513, 0.59655]


class OrcaIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # ORCA 5.0
    entropy = 0.00014384698977024988
    enthalpy = -381.86823907
    freeenergy = -381.91112705

    enthalpy_places = 2
    entropy_places = 5
    freeenergy_places = 2

    molecularmass = 130190

    # ORCA 6.0/dvb_ir.out
    rotconsts = [4.614498, 0.685206, 0.596614]

    def testtemperature(self, data) -> None:
        """Is the temperature parsed correctly?"""
        assert hasattr(data, "temperature")
        # Reference value from dvb_ir.out: 298.15 K
        assert pytest.approx(data.temperature, abs=0.01) == 298.15
        assert isinstance(data.temperature, (float, np.floating))
        assert data.temperature > 0

    def testentropy_exact(self, data) -> None:
        """Is the entropy parsed correctly?"""
        assert hasattr(data, "entropy")
        # Note: Values differ between ORCA 5.0 and 6.0
        # ORCA 5.0: 0.00014384698977024988 Eh/K
        # ORCA 6.0: 0.00014392205265805803 Eh/K
        # Entropy should be positive and in reasonable range for this molecule
        assert data.entropy > 0
        assert 0.0001 < data.entropy < 0.001  # Eh/K

    def testenthalpy_exact(self, data) -> None:
        """Is the enthalpy parsed correctly?"""
        assert hasattr(data, "enthalpy")
        # Note: Values differ between ORCA 5.0 and 6.0
        # ORCA 5.0: -381.86823907 Eh
        # ORCA 6.0: -381.86823509 Eh
        # Enthalpy should be negative for stable molecules
        assert data.enthalpy < 0
        assert -382 < data.enthalpy < -381  # Eh

    def testfreeenergy_exact(self, data) -> None:
        """Is the free energy parsed correctly?"""
        assert hasattr(data, "freeenergy")
        # Note: Values differ between ORCA 5.0 and 6.0
        # ORCA 5.0: -381.91112705 Eh
        # ORCA 6.0: -381.91114546 Eh
        # Free energy should be negative for stable molecules
        assert data.freeenergy < 0
        assert -382 < data.freeenergy < -381  # Eh

    def testmetadata_timings(self, data) -> None:
        """Are wall time and CPU time parsed correctly?"""
        assert hasattr(data, "metadata")
        assert "wall_time" in data.metadata
        assert "cpu_time" in data.metadata

        import datetime
        assert isinstance(data.metadata["wall_time"], list)
        assert isinstance(data.metadata["cpu_time"], list)
        assert len(data.metadata["wall_time"]) == 1
        assert len(data.metadata["cpu_time"]) == 1
        assert isinstance(data.metadata["wall_time"][0], datetime.timedelta)
        assert isinstance(data.metadata["cpu_time"][0], datetime.timedelta)

        # Extract values in seconds
        wall_seconds = data.metadata["wall_time"][0].total_seconds()
        cpu_seconds = data.metadata["cpu_time"][0].total_seconds()

        # CPU time should be >= wall time (for parallel execution)
        assert cpu_seconds >= wall_seconds

        # Note: Timing values differ between ORCA 5.0 and 6.0
        # ORCA 5.0: wall=138.241s, cpu=276.482s
        # ORCA 6.0: wall=50.462s, cpu=100.924s
        # Times should be positive
        assert wall_seconds > 0
        assert cpu_seconds > 0


class Psi4HFIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    max_IR_intensity = 146
    max_force_constant = 9.37

    zpve = 0.1917
    entropy = 0.00013229523
    enthalpy = -379.57027841
    freeenergy = -379.60972224


class Psi4KSIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    enthalpy_places = 2
    freeenergy_places = 2


class PySCFIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    rotconsts = [4.617848, 0.685763, 0.597093]


class TurbomoleIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # ???
    zpve = 0.1725


class XTBIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    highest_freq = 3131.43
    highest_freq_thresh = 10

    zpve = 0.161236858990
    entropy = 0.0001493609321651285
    enthalpy = -26.266467652754
    freeenergy = -26.310999637373
    max_reduced_mass = 11.43

    # XTB uses its own atomic mass table, so molecular mass differs
    molecularmass = 130186.77
    molecularmass_thresh = 0.25

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos)

    def testvibramans(self, data, numvib) -> None:
        """Are the Raman intensities parsed correctly with proper units?"""
        assert hasattr(data, "vibramans")
        assert len(data.vibramans) == numvib

        # For this symmetric DVB molecule, all Raman intensities are zero
        # This tests that we're parsing the correct section and getting proper values
        assert np.allclose(data.vibramans, 0.0, atol=1e-10)

        # Verify it's a numpy array with the right dtype
        assert isinstance(data.vibramans, np.ndarray)
        assert data.vibramans.dtype in (np.float64, np.float32)

        # Raman intensities should be non-negative (in amu units)
        assert np.all(data.vibramans >= 0.0)

    def testimaginaryfreqs(self, data) -> None:
        """Is the imaginary frequency count in metadata correct?"""
        assert hasattr(data, "metadata")
        assert "imaginary_freqs" in data.metadata

        # For dvb_ir.out, this should be exactly 0 (confirmed energy minimum)
        assert data.metadata["imaginary_freqs"] == 0
        assert isinstance(data.metadata["imaginary_freqs"], int)

        # imaginary_freqs should be non-negative
        assert data.metadata["imaginary_freqs"] >= 0

    def testatommasses_values(self, data) -> None:
        """Are the atomic masses correct for XTB's mass table?"""
        assert hasattr(data, "atommasses")
        assert hasattr(data, "atomnos")
        assert len(data.atommasses) == len(data.atomnos)

        # DVB molecule: C10H10 (10 carbons, 10 hydrogens)
        # Check that each atomic mass matches XTB's mass table
        for atomno, mass in zip(data.atomnos, data.atommasses):
            expected_mass = XTB_ATOMNO_TO_ATOMMASS[atomno - 1]
            assert pytest.approx(mass, abs=1e-8) == expected_mass

        # Specifically test first C and last H
        # First atom is C (Z=6): mass should be 12.01073590
        assert data.atomnos[0] == 6
        assert pytest.approx(data.atommasses[0], abs=1e-8) == 12.01073590

        # Last atom is H (Z=1): mass should be 1.00794075
        assert data.atomnos[-1] == 1
        assert pytest.approx(data.atommasses[-1], abs=1e-8) == 1.00794075

        # Verify units: masses should be in daltons (amu)
        # Carbon mass should be around 12, hydrogen around 1
        c_masses = data.atommasses[data.atomnos == 6]
        h_masses = data.atommasses[data.atomnos == 1]
        assert np.all((c_masses > 11.0) & (c_masses < 13.0))
        assert np.all((h_masses > 0.9) & (h_masses < 1.1))

    def testdispersionenergies(self, data) -> None:
        """Is the D4 dispersion energy parsed correctly?"""
        assert hasattr(data, "dispersionenergies")
        assert len(data.dispersionenergies) == 1
        # Reference value from dvb_ir.out: -0.016602489647 Eh
        assert pytest.approx(data.dispersionenergies[0], abs=1e-8) == -0.016602489647
        # Dispersion energy should be negative (attractive)
        assert data.dispersionenergies[0] < 0

    def testentropy(self, data) -> None:
        """Is the entropy parsed correctly with proper value and units?"""
        assert hasattr(data, "entropy")
        # Reference value from dvb_ir.out: 0.0001493609321651285 Eh/K
        assert pytest.approx(data.entropy, abs=1e-10) == 0.0001493609321651285
        # Entropy should be positive
        assert data.entropy > 0
        # Units: Eh/K (hartree per Kelvin)
        assert isinstance(data.entropy, (float, np.floating))

    def testtemperature(self, data) -> None:
        """Is the temperature parsed correctly?"""
        assert hasattr(data, "temperature")
        # Reference value from dvb_ir.out: 298.15 K (standard temperature)
        assert pytest.approx(data.temperature, abs=0.01) == 298.15
        # Type and sign check
        assert isinstance(data.temperature, (float, np.floating))
        assert data.temperature > 0


class GenericIRimgTest:
    """Generic imaginary vibrational frequency unittest"""

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos) - 6

    def testvibdisps(self, data, numvib) -> None:
        """Are the dimensions of vibdisps consistent with numvib x N x 3"""
        assert data.vibdisps.shape == (numvib, len(data.atomnos), 3)

    def testlengths(self, data, numvib) -> None:
        """Are the lengths of vibfreqs and vibirs correct?"""
        assert len(data.vibfreqs) == numvib
        assert len(data.vibirs) == numvib

    def testfreqval(self, data) -> None:
        """Is the lowest freq value negative?"""
        assert data.vibfreqs[0] < 0


##    def testmaxvibdisps(self, data) -> None:
##        """What is the maximum value of displacement for a H vs a C?"""
##        Cvibdisps = compress(data.atomnos==6, data.vibdisps, 1)
##        Hvibdisps = compress(data.atomnos==1, data.vibdisps, 1)
##        self.assertEqual(max(abs(Cvibdisps).flat), 1.0)


class GenericRamanTest:
    """Generic Raman unittest"""

    # This value is in amu.
    max_raman_intensity = 575

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos) - 6

    def testlengths(self, data, numvib) -> None:
        """Is the length of vibramans correct?"""
        assert len(data.vibramans) == numvib

    # The tolerance for this number has been increased, since ORCA
    # failed to make it inside +/-5, but it would be nice in the future
    # to determine is it's not too much work whether this is due to
    # algorithmic differences, or to differences in the input basis set
    # or coordinates. The first would be OK, but in the second case the
    # unit test jobs should be made more comparable. With cclib, we first
    # of all want to succeed in parsing, but would also like to remain
    # as comparable between programs as possible (for these tests).
    # Note also that this value is adjusted for Gaussian and DALTON - why?
    def testramanintens(self, data) -> None:
        """Is the maximum Raman intensity correct?"""
        assert abs(max(data.vibramans) - self.max_raman_intensity) < 8

        # We used to test this, but it seems to vary wildly between
        # programs... perhaps we could use it if we knew why...
        # self.assertInside(data.vibramans[1], 2.6872, 0.0001)

    def testvibdisps(self, data) -> None:
        """Is the length and value of vibdisps correct?"""
        assert hasattr(data, "vibdisps")
        assert len(data.vibdisps) == 54


class DALTONRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 745


class GaussianRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 1066


class OrcaRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 1045

    def testvibramans_exact(self, data, numvib) -> None:
        """Are Raman intensities parsed correctly with exact values?"""
        assert hasattr(data, "vibramans")
        assert len(data.vibramans) == numvib

        # Reference values from dvb_raman.out (ORCA 5.0)
        # Max intensity should be 1037.104624
        assert pytest.approx(data.vibramans.max(), abs=0.001) == 1037.104624

        # Check specific values for first few modes
        # Mode 0: 0.0, Mode 1: 2.498442, Mode 2: 0.0, Mode 3: 0.0, Mode 4: 7.649875
        expected_first_five = [0.0, 2.498442, 0.0, 0.0, 7.649875]
        for i in range(5):
            assert pytest.approx(data.vibramans[i], abs=0.001) == expected_first_five[i]

        # Raman intensities should be non-negative
        assert np.all(data.vibramans >= 0)

        # Type check
        assert isinstance(data.vibramans, np.ndarray)
        assert data.vibramans.dtype in (np.float64, np.float32)


class Orca6RamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    # This value has changed again in Orca 6 for some reason...
    max_raman_intensity = 1037


class QChemRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 588
