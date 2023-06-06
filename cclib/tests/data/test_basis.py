# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles related to basis sets"""

import os
import unittest

from cclib.constants import *
from cclib.tests.data.skip import skipForParser
from cclib.tests.data.suite import DataSuite


class GenericBasisTest(DataSuite):
    """Generic basis set unittest"""

    # The number of contraction per atom, by atom number.
    contractions = {1: 1, 6: 3}

    # Number of components in each contraction by subshell type,
    # so that we can infer nbasis from gbasis. Note how we assume
    # the basis set is not is spherical representation.
    names = ["S", "P", "D", "F", "G"]
    multiple = {"S": 1, "P": 3, "D": 6, "F": 10, "G": 15}
    multiple_spher = {"S": 1, "P": 3, "D": 5, "F": 7, "G": 9}
    spherical = False

    # These are the expected exponents and coefficients for the first
    # Gaussians in particular shells for hydrogen and carbon atoms.
    gbasis_H_1s_func0 = [3.42525, 0.15433]
    gbasis_C_2s_func0 = [2.9412, -0.1000]
    gbasis_C_2p_func0 = [2.9412, 0.1559]

    def test_basis(self):
        for self.data, self.logfile in zip(self.all_data, self.all_logfile):
            for feature in [
                "_gbasis",
                "_names",
                "_sizeofbasis",
                "_contractions",
                "_primitives",
                "_coeffs",
            ]:
                with self.subTest(logfile=self.logfile, feature=feature):
                    getattr(self, feature)()

    def _gbasis(self):
        """Is gbasis the right length?"""
        assert self.data.natom == len(self.data.gbasis)

    def _names(self):
        """Are the name of basis set functions acceptable?"""
        for atom in self.data.gbasis:
            for fns in atom:
                assert fns[0] in self.names, f"{fns[0]} not one of S or P"

    def _sizeofbasis(self):
        """Is the basis set the correct size?"""

        total = 0
        multiple = self.multiple_spher if self.spherical else self.multiple
        for atom in self.data.gbasis:
            for ftype, contraction in atom:
                total += multiple[ftype]

        assert self.data.nbasis == total

    def _contractions(self):
        """Are the number of contractions on all atoms correct?"""
        for iatom, atom in enumerate(self.data.gbasis):
            atomno = self.data.atomnos[iatom]
            assert len(atom) == self.contractions[atomno]

    def _primitives(self):
        """Are all primitives 2-tuples?"""
        for atom in self.data.gbasis:
            for ftype, contraction in atom:
                for primitive in contraction:
                    assert len(primitive) == 2

    def _coeffs(self):
        """Are the atomic basis set exponents and coefficients correct?"""

        for iatom, atom in enumerate(self.data.gbasis):
            if self.data.atomnos[iatom] == 1:
                coeffs = atom[0][1]
                assert round(abs(coeffs[0][0] - self.gbasis_H_1s_func0[0]), 4) == 0
                assert round(abs(coeffs[0][1] - self.gbasis_H_1s_func0[1]), 4) == 0
            else:
                s_coeffs = atom[1][1]
                p_coeffs = atom[2][1]
                assert round(abs(s_coeffs[0][0] - self.gbasis_C_2s_func0[0]), 4) == 0
                assert round(abs(p_coeffs[0][0] - self.gbasis_C_2p_func0[0]), 4) == 0
                assert round(abs(s_coeffs[0][1] - self.gbasis_C_2s_func0[1]), 4) == 0
                assert round(abs(p_coeffs[0][1] - self.gbasis_C_2p_func0[1]), 4) == 0


class JaguarBasisTest(GenericBasisTest):
    """Customized basis set unittest"""

    # For some reason, Jaguar seems to use slightly different coefficients for
    # contractions in the STO-3G basis set. Or perhaps we don't understand something.
    gbasis_H_1s_func0 = [3.42525, 0.24050]
    gbasis_C_2s_func0 = [2.941249, -0.29565]
    gbasis_C_2p_func0 = [2.941249, 0.22135]

    parsers = ["Jaguar"]
    modules = ["Basis"]


class GenericBigBasisTest(GenericBasisTest):
    """Generic big basis set unittest"""

    contractions = {6: 20}

    '''
    @unittest.skip('Write up a new test, and/or revise the one inherited.')
    def testcoeffs(self):
        """Are the basis set coefficients correct?"""
        assert True

    @unittest.skip('# of contractions is 20 for VQZ, but 29 for CVQZ; unify files first.')
    def testcontractions(self):
        """"""
        assert True
    '''


class DALTONBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""

    spherical = True


class GaussianBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""

    spherical = True

    parsers = ["Gaussian"]
    modules = ["Basis"]


class JaguarBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""

    spherical = True

    # Jaguar only goes up to F functions.
    names = ["S", "P", "D", "F"]

    parsers = ["Jaguar"]
    modules = ["Basis"]


class MolcasBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""

    spherical = True
    parsers = ["Molcas"]
    modules = ["Basis"]


class MolproBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""

    spherical = True
    parsers = ["Molpro"]
    modules = ["Basis"]


class Psi4BigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""

    spherical = True

    parsers = ["Psi4"]
    modules = ["Basis"]


class QChemBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""

    spherical = True

    parsers = ["QChem"]
    modules = ["Basis"]
