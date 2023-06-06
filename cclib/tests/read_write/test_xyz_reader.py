# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for XYZ reader."""

import os
import unittest

import numpy as np

import cclib
from cclib.data.data_handler import get_data_path


class XYZReaderTest(unittest.TestCase):
    def test_attributes_one(self):
        """Is an XYZ file with a single geometry read into a ccData properly?"""
        fpath = get_data_path("XYZ/uracil.xyz")
        xyz = cclib.io.xyzreader.XYZ(fpath)
        data = xyz.parse()

        assert data.natom == 12

        atomnos = np.array([7, 6, 1, 8, 7, 6, 8, 6, 6, 1, 1, 1], dtype=int)
        np.testing.assert_equal(data.atomnos, atomnos)

        assert data.atomcoords.shape == (1, 12, 3)

    def test_attributes_two(self):
        """Is an XYZ file with a two geometries read into a ccData properly?"""
        fpath = get_data_path("XYZ/uracil_two.xyz")
        xyz = cclib.io.xyzreader.XYZ(fpath)
        data = xyz.parse()

        assert data.natom == 12

        atomnos = np.array([7, 6, 1, 8, 7, 6, 8, 6, 6, 1, 1, 1], dtype=int)
        np.testing.assert_equal(data.atomnos, atomnos)

        assert data.atomcoords.shape == (2, 12, 3)
        assert data.metadata["comments"] == ["uracil", "	Energy:     -97.0646597"]


if __name__ == "__main__":
    unittest.main()
