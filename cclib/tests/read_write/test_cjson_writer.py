# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for the CJSON writer."""

import json
import os
import unittest
from math import sqrt

import cclib

from cclib.data.data_handler import read_data


class CJSONWriterTest(unittest.TestCase):
    """Unit tests for the CJSON writer."""

    def test_init(self):
        """Does the class initialize correctly?"""
        data = read_data("ADF/basicADF2007.01/dvb_gopt.adfout")
        cjson = cclib.io.cjsonwriter.CJSON(data)

        # The object should keep the ccData instance passed to its constructor.
        assert cjson.ccdata == data

    def test_cjson_generation(self):
        """Does the CJSON format get generated properly?"""
        data = read_data("ADF/basicADF2007.01/NH3.adfout")

        cjson = cclib.io.cjsonwriter.CJSON(data).generate_repr()

        # The data available in the cjson and ccdata objects should be equal.
        json_data = json.loads(cjson)
        number_of_atoms = json_data["properties"]["number of atoms"]
        assert number_of_atoms == data.natom

        dipole_moment = json_data["properties"]["total dipole moment"]
        assert round(abs(dipole_moment - sqrt(sum(data.moments[1] ** 2))), 7) == 0

        # Ensure the bond connectivity index starts from 0
        bonds = json_data.get("bonds", None)
        assert bonds is not None
        indices = bonds["connections"]["index"]
        assert min(indices) == 0
        assert max(indices) < number_of_atoms

    def test_zero_dipole_moment(self):
        """Does the CJSON writer handle zero dipole moment correctly?"""
        data = read_data("GAMESS/basicGAMESS-US2017/C_bigbasis.out")

        cjson = cclib.io.cjsonwriter.CJSON(data).generate_repr()

        json_data = json.loads(cjson)
        assert round(abs(json_data["properties"]["total dipole moment"]), 7) == 0

    def test_missing_dipole_moment(self):
        """Does the CJSON writer handle missing properties correctly?"""
        data = read_data("GAMESS/basicGAMESS-US2017/C_bigbasis.out")
        del data.moments

        cjson = cclib.io.cjsonwriter.CJSON(data).generate_repr()

        json_data = json.loads(cjson)
        assert "total dipole moment" not in json_data["properties"]


if __name__ == "__main__":
    unittest.main()