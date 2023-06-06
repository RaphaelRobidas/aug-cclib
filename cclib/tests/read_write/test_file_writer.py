# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for general file writer."""

import inspect
import os
import unittest

import cclib
import pytest

from cclib.data.data_handler import read_data


class FileWriterTest(unittest.TestCase):
    def test_init(self):
        """Does the class initialize properly?"""

        # You cannot instantiate a class with abstract methods.
        data = read_data("ADF/basicADF2007.01/dvb_gopt.adfout")
        with pytest.raises(TypeError):
            cclib.io.filewriter.Writer(data)


if __name__ == "__main__":
    unittest.main()
