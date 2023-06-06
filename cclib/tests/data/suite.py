import logging
import sys
import os
import unittest
from cclib.data.data_handler import get_test_data, get_data

from typing import List, Optional, Type, Union

def ccdata_getattribute_with_coverage(self, attr):
    """A bookkeeping version of __getattribute__ for ccData objects."""
    if attr != '_attrlist' and attr in self._attrlist:
        if not hasattr(self, 'coverage'):
            self.coverage = {}
        self.coverage[attr] = self.coverage.get(attr, 0) + 1
    return object.__getattribute__(self, attr)


class DataSuite(unittest.TestCase):
    """Suite containing data (logfile) tests in cclib.

    This is supposed to represent a single run of the entire data test suite in cclib or
    a subset of it. The main functions are to load data, run test cases in the data/
    subdirectory, and do some basic bookkeeping.
    """

    @classmethod
    def setUpClass(cls):

        # Load the test data and filter with parsers and modules.
        cls.testdata = get_test_data()
        cls.testdata = [td for td in cls.testdata if td['parser'] in cls.parsers]
        cls.testdata = [td for td in cls.testdata if td['module'] in cls.modules]

        # We want to gather the unit tests and results in several lists/dicts,
        # in order to easily generate summaries at the end.
        cls.errors = []
        cls.failures = []
        cls.alltests = []

        cls.all_data = []
        cls.all_logfile = []

        cls.perpackage = {p: [0, 0, 0, 0] for p in cls.parsers}

        for td in cls.testdata:
            description = ''
            data, logfile = get_data(
                td['parser'], td['subdir'], td['files']
            )
            cls.all_data.append(data)
            cls.all_logfile.append(os.path.join(td['subdir'], td['files'][0]))

    def summary(self) -> None: # DOOD?
        """Prints a summary of the suite after it has been run."""

        if self.errors:
            print("\n********* SUMMARY OF ERRORS *********\n", file=self.stream)
            print("\n".join(self.errors), file=self.stream)

        if self.failures:
            print("\n********* SUMMARY OF FAILURES *********\n", file=self.stream)
            print("\n".join(self.failures), file=self.stream)

        print("\n********* SUMMARY PER PACKAGE ****************", file=self.stream)
        names = sorted(self.perpackage.keys())
        total = [0, 0, 0, 0]
        print(" "*14, "\t".join(["Total", "Passed", "Failed", "Errors", "Skipped"]), file=self.stream)

        fmt = "%3d\t%3d\t%3d\t%3d\t%3d"
        for name in names:
            l = self.perpackage[name]
            args = (l[0], l[0]-l[1]-l[2]-l[3], l[2], l[1], l[3])
            print(name.ljust(15), fmt % args, file=self.stream)
            for i in range(4):
                total[i] += l[i]

        print("\n********* SUMMARY OF EVERYTHING **************", file=self.stream)
        print(
            f"TOTAL: {int(total[0])}\tPASSED: {int(total[0] - (total[1] + total[2] + total[3]))}\tFAILED: {int(total[2])}\tERRORS: {int(total[1])}\tSKIPPED: {int(total[3])}",
            file=self.stream,
        )

    def visualtests(self, stream=sys.stdout) -> None:
        """These are not formal tests -- but they should be eyeballed."""

        parsers_to_test = {
            'ADF2013.01' : getdatafile('ADF', "basicADF2013.01", ["dvb_gopt.adfout"])[0],
            'DALTON2015' : getdatafile('DALTON', "basicDALTON-2015", ["dvb_gopt_ks.out"])[0],
            'Firefly8.0' : getdatafile('GAMESS', "basicFirefly8.0", ["dvb_gopt_a.out"])[0],
            'Gaussian16' : getdatafile('Gaussian', "basicGaussian16", ["dvb_gopt.out"])[0],
            'GAMESS-US2018' : getdatafile('GAMESS', "basicGAMESS-US2018", ["dvb_gopt_a.out"])[0],
            'Jaguar8.0' : getdatafile('Jaguar', "basicJaguar8.3", ["dvb_gopt_ks.out"])[0],
            'Molpro2012' : getdatafile('Molpro', "basicMolpro2012", ["dvb_gopt.out", "dvb_gopt.log"])[0],
            # Note that it doesn't make sense to put MOPAC here, as it
            # is a semiempirical-only program.
            'NWChem6.5' : getdatafile('NWChem', "basicNWChem6.5", ["dvb_gopt_ks.out"])[0],
            'ORCA4.2' : getdatafile('ORCA', "basicORCA4.2", ["dvb_gopt.out"])[0],
            'Psi4-1.3.1' : getdatafile('Psi4', "basicPsi4-1.3.1", ["dvb_gopt_rks.out"])[0],
            'QChem5.4' : getdatafile('QChem', "basicQChem5.4", ["dvb_gopt.out"])[0],
        }
        parser_names = sorted(parsers_to_test.keys())
        output = [parsers_to_test[pn] for pn in parser_names]

        print("\n*** Visual tests ***", file=self.stream)
        print("MO energies of optimised dvb", file=self.stream)
        print(
            "      ", "".join([f"{pn:12s}" for pn in parser_names]), file=self.stream
        )
        print(
            "HOMO",
            "   ".join([f"{out.moenergies[0][out.homos[0]]:+9.4f}" for out in output]),
            file=self.stream,
        )
        print(
            "LUMO",
            "   ".join(
                [f"{out.moenergies[0][out.homos[0] + 1]:+9.4f}" for out in output]
            ),
            file=self.stream,
        )
        print(
            "H-L ",
            "   ".join(
                [
                    f"{out.moenergies[0][out.homos[0] + 1] - out.moenergies[0][out.homos[0]]:9.4f}"
                    for out in output
                ]
            ),
            file=self.stream,
        )


def test_all(parsers, modules, terse: bool, silent: bool, loglevel: int, summary: bool, visual_tests: bool) -> None:
    parsers = parsers or all_parsers
    modules = modules or all_modules
    data_suite = DataSuite(parsers, modules, terse=terse, silent=silent, loglevel=loglevel)
    errors_or_failures = data_suite.testall()
    if summary and not silent:
        data_suite.summary()
    if visual_tests and not silent:
        data_suite.visualtests()
    if errors_or_failures:
        sys.exit(1)


