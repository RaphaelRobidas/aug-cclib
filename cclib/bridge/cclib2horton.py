# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in horton (http://theochem.github.io/horton)."""
# Support for horton 2 will likely be dropped when Python 2 support is dropped from cclib.

import numpy
from cclib.parser.data import ccData
from cclib.parser.utils import find_package

# First check horton version
_old_horton = False
_found_horton = find_package("horton")
_found_iodata = find_package("iodata")

# Detect whether horton 2 is present or not
if _found_horton:
    # Older horton packages do not have __version__, causing exceptions.
    try:
        from horton import __version__
    except:
        _old_horton = True
    else:
        if __version__[0] == "2":
            from horton.io.iodata import IOData

# Detect whether iodata (part of horton 3) is present or not
# Horton 3 is divided into smaller (sub)packages that each take different functionalities.
if _found_iodata:
    from iodata import IOData
    from iodata.orbitals import MolecularOrbitals


def check_horton():
    if _old_horton:
        raise ImportError("You must have at least version 2 of `horton` to use this function.")
    elif not _found_horton and not _found_iodata:
        raise ImportError("You must install `horton` to use this function.")
    if _found_iodata:
        return 3
    elif _found_horton:
        return 2


def makehorton(data):
    """ Create horton IOData object from ccData object """

    hortonver = check_horton()
    attributes = {}

    if hortonver == 2:
        renameAttrs = {"atomnos": "numbers", "mult": "ms2", "polarizability": "polar"}
        inputattrs = data.__dict__

        attributes = dict(
            (renameAttrs[oldKey], val)
            for (oldKey, val) in inputattrs.items()
            if (oldKey in renameAttrs)
        )

        # Rest of attributes need some manipulation in data structure.
        if hasattr(data, "atomcoords"):
            # cclib parses the whole history of coordinates in the list, horton keeps the last one.
            attributes["coordinates"] = data.atomcoords[-1]
        if hasattr(data, "mocoeffs"):
            attributes["orb_alpha"] = data.mocoeffs[0]
            if len(data.mocoeffs) == 2:
                attributes["orb_beta"] = data.mocoeffs[1]
        if hasattr(data, "coreelectrons"):
            # cclib stores num of electrons screened out by pseudopotential
            # horton stores num of electrons after applying pseudopotential
            attributes["pseudo_numbers"] = data.atomnos - numpy.asanyarray(data.coreelectrons)
        if hasattr(data, "atomcharges") and ("mulliken" in data.atomcharges):
            attributes["mulliken_charges"] = data.atomcharges["mulliken"]
        if hasattr(data, "atomcharges") and ("natural" in data.atomcharges):
            attributes["npa_charges"] = data.atomcharges["natural"]

    elif hortonver == 3:
        pass  # Populate with bridge for horton 3 in later PR

    return IOData(**attributes)  # Pass collected attributes into IOData constructor


def makecclib(iodat):
    """ Create cclib ccData object from horton IOData object """

    hortonver = check_horton()
    attributes = {}

    if hortonver == 2:
        # For a few attributes, a simple renaming suffices
        renameAttrs = {
            "numbers": "atomnos",
            "ms2": "mult",
            "polar": "polarizability",
        }
        inputattrs = iodat.__dict__

        attributes = dict(
            (renameAttrs[oldKey], val)
            for (oldKey, val) in inputattrs.items()
            if (oldKey in renameAttrs)
        )

        # Rest of attributes need some manipulation in data structure.
        if hasattr(iodat, "coordinates"):
            # cclib parses the whole history of coordinates in the list, horton keeps the last one.
            attributes["atomcoords"] = [iodat.coordinates]
        if hasattr(iodat, "orb_alpha"):
            attributes["mocoeffs"] = [iodat.orb_alpha]
        if hasattr(iodat, "orb_beta"):
            attributes["mocoeffs"].append(iodat.orb_beta)
        if hasattr(iodat, "pseudo_numbers"):
            # cclib stores num of electrons screened out by pseudopotential
            # horton stores num of electrons after applying pseudopotential
            attributes["coreelectrons"] = iodat.numbers - iodat.pseudo_numbers
        if hasattr(iodat, "mulliken_charges"):
            attributes["atomcharges"] = {"mulliken": iodat.mulliken_charges}
            if hasattr(iodat, "npa_charges"):
                attributes["atomcharges"]["natural"] = iodat.npa_charges
        elif hasattr(iodat, "npa_charges"):
            attributes["atomcharges"] = {"natural": iodat.npa_charges}

    elif hortonver == 3:
        # Horton 3 IOData class uses attr and does not have __dict__.
        # In horton 3, some attributes have a default value of None.
        # Therefore, second hasattr statement is needed for mo attribute.
        if hasattr(iodat, "atcoords"):
            # cclib parses the whole history of coordinates in the list, horton keeps the last one.
            attributes["atomcoords"] = [iodat.atcoords]
        if hasattr(iodat, "mo") and hasattr(iodat.mo, "norba"):
            # MO coefficient should be transposed to match the dimensions.
            attributes["mocoeffs"] = [iodat.mo.coeffs[: iodat.mo.norba].T]
            if iodat.mo.kind == "unrestricted":
                attributes["mocoeffs"].append(iodat.mo.coeffs[iodat.mo.norba :].T)
        if hasattr(iodat, "spinpol"):
            # IOData stores 2S, ccData stores 2S+1.
            attributes["mult"] = iodat.spinpol + 1
        if hasattr(iodat, "atnums"):
            # cclib stores num of electrons screened out by pseudopotential
            # horton stores num of electrons after applying pseudopotential
            attributes["atnums"] = numpy.asanyarray(iodat.atnums)
        if hasattr(iodat, "atcorenums"):
            attributes["coreelectrons"] = numpy.asanyarray(iodat.atnums) - numpy.asanyarray(
                iodat.atcorenums
            )
        if hasattr(iodat, "atcharges"):
            attributes["atomcharges"] = iodat.atcharges
    return ccData(attributes)
