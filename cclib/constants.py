import cclib
import importlib

parser_names = [
    "ADF", "DALTON", "FChk", "GAMESS", "GAMESSUK", "Gaussian", "Jaguar",
     "Molpro", "Molcas", "MOPAC", "NWChem", "ORCA", "Psi4", "QChem",
    "Turbomole",
]
all_parsers = {name: getattr(cclib.parser, name) for name in parser_names}

# Not used currently, but keeping in a list to keep track of which parsers
# are in the legacy bin.
legacy_parser_names = ["Psi3"]


module_names = [
    "SP", "SPun", "GeoOpt", "Basis", "Core",    # Basic calculations.
    "MP", "CC", "CI", "TD", "TDun",             # Post-SCF calculations.
    "BOMD", "NMR", "Polar", "Scan", "vib"       # Other property calculations.
]

