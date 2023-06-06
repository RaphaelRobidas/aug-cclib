import os
import sys
import logging
import cclib
from cclib.constants import *

from typing import List, Optional, Type, Union

__filedir__ = os.path.dirname(__file__)

def get_data_path(short_path):
    # Might not be compatible with Windows if short_path contains forward slashes
    # Easy fix: split and reform a correct (OS-dependent) path
    return os.path.join(__filedir__, short_path)

def get_data(parser: str, subdir, files, stream=None, loglevel: int = logging.ERROR, datatype: Optional[Type[cclib.parser.data.ccData]] = None):
    """Returns a parsed logfile.

    Inputs:
        parser   - a logfile parser class (subclass of LogFile)
        subdir   - subdirectory containing data files (program version)
        files    - data filename(s)
        stream   - where to log to (sys.stdout by default)
        loglevel - what level to log at
        datatype - ccData or child class

    Outputs:
        data - the resulting data object
        logfile - the parser object used for parsing
    """

    if not isinstance(parser, str):
        print("Invalid call")
        import traceback
        traceback.print_stack()

    inputs = [get_data_path(f"{parser}/{subdir}/{fn}") for fn in files]

    # We should be able to pass a list of length one here, but for some reason
    # this does not work with some parsers and we get errors.
    if len(inputs) == 1:
        inputs = inputs[0]

    parser_obj = all_parsers[parser]

    stream = stream or sys.stdout
    logfile = parser_obj(inputs, logstream=stream, loglevel=loglevel,
                     datatype=datatype or cclib.parser.data.ccData)

    data = logfile.parse()
    return data, logfile

def read_data(short_path):
    return cclib.io.ccread(get_data_path(short_path))
