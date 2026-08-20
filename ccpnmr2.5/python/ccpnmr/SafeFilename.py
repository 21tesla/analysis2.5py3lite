"""
Functions to append a number to the end of a filename if it already exists
"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2019"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: CCPN $"
__dateModified__ = "$dateModified: 2017-07-07 16:32:59 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-06-16 16:28:31 +0000 (Fri, June 16, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import itertools
import errno
import os
import sys
from contextlib import contextmanager


def _iter_incrementing_file_names(path):
    """
    Iterate incrementing file names. Start with path and add '(n)' before the
    extension, where n starts at 1 and increases.

    :param path: Some path
    :return: An iterator.
    """
    yield path
    prefix, ext = os.path.splitext(path)
    for i in itertools.count(start=1, step=1):
        yield prefix + '({0})'.format(i) + ext


@contextmanager
def safeOpen(path, mode):
    """
    Open path, but if it already exists, add '(n)' before the extension,
    where n is the first number found such that the file does not already
    exist.
    Returns an open file handle.

    Usage:  with safeOpen(path, [options]) as (fd, safeFileName):
                ...

        fd is the file descriptor, to be used as with open, e.g., fd.read()
        safeFileName is the new safe filename.

    :param path: filepath and filename.
    :param mode: open flags
    :return: Open file handle and new fileName
    """
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY

    if 'b' in mode and sys.platform.system() == 'Windows' and hasattr(os, 'O_BINARY'):
        flags |= os.O_BINARY

    # repeat over filenames with iterating number
    for filename in _iter_incrementing_file_names(path):
        try:
            file_handle = os.open(filename, flags)
        except OSError as e:
            if e.errno == errno.EEXIST:
                # force repeat of the loop if file exists (file not opened)
                pass
            else:
                raise
        else:
            # yield the file descriptor and new, safe filename
            with os.fdopen(file_handle, mode) as fd:
                yield fd, filename

            # ...and exit
            return

def getSafeFilename(path, mode='w'):
    """Get the first safe filename from the given path

    :param path: filepath and filename.
    :param mode: open flags
    :return: Open file handle and new fileName
    """
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY

    if 'b' in mode and sys.platform.system() == 'Windows' and hasattr(os, 'O_BINARY'):
        flags |= os.O_BINARY

    # repeat over filenames with iterating number
    for filename in _iter_incrementing_file_names(path):
        try:
            file_handle = os.open(filename, flags)
        except OSError as e:
            if e.errno == errno.EEXIST:
                # force repeat of the loop if file exists (file not opened)
                pass
            else:
                raise
        else:
            # return the new filename
            return filename
