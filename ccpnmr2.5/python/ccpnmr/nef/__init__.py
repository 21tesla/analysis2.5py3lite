"""Star, Cif, and NEF I/O.

Files are converted to and from a tree of nested objects (essentially enhanced, ordered dictionaries)

The GenericStarParser.py will work for any valid Star file; see this file
for the precise reading and writing behaviour. Parsing comes in several modes, depending on how
strictly the STAR standard is enforced.

NmrStar and NEF files have a more restricted syntax, with all items and loops contained within
saveframes, and with saveframe and loop tags beginning with a prefix matching the saveframe and loop
category. The StarIo.py code is used for these files; it converts the output of the GenericStarParser
into a simpler object tree that makes use of these restrictions, strips mandatory tag prefixes,
and converts strings to numerical values with some limited heuristics.

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (https://www.ccpn.ac.uk) 2014 - 2022"
__credits__ = ("Ed Brooksbank, Joanna Fox, Victoria A Higman, Luca Mureddu, Eliza Płoskoń",
               "Timothy J Ragan, Brian O Smith, Gary S Thompson & Geerten W Vuister")
__licence__ = ("CCPN licence. See https://ccpn.ac.uk/software/licensing/")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Geerten Vuister $"
__dateModified__ = "$dateModified: 2022-02-07 17:13:54 +0000 (Mon, February 07, 2022) $"
__version__ = "$Revision: 3.1.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

NEF_ROOT_PATH = __path__[0]


