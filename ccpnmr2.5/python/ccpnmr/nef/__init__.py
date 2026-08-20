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
__dateModified__ = "$dateModified: 2017-07-07 16:33:01 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================



