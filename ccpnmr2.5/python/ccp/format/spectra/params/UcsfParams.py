
"""
======================COPYRIGHT/LICENSE START==========================

UcsfParams.py: Part of the CcpNmr Analysis program

Copyright (C) 2005 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
import array

from memops.universal.Util import isBigEndian

from memops.general.Implementation import ApiError

from ccp.format.spectra.params.ExternalParams import ExternalParams

ucsf_file_header = 180
ucsf_dim_header = 128

class UcsfParams(ExternalParams):

  format = 'UCSF'

  def __init__(self, file, **kw):

    self.dataFile = file
    ExternalParams.__init__(self, **kw)

  # ExternalParams requires this to be defined
  def parseFile(self):

    try:
      fp = open(self.dataFile, 'rb')
    except IOError, e:
      raise ApiError(str(e))

    s = fp.read(ucsf_file_header)
    if (len(s) < ucsf_file_header):
      raise ApiError('file shorter than expected length (%d bytes) of first part of header (never mind data)' % ucsf_file_header)

    if (s[:8] != 'UCSF NMR'):
      raise ApiError('first eight bytes of header = "%s", should be "UCSF NMR"' % s[:8])

    ndim = self.ndim = ord(s[10])

    self.head = ucsf_file_header + ndim*ucsf_dim_header
    self.big_endian = True # UCSF files always big endian
    self.swap = not isBigEndian()

    s = fp.read(ndim*ucsf_dim_header)
    if (len(s) < ndim*ucsf_dim_header):
      raise ApiError('file shorter than expected length (%d bytes) of header (never mind data)' % (ucsf_file_header+ndim*ucsf_dim_header))

    fp.close()

    x = array.array('i')
    y = array.array('f')

    x.fromstring(s)
    y.fromstring(s)

    if (self.swap):
      x.byteswap()
      y.byteswap()

    self.initDims()

    for i in range(ndim):
      j = ndim - i - 1 # UCSF does dims backwards
      base = (j * ucsf_dim_header) / 4
      self.npts[i] = x[base+2]
      self.block[i] = x[base+4]
      self.sf[i] = y[base+5]
      self.sw[i] = y[base+6]
      self.refppm[i] = y[base+7]
      # 27 Oct 2006: added 1.0, which seems to be the correct thing to do
      self.refpt[i] = 1.0 + 0.5 * float(self.npts[i])
      nuc = s[4*base:4*base+6]
      n = nuc.find(chr(0))
      if (n >= 0):
        nuc = nuc[:n]
      self.nuc[i] = self.standardNucleusName(nuc)

if (__name__ == '__main__'):

  import sys
  if (len(sys.argv) != 2):
    print 'Error: correct syntax: <script> <UCSF_file>'
    sys.exit(1)

  ucsf_file = sys.argv[1]
  params = UcsfParams(ucsf_file)
