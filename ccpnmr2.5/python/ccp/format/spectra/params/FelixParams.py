
"""
======================COPYRIGHT/LICENSE START==========================

FelixParams.py: Part of the CcpNmr Analysis program

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

from memops.general.Implementation import ApiError

from ccp.format.spectra.params.ExternalParams import ExternalParams

class FelixParams(ExternalParams):

  format = 'Felix'

  def __init__(self, file, **kw):

    self.dataFile = file
    self.head = 4*4096
    ExternalParams.__init__(self, **kw)

  # ExternalParams requires this to be defined
  def parseFile(self):

    try:
      fp = open(self.dataFile, 'rb')
    except IOError, e:
      raise ApiError(str(e))

    s = fp.read(self.head)
    if (len(s) < self.head):
      raise ApiError('file shorter than expected length (%d bytes)of header (never mind data)' % self.head)

    fp.close()

    x = array.array('i')
    y = array.array('f')

    x.fromstring(s)
    y.fromstring(s)

    matrix_type = x[1]
    if (matrix_type != 1):
      self.swap = True
      x.byteswap()
      y.byteswap()
      if (x[1] != 1):
        raise ApiError('normal byte ordered matrix type (word 1) = %d and swapped byte ordered matrix type = %d, one of these should be 1' % (matrix_type, x[1]))

    ndim = self.ndim = x[0]
    self.initDims()

    for i in range(ndim):
      self.npts[i] = x[20+1*ndim+i]
      self.block[i] = x[20+4*ndim+i]
      self.sf[i] = y[20+6*ndim+i]
      self.sw[i] = y[20+7*ndim+i]
      self.refpt[i] = y[20+8*ndim+i]
      self.refppm[i] = y[20+9*ndim+i] / self.sf[i]
      nuc = ''
      j = 0
      while ((j < 8) and x[220+8*i+j]):
        nuc = nuc + chr(x[220+8*i+j])
        j = j + 1
      self.nuc[i] = self.standardNucleusName(nuc)

if (__name__ == '__main__'):

  import sys
  if (len(sys.argv) != 2):
    print 'Error: correct syntax: <script> <matrix_file>'
    sys.exit(1)

  matrix_file = sys.argv[1]
  params = FelixParams(matrix_file)
