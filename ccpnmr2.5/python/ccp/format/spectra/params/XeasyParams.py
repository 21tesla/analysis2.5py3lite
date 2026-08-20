
"""
======================COPYRIGHT/LICENSE START==========================

XeasyParams.py: Part of the CcpNmr Analysis program

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

from ccp.format.spectra.params.ExternalParams import ExternalParams

XeasyParamDict = {
  'ndim': 'Number of dimensions', 
  'nbits': '16 or 8 bit file type', 
  'sf': 'Spectrometer frequency in w',
  'sw': 'Spectral sweep width in w',
  'maxppm': 'Maximum chemical shift in w',
  'npts': 'Size of spectrum in w',
  'block': 'Submatrix size in w',
  'order': 'Permutation for w',
  'fold': 'Folding in w',  # not used
  'type': 'Type of spectrum',  # not used
  'nuc': 'Identifier for dimension w',
}

class XeasyParams(ExternalParams):

  format = 'XEASY'

  def __init__(self, file, **kw):

    self.paramFile = file
    ExternalParams.__init__(self, **kw)

  # ExternalParams requires this to be defined
  def parseFile(self):

    fp = open(self.paramFile, 'rU')

    firstLine = 'Version ....................... '
    line = fp.readline().strip()
    if line[:32] != firstLine:
      raise IOError('The file %s does not look like an XEASY param file because the first line does not start "%s"' % (self.paramFile, firstLine))

    if line[-1] != '1':
      print 'Warning: this XEASY param file is version != 1 so might not be interpreted correctly'

    lines = fp.readlines()
    fp.close()

    dd = {}
    for line in lines:
      key = line[:32].replace('.', '').strip()
      value = line[32:].strip()
      dd[key] = value

    ndim = self.ndim = int(dd[XeasyParamDict['ndim']])

    nbits = int(dd[XeasyParamDict['nbits']])
    self.nbytes = nbits / 8

    self.dataFile = self.paramFile[:-5] + str(nbits)

    self.initDims()

    for i in range(ndim):
      ss = str(i+1)
      j = int(dd[XeasyParamDict['order']+ss]) - 1
      self.npts[j] = int(dd[XeasyParamDict['npts']+ss])
      self.block[j] = int(dd[XeasyParamDict['block']+ss])
      self.sf[j] = float(dd[XeasyParamDict['sf']+ss])
      self.sw[j] = float(dd[XeasyParamDict['sw']+ss])
      self.sw[j] *= self.sf[j]  # convert from ppm to Hz
      self.refppm[j] = float(dd[XeasyParamDict['maxppm']+ss])
      self.refpt[j] = 1.0
      nuc = dd[XeasyParamDict['nuc']+ss]
      self.nuc[j] = self.standardNucleusName(nuc)

if __name__ == '__main__':

  import sys
  if len(sys.argv) != 2:
    print 'Error: correct syntax: <script> <XEASY_file>'
    sys.exit(1)

  xeasy_file = sys.argv[1]
  params = XeasyParams(xeasy_file)
