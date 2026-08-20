
"""
======================COPYRIGHT/LICENSE START==========================

VarianParams.py

Copyright (C) 2009 Wayne Boucher (University of Cambridge)

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
import os

from memops.general.Implementation import ApiError

from memops.universal.Io import joinPath
from ccp.format.varian.varianFile import parseProcparFile, readDataFileHeader, VNMR_FILE_HEADER_SIZE, VNMR_BLOCK_HEADER_SIZE

from ccp.format.spectra.params.ExternalParams import ExternalParams

class VarianParams(ExternalParams):

  format = 'Varian'

  def __init__(self, procpar_file, data_file=None, **kw):

    if not data_file:
      data_file = joinPath(os.path.dirname(procpar_file), 'datdir', 'phasefile')

    self.procpar_file = procpar_file
    self.data_file = data_file

    ExternalParams.__init__(self, **kw)

  # ExternalParams requires this to be defined
  def parseFile(self):

    try:
      procparParams = parseProcparFile(self.procpar_file)
      dataFileParams = readDataFileHeader(self.data_file)
    except IOError, e:
      raise ApiError(str(e))

    ccpnParams = procparParams['ccpnParams']
    self.ndim = ccpnParams['ndim']
    self.initDims()

    for key in ccpnParams.keys():
      if key != 'ndim': # already done (could do again...)
        value = ccpnParams[key]
        if type(value) == type([]) and (dataFileParams['firstBlock']==0) and self.ndim == 2:
          # in this case are using transposed data so everything is backwards
          # TBD: not sure what to do in 3D, etc.
          value.reverse()
        setattr(self, key, value)

    blockSize = dataFileParams['np'] * dataFileParams['ntraces']
    blockHead = dataFileParams['nbheaders'] * VNMR_BLOCK_HEADER_SIZE
    self.head = VNMR_FILE_HEADER_SIZE + dataFileParams['firstBlock']*(4*blockSize+blockHead)
    self.blockHead = blockHead
    if self.ndim ==  1:
      self.block = [blockSize]
    else:
      self.block = self.ndim*[1]
      #self.block[0] = dataFileParams['ntraces']
      #self.block[1] = dataFileParams['np']
      self.block[0] = self.npts[0]
      self.block[1] = blockSize // self.npts[0]  # TBD: is it always a multiple??
    self.swap = dataFileParams['swapped']
    self.dataFile = self.data_file

if (__name__ == '__main__'):

  import sys
  if (len(sys.argv) != 2):
    print 'Error: required argument: <procpar_file>'
    sys.exit(1)

  procs_file = sys.argv[1]
  params = VarianParams(procs_file)
