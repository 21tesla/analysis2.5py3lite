
"""
======================COPYRIGHT/LICENSE START==========================

NmrPipeParams.py: Part of the CcpNmr Analysis program

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
import os

from memops.universal.Util import isBigEndian

from memops.general.Implementation import ApiError

from ccp.format.spectra.params.ExternalParams import ExternalParams

ndim_index = 9
npts_index = (99, 219, 15, 32)
complex_index = (55, 56, 51, 54)
order_index = (24, 25, 26, 27)
sw_index = (229, 100, 11, 29)
sf_index = (218, 119, 10, 28)
origin_index = (249, 101, 12, 30)
nuc_index = (18, 16, 20, 22)

value_index = 199

head = 4*512

class NmrPipeParams(ExternalParams):

  format = 'NmrPipe'

  def __init__(self, file, **kw):

    self.dataFile = file
    self.head = head
    ExternalParams.__init__(self, **kw)

  # ExternalParams requires this to be defined
  def parseFile(self):

    (s, x, self.big_endian, self.swap) = getHeader(self.dataFile)

    ndim = self.ndim = int(x[ndim_index])

    if ((ndim < 1) or (ndim > 4)):
      raise ApiError('ndim = %d, should be between 1 and 4' % ndim)

    self.initDims()

    for i in range(ndim):
      j = int(x[order_index[i]]) - 1
      #c = int(x[complex_index[j]])
      c = int(x[complex_index[i]])
      if (c == 0):
        raise ApiError('data is complex in dim %d, can only cope with real data so far' % (i+1))
      #self.npts[i] = int(x[npts_index[j]])
      self.npts[i] = int(x[npts_index[i]])
      if (i == 0):
        self.block[i] = self.npts[i]
      else:
        self.block[i] = 1
      self.sw[i] = x[sw_index[j]]
      if (self.sw[i] == 0): self.sw[i] = 1000 # TBD: TEMP hack
      self.sf[i] = x[sf_index[j]]
      o = x[origin_index[j]]
      self.refppm[i] = (self.sw[i] + o) / self.sf[i]
      self.refpt[i] = 0
      n = 4 * nuc_index[j]
      nuc = s[n:n+4].strip()
      # get rid of null termination
      m = nuc.find(chr(0))
      if m >= 0:
        nuc = nuc[:m]
      if nuc == 'ID' or nuc == 'TAU': # TBD: do not know if this is necessary or sufficient
        self.nuc[i] = None
      else:
        self.nuc[i] = standardNucleusNameFromPipe(nuc)

def getHeader(fileName):

    try:
      fp = open(fileName, 'rb')
    except IOError, e:
      raise ApiError(str(e))

    s = fp.read(head)
    fp.close()

    if (len(s) < head):
      raise ApiError('file shorter than expected length (%d bytes) of header (never mind data)' % head)

    x = array.array('f')
    x.fromstring(s)

    if (x[0] != 0):
      raise ApiError('first word of header = %s, should be 0' % x[0])

    byte_order = [ 0x40, 0x16, 0x14, 0x7b ]
    t = [ ord(c) for c in s[8:12] ]
    if (t == byte_order):
      big_endian = True
    else:
      t.reverse()
      if (t == byte_order):
        big_endian = False
      else:
        raise ApiError('bytes 8 through 11 should be [ 0x40, 0x16, 0x14, 0x7b ] or reverse')

    big_endian = big_endian
    swap = not (big_endian == isBigEndian())

    if (swap):
      x.byteswap()

    return (s, x, big_endian, swap)

def getDataFileName(template, z, a = None):
  """Get NMRPipe fileName given template and z and (optionally) a.
     Note that z and a start counting at 0, not 1.
  """

  if a is None:
    name = template % (z+1)
  else:
    name = template % (a+1, z+1)

  return name

def getSampledValue(fileName):
  """Get sampled value from NMRPipe header for file specified by template in directory
  """

  (s, x, big_endian, swap) = getHeader(fileName)
  value = float(x[value_index])

  return value

def getSampledValues(directory, template):
  """Get sampled values from NMRPipe 2D headers for files specified by template in directory
  """

  values = []

  z = 0
  while 1:
    try:
      file = getDataFileName(template, z)
    except:
      raise ApiError('template not valid')
    fullfile = os.path.join(directory, file)
    if os.path.exists(fullfile) and os.path.isfile(fullfile):
      value = getSampledValue(fullfile)
      values.append(value)
    else:
      break
    z = z + 1

  return values

def standardNucleusNameFromPipe(name):
  """ Given Pipe nucleus name, return standard nucleus name
      Replaces standardNucleusName in format.general.Util 
  """

  full_nucs = ('1H',    '13C',    '15N',    '31P')
  nucs      = {'H':'1H','C':'13C','N':'15N','P':'31P'}

  for f in full_nucs:
    if name.count( f ) > 0:
      return f

  if name and name[0].upper() in nucs.keys():
    return nucs[ name[0].upper() ]

  else:
    # default name (as in standardNucleusName)
    return '1H'

if (__name__ == '__main__'):

  import sys
  if (len(sys.argv) != 2):
    print 'Error: correct syntax: <script> <NMRPipe_file>'
    sys.exit(1)

  nmrpipe_file = sys.argv[1]
  params = NmrPipeParams(nmrpipe_file)
