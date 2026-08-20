"""
======================COPYRIGHT/LICENSE START==========================

ContourStore.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

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

from memops.universal.Io import splitPath
from memops.universal.Util import isBigEndian
from ccpnmr.analysis.core.Util import getAnalysisSpectrum

# TODO: view --> spectrum, contourFile, xdim, ydim
# indeed could just use block_file and get xdim, ydim, npoints, ndim from there
# (and also contourFile)
# saveViewContours for use when you want to contour existing view
def saveViewContours(view, fileName, levels, firstInt, lastInt, isBigEndian = None):

  from ccpnmr.c import ContourFile
  from memops.c import StoreHandler

  if not hasattr(view, 'contourFile'):
    return

  spectrum = view.analysisSpectrum.dataSource
  if not hasattr(spectrum, 'block_file'):
    return

  if not spectrum.block_file:
    return

  try:
    saveContours(spectrum, view.contourFile, fileName, levels, firstInt, lastInt, isBigEndian)
  except ContourFile.error, e:
    raise Exception(str(e))
  except StoreHandler.error, e:
    raise Exception(str(e))

# saveSpectrumContours for use when you want to contour specified region
def saveSpectrumContours(spectrum, fileName, xdim, ydim, levels,
                         firstInt, lastInt, isBigEndian = None, mem_cache = None):

  from ccpnmr.analysis.core.BlockUtil import getBlockFile
  from ccpnmr.c import ContourFile
  from memops.c import StoreHandler
  from memops.c.MemCache import MemCache
  from memops.general import Implementation

  if not mem_cache:
    cache_size = 64 * 1024 * 1024
    mem_cache = MemCache(cache_size)

  if hasattr(spectrum, 'block_file') and spectrum.block_file:
    block_file = spectrum.block_file
  else:
    block_file = getBlockFile(spectrum, mem_cache)
  if not block_file:
    raise Exception('could not get block file')

  try:
    # -1 because of dim convention
    contourFile = ContourFile.ContourFile(xdim-1, ydim-1, block_file, mem_cache)

    saveContours(spectrum, contourFile, fileName, levels, firstInt, lastInt, isBigEndian)
  except ContourFile.error, e:
    raise Exception(str(e))
  except StoreHandler.error, e:
    raise Exception(str(e))

# saveContours internal, called by both saveViewContours and saveSpectrumContours
def saveContours(spectrum, contourFile, fileName, levels, firstInt, lastInt, isBigEndian = None):

  from memops.c.StoreHandler import StoreHandler
  import ccpnmr.c.ContourStyle as ContourStyle
  import ccpnmr.c.ContourLevels as ContourLevels

  # set things up

  if isBigEndian is None or isBigEndian == isBigEndian():
    swap = 0
  else:
    swap = 1

  handler = StoreHandler(fileName, swap)

  # below irrelevant but mandatory argument
  pos_colors = neg_colors = [(0.0, 0.0, 0.0)]
  contourStyle = ContourStyle.ContourStyle(pos_colors, neg_colors, 0, 0)

  contourLevels = ContourLevels.ContourLevels(levels)

  contourFile.draw(handler, firstInt, lastInt, contourLevels, contourStyle)

# returns dict with keys: ndim, xdim, ydim, npoints, blockSize, levels
def getStoredContourHeader(fileName):

  # this needs to be kept consistent with memops/global/store_file.c and store_handler.c
  n = 4
  fp = open(fileName, 'rb')
  s = fp.read(6*n)
  x = array.array('i')
  x.fromstring(s)

  magic = x[0]
  if magic == 1789:
    swap = False
  else:
    swap = True
    x.byteswap()
    magic = x[0]
    if magic != 1789:
      raise Exception('magic number invalid')

  header = {}
  header['ndim'] = ndim = x[3]
  header['xdim'] = x[4] + 1  # + 1 because of dim convention
  header['ydim'] = x[5] + 1

  s = fp.read((4*ndim+1)*n)
  x = array.array('i')
  x.fromstring(s)
  if swap:
    x.byteswap()

  # x is an array object, so convert into ordinary lists
  header['npoints'] = list(x[:ndim])
  header['blockSize'] = list(x[3*ndim:4*ndim])

  nlevels = x[4*ndim]

  s = fp.read(nlevels*n)
  x = array.array('f')
  x.fromstring(s)
  if swap:
    x.byteswap()

  header['levels'] = list(x)

  fp.close()

  return header

# Note: in v1 fileName was absolute, now it is relative
def createStoredContour(spectrum, fileName, xdim, ydim):

  analysisSpectrum = spectrum.analysisSpectrum
  if not analysisSpectrum:
    print 'Warning: analysisSpectrum not set for', spectrum
    return

  contourDir = analysisSpectrum.contourDir
  if not contourDir:
    print 'Warning: contourDir not set for', analysisSpectrum
    return

  path = contourDir.dataLocation
  if fileName.startswith(path):  # dangerous
    fileName = fileName[len(path)+1:]
  analysisSpectrum.newStoredContour(dims=(xdim, ydim), path=fileName)

