
"""
======================COPYRIGHT/LICENSE START==========================

BlockUtil.py: Part of the CcpNmr Analysis program

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

import os

from memops.universal.BlockData import determineBlockSizes

from ccp.general.Io import getDataSourceFileName
from ccp.util.ShapeUtil import getShapeFile, get1dShapeFile

from ccp.api.general.DataLocation import BlockedBinaryMatrix

from ccpnmr.analysis.core.Util import getDimWrapped

def getBlockFile(spectrum, mem_cache, writeable = False):

  from memops.c import BlockFile

  fileName = getDataSourceFileName(spectrum)
  if (not fileName or (not writeable and not os.path.exists(fileName))):
    msg = 'Warning: spectrum (%s, %s): data file %s not accessible'
    print msg % (spectrum.experiment.name, spectrum.name, fileName)
    block_file = None

  else:

    points = [ dataDim.numPoints for dataDim in spectrum.sortedDataDims() ]
    
    dataStore = spectrum.dataStore
    
    if isinstance(dataStore, BlockedBinaryMatrix):
      blockSize = dataStore.blockSizes

      dimWrapped = getDimWrapped(spectrum)

      if dataStore.isBigEndian:
        isBigEndian = 1
      else:
        isBigEndian = 0

      isPadded = 1 # TBD: allow something else???

      if dataStore.numberType == 'int':
        isInteger = 1
      else:
        isInteger = 0

      blockHeaderSize = 0
      if hasattr(dataStore.root, 'application'):
        # TBD: was temporary hack, now use just in case it was set this way
        blockHeaderSize = dataStore.root.application.getValue(dataStore, 
           'blockHeaderSize', defaultValue=0, deleteAppData=True)
      if blockHeaderSize == 0:
        blockHeaderSize = dataStore.blockHeaderSize
      else:
        dataStore.blockHeaderSize = blockHeaderSize

      try:
        block_file = BlockFile.BlockFile(fileName, spectrum.numDim,
                             points, blockSize, dimWrapped, mem_cache,
                             dataStore.nByte, isBigEndian, isPadded,
                             dataStore.headerSize, isInteger, writeable, 
                             blockHeaderSize)
      except:  # TBD: temporary until blockHeaderSize C code is more widely available
        try:
          block_file = BlockFile.BlockFile(fileName, spectrum.numDim,
                             points, blockSize, dimWrapped, mem_cache,
                             dataStore.nByte, isBigEndian, isPadded,
                             dataStore.headerSize, isInteger, writeable)
        except BlockFile.error, e:
          print 'Warning, BlockFile error:', e
          block_file = None
    
    else:
      # probably a SHapeMatrix NBNB TBD how to handle this?
      block_file = None

  return block_file

def getShapeBlockFile(spectrum):

  from memops.c import BlockFile

  block_file = None

  # below is short-term hack
  if hasattr(spectrum, 'valuesList'):
    valuesList = spectrum.valuesList
  else:
    valuesList = None

  # TBD: why should we care about fileName: it's not used, is it?
  fileName = getDataSourceFileName(spectrum)
  if ((not fileName or not os.path.exists(fileName)) and 
      (spectrum.numDim > 1 or not valuesList)):
    msg = 'Warning: spectrum (%s, %s): data file %s not accessible'
    print msg % (spectrum.experiment.name, spectrum.name, fileName)

  else:

    points = [ dataDim.numPoints for dataDim in spectrum.sortedDataDims() ]
    
    ###blockSize = determineBlockSizes(points, totalBlockSize=64*1024)
    ###blockSize = determineBlockSizes(points, totalBlockSize=1)
    # below is short-term hack
    # want big block size in dimensions that being contoured
    # and small block size in orthogonal dimensions
    # so assume that first two dims are contoured (false in general)
    ndim = spectrum.numDim
    blockSize = ndim * [1]
    if ndim > 1:
      for i in range(2):
        blockSize[i] = min(8, points[i])
    else:
      blockSize[0] = points[0]

    dimWrapped = getDimWrapped(spectrum)

    if valuesList:
      shapeFile = get1dShapeFile(spectrum, valuesList)
    else:
      shapeFile = getShapeFile(spectrum)

    if shapeFile:
      if not fileName:
        fileName = ""
      try:
        block_file = BlockFile.ShapeBlockFile(fileName, ndim,
                             points, blockSize, dimWrapped, shapeFile)
      except BlockFile.error, e:
        print 'Warning, BlockFile error:', e

  return block_file
