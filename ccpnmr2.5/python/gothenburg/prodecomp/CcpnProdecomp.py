"""Code to prepare, analyse, wrap, etc. Prodecomp

WARNING!  TOXIC CODE!

This code contains numpy objects that may break normal Python operations.
Please make sure that numpy objects are never passed to code in other files,
unless the other file also imports numpy.
"""

import sys, os, struct

try:
  from gothenburg.prodecomp import prodecomp
except ImportError:
  class prodecomp:
    publicDocumentation = '*WARNING* PRODECOMP is not installed'

from memops.general.Implementation import ApiError
from ccp.general.Io import getDataSourceFileName

def runProdecomp(dataSources, pythonDefsMatrix, intl, cmps, rglf, itrs):

  from numpy import zeros, array

  # get DataMatrix
  dataMatrix = getDataMatrix(dataSources)

  # convert defs matrix to numpy
  defs = zeros((len(pythonDefsMatrix), len(pythonDefsMatrix[0])))
  for ii, factors in enumerate(pythonDefsMatrix):
    defs[ii,:] = array(factors)

  # run prodecomp
  fdir, f = prodecomp.prodecomp(dataMatrix, defs, intl, cmps, rglf, itrs)

  return (fdir, f)


def getDataMatrix(dataSources):

  from numpy import zeros


  # get data matrix
  dataPlane = getDataPlane(dataSources[0])
  size1, size2 = dataPlane.shape

  dataMatrix = zeros((size1, size2, len(dataSources)))
  dataMatrix[:,:,0] = dataPlane

  for ii, dataSource in enumerate(dataSources[1:]):
    dataMatrix[:,:,ii+1] = getDataPlane(dataSource)
  #
  return dataMatrix


def getDataPlane(dataSource):
  """ get plane of data.
  Currently returns acquisition dimension as fastest varying
  Currently always returns entire plane
  """

  from numpy import array

  # data type in file - may vary between formats
  btype = 'i'

  dataStore = dataSource.dataStore

  # tests for not-yet-implemented:
  # may be handled in the future, using new optional function parameters
  if dataSource.numDim != 2:
    # Assumes 2D
    raise ApiError('%s: num dims %s not supported'
                   % (dataSource, dataSource.numDim))

  if dataStore.blockSizes!= dataStore.numPoints:
    # assumes unblocked data
    raise ApiError('%s: file is blocked. numPoints:%s, blockSizes:%s'
                   % (dataSource, dataStore.numPoints, dataStore.blockSizes))

  if True in dataStore.isComplex:
    # assumes all-real data
    raise ApiError('%s: data have complex dimension - isComplex:'
                   % (dataSource, dataStore.isComplex))

  # active code starts

  # check number of points
  points1 = dataStore.numPoints
  points2 = tuple(x.numPoints for x in dataSource.sortedDataDims())
  if points1 != points2:
    raise ApiError("Number of points differ between file (%s) and spectrum (%s)"
                   % (points1, points2))

  # get and check file
  filePath = getDataSourceFileName(dataSource)
  if filePath is None:
    raise ApiError('No data file found for %s' % dataSource)

  # read data from file
  buf = open(filePath).read()
  nPoints = len(buf)/struct.calcsize(btype)
  data = array(struct.unpack(str(nPoints)+btype, buf))

  # convert data to array
  data = data.reshape(points1[1],points1[0])

  #
  return data
