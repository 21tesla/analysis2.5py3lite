"""
======================COPYRIGHT/LICENSE START==========================

BoxIntegral.py: Part of the CcpNmr Analysis program

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
import operator

# code based on rewriting and modification of integrator code of
# Felician Dancea and Ulrich Gunther, Birmingham, copyright (c) 2005

from memops.universal.BlockData import cumulativeArray, arrayOfIndex, indexOfArray

def distance2Center(posn, center, scale):

  d2 = 0
  for i in range(len(posn)):
    dx = posn[i] - center[i]
    d2 = d2 + scale[i] * dx * dx

  return d2

def cmpDist(posn0, posn1, dists):

  d0 = dists[posn0]
  d1 = dists[posn1]

  if d0 < d1:
    return -1
  elif d0 == d1:
    return 0
  else:
    return 1

def firstOrderPosns(posn, dists, boxSize, center, neighbors):

  ndim = len(boxSize)
  d = dists[posn]
  p = ndim * [0]

  positions = []
  for neighbor in neighbors:
    for i in range(ndim):
      # first check that neighbor going in direction towards center
      r = center[i] - posn[i]
      if (r > 0):
        if (neighbor[i] < 0):
          break
      elif (r < 0):
        if (neighbor[i] > 0):
          break
      else: # r = 0
        if (neighbor[i] != 0):
          break
      # now check that neighbor within box
      p[i] = posn[i] + neighbor[i]
      if (p[i] < 0 or p[i] >= boxSize[i]):
        break
    else:
      q = tuple(p)
      # now check that neighbor closer to center
      if (dists[q] < d):
        positions.append(q)

  return positions

def calcNeighbors(ndim):

  nearSize = ndim * [3]
  (n, cumulNear) = cumulativeArray(nearSize)

  center = tuple(ndim * [0])
  neighbors = []
  for i in range(n):
    neighbor = arrayOfIndex(i, cumulNear)
    neighbor = tuple([n-1 for n in neighbor])
    if (neighbor != center):
      neighbors.append(neighbor)

  return neighbors

def truncatedBoxIntegral(values, boxSize, center = None, scale = None, noiseLevel = 0):
  """ Integrate a peak represented as a list of values.
      values = the peak values as a 1D tuple or list (so not an ND tuple or list).
      boxSize = a tuple or list specifying the dimensions of values.
      center = a tuple or list specifying what is the peak center point.
        if center = None then it is set to center of box.
      scale = a tuple or list specifying how (squared) distances should be weighted in different dimensions.
        if scale = None then it is set to 1 in all dimensions.
      noiseLevel = number specifying what the noiseLevel of the spectrum is.
        if noiseLevel = 0 or if it is greater than the central peak value then it is set to 0.05 * central peak value.
      The function zeroes those elements of values it considers to be outside the peak.
      Function returns the sum of values (but as said, values is also updated).
  """

  ndim = len(boxSize) # = len(center)

  if scale is None:
    scale = ndim * [1]

  if (len(center) != ndim):
    raise Exception('len(center) != len(boxSize)')

  if (len(scale) != ndim):
    raise Exception('len(scale) != len(boxSize)')

  (n, cumulBox) = cumulativeArray(boxSize)
  if (len(values) != n):
    raise Exception('len(values) != product(boxSize)')

  centerInd = indexOfArray(center, cumulBox)

  neighbors = calcNeighbors(ndim)

  v = abs(values[centerInd])
  if noiseLevel == 0 or noiseLevel > v:
    noiseLevel = 0.05 * v

  dists = {}
  posns = []
  for i in range(n):
    posn = arrayOfIndex(i, cumulBox)
    posns.append(posn)
    dists[posn] = distance2Center(posn, center, scale)

  func = lambda posn1, posn2: cmpDist(posn1, posn2, dists)
  posns.sort(func)
  posns = posns[1:] # exclude center itself

  for posn in posns:
    index = indexOfArray(posn, cumulBox)
    v = abs(values[index])
    if v <= noiseLevel:
      values[index] = 0
    else:
      nearest = firstOrderPosns(posn, dists, boxSize, center, neighbors)
      for near in nearest:
        ind = indexOfArray(near, cumulBox)
        # if nearer neighbor has lower value then index is not in peak so zero
        if abs(values[ind]) < v:
          values[index] = 0
          break

  s = reduce(operator.add, values)

  return s

if __name__ == '__main__':

  import sys

  def trans2dmat(mat,boxSize):
    
    matout = len(mat)*[0]
    (n, cumulBox)  = cumulativeArray(boxSize)
    (n, cumulBoxt) = cumulativeArray([boxSize[1],boxSize[0]])
    
    for i in range(n):
      posn = arrayOfIndex(i, cumulBox)
      ind  = indexOfArray((posn[1], posn[0]), cumulBoxt)
      matout[ind] = mat[i]
        
    return matout

  def printMat(mat, dims):

    ndim = len(dims)
    if (ndim == 1):
      print mat
    elif (ndim == 2):
      for i in range(dims[1]):
        for j in range(dims[0]):
          sys.stdout.write('%5.1f' % mat[i*dims[0]+j])
        sys.stdout.write('\n')
    elif (ndim == 3):
      n = dims[0] * dims[1]
      d = dims[:2]
      for i in range(dims[2]):
        printMat(mat[i*n:(i+1)*n], d)
        sys.stdout.write('\n')
    else:
      raise 'only does ndim <= 3'

  values = [ \
    1, 1, 1, 1, 1,
    2, 2, 2, 2, 1,
    1, 2, 3, 2, 1,
    1, 2, 2, 2, 1,
    1, 1, 1, 1, 1,
  ]
  boxSize = [5, 5]
  center = [2, 2]

  print 'values before:'
  printMat(values, boxSize)
  s = truncatedBoxIntegral(values, boxSize, center)
  print 'values after:'
  printMat(values, boxSize)
  print 'integral = ', s
  print

  values = [ \
    2, 3, 1,
    3, 2, 2,
    2, 3, 2,
    2, 2, 2,
    1, 1, 1,
  ]
  boxSize = [3, 5]
  center = [1, 2]

  print 'values before:'
  printMat(values, boxSize)
  s = truncatedBoxIntegral(values, boxSize, center)
  print 'values after:'
  printMat(values, boxSize)
  print 'integral = ', s
  print

  mat1d = [ 1, 2, 3, 2, 1 ]
  center1d = [ 2 ]
  boxSize1d = [ 5 ]

  mat2d = [ 1, 1, 1, 1, 1, 1, 2, 2, 2, 1 ] + mat1d + [ 1, 2, 2, 2, 1, 1, 1, 1, 1, 1 ]
  center2d = [ 2, 2 ]
  boxSize2d = [ 5, 5 ]

  mat3d = [ 0.5 * x for x in mat2d ] + mat2d + [ 0.5 * x for x in mat2d ]
  center3d = [ 2, 2, 1 ]
  boxSize3d = [ 5, 5, 3 ]

  for (mat, center, boxSize) in ((mat1d, center1d, boxSize1d), (mat2d, center2d, boxSize2d), (mat3d, center3d, boxSize3d)):
    dims = boxSize
    print 'values before:'
    printMat(mat, dims)
    s = truncatedBoxIntegral(mat, boxSize, center)
    print 'values after:'
    printMat(mat, dims)
    print 'integral = ', s
    print

  values = [ \
    1, 1, 1, 1, 1, 1, 1,
    1, 2, 2, 2, 2, 2, 1,
    1, 2, 4, 4, 4, 2, 1,
    2, 2, 4, 5, 4, 2, 1,
    5, 4, 4, 4, 4, 3, 1,
    5, 2, 2, 2, 3, 3, 1,
    1, 1, 1, 1, 1, 1, 1
  ]
  boxSize = [7, 7]
  center  = [3, 3]
    
  values2  = trans2dmat(values,boxSize)
  boxSize2 = [7, 7]
  center2  = [3, 3]    
  
  print 'values before:'
  printMat(values, boxSize) 
  s = truncatedBoxIntegral(values, boxSize, center)
  print 'values after:'
  printMat(values, boxSize)    
  print 'integral = ', s
  print
  
  print 'values before:'
  printMat(values2, boxSize2) 
  s = truncatedBoxIntegral(values2, boxSize2, center2)
  print 'values after:'
  printMat(values2, boxSize2)    
  print 'integral = ', s
  print  
    
  values = [ \
    1, 1, 1, 1, 1,
    2, 2, 2, 2, 2,
    2, 4, 4, 4, 2,
    2, 4, 5, 4, 2,
    2, 4, 4, 2, 3,
    2, 0, 2, 3, 3,
    1, 1, 1, 1, 1,
  ]
  boxSize = [5, 7]
  center  = [2, 3]
  
  values2 = trans2dmat(values,boxSize)
  boxSize2 = [7,5]
  center2 = [3,2]
  
  print 'values before:'
  printMat(values, boxSize) 
  s = truncatedBoxIntegral(values, boxSize, center)
  print 'values after:'
  printMat(values, boxSize)
  print 'integral = ', s
  print
  
  print 'values before:'
  printMat(values2, boxSize2) 
  s2 = truncatedBoxIntegral(values2, boxSize2, center2)
  print 'values after:'
  printMat(values2, boxSize2)
  print 'integral = ', s2
  print  
  
  values = [1, 2, 3, 2, 3, 4, 5, 4, 3, 2, 6, 5, 4, 7, 3, 2, 1]
  boxSize = [17]
  center  = [6]
  print 'values before:'
  printMat(values, boxSize)
  s = truncatedBoxIntegral(values, boxSize, center)
  print 'values after:'
  printMat(values, boxSize)     
  print 'integral = ', s
  print
