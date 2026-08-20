"""
======================COPYRIGHT/LICENSE START==========================

Geometry.py: Utility code for CCPN code generation framework

Copyright (C) 2005  (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""
import math

def matrixMultiply(matrix1, matrix2):
  """ Multiply matrix1 by matrix2.
  """

  m = len(matrix1)
  for i in range(m):
    if i:
      if n != len(matrix1[i]):
        raise Exception('matrix1 not rectangular')
    else:
      n = len(matrix1[i])

  if n != len(matrix2):
    raise Exception('number of matrix1 columns != number of matrix2 rows')

  for j in range(n):
    if j:
      if p != len(matrix2[j]):
        raise Exception('matrix2 not rectangular')
    else:
      p = len(matrix2[j])

  matrix = []
  for i in range(m):
    matrix.append([])
    for j in range(p):
      matrix[i].append(0)
      for k in range(n):
        matrix[i][j] += matrix1[i][k] * matrix2[k][j]

  return matrix

def applyMatrix(matrix, v):
  """ Multiply matrix by vector v.
  """

  m = len(matrix)
  for i in range(m):
    if i:
      if n != len(matrix[i]):
        raise Exception('matrix not rectangular')
    else:
      n = len(matrix[i])

  if n != len(v):
    raise Exception('number of matrix columns != length of v')

  v2 = []
  for i in range(m):
    v2.append(0)
    for j in range(n):
      v2[i] += matrix[i][j] * v[j]

  return v2

def matrixRotX(angle):
  """ Return rotation matrix about x axis of angle (in radians).
  """

  c = math.cos(angle)
  s = math.sin(angle)

  R = [[ 1, 0, 0 ],
       [ 0, c,-s ],
       [ 0, s, c ]]

  return R

def matrixRotY(angle):
  """ Return rotation matrix about y axis of angle (in radians).
  """

  c = math.cos(angle)
  s = math.sin(angle)

  R = [[ c, 0, s ],
       [ 0, 1, 0 ],
       [-s, 0, c ]]

  return R

def matrixRotZ(angle):
  """ Return rotation matrix about z axis of angle (in radians).
  """

  c = math.cos(angle)
  s = math.sin(angle)

  R = [[ c,-s, 0 ],
       [ s, c, 0 ],
       [ 0, 0, 1 ]]

  return R

def matrixRot(angle, axis):
  """ Return rotation matrix about axis of angle (in radians).
  """

  axis = normaliseVector(axis)

  c = math.cos(angle)
  s = math.sin(angle)

  (n0, n1, n2) = axis

  R = [[ c + (1-c)*n0*n0,    (1-c)*n0*n1 - s*n2, (1-c)*n0*n2 + s*n1 ],
       [ (1-c)*n1*n0 + s*n2, c + (1-c)*n1*n1,    (1-c)*n1*n2 - s*n0 ],
       [ (1-c)*n2*n0 - s*n1, (1-c)*n2*n1 + s*n0, c + (1-c)*n2*n2 ]]

  return R

def vectorRotX(v, angle):
  """ Rotate the vector v by angle (in radians) about the x axis.
  """

  R = matrixRotX(angle)
  v = applyMatrix(R, v)

  return v

def vectorRotY(v, angle):
  """ Rotate the vector v by angle (in radians) about the y axis.
  """

  R = matrixRotY(angle)
  v = applyMatrix(R, v)

  return v

def vectorRotZ(v, angle):
  """ Rotate the vector v by angle (in radians) about the z axis.
  """

  R = matrixRotZ(angle)
  v = applyMatrix(R, v)

  return v

def vectorRot(v, angle, axis):
  """ Rotate the vector v by angle (in radians) about the axis.
  """

  R = matrixRot(angle, axis)
  v = applyMatrix(R, v)

  return v

def vectorLength(v):
  """ The length of v.
  """

  return math.sqrt(innerProduct(v, v))

def innerProduct(v1, v2):
  """ The inner product between v1 and v2.
  """

  n = len(v1)
  if (n != len(v2)):
    raise Exception('v1 and v2 must be same length')

  d = 0
  for i in range(n):
    d = d + v1[i] * v2[i]

  return d

def crossProduct(v1, v2):
  """ Returns the cross product of v1 and v2.
  Both must be 3-dimensional vectors.
  """

  if (len(v1) != 3 or len(v2) != 3):
    raise Exception('v1 and v2 must be of length 3')

  return [ v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0] ]

def normaliseVector(v):
  """ Normalises the vector v to be of unit length (not in place)
  """

  d = vectorLength(v)
  if (d == 0):
    raise Exception('v must be of non-zero length')

  v = [ x/d for x in v ]

  return v

def scaleVector(v, s):
  """ Scales the vector v by s
  """

  return [ s*x for x in v ]

def vectorsAngle(v1, v2):
  """ The angle (in radians) between v1 and v2.
  """

  l1 = vectorLength(v1)
  l2 = vectorLength(v2)
  if l1 == 0 or l2 == 0:
    return 0 # arbitrary

  d = innerProduct(v1, v2) / (l1*l2)

  return math.acos(d)

def vectorsAdd(v1, v2):
  """ Add vectors v1 and v2.
  """

  n = len(v1)
  if (n != len(v2)):
    raise Exception('length of v1 != length of v2')

  v = n * [0]
  for i in range(n):
    v[i] = v1[i] + v2[i]

  return v

def vectorsSubtract(v1, v2):
  """ Subtract vectors v1 and v2.
  """

  n = len(v1)
  if n != len(v2):
    raise Exception('length of v1 != length of v2')

  v = n * [0]
  for i in range(n):
    v[i] = v1[i] - v2[i]

  return v

def rotationOfVectorTowardsVector(v1, v2):
  """ Returns (angle, axis) for rotating the
  vector v1 to point towards the vector v2.
  Both vectors must be 3-dimensional.
  """

  defaultRotationAxis = (0, (1, 0, 0)) # arbitrary

  l1 = vectorLength(v1)
  if l1 == 0:
    return defaultRotationAxis

  l2 = vectorLength(v2)
  if l2 == 0:
    return defaultRotationAxis

  angle = vectorsAngle(v1, v2)
  try:
    axis = normaliseVector(crossProduct(v1, v2))
  except: # v1 parallel to v2
    if angle == 0: # parallel, no rotation needed
      return defaultRotationAxis

    # anti-parallel, need to find orthogonal direction
    (x, y, z) = v1
    if x != 0 or y != 0:
      axis = normaliseVector((y, -x, 0))
    else:
      axis = (1, 0, 0)

  return (angle, axis)

def rotationMatrixOfVectorTowardsVector(v1, v2):
  """ Returns rotation matrix for rotating the
  vector v1 to point towards the vector v2.
  Both vectors must be 3-dimensional.
  """

  (angle, axis) = rotationOfVectorTowardsVector(v1, v2)

  return matrixRot(angle, axis)

def determinantMatrix(matrix):

  if (len(matrix) != 3):
    raise Exception('matrix must be of length 3')

  for i in range(3):
    if (len(matrix[i]) != 3):
      raise Exception('matrix must be 3 x 3')

  return matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1]) \
       - matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[1][2]*matrix[2][0]) \
       + matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0])
  
def inverseMatrix(matrix):

  if (len(matrix) != 3):
    raise Exception('matrix must be of length 3')

  for i in range(3):
    if (len(matrix[i]) != 3):
      raise Exception('matrix must be 3 x 3')

  det = float(determinantMatrix(matrix))

  m  = [ [  matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1],
           -matrix[0][1]*matrix[2][2]+matrix[0][2]*matrix[2][1],
            matrix[0][1]*matrix[1][2]-matrix[0][2]*matrix[1][1] ],
         [ -matrix[1][0]*matrix[2][2]+matrix[1][2]*matrix[2][0],
            matrix[0][0]*matrix[2][2]-matrix[0][2]*matrix[2][0],
           -matrix[0][0]*matrix[1][2]+matrix[0][2]*matrix[1][0] ],
         [  matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0],
           -matrix[0][0]*matrix[2][1]+matrix[0][1]*matrix[2][0],
            matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0] ] ]

  for i in range(3):
    for j in range(3):
      m[i][j] = m[i][j] / det

  return m

def angleNormaliseAroundZero(angleDegrees):
  
  """
  Returns angle in degrees within 180 to -180 range
  """
  
  if angleDegrees > 180.0:
    angleDegrees -= (360.0 * (int( (angleDegrees-180) / 360) + 1))
  elif angleDegrees < -180.0:
    angleDegrees += (360.0 * (int( -(angleDegrees+180) / 360) + 1))

  return angleDegrees

def angleDifference(angle1, angle2, units = 'radians'):
  """ Returns difference in angles, taking into account
  that can have wrap-around.  Result is between 0 and pi
  for units of radians, and 0 and 180 for units of degrees.
  """

  if units == 'radians':
    w = 2 * math.pi
  else:
    w = 360

  angle = abs(angle1 - angle2) % w
  if angle > 0.5 * w:
    angle = w - angle

  return angle

def angleStats(angles, units = 'radians'):
  """ Returns approximate average and standard deviation
  of angles, using the largest gap between the sorted
  angles to determine where the angle boundary is.
  """

  pi = math.pi
  twopi = 2 * pi

  # convert to radians
  if units != 'radians':
    angles = [ pi * angle / 180.0 for angle in angles ]

  # put in range 0 to 2 pi
  angles = [ angle % twopi for angle in angles ]

  # find biggest gap between angles and consider the midpoint of that
  # to be the boundary of the angles
  # (an alternative was to use cosines to transform the problem,
  #  but it gave not great answers on simple tests)
  n = len(angles)
  angles.sort()
  dmax = 0
  imax = 0
  for i in range(n):
    d = (angles[(i+1)%n] - angles[i]) % twopi
    if d > dmax:
      imax = i
      dmax = d

  if imax == n-1: #  largest gap is between last and first angles
    boundary = (pi + 0.5 * (angles[0] + angles[n-1])) % twopi
  else:
    boundary = 0.5 * (angles[imax] + angles[imax+1])

  # if so then use standard (adjusted average)
  # average and standard deviation for cosines
  sum = 0
  sum2 = 0
  for angle in angles:
    d = (angle - boundary) % twopi
    sum = sum + d
    sum2 = sum2 + d*d
  
  n = float(n)
  avg = sum / n
  std = math.sqrt(max(0, sum2/n - avg*avg))
  avg = (avg + boundary) % twopi

  if units != 'radians':
    avg = 180.0 * avg / pi
    std = 180.0 * std / pi

  return (avg, std)

def superposeNewVectorsOnOld(oldVectors, newVectors, phi = 0):
  """  Superposes new vectors on old.  The transformed vectors are returned.
  Specifically, oldVectors should be of length 2, and newVectors of length >= 2.
  newVectors is translated and rotated so that newVectors[0] ends up at oldVectors[0]
  and so that newVectors[1]-newVectors[0] ends up in same direction as
  vector oldVectors[1]-oldVectors[1].  Finally, a rotation of phi (which should
  be in radians) is done about the (transformed) 01 vector.
  """

  # first translate so that vectors relative to first point of old and new
  old01 = vectorsSubtract(oldVectors[1], oldVectors[0])

  vectors = []
  for vector in newVectors[1:]:
    vectors.append(vectorsSubtract(vector, newVectors[0]))

  # now rotate new 01 vector to lie along old 01 vector
  (angle, axis) = rotationOfVectorTowardsVector(vectors[0], old01)
  for i in range(len(vectors)):
    vectors[i] = vectorRot(vectors[i], angle, axis)

  # now rotate around new 01 vector by phi
  for i in range(1, len(vectors)):
    vectors[i] = vectorRot(vectors[i], phi, vectors[0])

  # finally translate so that first point of new is where first point of old was
  for i in range(len(vectors)):
    vectors[i] = vectorsAdd(vectors[i], oldVectors[0])
  
  vectors.insert(0, oldVectors[0])

  return vectors
