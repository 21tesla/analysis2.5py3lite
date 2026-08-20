"""
======================COPYRIGHT/LICENSE START==========================

Geometry.py: Geometry code relevant to molecular systems

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""
import math

import memops.universal.Geometry as uniGeometry

EPSILON = 1.0e-6

def convertDegreeRadian(angle):

  return math.pi * angle / 180

def calcPlanarCoord(coord1, coord2, coord3, dist24, angle124 = None):
  """ Calculate the coordinates of the fourth point of a planar set
      given the coordinates of the first three points, the distance
      from the central point (number 2) and optionally the angle 124
      (if latter not given then splits the difference).
  """

  x12 = uniGeometry.vectorsSubtract(coord2, coord1)
  x32 = uniGeometry.vectorsSubtract(coord2, coord3)

  if angle124 is None:
    angle123 = uniGeometry.vectorsAngle(x12, x32)
    angle124 = math.pi - 0.5*angle123

  n12 = uniGeometry.normaliseVector(x12)
  d = uniGeometry.innerProduct(n12, x32)
  x32par = uniGeometry.scaleVector(n12, d) # component of x32 parallel to x12
  x32perp = uniGeometry.vectorsSubtract(x32, x32par) # component of x32 perpendicular to x12
  n32perp = uniGeometry.normaliseVector(x32perp)

  d24par = - dist24 * math.cos(angle124)
  d24perp = dist24 * math.sin(angle124)

  x24par = uniGeometry.scaleVector(n12, d24par)
  x24perp = uniGeometry.scaleVector(n32perp, d24perp)
  x24 = uniGeometry.vectorsAdd(x24par, x24perp)

  coord4 = uniGeometry.vectorsAdd(coord2, x24)

  return coord4

def calcTorsionAngleRadians(coord1, coord2, coord3, coord4):
  """ Calculate the torsion angle between four points in radians.
  """

  v12 = uniGeometry.vectorsSubtract(coord1, coord2)
  v32 = uniGeometry.vectorsSubtract(coord3, coord2)
  v43 = uniGeometry.vectorsSubtract(coord4, coord3)

  vn13 = uniGeometry.crossProduct(v12, v32)
  vn24 = uniGeometry.crossProduct(v43, v32)

  d12 = uniGeometry.innerProduct(vn13, vn24)
  d11 = uniGeometry.innerProduct(vn13, vn13)
  d22 = uniGeometry.innerProduct(vn24, vn24)

  if d11 < EPSILON or d22 < EPSILON:  # points collinear
    ang = 0.0  # arbitrary
  else:
    ang = d12 / math.sqrt(d11*d22)
    if ang >= 1.0:
      ang = 0.0
    elif ang <= -1.0:
      ang = - math.pi
    else:
      ang = math.acos(ang)
      d = uniGeometry.innerProduct(vn13, uniGeometry.crossProduct(vn24, v32))
      if d < 0:
        ang = -ang

  return ang

def calcTorsionAngleDegrees(coord1, coord2, coord3, coord4):
  """ Calculate the torsion angle between four points in degrees.
  """

  ang = calcTorsionAngleRadians(coord1, coord2, coord3, coord4) * 180.0 / math.pi
  ang = max(min(ang, 180.0), -180.0)

  return ang


def calcAngleViolation(angleDegrees,upperLimit,lowerLimit):
  
  """ 
  Calculates (in degrees) the violation for an angle given an upperLimit and lowerLimit, where the allowed range
  goes from the upperLimit in decreasing degrees to the lowerLimit
  """

  violationDegrees = 0.0
  
  if angleDegrees > upperLimit:
    # Can be either violated, or still within lower limit (but only if lower limit *higher*!)
    if lowerLimit > upperLimit:
      if angleDegrees > lowerLimit:
        pass
      else:
        violationDegrees = min(abs(angleDegrees - upperLimit), abs(lowerLimit - angleDegrees))
    else:
      violationDegrees = min(abs(angleDegrees - upperLimit), abs(lowerLimit + 360 - angleDegrees))

  elif angleDegrees < lowerLimit:
    # Can be either violated, or still within upper limit (but only if upper limit *lower*!)
    if upperLimit < lowerLimit :
      if angleDegrees < upperLimit:
        pass
      else:
        violationDegrees = min(abs(angleDegrees - upperLimit), abs(lowerLimit - 360 - angleDegrees))
    else:
      violationDegrees = min(abs(angleDegrees - upperLimit), abs(angleDegrees - lowerLimit))
      
  return violationDegrees
    
def getDistanceFromCoordinates(coord1,coord2):
  
  """   
  Function to calculate distance between two Structure.Coordinates
  """
  
  return math.sqrt(math.pow(coord1.x-coord2.x,2) + math.pow(coord1.y-coord2.y,2) + math.pow(coord1.z-coord2.z,2))

if __name__ == '__main__':

  coord1 = (-0.1858, -0.5331, 0.4941)
  coord2 = (-1.4475, 0.273, 0.6661)
  coord3 = (-2.3181, 0.2173, -0.1689)

  x12 = uniGeometry.vectorsSubtract(coord2, coord1)
  x32 = uniGeometry.vectorsSubtract(coord2, coord3)

  print 'dist12 = %3.2f' % uniGeometry.vectorLength(x12)
  print 'dist32 = %3.2f' % uniGeometry.vectorLength(x32)
  print 'angle123 = %3.2f' % (uniGeometry.vectorsAngle(x12, x32)*180/math.pi)

  dist24 = 1.33
  #angle124 = None
  angle124 = 115.6 * math.pi / 180
  coord4 = calcPlanarCoord(coord1, coord2, coord3, dist24, angle124=angle124)

  x42 = uniGeometry.vectorsSubtract(coord2, coord4)

  print 'coord4 = %s' % (['%4.3f' % x for x in coord4])
  print 'dist24 = %3.2f' % uniGeometry.vectorLength(x42)
  print 'angle124 = %3.2f' % (uniGeometry.vectorsAngle(x12, x42)*180/math.pi)
  print 'angle324 = %3.2f' % (uniGeometry.vectorsAngle(x32, x42)*180/math.pi)

  print 'torsion_angle1234 = %3.2f' % calcTorsionAngleDegrees(coord1, coord2, coord3, coord4)

  coord4 = coord3
  print 'torsion_angle1234 = %3.2f' % calcTorsionAngleDegrees(coord1, coord2, coord3, coord4)
