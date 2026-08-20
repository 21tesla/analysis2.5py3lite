"""
======================COPYRIGHT/LICENSE START==========================

Structure.py: Utility functions

Copyright (C) 2009 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

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

# Entire File obsolete and incorrect Deprecated
pass

'''
def alignCoordinates(coords0, coords1, weight=None, translate0=True, translate1=True):
  """Descrn: Align coords1 to coords0, using weight and translating coords0 if translate0 and translating coords1 if translate1
     Inputs: List of List of floats, List of List of floats (both of size [ncoords][3]), List of floats (of size ncoords), boolean, boolean
     Output: (coords: List of List of floats, err: float, translation0: List of floats, translation1: List of floats, rotation: List of List of floats)
  """

  try:
    from numpy import array, cross, dot, outer, sqrt, sum, transpose
    from numpy.linalg import eig
  except:
    print 'ERROR: cannot use alignCoordinates without numpy'
    raise

  coords0 = array(coords0)
  coords1 = array(coords1)

  if weight:
    weight = array(weight)
    w = float(sum(weight))
    if translate0:
      translation0 = - sum(transpose(weight*transpose(coords0))) / w
      coords0 = coords0 + translation0
    else:
      translation0 = None

    if translate1:
      translation1 = - sum(transpose(weight*transpose(coords1))) / w
      coords1 = coords1 + translation1
    else:
      translation1 = None

    e = sum(weight*transpose(coords0*coords0)) + sum(weight*transpose(coords1*coords1))
    r = dot(weight*transpose(coords0), coords1)

  else:
    ncoords = float(len(coords0))
    if translate0:
      translation0 = - sum(coords0, 0) / ncoords
      coords0 = coords0 + translation0
    else:
      translation0 = None

    if translate1:
      translation1 = - sum(coords1, 0) / ncoords
      coords1 = coords1 + translation1
    else:
      translation1 = None

    e = sum(coords0*coords0) + sum(coords1*coords1)
    r = dot(transpose(coords0), coords1)

  e = 0.5 * e

  rtr = dot(transpose(r), r)

  (eigenvalues, eigenvectors) = eig(rtr)
  eigenvectors = transpose(eigenvectors)

  # sort into decreasing order
  (lambda0, lambda1, lambda2) = eigenvalues
  lambdaMax = max(eigenvalues)
  if lambda0 == lambdaMax:
    if lambda1 >= lambda2:
      v = (0, 1, 2)
    else:
      v = (0, 2, 1)
  elif lambda1 == lambdaMax:
    if lambda0 >= lambda2:
      v = (1, 0, 2)
    else:
      v = (1, 2, 0)
  else:
    if lambda0 >= lambda1:
      v = (2, 0, 1)
    else:
      v = (2, 1, 0)

  # the max with 0 is in case of rounding errors which can give small negative number
  eigenvalues = array((max(eigenvalues[v[0]], 0), max(eigenvalues[v[1]], 0), max(eigenvalues[v[2]], 0)))
  eigenvectors = array((eigenvectors[v[0]], eigenvectors[v[1]], eigenvectors[v[2]]))

  (s0, s1, s2) = sqrt(eigenvalues)

  e0 = eigenvectors[0]
  e1 = eigenvectors[1]
  e2 = cross(e0, e1)

  v0 = dot(r, e0)
  v1 = dot(r, e1)
  v2 = dot(r, e2)

  b0 = v0 / s0
  b1 = v1 / s1
  b2 = cross(b0, b1)

  rotation = outer(b0, e0) + outer(b1, e1) + outer(b2, e2)

  if sum(b2*v2) < 0:
    s = -1
  else:
    s = 1

  err = e - (s0 + s1 + s*s2)

  coords1 = transpose(dot(rotation,transpose(coords1)))

  return (coords1.tolist(), err, translation0.tolist(), translation1.tolist(), rotation.tolist())
'''
