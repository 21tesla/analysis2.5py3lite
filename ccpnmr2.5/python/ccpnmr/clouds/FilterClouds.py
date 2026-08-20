"""
  These functions provide quality control filtering of CLOUDS

======================COPYRIGHT/LICENSE START==========================

FilterClouds.py: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
from ccpnmr.clouds.CloudBasic import writeTypedPdbCloud

from math import sqrt

def alignClouds(clouds, names):

  resonances = clouds[0].keys()

  cloudsList     = []
  for cloud in clouds:
    orderCloud = []
    for resonance in resonances:
      orderCloud.append(cloud.get(resonance) or (0.0,0.0,0.0))
    cloudsList.append(orderCloud)

  (meanCloud,cloudsList) = alignToMeanCloud(cloudsList)
  alignCloudsToRef(meanCloud,cloudsList)
  
  for i in range(len(clouds)):
    pdbFileName = names[i]
    print pdbFileName
    writeTypedPdbCloud(cloudsList[i], pdbFileName, resonances)

def filterClouds(clouds, atomTypes=None):

  cloudsList     = []
  meanGroupRmsds = []
  numClouds      = len(clouds)
      
  resonances = []
  for resonance in clouds[0].keys():
    if atomTypes:
      for assignName in resonance.assignNames:
        if assignName in atomTypes:
          resonances.append(resonance)
          break
      
    else:
      resonances.append(resonance)
  
  for cloud in clouds:
  
    orderCloud = []
    for resonance in resonances:
      orderCloud.append(cloud.get(resonance) or (0.0,0.0,0.0))
    
    cloudsList.append(orderCloud)
      
  print "Generating mean and aligning"
  
  (meanCloud,cloudsList) = alignToMeanCloud(cloudsList)
  #minToMeanRmsd          = getMeanPairRmsd(meanCloud,[cloudsList[0],])
  #bestCloud = cloudsList[0]
  #for i in range(1, numClouds):
  #  meanPairRmsd = getMeanPairRmsd(meanCloud,[cloudsList[i],])
  #  if meanPairRmsd < minToMeanRmsd:
  #    minToMeanRmsd = meanPairRmsd
  #    bestCloud     = cloudsList[i]
 
  scores = []
  alignCloudsToRef(meanCloud,cloudsList)
  for i in range(numClouds):
    
    #(newMean,newClouds) = alignCloudsToRef(bestCloud,[cloudsList[i],])
    meanPairRmsd         = getMeanPairRmsd(meanCloud,[cloudsList[i],])
    meanGroupRmsds.append(meanPairRmsd)
            
  return meanGroupRmsds

def alignCloudsToRef(refCoords,coordsList):

  assert refCoords and coordsList
  
  N = len(coordsList)
  
  assert len(coordsList[0]) == len(refCoords)

  # align centres of masses of clouds to the origin
  centerCoords(refCoords)
  for i in range(N):
    cloud = list(coordsList[i])
    coordsList[i] = centerCoords(cloud)

  # calculate weights: generically not always 1, e.g. could be according to atomic mass
  weights = []
  for atom in refCoords:
    weights.append(1)
    
  # align each cloud and its inverse with the refCoords
  for i in range(N):
  
    coords  = coordsList[i]
    inverse = []
    for (x,y,z) in coords:
      inverse.append([-x,-y,-z])
 
    (tryCoords1,error1,rotMat1) = alignCoordinates(refCoords, coords,  weights)
    (tryCoords2,error2,rotMat2) = alignCoordinates(refCoords, inverse, weights)


    # update cloud coords given optimum rotation
    if error1 <= error2:
      coordsList[i] = tryCoords1
    else:
      coordsList[i] = tryCoords2
    
  # calc mean pairwise rmsd
    
  return (refCoords,coordsList)

def matrixDeterminant3d(M):

  d  = M[0][0]*((M[1][1]*M[2][2])-(M[1][2]*M[2][1]))
  d -= M[0][1]*((M[1][0]*M[2][2])-(M[1][2]*M[2][0]))
  d += M[0][2]*((M[1][0]*M[2][1])-(M[1][1]*M[2][0]))  
  
  return d
  
def scalarProduct3d(V,W):

  N = len(V)
  assert N == len(W)
  
  p = 0
  for i in range(N):
    p += (V[i]*W[i])

  return p 

def vectorScale(V, f):

  W = []
  for i in range(len(V)):
    W.append(V[i]*f)
    
  return W  
   
def vectorModulus(V):

  s = 0
  for v in V:
    s += v*v
  
  return sqrt(s)

def unitVector(V):
 
  m = vectorModulus(V)
  
  if m != 0:
    U = []
    for i in range(len(V)):
      U.append(V[i]/m)
 
    return U

def matrixColumnVectors(M):

  C = []
  for j in range(len(M[0])):
    C.append([])    
    for i in range(len(M)):
      C[j].append(M[i][j]) 

  return C

def vectorProduct3d(v,w):

  assert len(v) == len(w) == 3
  
  return [v[1]*w[2]-v[2]*w[1],v[2]*w[0]-v[0]*w[2],v[0]*w[1]-v[1]*w[0]]

def matrixIdentity(n):

  I = []
  for p in range(n):
    I.append([])
    for q in range(n):
      I[p].append(0.0)
    I[p][p] = 1.0
  
  return I

def matrixTranspose(M):

  Mt = []
  for j in range(len(M[0])):
    Mt.append([])
    for i in range(len(M)):
      Mt[j].append(M[i][j])
      
  return Mt

def matrixVecMultiply(M, V):

  assert len(V) == len(M[0])
  
  L = []
  for i in range(len(M)):
    l = 0
    for j in range( len(M[0]) ):
      l += M[i][j]*V[j]
    L.append(l)
  
  return L

def matrixMultiply(M, N):

  assert len(N) == len(M[0])

  L = []
  for i in range(len(M)):
    L.append([])
    for j in range(len(N)):
      l = 0
      for k in range(len(M[i])):
        l += N[k][i]*M[j][k]
      L[i].append(l)
 
  return L

def alignToMeanCloud(coordsList):

  meanCoords = []
  N = len(coordsList[0])
  M = len(coordsList)
  alignCloudsToRef(coordsList[0],coordsList[1:])
  
  for i in range(N):
    x = 0
    y = 0
    z = 0
    for j in range(M):
      x += coordsList[j][i][0]
      y += coordsList[j][i][1]
      z += coordsList[j][i][2]
    x /= M
    y /= M
    z /= M
    meanCoords.append([x,y,z])
      
  return (meanCoords,coordsList)

def getMeanPairRmsd(refCoords,coordsList):
 
  numAtoms  = len(refCoords)
  numClouds = float(len(coordsList))
  rmsdSum   = 0.0

  for cloud in coordsList:
    d2 = 0.0
    for i in range(numAtoms):
      dx  = cloud[i][0] - refCoords[i][0]
      dy  = cloud[i][1] - refCoords[i][1]
      dz  = cloud[i][2] - refCoords[i][2]
      d2 += (dx*dx) + (dy*dy) + (dz*dz)
    rmsdSum += sqrt(d2/numAtoms)

  return rmsdSum/numClouds

def centerCoords(coords):

  if coords:
    centreOfMass = getMeanCoords(coords)
    vector = [-c for c in centreOfMass]
    moveCoords(coords, vector)
    return coords
 
def moveCoords(coords, vector):

  (dx,dy,dz) = vector
  
  for coord in coords:
    coord[0] += dx
    coord[1] += dy
    coord[2] += dz

  return coords

def getMeanCoords(coords):

  N = len(coords)
  x = 0
  y = 0
  z = 0
  for coord in coords:
    x += coord[0]
    y += coord[1]
    z += coord[2]
    
  if N > 0:
    x /= N
    y /= N
    z /= N
   
  return (x,y,z) 

def removeDisconnectedAtoms(atomCoordList):

  coords = []
  newCoords = []
  for atomCoord in atomCoordList:
    coords.append( [atomCoord[0],atomCoord[1],atomCoord[2]] )

  (cx,cy,cz) = getMeanCoords(coords)
  
  m = 0
  n = 0
  D = []
  for coord in coords:
    dx2 = ( coord[0] - cx ) * ( coord[0] - cx )
    dy2 = ( coord[1] - cy ) * ( coord[1] - cy )
    dz2 = ( coord[2] - cz ) * ( coord[2] - cz )
    d = math.sqrt(dx2+dy2+dz2)
    D.append(d)
    m += d
    n += 1

  m /= n
  s = 0
  for i in range(n):
    s += ( D[i] - m ) * ( D[i] - m )
  s /= n
  s2 = sqrt(s)

  for i in range(n):
    if abs(D[i] - m) <= 2*s2:
      newCoords.append(coords[i])
     
  return newCoords

def eigenSort(d,V):

  n = len(d)

  if (n != len(V)) or (n != len(V[0])):
    return

  for i in range(n-1):
    k = i
    p = d[k]
    for j in range(i+1, n):
      if d[j] >= p:
        k = j
        p = d[k]
    if k != i:
      d[k] = d[i]
      d[i] = p
      for j in range(n):
        p = V[j][i]
        V[j][i] = V[j][k]
        V[j][k] = p

  return (d,V)

def alignCoordinates(coords0, coords1, W):

  N = len(coords1)
  assert N == len(coords0) == len(W)
  
  X = coords1
  Y = coords0

  E0 = 0
  for i in range(N):
    E0 += W[i] * ( scalarProduct3d(X[i],X[i]) + scalarProduct3d(Y[i],Y[i]) )

  E0 = E0/2

  R = []
  for i in range(3):
    R.append([])
    for j in range(3):
      r = 0
      for n in range(N):
        r += W[n]*Y[n][i]*X[n][j]
      R[i].append(r)

  RtR = matrixMultiply(matrixTranspose(R),R)
  
  (eigenVals,eigenVecsMatrix,nrot) = eigenJacobi(RtR)

  (eigenVals,eigenVecsMatrix)      = eigenSort(eigenVals,eigenVecsMatrix)
  
  eigenVecs = matrixTranspose(eigenVecsMatrix)
    
  eigenVecs[2] = vectorProduct3d(eigenVecs[0],eigenVecs[1])

  Ra0 = matrixVecMultiply(R,eigenVecs[0])
  Ra1 = matrixVecMultiply(R,eigenVecs[1])
  Ra2 = matrixVecMultiply(R,eigenVecs[2])
  
  b0 = vectorScale(Ra0,1/sqrt(eigenVals[0]))
  b1 = vectorScale(Ra1,1/sqrt(eigenVals[1]))
  
  b  = [b0,b1,vectorProduct3d(b0,b1)]
  
  U = []
  for i in range(3):
    U.append([])
    for j in range(3):
      u = 0
      for k in range(3):
        u += b[k][i]*eigenVecs[k][j]
      U[i].append(u)
  
  s2 = 1
  if scalarProduct3d(b[2],Ra2) < 0:
    s2 = -1
    
  E = E0 - sqrt(eigenVals[0]) - sqrt(eigenVals[1]) - (s2*sqrt(eigenVals[2]))
  
  for n in range(N):
    [coords1[n][0],coords1[n][1],coords1[n][2]] = matrixVecMultiply(U,X[n])
    
  return coords1, E, U

def rotateForJacobi(a,i,j,k,l,s,tau):

  g=a[i][j]
  h=a[k][l]
  a[i][j]=g-s*(h+g*tau)
  a[k][l]=h+s*(g-h*tau)

def eigenJacobi(M):

  # translated from Numerical Recipes in C: 11.1 Jacobi Transformations of a Symmetric Matrix 
  n = len(M)
  if n != len(M[0]):
    return 

  thresh = 0.0
  d      = []
  b      = []
  z      = []

  V = matrixIdentity(n)
    
  for p in range(n):
    b.append(M[p][p])
    d.append(M[p][p])
    z.append(0.0)

  nrot = 0
  for i in range(50):
    sm = 0.0
    for p in range(n-1):
      for q in range(p+1,n):
        sm += abs(M[p][q])

    if sm == 0.0:
      return (d, V, nrot)
    elif i < 3:
      thresh = 0.2 * sm/(n*n)
    else:
      thresh = 0.0
  
    for p in range(n-1):
      for q in range(p+1,n):
        g = 100.0 * abs(M[p][q])
        
        if (i > 3) and ( (abs(d[p])+g) == abs(d[p])) and ( (abs(d[q])+g) == abs(d[q])):
          M[p][q] = 0.0
        elif abs(M[p][q]) > thresh:
          h = d[q] - d[p]
          
          if (abs(h)+g) == abs(h):
            t = M[p][q]/h
          else:
            theta = 0.5*h/M[p][q]
            t = 1.0/(abs(theta)+sqrt(1.0+theta*theta))
            if theta < 0.0:
              t = -t
              
          c = 1.0/sqrt(1+t*t)
          s = t*c
          tau = s/(1.0+c)
          h = t*M[p][q]
          z[p] -= h
          z[q] += h
          d[p] -= h
          d[q] += h
          M[p][q] = 0.0
          for j in range(p):
            rotateForJacobi(M,j,p,j,q,s,tau)
          for j in range(p+1,q):
            rotateForJacobi(M,p,j,j,q,s,tau)
          for j in range(q+1,n):
            rotateForJacobi(M,p,j,q,j,s,tau)
          for j in range(n):
            rotateForJacobi(V,j,p,j,q,s,tau)
          nrot +=1
    for p in range(n):
      b[p] += z[p]
      d[p]  = b[p]
      z[p]  = 0.0
    
  raise 'Too many rotations:', nrot
    
  return None


