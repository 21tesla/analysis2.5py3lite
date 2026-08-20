"""
======================COPYRIGHT/LICENSE START==========================

HydrogenDynamics.py: Part of the CcpNmr Clouds program

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
import random, os

from ccpnmr.c.AtomCoordList import AtomCoordList
from ccpnmr.c.DistConstraintList import DistConstraintList
from ccpnmr.c.DistForce import DistForce
from ccpnmr.c.Dynamics import Dynamics

from ccpnmr.clouds.CloudBasic   import writeTypedPdbCloud
from ccpnmr.clouds.FilterClouds import removeDisconnectedAtoms

"""
  These are functions to perform CLOUDS molecular dynamicsEngine
"""

def startMdProcess(numClouds, constraints, resonances, coolingScheme, filePrefix):

  pid = os.fork()
  if pid:
    #print "Child PID:", pid
    return 
  else:
    generateClouds(numClouds, constraints, resonances, coolingScheme, filePrefix)
    os._exit(0)

def generateClouds(numClouds, constraints, resonances, coolingScheme, filePrefix):

  resonances, noeConstrList = makeCloudsDistConstraints(constraints,resonances)
    
  # make follwing params user configureable?
  # - maybe not should be set according to the noe type and whether anti-distance 
  (dynamicsEngine,noeForceField) = makeDynamicsProtocol(noe_force_const=25,noe_exponent=2,noe_soft_exponent=1,
                                                        noe_r_switch=0.5,noe_asymptote=1,beta=10,rmin=2.25,drzap=2)
  
  numAtoms = len(resonances)
  
  print "Number of atoms: %d" % (numAtoms)
  
  for i in range(numClouds):
    pdbFileName = '%s%3.3d.pdb' % (filePrefix,i)
    atomCoordList = makeRandomAtomCoords(numAtoms)
    newAtomCoords = runDynamicsProtocol(atomCoordList, noeConstrList, dynamicsEngine,
                                        noeForceField, coolingScheme)
    
    # newAtomCoords = removeDisconnectedAtoms(newAtomCoords)
    # need violations check here
    writeTypedPdbCloud(newAtomCoords, pdbFileName, resonances)

def makeCloudsDistConstraints(constraints,resonances):

  noeConstrList = DistConstraintList()

  # need to cluster the constraints to make sure there are no bits that will fall off
  cluster = {}
  for constraint in constraints:
  
    r1, r2 = constraint.items[0].resonances
    cluster0 = cluster.get(r1)
    cluster1 = cluster.get(r2)
    
    if (cluster0 is None) and (cluster1 is None):
      cluster[r1] = [r1,r2]
      cluster[r2] = cluster[r1]
    
    elif cluster0 is None:
      cluster1.append(r1)
      cluster[r1] = cluster1
      
    elif cluster1 is None:
      cluster0.append(r2)
      cluster[r2] = cluster0   

    elif cluster1 != cluster0:
      cluster0.extend(cluster1)
      for r3 in cluster1:
        cluster[r3] = cluster0

  clusters = cluster.values()
  core = clusters[0]
  for resonances0 in clusters:
    if len(resonances0) > len(core):
      core = resonances0

  coreDict = {}
  for resonance in core:
    coreDict[resonance] = 1

  constraints0 = []
  for constraint in constraints:
    r1, r2 = constraint.items[0].resonances
    if coreDict.get(r1) is not None:
      constraints0.append(constraint)
  
  resDict = {}
  for resonance in resonances:
    resDict[resonance.serial] = resonance

  dict = {}
  resonances0 = []
  for i in range(len(core)):
    dict[core[i].resonanceSerial] = i
    coreResonance = resDict.get(core[i].resonanceSerial)
    if coreResonance:
      resonances0.append(coreResonance)
    else:
      print "Core resonance missing", core[i].resonanceSerial, core[i]
    
    
  for constraint in constraints0:
    item = constraint.items[0]
    if dict.get(item.resonances[0].resonanceSerial) is None:
      print 'Warning: Missing resonance serial %d' % item.resonances[0].resonanceSerial
      continue 

    if dict.get(item.resonances[1].resonanceSerial) is None:
      print 'Warning: Missing resonance serial %d' % item.resonances[1].resonanceSerial
      continue 
    
    atom0 = dict[item.resonances[0].resonanceSerial]
    atom1 = dict[item.resonances[1].resonanceSerial]
    noeConstrList.add(atom0, atom1, constraint.lowerLimit, constraint.upperLimit)

  print 'Total  C:%d  R:%d' % (len(constraints),  len(resonances) )
  print 'Core   C:%d  R:%d' % (len(constraints0), len(resonances0))
    
  return resonances0, noeConstrList

def makeRandomAtomCoords(natoms):

  atomCoordList = AtomCoordList()
  hmass = 25
  r = 30
  for n in range(natoms):
    x = r * (2*random.uniform(0,1) - 1)
    y = r * (2*random.uniform(0,1) - 1)
    z = r * (2*random.uniform(0,1) - 1)
    atomCoordList.add(hmass, x, y, z)

  return atomCoordList

def makeDynamicsProtocol(noe_force_const=25,noe_exponent=2 ,noe_soft_exponent=1,noe_r_switch=0.5,
                         noe_asymptote=1,beta=10,rmin=2.25,drzap=2,nprint = 3000):

  noeForceField  = DistForce(noe_force_const,noe_exponent,noe_soft_exponent,noe_r_switch,noe_asymptote)
  dynamicsEngine = Dynamics(beta=beta, rmin=rmin, drzap=drzap, elapsed_time=0, nprint=nprint)
  
  return (dynamicsEngine,noeForceField)

def runDynamicsProtocol(atomCoordList, noeConstrList, dynamicsEngine,
                        noeForceField, coolingScheme, rp_force_const=1):

  for (i,temp_i,temp_f,ncooling,nsteps,tau,rp_scale) in coolingScheme:

    tref = temp_i
    rpf  = rp_force_const * rp_scale

    if (ncooling > 1):
      dtemp = (temp_f - temp_i) / float(ncooling - 1)
    else:
      dtemp = 0

    dynamicsEngine.rp_force_const = rpf
    dynamicsEngine.tref           = tref
    dynamicsEngine.tau            = tau
    dynamicsEngine.nsteps         = nsteps

    for j in range(ncooling):
      dynamicsEngine.run(atomCoordList, noeConstrList, noeForceField)
      dynamicsEngine.tref = dynamicsEngine.tref + dtemp

  return atomCoordList


