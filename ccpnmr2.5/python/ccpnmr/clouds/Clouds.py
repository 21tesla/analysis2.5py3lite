"""
  These functions control the entire CLOUDS protocol

======================COPYRIGHT/LICENSE START==========================

Clouds.py: Part of the CcpNmr Clouds program

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

from math import sqrt

from FileIO import writePdbCloud, readPdbCloud
from HydrogenDynamics import *
from FilterClouds import alignCloudsToRef, alignToMeanCloud, getMeanPairRmsd, filterClouds, alignClouds
from NoeRelaxation import optimiseRelaxation
from ResonanceIdentification import getCloudsResonanceList

from ccpnmr.analysis.core.StructureBasic import getAtomSetsDistance

def midge(argServer=None):
  
  assert argServer
  
  (resonances,noesyPeaks,intensityFactors) = getCloudsResonanceList(argServer)
  
  # now make the noesy F2 assignments!
  print len(resonances), 'Resonances'
  print len(noesyPeaks), 'Noesy peaks'
  
  constraintList = optimiseRelaxation(resonances,noesyPeaks,intensityMax=36000000,intensityFactors=intensityFactors,tmix=60,sf=500,tcor=3,rleak=2)
  
  structure = argServer.getStructure()
  
  if structure:
    for constraint in constraintList.constraints:
      resonances = list(constraint,findFirstItem().resonances)
    
      atomSets1 = list(resonances[0].resonanceSet.atomSets)
      atomSets2 = list(resonances[1].resonanceSet.atomSets)
      distance = getAtomSetsDistance(atomSets1, atomSets2, structure)
      constraint.setDetails('Known Dist: %4.3f' % (distance))
  
  return constraintList

def hydrogenCloudsDynamics(numClouds,constraintList, resonances):

  
  coolingScheme = []
  # initial temp, final temp, cooling steps, MD steps, MD tau, rep scale
  coolingScheme.append([    1,    1,  3,  500, 0.001, 0])
  coolingScheme.append([80000, 4000, 19, 1000, 0.001, 0])
  coolingScheme.append([ 4000,    1,  5,  500, 0.001, 0])
  coolingScheme.append([15000,    1,  3, 1000, 0.001, 0])
  coolingScheme.append([    1,    1,  5,  500, 0.001, 0])
  coolingScheme.append([ 8000,    1,  3, 1000, 0.001, 0])
  coolingScheme.append([    1,    1,  5,  500, 0.001, 0])
  coolingScheme.append([ 3000,   25, 60, 2500, 0.001, 1])
  coolingScheme.append([   25,   25,  1, 7500, 0.001, 1])
  coolingScheme.append([   10,   10,  1, 7500, 0.001, 1])
  coolingScheme.append([ 0.01, 0.01,  1, 7500,0.0005, 1])
   
  cloudsFiles = generateClouds(numClouds, constraintList, resonances, coolingScheme)
  return cloudsFiles

   
if __name__ == '__main__':
  
  pass
  #cloudsFiles = hydrogenCloudsDynamics(1000,None,None)
  #cloudsFiles = []
  #for i in range(100):
  #  fileName = '/ccpn/clouds/coreCloud/core_cloud_%3.3d.pdb' % (i)
  #  cloudsFiles.append(fileName)
  #
  #alignClouds(cloudsFiles,'/ccpn/clouds/coreCloud/core_align_')
