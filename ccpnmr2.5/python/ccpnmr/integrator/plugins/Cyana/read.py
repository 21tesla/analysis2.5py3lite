
"""
======================COPYRIGHT/LICENSE START==========================

read.py: code for CCPN data model and code generation framework

Copyright (C) 2012  (CCPN Project)

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
import sys, os, operator, json

from ccpnmr.integrator.core import Io as intIo
from memops.universal import Io as uniIo
from ccp.lib import StructureIo

"""
Cyana .tgz expands in top - there is no top directory.

Overview file: final.ovw
cyanatable.txt per-cycle data.

final.pdb. Structures:
No pseudoatoms

final.upl. Singly assigned restraints. 
All lists mixed - have peak no but not parent list
INcludes pseudoatoms

Peaklists:
For each list input.xpk and ...-cycle7.xpk.
  NBNB cycle7 ones are EMPTY!.
 Skip for now.

ShiftLists:
Starting and final. Identical except where prochirals are switched.
Includes pseudoatoms
Skip

finalstereo.cya. List of stereospecific assignments

log file (std.out) and error file (std.err) 

cycle,aco and final.aco.
Looks like initial and final backbone dihedral constraints.
Looks like default (in my case), differe in interval widths.
Skip

So:
1) read final.ovw and cyanatable.txt into output:Summary
   read final.pdb in structure ensemble.

2) consider .upl file, and shiftlist / stereospecifics handling

3) check what is wrong with returned peak lists.

 
"""

stdFileNames = {
 'pdb' : 'final.pdb',
 'overview1' : 'final.ovw',
 'overview2' : 'cyanatable.txt',
 #'log' : 'std.out',
}

def read(nmrCalcRun, dataDir):
  """ Read Output files for Cyana run
    Input:
      nmrCalcRun: NmrCalc.Run 
      dataDir: directory directly containing program output.
  """
  
  fileNames = stdFileNames.copy()
  
  # make summary text element
  summary = ""
  for tag in ('overview1', 'overview2'):
    path = uniIo.joinPath(dataDir, fileNames[tag])
    if os.path.isfile(path):
      txt = open(path).read()
    else:
      txt = 'WARNING, file %s not found' % fileNames[tag]
    summary += ('-'*80 + "\nFROM FILE :      %s\n----------------\n%s\n\n\n"  
                % (fileNames[tag], txt))
  
  nmrCalcRun.newRunParameter(name='summary', ioRole='output', 
                             textValue=summary)
  
  
  # Get MolSystem
  molResidueData = (nmrCalcRun.findFirstData(className='MolResidueData', 
                                             ioRole='input')
                    or
                    nmrCalcRun.findFirstData(className='MolSystemData', 
                                             ioRole='input')
                   )
  molSystem = molResidueData.molSystem
  
  # read pdb structures
  path = uniIo.joinPath(dataDir, fileNames['pdb'])
  if os.path.isfile(path):
    ensemble = StructureIo.getStructureFromFile(molSystem, path)
  
  # make EnsembleData
  nmrCalcRun.newStructureEnsembleData(name='result', ioRole='output',
                                structureEnsemble=ensemble,
                                details="Cyana calculated structures")

  
  # Make Nmr.StructureCalculation
  nmrProject = nmrCalcRun.nmrCalcStore.nmrProject
  nmrCalcRun.structureGeneration = nmrProject.newStructureGeneration(
      generationType='denovo', name='CYANA', structureEnsemble=ensemble)
  
#

if __name__ == '__main__':
  """ Run read function from command line.
  Input is projectDir NmrCalcRun.IDstring, Rosetta.outputdir
  projectDir must contain the desired project (and no others)
  NmrCalcRun.IDstring is of the form 
  '%s+%s' % (NmrCalcStore.guid, Run.serial)
  """
  
  if len(sys.argv) == 4:
    
    # set up input
    junk, projectDir, nmrCalcRunId, dataDir = sys.argv
    
    nmrCalcRun = intIo.getNmrCalcRun(projectDir, nmrCalcRunId)
    if nmrCalcRun is None:
      print "No NmrCalcRun found. Aborting"
    else:
      read(nmrCalcRun, dataDir)
      nmrCalcRun.root.saveModified()
    
  else:
    print "Usage: read projectDir NmrCalcRun.IDstring CyanaDataDir"
  
  
