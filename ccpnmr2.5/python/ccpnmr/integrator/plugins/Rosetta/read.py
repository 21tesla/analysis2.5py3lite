
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
from ccpnmr.integrator.core import Util as intUtil
from memops.universal import Io as uniIo
from ccp.lib import StructureIo

from ccpnmr.integrator.plugins.Talos import Io as talosIo


##########################################################################
#
# Old code from here on
#

stdFileNames = {
 'scores' : 'name.rescore.txt',
 'rawRms' : 'rms2LowRawScore.txt',
 'rms' : 'rms2LowReScore.txt',
 'input' : 'talos.tab',
}


def read(nmrCalcRun, dataDir):
  """ Read Output files for Rosetta run
    Input:
      nmrCalcRun: NmrCalc.Run 
      dataDir: directory directly containing program output.
  """
  
  # set up
  files = {}
  data = {}
  
  fileNames = stdFileNames.copy()
  #fileNames['properties'] = intIo.propFileName
  
  # check presence of files:
  if not os.path.isdir(dataDir):
    raise IOError("Data directory %s not found" % dataDir)
  
  for tag in ('scores', 'rms', 'rawRms', ):
    ss = uniIo.joinPath(dataDir, fileNames[tag])
    files[tag] = ss
    if not os.path.isfile(ss):
      raise IOError("%s file %s not found in directory %s" % (tag, ss, dataDir))
  
  # check NmrCalcRun:
  if not intUtil.runIsWritable(nmrCalcRun):
    raise Exception("Cannot read into %s - has already been used"
                    % nmrCalcRun)
  
  # get simple data records
  #ss = fileNames.get('input')
  #if ss:
  #  path = uniIo.joinPath(dataDir, ss)
  #  if os.path.isfile(path):
  #    data = open(path).read()
  #    nmrCalcRun.newRunParameter(name='inputUsed', code=stdFileNames['input'],
  #                               ioRole='output', textValue=data,
  #                               )
                                 #details='Input data, as used by Rosetta')
  
  
  #NBNB TBD remove. Not there. Properties is an input thing.
  #ss = fileNames.get('properties')
  #if ss:
  #  path = uniIo.joinPath(dataDir, ss)
  #  if os.path.isfile(path):
  #    data = open(ss).read()
  #    nmrCalcRun.newRunParameter(name='propertiesUsed', 
  #                               code=stdFileNames['properties'],
  #                               ioRole='output', textValue=data,
  #                               details='Properties file, as used by Rosetta')
  
  # read data arrays and sort by Rosetta score
  scoreFile = open(files['scores'])
  rmsFile = open(files['rms'])
  rawRmsFile = open(files['rawRms'])
  data = []
  for ss in scoreFile:
    rms = float(rmsFile.next().split()[1])
    rawRms = float(rawRmsFile.next().split(None,1)[0])
    tt = ss.split()
    # Order is: score, rawscore shiftChi2, rms, rawrms, name
    data.append((float(tt[3]), float(tt[1]), float(tt[2]), rms, rawRms, tt[0]))
  data.sort()
  scores, rawScores, shiftChi2, rms, rawRms, names = zip(*data)
  
  # Set scores data in NmrCalc matrices
  docTemplate = (
   "Rosetta %s for all calculated structures, sorted by corrected score"
  )
  ss = 'corrected scores'
  nmrCalcRun.newFloatMatrixData(code='scores', name = ss, ioRole='output',
                                shape=(0,), data=scores, 
                                details=docTemplate % ss)
  ss = 'raw scores'
  nmrCalcRun.newFloatMatrixData(code='rawScores', name = ss, ioRole='output',  
                                shape=(0,), data=rawScores, 
                                details=docTemplate % ss)
  ss = 'shift chi2 scores'
  nmrCalcRun.newFloatMatrixData(code='shiftChi2', name = ss, ioRole='output',  
                                shape=(0,), data=shiftChi2, 
                                details=docTemplate 
                                 % 'Chemical' + ss)
  ss = 'RMS to best struct'
  nmrCalcRun.newFloatMatrixData(code='rms', name = ss, ioRole='output',  
                                shape=(0,), data=rms, 
                                details=docTemplate % 
                                 'RMS to structure with best corrected score')
  ss = 'RMS to best raw struct'
  nmrCalcRun.newFloatMatrixData(code='rawRms', name = ss, ioRole='output',  
                                shape=(0,), data=rawRms, 
                                details=docTemplate % 
                                 'RMS to structure with best raw score')
  
      
  # Get MolSystem
  molResidueData = nmrCalcRun.findFirstData(className='MolResidueData', 
                                            ioRole='input')
  molSystem = molResidueData.molSystem
  
  
  # read structures in order of lowest corrected score
  paths = []
  for modelNum,fileName in enumerate(names): 
    path = uniIo.joinPath(dataDir, fileName + '.pdb')
    if os.path.isfile(path):
      paths.append(path)
    else:
      break
  
  ensemble = StructureIo.getStructureFromFiles(molSystem, paths,
                                               fileType='rough')
  
  ss =  "Rosetta calculated structures, sorted by corrected score"
  nmrCalcRun.newStructureEnsembleData(name='result', ioRole='output',
                                structureEnsemble=ensemble,
                                details=docTemplate % ss)
  
  # Make Nmr.StructureCalculation
  nmrProject = nmrCalcRun.nmrCalcStore.nmrProject
  nmrCalcRun.structureGeneration = nmrProject.newStructureGeneration(
      generationType='denovo', name='CSRosetta', structureEnsemble=ensemble)
  
  

if __name__ == '__main__':
  """ Run read function from command line.
  Input is projectDir NmrCalcRun.IDstring, Rosetta.outputdir
  projectDir must contain the desired project (and no others)
  NmrCalcRun.IDstring is of the form 
  '%s+%s' % (NmrCalcStore.guid, Run.serial)
  """
  
  if len(sys.argv) == 4:
    
    # set up input
    junk, projectDir, nmrCalcRunId, rosettaDir = sys.argv
    
    nmrCalcRun = intIo.getNmrCalcRun(projectDir, nmrCalcRunId)
    if nmrCalcRun is None:
      print "No NmrCalcRun found. Aborting"
    else:
      read(nmrCalcRun, rosettaDir)
      nmrCalcRun.root.saveModified()
    
  else:
    print "Usage: read projectDir NmrCalcRun.IDstring RosettaDataDir"
  
  
