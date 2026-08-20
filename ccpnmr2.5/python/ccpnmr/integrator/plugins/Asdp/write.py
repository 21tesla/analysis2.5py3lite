
"""
======================COPYRIGHT/LICENSE START==========================

write.py: code for CCPN data model and code generation framework

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
#from ccpnmr.integrator.core import Util as intUtil
#from ccpnmr.integrator.plugins.Unio import Util as unioUtil
from memops.universal import Io as uniIo

pluginName = 'ccpnmr.integrator.plugins.Asdp'

programBin = '${ASBIN}' 

def write(nmrCalcRun, targetDir):
  """ Write input files for Program run
    Input:
      nmrCalcRun: NmrCalc.Run 
      targetDir: destination directory.
  """
  
  intIo.writeDataFiles(nmrCalcRun, targetDir)
  
  jsonDict = intIo.makeJsonDict(nmrCalcRun)
  
  # write properties file 
  propFile = uniIo.joinPath(targetDir, intIo.propFileName)
  open(propFile,'w').write(json.dumps(jsonDict, sort_keys=True, 
                                      indent=intIo.propIndent))
  
  # Write program setup file
  fileNameObj = nmrCalcRun.findFirstRunParameter(name='fileNameSetup')
  if fileNameObj is not None:
    filePath = uniIo.joinPath(targetDir, fileNameObj.textValue)
    writeSetupFile(filePath, jsonDict)
  
def writeSetupFile(filePath, jsonDict):
  """ Write ASDP-style setup file with parameters in jsonDict
  """
  
  lineFormat = "%s=%s\n"
  headingFormat = "[%s]\n"
  
  fp = open(filePath, 'w')
  
  runParameters = jsonDict['RunParameter']
  
  try:
  
    fp.write(headingFormat % 'General')
    for tag in ('ACO', 'nCycles', 'proteinName', 'center'):
      val = runParameters.get(tag)
      if val is not None:
        fp.write(lineFormat % (tag,val))
    shiftListDict = jsonDict['MeasurementListData'][0]
    tag = 'chemicalShiftFile'
    fp.write(lineFormat % (tag,shiftListDict.get(tag)))
    
    
    calcHeader = jsonDict.get('calcHeader')
    if calcHeader:
      fp.write(headingFormat % calcHeader)
      ll = [jsonDict['calcCommand']]
      for tag in ('np', 'nb', 'que', 'av', 'nstr', 'hotcy', 'hotad'):
        ll.append('-'+tag)
        ll.append(str(runParameters[tag]))
      fp.write(' '.join(ll) + '\n')
    
    refineHeader = jsonDict.get('refineHeader')
    if refineHeader:
      fp.write(headingFormat % refineHeader)
      ll = [jsonDict['refineCommand']]
      for tag in ('np', 'nb'):
        ll.append('-'+tag)
        ll.append(str(runParameters[tag]))
      fp.write(' '.join(ll) + '\n')
        
    
    for dd in jsonDict['PeakListData']:
      fp.write('\n')
      fileName = dd.pop('fileNamePeak')
      fp.write(headingFormat % fileName)
      for tag, val in sorted(dd.items()):
        if tag != 'serial':
          fp.write(lineFormat % (tag, val))

  finally:
    fp.close()


def prepareLocalExecution(nmrCalcRun, targetDir):
  """ return [procArgs, targetDir, logFile list for later, local execution
  And carry out any preliminary commands
  """

  # Set up parameters for program call
  shellParObj = nmrCalcRun.findFirstRunParameter(name='programCall')
  if shellParObj is None:
    raise Exception("Parameter name=programCall not found")
  else:
    shellCall = shellParObj.textValue
    
  executeScript = os.path.join(os.path.expandvars(programBin), shellCall)
  
  fileNameObj = nmrCalcRun.findFirstRunParameter(name='fileNameSetup')
  if fileNameObj is None:
    raise Exception("no 'fileNameSetup' parameter found")
  filePath = uniIo.joinPath(targetDir, fileNameObj.textValue)
  
  runDir =  os.path.join(targetDir,'asdpRun')
  procargs = [executeScript, '-c', filePath, '-o', runDir]
  
  #
  return (procargs, targetDir)

if __name__ == '__main__':
  """ Run write function from command line.
  Input is projectDir NmrCalcRun.IDstring
  projectDir must contain the desired project (and no others)
  NmrCalcRun.IDstring is of the form 
  '%s+%s' % (NmrCalcStore.guid, Run.serial)
  """
  
  from memops.general import Io as genIo
  
  if len(sys.argv) == 3:
    
    # set up input
    junk, projectDir, nmrCalcRunId = sys.argv
    
    # NB necessary. 
    # As side effect prints message that passes newCalcId to WMS Java
    nmrCalcRun = intIo.getNmrCalcRun(projectDir, nmrCalcRunId, pluginName)
      
  else:
    print "Usage: write projectDir NmrCalcRun.IDstring"
  
  
