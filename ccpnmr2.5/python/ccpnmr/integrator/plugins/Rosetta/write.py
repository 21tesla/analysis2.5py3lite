
"""
======================COPYRIGHT/LICENSE START==========================

write.py: code for CCPN data model and code generation framework

Copyright (C) 2011  (CCPN Project)

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

pluginName = 'ccpnmr.integrator.plugins.Rosetta'

#from ccpnmr.integrator.plugins.Talos import Io as talosIo


def write(nmrCalcRun, targetDir):
  """ Write input files for Program run
    Input:
      nmrCalcRun: NmrCalc.Run 
      targetDir: destination directory.
  """
  
  intIo.writeDataFiles(nmrCalcRun, targetDir)
  
  jsonDict = intIo.makeJsonDict(nmrCalcRun)
  
  
  # write properties file (must be done at the end
  propFile = uniIo.joinPath(targetDir, intIo.propFileName)
  print 'About to write', propFile
  open(propFile,'w').write(json.dumps(jsonDict, sort_keys=True, 
                                      indent=intIo.propIndent))



##########################################################################
#
# Old code from here on
#

def oldwrite(nmrCalcRun, topDir=None):
  """ Write input files for Program run
    Input:
      nmrCalcRun: NmrCalc.Run 
      topDir: optional destination directory.
  """
  targetDir = intIo.createTargetDir(nmrCalcRun, topDir=topDir)
  
  # write Properties file
  propDict, dummy = intIo.makeParameterDict(nmrCalcRun)
  propFile = uniIo.joinPath(targetDir, intIo.propFileName)
  open(propFile,'w').write(json.dumps(propDict, sort_keys=True, 
                                      indent=intIo.propIndent))
  
  # Write data file
  
  # get file
  inpFile = propDict.get('INPUT_FILE')
  if inpFile:
    talosInputFile = uniIo.joinPath(targetDir, inpFile)
  else:
    raise Exception("No code=INPUT_FILE RunParameter in nmrCalcRun")
  
  # Get shiftList
  shiftList = intUtil.getMeasurementList(nmrCalcRun)
  if shiftList is None:
    raise Exception("Run must have exactly one shift list, %s found" % len(ll))
  
  # Get residues
  ll = list(nmrCalcRun.findAllData(className='MolResidueData'))
  if len(ll) == 1:
    obj =  ll.pop()
    residues = list(obj.residues)
    residues.sort(key=operator.attrgetter('seqId'))
    if len(set(obj.chainCodes)) != 1:
      raise Exception("Run MolResidueData did not have a (single) chain" )
  else:
    raise Exception("Run must have excactly one MolResidueData, %s found" % len(ll))
  
  #  
  talosIo.writeShiftFile(open(talosInputFile,'w'), residues, shiftList,
                         propDict.get('minShiftQuality'))
  


if __name__ == '__main__':
  """ Run write function from command line.
  Input is projectDir NmrCalcRun.IDstring
  projectDir must contain the desired project (and no others)
  NmrCalcRun.IDstring is of the form 
  '%s+%s' % (NmrCalcStore.guid, Run.serial)
  """
 
  if len(sys.argv) >= 3:
 
    # set up input
    junk, projectDir, nmrCalcRunId = sys.argv[:3]
  
    if len(sys.argv) >= 4:
      targetDir = sys.argv[3]
    else:
      targetDir=None
    
    #intIo.prepareStdWmsRun(pluginName, projectDir, nmrCalcRunId, targetDir)
    
    # NB necessary. 
    # As side effect prints message that passes newCalcId to WMS Java
    nmrCalcRun = intIo.getNmrCalcRun(projectDir, nmrCalcRunId, pluginName)
 
  else:
    print "Usage: write projectDir NmrCalcRun.IDstring targetDir(optional)"

