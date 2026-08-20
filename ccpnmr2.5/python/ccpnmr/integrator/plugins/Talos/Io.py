
"""
======================COPYRIGHT/LICENSE START==========================

Io.py: code for CCPN data model and code generation framework

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

import operator

atomNameMap = {
 'H':'HN',
}

from ccpnmr.analysis.core import AssignmentBasic
from ccp.util import Assignment as AssignmentUtil
from ccpnmr.integrator.plugins.Talos.TalosShiftFormat import TalosShiftFormat

shiftColumns = ['RESID', 'RESNAME', 'ATOMNAME', 'SHIFT',]
shiftFormats = ['%4d', '%1s', '%4s', '%8.3f']

class TalosFormat:
  """ Adapter class to allow FormatConverter-like interface
  """
  
  def __init__(self, project=None, residues=None):
    self.project = project
    self.IOkeywords = {}
  
  def writeShifts(self, filePath, measurementList, **kw):
    """ Write talos format shift list
    """
    minShiftQuality = self.IOkeywords.get('minShiftQuality', 0.0)
    residues = self.residues
    if not residues:
      residues = list(set(AssignmentUtil.getResonanceResidue(x.reaonance) 
                          for x in measurementList.measurements))
    if len(set(x.chain for x in residues)) != 1:
      raise Exception("Talos shifts must all be from same chain")
    
    residues.sort(key=operator.attrgetter('seqId'))
    stream = open(filePath,'w')
    try:
      if 'atomNames' in self.IOkeywords:
        writeShiftFile(stream, residues, measurementList, minShiftQuality,
                       atomNames=atomNames)
      else:
        writeShiftFile(stream, residues, measurementList, minShiftQuality)
    finally:
      stream.close()
    
def writeShiftFile(stream, residues, shiftList, minShiftQuality=0.0,
                    atomNames = ('H','N','C','CA','CB','HA','HA2','HA3')):
  """ Write Talos type shifts file to stream
  """
  formatObj = TalosShiftFormat()
  
  formatObj.addSequence(residues)
  formatObj.startTable(colNames=shiftColumns, formats=shiftFormats)
  
  for ii,res in enumerate(residues):
    resId = ii + 1
    resName = formatObj.getResidueTag(res)
    usedShifts = set()
    for name in atomNames:
      atom = res.findFirstAtom(name=name)
      if atom is not None:
        atomSet = atom.atomSet
        if atomSet is not None:
          atName = atomNameMap.get(name,name)
          shifts = AssignmentBasic.getAtomSetShifts(atomSet,shiftList)
          if len(shifts) in (1,2):
            for shift in shifts:
              if shift not in usedShifts:
                usedShifts.add(shift)
                if shift.figOfMerit < minShiftQuality:
                  continue
              value = shift.value
              break
            else:
              continue
            
            # We have all we need - write the line
            formatObj.addTableLine((resId, resName, atName, value))
  
  # finished. write it out
  formatObj.endTable()
  formatObj.export(stream)
  
  
