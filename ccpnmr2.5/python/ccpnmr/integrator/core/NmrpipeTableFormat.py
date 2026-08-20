
"""
======================COPYRIGHT/LICENSE START==========================

NmrpipeTableFormat.py: code for CCPN data model and code generation framework

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

from ccpnmr.integrator.core.TabularFormat import TabularFormat

class NmrpipeTableFormat(TabularFormat):
  """ General Nmrpipe table format. Used also by Talos, Rosetta, PALES, ...
  """
  
  # format settings. 
  # NB should be called as self.xyz, so they can be overridden
  format = 'NmrpipeTable'
  
  dummiesAllowed=False
  startDataLine = 'DATA '
  startComment = None
  startCommentLine = 'REMARK '
  
  seqLineStart = 'DATA SEQUENCE '
  seqGroupsPerLine = 50
  seqSeparator = ''

  def getSequenceList(self, sequence):
    """ Get sequence list 
    Input: a list of MolResidues or Residues in correct order
    """
    
    # Get sequence string
    return [self.getResidueTag(x) for x in sequence]
    

  def addSequence(self, sequence):
    """ Add sequence record, breaking
    Input: a list of MolResidues or Residues in correct order
    """
    sequenceList = self.getSequenceList(sequence)
    
    groupSize = self.seqGroupsPerLine or len(sequenceList)
    for ii in range(0,len(sequenceList), groupSize):
      self.lines.append(self.seqLineStart + 
                     self.seqSeparator.join(sequenceList[ii:ii+groupSize]))
    
  
  def getResidueTag(self, residue):
    """ get residue type idnetifying string
    """
    code1Letter = residue.chemCompVar.chemComp.code1Letter
    if code1Letter is None:
      raise ValueError("Residue %s has no one-letter code" % res)
    #
    return code1Letter
  
  def startTable(self, colNames, formats=None):
    """ Add table start record 
    Input:
    colNames: Column names
    formats (optional) Format expressions, in C/Python '%' syntax 
    """
    self._initTable(colNames, formats)
    self.addLine()
    self.lines.append(self.fieldSep.join(['VARS  '] + colNames))
    self.lines.append(self.fieldSep.join(['FORMAT'] + formats))
  
  
