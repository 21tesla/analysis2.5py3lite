
"""
======================COPYRIGHT/LICENSE START==========================

ViewResidueFrame.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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

===========================REFERENCE END===============================
"""

import Tkinter, math

from memops.gui.CheckButton          import CheckButton
from memops.gui.Frame                import Frame
from memops.gui.Label                import Label
from ccp.gui.ViewChemCompVarFrame    import ViewChemCompVarFrame
from ccpnmr.analysis.core.AssignmentBasic import getAtomSetShifts
from ccpnmr.analysis.core.MoleculeBasic   import getResidueCode


#colorDict = {'C':'#404040','N':'#0000A0','O':'#A00000','H':'#A0A0A0','P':'#00A000','S':'#A0A000'}

class ViewResidueFrame(Frame):

  showAssign   = True
  chemCompVar  = None
  residue      = None

  def __init__(self, parent, residue=None, resizeCallback=None,
               project=None, tipText=None, shiftList=None, *args, **kw):
    
    self.shiftList = shiftList

    Frame.__init__(self, parent, *args, **kw)
 
    self.grid_columnconfigure(2, weight=1)
 
    row = 0
    self.label = Label(self, text='', grid=(row,0))
    
    self.assignSelect = CheckButton(self, callback=self.setDisplayAssign, grid=(row,1),
                                    tipText='Whether to show chemical shifts of assigned atoms')
    
    label0 = Label(self, text='Show Assignments', grid=(row,2))
    
    row +=1 
    self.grid_rowconfigure(row, weight=1)
    self.varFrame = ViewChemCompVarFrame(self, chemCompVar=self.chemCompVar,
                                         project=project, tipText=tipText,
                                         grid=(row,0), gridSpan=(1,3))

    self.assignSelect.set(True)
    

  def setDisplayAssign(self, trueOrFalse):
    
    self.showAssign = trueOrFalse

    self.update(self.residue)

  def setResidue(self, residue):

    self.residue = residue
    shiftList = self.shiftList
    
    if residue:
      self.varFrame.update(self.residue.chemCompVar)
            
      cAtomDict = self.varFrame.cAtomDict
      for atom in self.residue.atoms:
        chemAtom = atom.chemAtom
        cAtom = cAtomDict.get(chemAtom)
          
        if cAtom and self.showAssign:
          shifts = getAtomSetShifts(atom.atomSet, shiftList)
          label = '/'.join(['%3.3f' % (shift.value) for shift in shifts])
          cAtom.setAnnotation(chemAtom.name + ' ' + label)
      
      self.varFrame.drawStructure()     
      chain = self.residue.chain
      self.label.set('Residue: %d%s ( %s %s )' % (self.residue.seqCode,getResidueCode(self.residue),chain.molSystem.code,chain.code))
        
    else:
      self.varFrame.update(None)
      self.label.set('Residue: <None>')

  def update(self, residue=None, shiftList=None):
  
    self.shiftList = shiftList
    self.setResidue(residue)
 
  def printStructure(self):

    self.varFrame.printStructure()
