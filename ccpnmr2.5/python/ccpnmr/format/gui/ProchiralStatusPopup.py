
"""
======================COPYRIGHT/LICENSE START==========================

ProchiralStatusPopup.py: GUI window for setting the prochiral status of resonances

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

- contact Wim Vranken (wim@ebi.ac.uk)
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
import Tkinter

import string

from memops.universal.Io import joinPath

from memops.universal.Util import returnStrings
from memops.universal.Util import returnInt

from memops.gui.Label import Label
from memops.gui.CheckButton import CheckButton
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledFrame import ScrolledFrame
from memops.gui.Separator import Separator

from ccp.api.nmr import Nmr, NmrConstraint

from ccpnmr.format.gui.CrossLine import CrossLine
from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

from ccpnmr.format.general.Io import getHelpUrlDir

class ProchiralStatusPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ProchiralStatus.html')

  def __init__(self, parent, chainList, chainAtomSetsDict):

    self.chainList = chainList
    self.chainAtomSetsDict = chainAtomSetsDict
    
    self.resParentDict = {}
    self.resParentList = []
    
    self.resParentChainLabelDict = {}
    self.resParentChainLabelList = {}
    
    self.chainResonancesDict = {}
    
    for resParent in self.chainAtomSetsDict.keys():
      if isinstance(resParent,NmrConstraint.NmrConstraintStore):
        strucGens = resParent.structureGenerations
        strucGenSerials = [str(sg.serial) for sg in strucGens]
        strucGenLabel = string.join(strucGenSerials,',')
        resParentLabel = "Constraint set %d (strucGens %s)" % (resParent.serial,strucGenSerials)
      else:
        resParentLabel = "Nmr project %s" % resParent.name
        
      self.resParentDict[resParentLabel] = resParent
      self.resParentList.append(resParentLabel)
      
      self.chainResonancesDict[resParent] = {}
    
      self.resParentChainLabelDict[resParent] = {}
      self.resParentChainLabelList[resParent] = []

      for chain in self.chainList:
        self.project = chain.root
        if self.chainAtomSetsDict[resParent].has_key(chain) and self.chainAtomSetsDict[resParent][chain]:
          chainLabel = "Molecular system '%s':Chain '%s' (molecule '%s')" % (chain.molSystem.code,chain.code,chain.molecule.name)
          self.resParentChainLabelList[resParent].append(chainLabel)
          self.resParentChainLabelDict[resParent][chainLabel] = chain
      
          self.chainResonancesDict[resParent][chain] = []
      
    # modal = true means that it won't continue unless this one returns value
    TemporaryBasePopup.__init__(self, parent=parent, title="Project '%s': " % self.project.name + 'Set status prochiral atoms/resonances', modal=False, transient=True)

  def body(self, master):

    #
    # Setup header 
    #

    row = 0

    label = Label(master, text = 'Select the resonance list:')
    label.grid(row=row, column=0, sticky=Tkinter.EW)

    label = Label(master, text = 'Select the chain:')
    label.grid(row=row, column=1, sticky=Tkinter.EW)

    row += 1

    self.resParentSelect = PulldownMenu(master, entries = self.resParentList, callback = self.setupResParentList, do_initial_callback = False)
    self.resParentSelect.grid(row=row, column=0, sticky=Tkinter.EW)

    self.chainSelect = PulldownMenu(master, entries = [], callback = self.setupProchiralList, do_initial_callback = False)
    self.chainSelect.grid(row=row, column=1, sticky=Tkinter.EW)

    row += 1

    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, columnspan = 2, sticky = Tkinter.EW)

    row += 1

    # THIS BIT TELLS MASTER TO CONFIGURE WINDOW ON INSIDE WIDGET!!
    master.grid_rowconfigure(row,weight = 1)
    for i in range(2):
      master.grid_columnconfigure(i,weight = 1)

    self.prochiralListFrame = ScrolledFrame(master, width = 200, doExtraConfig = False)
    self.prochiralListFrame.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.NSEW)
    self.prochiralListFrameRow = row
    self.prochiralListFrameMaster = master
    
    row = row + 1

    texts = [ 'Set all', 'Set chain']
    commands = [ self.ok, self.setStereo ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = 'Exit', help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.EW)

    self.resParentSelect.callback(0,self.resParentList[0])

  def setupResParentList(self,listIndex,resParentLabel):
  
    resParent = self.resParentDict[resParentLabel]
    self.resParent = resParent
     
    self.chainSelect.setup(self.resParentChainLabelList[self.resParent], 0)
    self.chainSelect.callback(0,self.resParentChainLabelList[self.resParent][0])

  def setupProchiralList(self,listIndex,chainLabel):
    
    chain = self.resParentChainLabelDict[self.resParent][chainLabel]
    self.chain = chain
    
    self.prochiralListFrame.destroy()
    self.prochiralListFrame = ScrolledFrame(self.prochiralListFrameMaster, width = 200, doExtraConfig = False)
    self.prochiralListFrame.grid(row=self.prochiralListFrameRow, column=0, columnspan = 2, sticky=Tkinter.NSEW)
    
    #
    # TODO: also handle case one resonance, two atomSets!?!?
    #
    
    frameRow = 0
    x = y = 0
    
    frame = self.prochiralListFrame.frame
    
    self.resonanceObjects = {}
    self.resonanceLabelDict = {}
    
    #
    # Reuse previous settings...
    #

    if len(self.chainResonancesDict[self.resParent][chain]) == len(self.chainAtomSetsDict[self.resParent][chain]):
      usePreviousSettings = True
    else:
      usePreviousSettings = False
    
    for i in range(0,len(self.chainAtomSetsDict[self.resParent][chain])):
      atomSetGroup = self.chainAtomSetsDict[self.resParent][chain][i]
      
      isStereo = False
      hasResonances = False
      isValidCase = False
      
      resonances = []
      
      if usePreviousSettings and self.chainResonancesDict[self.resParent][chain][i]:
      
        resonances = self.chainResonancesDict[self.resParent][chain][i][3]
        isStereo = self.chainResonancesDict[self.resParent][chain][i][0]
        hasResonances = True
        isValidCase = True
               
      else:
      
        for atomSet in atomSetGroup:
          resonanceSets = atomSet.sortedResonanceSets()

          if resonanceSets:
            hasResonances = True

            for resonanceSet in resonanceSets:
              for resonance in resonanceSet.sortedResonances():
                if resonance not in resonances:
                  resonances.append(resonance)

            if len(resonanceSets) == 1:
              isValidCase = True
              if len(resonanceSet.atomSets) == 1:
                isStereo = True

          else:
            # Case only one resonance but stereo assigned
            resonances.append(0)
        
        if len(resonances) == 1:
          # Case only one resonance, not stereo assigned
          resonances.append(None)

      #
      # Set up only for valid cases
      #

      if hasResonances and isValidCase:
      
        if not usePreviousSettings:
          self.chainResonancesDict[self.resParent][chain].append([isStereo,False,False,resonances])
      
        residue = atomSetGroup[0].findFirstAtom().residue
        ccpCode = residue.molResidue.ccpCode
        seqCode = residue.seqCode

        rowspan = 2

        checkButton = CheckButton(frame, selected = isStereo, callback = lambda status = isStereo, fr = frameRow, item = i, setup = False: self.toggleResonances(status,fr,item,setup))
        checkButton.grid(row = frameRow, rowspan = rowspan, column = 0, sticky = Tkinter.EW)
        
        label = Label(frame, text = ccpCode, width = 5)
        label.grid(row = frameRow, rowspan = rowspan, column = 1, sticky = Tkinter.E)
        label = Label(frame, text = str(seqCode), width = 5)
        label.grid(row = frameRow, rowspan = rowspan, column = 2, sticky = Tkinter.W)

        #
        # Set up the atom names
        #

        for j in range(0,len(atomSetGroup)):
          atomSet = atomSetGroup[j]
          label = Label(frame, text = atomSet.name)
          label.grid(row = frameRow + j, column = 3, sticky = Tkinter.E, ipadx = 4)
                  
        #
        # Set up the resonances
        #
        
        self.resonanceObjects[frameRow] = {True: [[],[]],
                                           False: [[],[]]}

        resonanceLabelList = []

        for resonance in resonances:
        
          if not resonance:
          
            resonanceLabel = 'unknown'
            
          else:

            if hasattr(resonance,'shifts') and resonance.shifts:
              shiftLabel = " (shifts: "
            
              for shift in resonance.sortedShifts():
                shiftLabel += "%.2f, " % shift.value
              
              shiftLabel = shiftLabel[:-2] + ")"
            
            else:
              shiftLabel = ""
            
            resonanceLabel = "'%d.%s'%s" % (seqCode,resonance.name,shiftLabel)
            self.resonanceLabelDict[resonanceLabel] = resonance        
  
          resonanceLabelList.append(resonanceLabel)

        separator = Separator(frame, height = 1, width = 30)
        separator.setColor('black', bgColor = 'black')
        separator.grid(row = frameRow, column = 4, sticky = Tkinter.EW)
        self.resonanceObjects[frameRow][True][0].append(separator)

        separator = CrossLine(frame, height = 20, width = 30, canvas_bg = frame['bg'], color = 'dimgray')
        self.resonanceObjects[frameRow][False][0].append(separator)

        resonanceList = PulldownMenu(frame,
                                     entries = resonanceLabelList,
                                     callback = lambda selInd = x, selItem = y, fr = frameRow, item = i: self.setProchiralLabels(selInd,selItem,fr,item),
                                     do_initial_callback = False, label_color = 'red')
        resonanceList.grid(row=frameRow, column = 5, sticky=Tkinter.EW)
        self.resonanceObjects[frameRow][True][0].append(resonanceList)

        resonanceLabel = Label(frame,text = "%s" % resonanceLabelList[0], fg = 'dimgray')
        self.resonanceObjects[frameRow][False][0].append(resonanceLabel)

        separator = Separator(frame, height = 1, width = 30)
        separator.setColor('black', bgColor = 'black')
        separator.grid(row = frameRow + 1, column = 4)
        self.resonanceObjects[frameRow][True][1].append(separator)

        self.resonanceObjects[frameRow][False][1].append(None)

        resonanceLabel = Label(frame,text = "%s" % resonanceLabelList[1], fg = 'red')
        resonanceLabel.grid(row=frameRow + 1, column = 5, sticky=Tkinter.EW)
        self.resonanceObjects[frameRow][True][1].append(resonanceLabel)

        resonanceLabel = Label(frame,text = "%s" % resonanceLabelList[1], fg = 'dimgray')
        self.resonanceObjects[frameRow][False][1].append(resonanceLabel)

        if not isStereo:
          checkButton.callback(isStereo,setup = True)    

        separator = Separator(frame, height = 1)
        separator.setColor('black', bgColor = 'black')
        separator.grid(row = frameRow + rowspan, columnspan = 6, sticky = Tkinter.EW)

        frameRow += rowspan + 1

      elif not usePreviousSettings:
          self.chainResonancesDict[self.resParent][chain].append(None)
     
    return True

  def setProchiralLabels(self,selInd,selItem,fr,item):

    resonanceLabel = self.resonanceObjects[fr][True][0][1].entries[selInd]
    self.resonanceObjects[fr][False][0][1].set(resonanceLabel)
  
    selInd = not selInd
    resonanceLabel = self.resonanceObjects[fr][True][0][1].entries[selInd]
    self.resonanceObjects[fr][True][1][1].set(resonanceLabel)
    self.resonanceObjects[fr][False][1][1].set(resonanceLabel)
    
    self.chainResonancesDict[self.resParent][self.chain][item][2] = not self.chainResonancesDict[self.resParent][self.chain][item][2]
    self.chainResonancesDict[self.resParent][self.chain][item][3].reverse()

  def toggleResonances(self, status, frameRow, item, setup):
  
    if not setup:
      self.chainResonancesDict[self.resParent][self.chain][item][0] = status
      self.chainResonancesDict[self.resParent][self.chain][item][1] = not self.chainResonancesDict[self.resParent][self.chain][item][1]

    for i in range(0,len(self.resonanceObjects[frameRow])):
      
      #
      # Remove other widgets from grid
      #
      
      rpds = self.resonanceObjects[frameRow][not status][i]
      for rpd in rpds:
        if rpd:
          rpd.grid_forget()

      #
      # Set current widgets
      #

      rpds = self.resonanceObjects[frameRow][status][i]
      
      for j in range(0,len(rpds)):
        rpd = rpds[j]
        
        if rpd:
          keywds = {}

          if j == 1:
            keywds['sticky'] = Tkinter.W

          if j == 0:
            if not status:
              keywds['rowspan'] = 2
              
          rpd.grid(row=frameRow + i, column=4 + j, **keywds)
      
  def setStereo(self,chain = None):
  
    if not chain:
      chain = self.chain
  
    if self.chainResonancesDict[self.resParent][chain]:
      for i in range(0,len(self.chainResonancesDict[self.resParent][chain])):

        atomSetGroup = self.chainAtomSetsDict[self.resParent][chain][i]
        (stereoStatus,stereoStatusChanged,resChanged,resonances) = self.chainResonancesDict[self.resParent][chain][i]
        
        if stereoStatusChanged or (stereoStatus and resChanged):
        
          residue = atomSetGroup[0].findFirstAtom().residue
          ccpCode = residue.molResidue.ccpCode
          seqCode = residue.seqCode
          
          atomText = "%s/%s" % (atomSetGroup[0].name,atomSetGroup[1].name)

          if stereoStatusChanged:
            if stereoStatus:
              stereoText = "non-stereospecific to stereospecific."
            else:
              stereoText = "stereospecific to non-stereospecific."
              
            print "  Warning: Changed '%s.%d' prochiral resonances for atoms %s from %s" % (ccpCode,seqCode,atomText,stereoText)
            self.chainResonancesDict[self.resParent][chain][i][1] = False
            self.chainResonancesDict[self.resParent][chain][i][2] = False
            
            if stereoStatus:
              resonanceSet = atomSetGroup[0].findFirstResonanceSet()
              
              if resonances[0]:
                resonanceSet.removeAtomSet(atomSetGroup[1])
                
              if resonances[1]:
                
                if not resonances[0]:
                  # Remove right atom set if only one resonance...
                  resonanceSet.removeAtomSet(atomSetGroup[0])
                else:
                  # Only create a new resonance set if two resonances present!
                  resonanceSet.removeResonance(resonances[1])
                  
                  if isinstance(self.resParent,NmrConstraint.NmrConstraintStore):
                    newResonanceSet = self.resParent.newFixedResonanceSet(resonances = [resonances[1]], atomSets = [atomSetGroup[1]])
                  else:
                    newResonanceSet = self.resParent.newResonanceSet(resonances = [resonances[1]], atomSets = [atomSetGroup[1]])
              
            else:
            
              if resonances[0]:
                resonanceSet = atomSetGroup[0].findFirstResonanceSet()
                if resonances[1]:
                  atomSetGroup[1].findFirstResonanceSet().delete()
                  resonanceSet.addResonance(resonances[1])
                resonanceSet.addAtomSet(atomSetGroup[1])
                
              elif resonances[1]:
                resonanceSet = atomSetGroup[1].resonanceSets[0]
                resonanceSet.addAtomSet(atomSetGroup[0])
              
              
          elif stereoStatus and resChanged:
            print "  Warning: Changed resonance assignment for '%s.%d' stereospecifically assigned atoms %s." % (ccpCode,seqCode,atomText)
            self.chainResonancesDict[self.resParent][chain][i][2] = False
            
            resonanceSets = []
            
            for i in range(0,2):
              if atomSetGroup[i].resonanceSets:
                resonanceSets.append(atomSetGroup[i].findFirstResonanceSet())
              else:
                resonanceSets.append(None)
              
            for i in range(0,2):
              resonanceSet = resonanceSets[i]
              
              if resonanceSet:
                resonanceSet.addAtomSet(atomSetGroup[not i])
                resonanceSet.removeAtomSet(atomSetGroup[i])
            
  def apply(self):
    
    for chain in self.chainList:
      self.setStereo(chain = chain)
    
    return True
