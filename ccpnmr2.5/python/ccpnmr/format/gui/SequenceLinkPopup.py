
"""
======================COPYRIGHT/LICENSE START==========================

SequenceLinkPopup.py: GUI window for determining sequence information from original format file

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

import string, copy

from memops.universal.Io import joinPath
from memops.universal.Util import returnInt

from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.MessageReporter import showError
from memops.gui.Entry import Entry
from memops.gui.DataEntry import askString

from ccpnmr.format.general.Io import getHelpUrlDir

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup
from ccpnmr.format.gui.SequenceChemCompSelect import SequenceChemCompSelect

class SequenceLinkPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'SequenceLink.html')

  def __init__(self, parent, project, createMoleculeDict):
   
    self.project = project
    self.localCreateMoleculeDict = createMoleculeDict
    self.createMoleculeDict = createMoleculeDict
    
    self.resetInfo()
    
    TemporaryBasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + 'Import original sequence', modal=False, transient=True)

  def body(self, master):
  
    self.frameMaster = master
    
    #
    # Setup header 
    #

    row = 0


    label = Label(master, text= "Original import")
    label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.EW)

    label = Label(master, text= "Data model", fg = 'blue')
    label.grid(row=row, column=1, columnspan = 3, sticky=Tkinter.EW)

    row += 1
    
    label = Label(master, text= "Sequence")
    label.grid(row=row, column=0, sticky=Tkinter.EW)

    label = Label(master, text= "Molecular type", fg = 'blue')
    label.grid(row=row, column=1, sticky=Tkinter.EW)

    label = Label(master, text= "Molecule name", fg = 'blue')
    label.grid(row=row, column=2, sticky=Tkinter.EW)

    label = Label(master, text= "Number of chains", fg = 'blue')
    label.grid(row=row, column=3, sticky=Tkinter.EW)
    
    row = self.setupWidgets(row = row)

    
  def makeNewInfo(self):

    newCreateMoleculeDict = {}
    
    #
    # Then get infor from self.localCreateMoleculeDict
    #
    
    self.localCreateMoleculeDict = newCreateMoleculeDict
    
  def setupWidgets(self,row = None):

    if row:
      self.widgetStartRow = row
    else:
      row = self.widgetStartRow
      
    #
    # Clear previous stuff...
    #
    
    if hasattr(self,'origSeqButton'):
      for molName in self.origSeqButton.keys():      
        row = self.molNameRow[molName]
        self.origSeqButton[molName].grid_forget()
        self.molType[molName].grid_forget()
        self.molNameEntry[molName].grid_forget()
        self.numChains[molName].grid_forget()
        self.warningCode[molName].grid_forget()
      self.buttons.grid_forget()

    #
    # Setup new widget dicts and lists
    #
    
    self.molNameEntry = {}
    self.origSeqButton = {}
    self.molType = {}
    self.numChains = {}
    self.warningCode = {}
    self.molNameRow = {}
    
    #
    # Create the specific information...
    #

    for molName in self.molNames:

      row = row + 1
      
      self.molNameRow[molName] = row
      
      self.origSeqButton[molName] = Tkinter.Button(self.frameMaster, text = self.origSeqLabel[molName], command = lambda molName = molName: self.popupSequenceChemCompSelect(molName))
      self.origSeqButton[molName].grid(row=row, column=0, sticky=Tkinter.EW)

      self.molType[molName] = Label(self.frameMaster, text = self.currentMolType[molName], fg = 'blue')
      self.molType[molName].grid(row=row, column=1, sticky=Tkinter.EW)

      self.molNameEntry[molName] = Tkinter.Button(self.frameMaster, text = molName, command = lambda molName = molName: self.changeMolName(molName))
      self.molNameEntry[molName].grid(row=row, column=2, sticky=Tkinter.EW)

      self.numChains[molName] = Entry(self.frameMaster, text = str(self.ccpNumChains[molName]), width = 2)
      self.numChains[molName].grid(row=row, column=3)

      self.warningCode[molName] = Tkinter.StringVar()
      self.warningCode[molName] = Label(self.frameMaster, textvariable = self.warningCode[molName], fg = 'red')
      self.warningCode[molName].grid(row=row, column=4, sticky=Tkinter.N)

    self.updateWarnings()

    row = row + 1

    texts = [ 'OK', 'Reset to original' ]
    commands = [ self.ok, self.resetOriginal ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    self.buttons = createDismissHelpButtonList(self.frameMaster, texts=texts, commands=commands, dismiss_text = 'Cancel', help_url=self.help_url)
    self.buttons.grid(row=row, column=0, columnspan = 5)
    
    return row

  def popupSequenceChemCompSelect(self,molName):
  
    popup = SequenceChemCompSelect(self,self.project,molName,self.createMoleculeDict,self.localCreateMoleculeDict)
    
    if popup and popup.hasChanged:
      self.localCreateMoleculeDict = popup.localCreateMoleculeDict
      self.resetInfo()
      self.setupWidgets()
    
  def changeMolName(self,molName):
    
    newMolName = askString('Enter molecule name','New name for molecule',molName,self)
    
    if newMolName and newMolName != molName:
    
      if newMolName in self.molNames:
        showError("Molecule name present","This molecule name is already present!")
      
      elif self.project.findFirstMolecule(name = newMolName):
        # TODO TODO: allow user to select an existing molecule?!?!
        showError("Molecule exists","A molecule with this name already exists in the data model.")
      
      else:    
        moleculeInfo = self.localCreateMoleculeDict[molName]
        self.localCreateMoleculeDict[newMolName] = moleculeInfo
        del self.localCreateMoleculeDict[molName]
        
        #
        # Clean out text...
        #
        
        for molName in self.molNames:
          self.warningCode[molName].set("")
        
        self.resetInfo()
        self.setupWidgets()

  def updateWarnings(self):
  
    for molName in self.molNames:
            
      #
      # Just use the first one as reference!
      #
      
      (origSequence,origSequenceList) = self.localCreateMoleculeDict[molName][1][0]
      
      texts = {}

      for seqIndex in range(0,len(origSequenceList)):
        (seqElement,chemComp,chemCompVar) = origSequenceList[seqIndex]
        
        text = None

        if not chemComp:
          text = 'chemComp(s) missing'

        elif not chemCompVar:
          text = 'chemCompVar(s) missing'
          
        if text:
          if not texts.has_key(text):
            texts[text] = 0
          texts[text] += 1

      text = ""
      
      if texts:
        
        textSorted = texts.keys()
        textSorted.sort()
        
        
        for textItem in textSorted:
          text += "Warning: %d %s in sequence.\n" % (texts[textItem],textItem)

      self.warningCode[molName].set(text)

  def resetOriginal(self):
  
    self.localCreateMoleculeDict = self.createMoleculeDict
    self.resetInfo()
    self.setupWidgets()

  def resetInfo(self):
  
    self.origSeqLabel = {}
    self.origSeqDict = {}
    
    self.currentMolType = {}
    
    self.ccpNumChains = {}
    
    self.molNames = self.localCreateMoleculeDict.keys()
    self.molNames.sort()
    
    for molName in self.molNames:
      
      molecule = self.localCreateMoleculeDict[molName][0]
      sequences = self.localCreateMoleculeDict[molName][1]
        
      self.ccpNumChains[molName] = len(sequences)
      
      #
      # Just use the first one as reference!
      #
      
      (origSequence,origSequenceList) = sequences[0]

      self.origSeqDict[molName] = {}

      #
      # Set the text to appear
      #

      seqTexts = []

      if len(origSequenceList) > 9:
        seqIndexes = [0,1,2,'...',-3,-2,-1]

      else:
        seqIndexes = range(0,len(origSequenceList))

      for seqIndex in seqIndexes:
        if type(seqIndex) == type(''):
          seqTexts.append(seqIndex)
        else:
          ccHead = origSequenceList[seqIndex][1]
          
          if not ccHead:
            seqEl = origSequenceList[seqIndex][0]
            if hasattr(seqEl,'code3Letter') and seqEl.code3Letter:
              seqText = seqEl.code3Letter
            else:
              seqText = seqEl.code1Letter
          else:
            seqText = ccHead.ccpCode
          
          seqTexts.append(seqText)

      self.origSeqLabel[molName] = string.join(seqTexts,'-')
      self.origSeqDict[molName][self.origSeqLabel[molName]] = origSequence

      #
      # Check for allowed molTypes
      #
      
      molTypes = []

      for seqIndex in range(0,len(origSequenceList)):
        (seqElement,chemComp,chemCompVar) = origSequenceList[seqIndex]

        if chemComp:
          if chemComp.molType not in molTypes:
            molTypes.append(chemComp.molType)
      
      if not molTypes or len(molTypes) > 1:
        molType = 'other'
      else:
        molType = molTypes[0]

      self.currentMolType[molName] = molType
  
  def apply(self):
  
    returnStatus = True
  
    self.createMoleculeDict = self.localCreateMoleculeDict
    
    #
    # Have to check the number of chains and check if any info available...
    #
    
    self.molNames = self.createMoleculeDict.keys()
    self.molNames.sort()
    
    for molName in self.molNames:
      
      molecule = self.createMoleculeDict[molName][0]
      sequences = self.createMoleculeDict[molName][1]
      
      #
      # Remove altogether if there's no info..
      #
      
      if not sequences or not sequences[0] or not sequences[0][1]:
        del self.createMoleculeDict[molName]
        continue
        
      #
      # Check number of chains
      #
      
      numChains = returnInt(self.numChains[molName].get(), default = None)
      
      if numChains != None:
        numSequences = len(sequences)
        
        if numChains == 0:
          for sequence in sequences:
            sequence[0] = None
            
        elif numChains > numSequences:
          (sequence,sequenceList) = sequences[0]
          for i in range(numSequences,numChains):
            self.createMoleculeDict[molName][1].append((sequence,sequenceList[:]))
        
        elif numChains < numSequences:
          for i in range(numSequences,numChains,-1):
            self.createMoleculeDict[molName][1].pop(i-1)
    
    return returnStatus
