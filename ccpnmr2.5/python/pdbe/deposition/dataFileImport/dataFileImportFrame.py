"""
======================COPYRIGHT/LICENSE START==========================

dataFileImportFrame.py: Frame code to handle import of files for PDBe deposition

Copyright (C) 2010 Wim Vranken (PDBe, EBI)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.

A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wim@ebi.ac.uk
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

import os

from memops.gui.Base import getPopup
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
#from memops.gui.DataEntry import askString
from memops.gui.Frame import Frame
from memops.gui.HelpPopup import showHelpText
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.MessageReporter import showYesNo, showOkCancel, showWarning, showInfo
from memops.gui.PulldownList import PulldownList
from memops.gui.Separator import Separator
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.Text import Text

from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FileSelect import FileType

from ccpnmr.format.gui.SelectionListPopup import SelectionListPopup
from ccpnmr.format.general.Constants import fileTypeDict

from ccpnmr.eci.EntryCompletionGui import EntryCompletionGui

class DataFileImportFrame(Frame):

  """
  TODO: should I also ask for sample conditions info with chemical shifts? Or do in ECI - possible? Check!
  """

  # TODO should go where?
  email = 'wim@ebi.ac.uk'

  defaultSelectText = "Select file."
  fontOkColor = 'green3'
  fontCheckColor = 'orange3'
  fontBadColor = 'red3'
  fontDefaultColor = 'black'

  def __init__(self, parent, basePopup, *args, **kw):

    #
    # Variable initialisation
    #
    
    self.fcWrapper = basePopup.fcWrapper
    self.project = basePopup.project
    self.basePopup = basePopup

    # TODO necessary?
    if self.project:
      self.nmrProject = self.project.currentNmrProject
      self.entry = self.project.currentNmrEntryStore.findFirstEntry()
      if not self.entry:
        self.entry = self.project.currentNmrEntryStore.newEntry(name=self.project.name)
    else:
      self.nmrProject = None  

    self.sequenceCoordinatesLoaded = False
    self.shiftsLoaded = False
    self.linkResDone = False
    
    self.currentShiftList = None
    self.shiftListChainPairs = []

    self.moleculeList = []
    self.moleculeDict = {}
      
    #
    # Frame setup
    #
    
    Frame.__init__(self, parent, **kw)

    self.grid_columnconfigure(1, weight=1)
    self.grid_rowconfigure(1, weight=1)

    options = ['Import', 'Deposition']

    tabbedFrame = TabbedFrame(self, options=options, callback=self.selectTab)
    tabbedFrame.grid(row=1, column=0, columnspan=2, sticky='nsew')

    self.tabbedFrame = tabbedFrame
    frameA, frameD = tabbedFrame.frames

    #
    # Main
    #

    frameA.grid_columnconfigure(0, weight=1)  
    #frameA.grid_columnconfigure(1, weight=1) # Change to 2 if want 2 columns
    frameA.grid_rowconfigure(12, weight=1)
    #frameA.grid_rowconfigure(12, weight=1)

    row = 0

    div = LabelDivider(frameA, text='Select the full coordinate file or a sequence file (with info for a single molecule).', justify = 'center', grid=(row,0), gridSpan=(1,1))
    
    row += 1
    
    self.sequenceCoordinatesImport = Button(frameA, text = self.defaultSelectText , command = self.importSequenceOrCoords, foreground = self.fontDefaultColor)
    self.sequenceCoordinatesImport.grid(row=row, column=0, sticky='ew')

    row += 1

    label =Label(frameA,text="")
    label.grid(row=row,column=0,sticky='ew')
    
    row += 1

    div = LabelDivider(frameA, text='Select the molecule relevant for your chemical shift file.', justify = 'center', grid=(row,0), gridSpan=(1,1))
    
    row += 1
    
    self.moleculeSelect = Label(frameA, text = "None available yet - import valid file first", foreground = self.fontBadColor)
    self.moleculeSelect.grid(row=row, column=0, sticky='ew')
    self.moleculeSelectRow = row
    
    row += 1

    label =Label(frameA,text="")
    label.grid(row=row,column=0,sticky='ew')
    
    row += 1

    div = LabelDivider(frameA, text='Select a chemical shift file with values only for the above molecule.', justify = 'center', grid=(row,0), gridSpan=(1,1))
    
    row += 1

    self.shiftImport = Button(frameA, text = self.defaultSelectText , command = self.importShifts, foreground = self.fontDefaultColor)
    self.shiftImport.grid(row=row, column=0, sticky='ew')
    
    row += 1

    label =Label(frameA,text="")
    label.grid(row=row,column=0,sticky='ew')

    row += 1              

    div = LabelDivider(frameA, text='Consistency check between molecule and shift information.', justify = 'center', grid=(row,0), gridSpan=(2,1))

    row += 1
    
    self.linkResCheckInfo = Label(frameA, text='')
    self.linkResCheckInfo.grid(row=row, column=0, sticky='ew')

    row += 1

    div = Separator(frameA, grid=(row,0), gridSpan=(1,1))

    row += 1

    texts = ['Import new sequence','Import new set of shifts']
    commands = [self.resetSequenceImport,self.resetShiftImport]

    self.mainButtons = ButtonList(frameA, texts=texts, commands=commands)
    self.mainButtons.grid(row=row, column=0, columnspan=2, sticky='ew')
    self.mainButtons.buttons[0].config(foreground = self.fontDefaultColor)
    
    #print row
    
    self.frameA = frameA

    #
    # Not in use...
    #

    #frameX.grid_columnconfigure(0, weight=1)
    #frameX.grid_columnconfigure(1, weight=1)
    #frameX.grid_rowconfigure(1, weight=1)
    #frameX.grid_rowconfigure(3, weight=1)


    #
    # Deposition, is updated after each successful import run.
    #

    frameD.grid_columnconfigure(0, weight=1)
    frameD.grid_rowconfigure(5, weight=1)
    
    self.frameD = frameD

    row = 0

    div = LabelDivider(frameD, text='Imported data.', justify = 'center', grid=(row,0), gridSpan=(1,2))
    
    row += 1
    
    self.depositionImportText = "\nImported %d shift list(s) for a total of %d shifts.\n\nImported %d molecule(s) and %d chain(s).\n\nImported %d model(s) for a total of %d atom coordinates.\n\nLinked %.2f%% of imported NMR information to %d chain(s)."
    self.depositionImportLoc = (row,0)
    
    # These used for setting text above...
    self.depositionImportNums = [0,0,0,0,0,0,0.0,0]
    self.importedShiftLists = []
    self.connectedChains = []
    
    self.depositionImportLabel = Label(frameD, text = self.depositionImportText % tuple(self.depositionImportNums), foreground = self.fontBadColor)
    self.depositionImportLabel.grid(row=row, column=0, sticky='ew')

    row += 1

    label =Label(frameD,text="")
    label.grid(row=row,column=0,sticky='ew')

    #
    # Finalize the import part, proceed to ECI.
    #
    
    row += 1

    div = LabelDivider(frameD, text='Import completed, save project and start Entry Completion Interface.', justify = 'center', grid=(row,0), gridSpan=(1,1))
    
    row += 1

    self.eciStart = Button(frameD, text = "Finalise import" , command = self.finaliseImport, foreground = self.fontBadColor)
    self.eciStart.grid(row=row, column=0, sticky='ew')
      
  def finaliseImport(self):
  
    if not self.depositionImportNums[6]:
      showWarning('Failure','Need connected molecule and shift information first')
      return  
  
    if showYesNo("Save project and continue annotation","Are you sure you want to save this project and continue?",parent=self):
      
      if self.depositionImportNums[5] or showYesNo("No coordinates","No coordinates are available - are you sure you want to continue?",parent=self):

        showInfo("Project name","Your project will be saved in the %s directory." % self.project.name,parent=self)

        self.project.saveModified()

        #userData = self.project.findFirstRepository(name='userData')
        #currentPath = userData.url.path
        #currentProjectName = self.project.name

        #projectDir = os.path.join(currentPath,currentProjectName)

        eci  = EntryCompletionGui(self.basePopup.root)
        eci.initProject(self.project)
  
  #
  # File type determination and import
  #
  
  def createFileTypes(self,dataTypes):

    fileTypes = [FileType('all', ['*'])]
    
    for dataType in dataTypes:
      formatNames = fileTypeDict[dataType].keys()
      formatNames.sort()
      for formatName in formatNames:
        if formatName in self.fcWrapper.formatNameLists[dataType]:
          fileTypes.append( FileType(formatName, fileTypeDict[dataType][formatName]))
          
    return fileTypes
  
  def getFileName(self,title,fileTypes):
  
    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                        title=title, dismiss_text='Cancel',
                        selected_file_must_exist=True, multiSelect=False,)

    self.fileName = fileSelectPopup.getFile()

    if not self.fileName:
      showWarning('Failure','Please select an existing file.', parent=self)
      return False
    
    return True
  
  def importSequenceOrCoords(self):
    
    dataTypes = ['sequence','coordinates']
    fileTypes = self.createFileTypes(dataTypes)

    if self.getFileName('Import sequence or coordinate file',fileTypes):
    
      formatNameSuggestions = {}
      
      for dataType in dataTypes:
      
        tmpList = self.fcWrapper.determineFormatNamesForFile(dataType,self.fileName)
        
        if tmpList:
          formatNameSuggestions[dataType] = tmpList

      if not formatNameSuggestions:       
        showWarning('Failure','This file cannot be read by this software.\nPlease send the file with an explanation to %s.' % self.email, parent=self)
        return False
      
      #
      # Let user select if multiple options, otherwise take default
      #
            
      if len(formatNameSuggestions) == 1 and len(formatNameSuggestions[formatNameSuggestions.keys()[0]]) == 1:
        
        dataType = formatNameSuggestions.keys()[0]
        formatName = formatNameSuggestions[dataType][0]

        if not showYesNo('File type detected', 'Reading as %s file in %s format. Is this correct?' % (dataType,formatName), parent=self):
          showWarning('Failure','This file cannot be read by this software.\nPlease send the file with an explanation to %s.' % self.email, parent=self)
          return False
       
      else:
      
        #
        # Create a selection (hopefully user-understandable)
        #

        selectionList = []
        selectionDict = {}

        for dataType in dataTypes:

          dataTypeString = dataType

          if formatNameSuggestions.has_key(dataType):
            formatNames = formatNameSuggestions[dataType]
            formatNames.sort()

            for formatName in formatNames:
              selectionString = "%s file in %s format." % (dataTypeString,formatName)
              selectionList.append(selectionString)
              selectionDict[selectionString] = (dataType,formatName)
                  
        interaction = SelectionListPopup(self, selectionList, title = 'File format selection', text = 'This is a:', selectionDict = selectionDict, dismissButton = True, modal = True)

        #
        # Check if anything was selected...
        #
        
        dataType = formatName = None
          
        if interaction.isSelected:
          (dataType,formatName) = interaction.selection
        else:
          showWarning('Failure','This file cannot by read without a format selection.\nIf the correct format is not available, please send the file with an explanation to %s' % self.email, parent=self)
          return False
      
      #
      # Now read the file, need to do some field updates!
      #
      
      (fileRead,fileInformation) = self.fcWrapper.readFile(dataType,formatName,self.fileName)
                      
      if not fileRead:
        showWarning('Failure','This file cannot be read by this software:%s\nPlease send the file with an explanation to %s.' % (fileInformation,self.email), parent=self)
        return False
      
      (conversionInfo,conversionSuccess) = (self.fcWrapper.formatConversion.conversionInfo,self.fcWrapper.formatConversion.conversionSuccess)
      
      if not conversionSuccess:
        showWarning('Failure','This file was read by the software but contains invalid information.\nPlease send the file with an explanation to %s.' % self.email, parent=self)
        return False
      
      #
      # Set info if import worked OK
      #
      
      conversionLines = conversionInfo.split(": ")
      
      showInfo("Import coordinate and/or sequence information",":\n".join(conversionLines),parent=self)
      
      if dataType == 'sequence':
        chains = self.fcWrapper.importReturns[dataType]
        models = []
      elif dataType == 'coordinates':
        models = self.fcWrapper.importReturns[dataType]
        chains = [cChain.chain for cChain in models[0].structureEnsemble.coordChains]
      
      self.sequenceCoordinatesImport.setText(self.fileName) # TODO change color or something?
      self.sequenceCoordinatesImport.configure(foreground = self.fontOkColor)

      self.sequenceCoordinatesLoaded = True
      
      #
      # Reset to list selector for further use
      #
      
      moleculeName = None
      
      for chain in chains:
      
        if chain in self.moleculeDict.values():
          continue
      
        numResidues = len(chain.residues)
        if numResidues == 1:
          residueText = "%s residue" % chain.findFirstResidue().ccpCode
        else:
          residueText = "%d residues" % numResidues
      
        moleculeName = "%s (chain '%s', %s)" % (chain.molecule.name,chain.code,residueText)
        
        self.moleculeList.append(moleculeName)
        self.moleculeDict[moleculeName] = chain
        
      self.moleculeList.sort()

      self.moleculeSelect.destroy()
      
      if len(chains) == 1 and moleculeName:
        selectedIndex = self.moleculeList.index(moleculeName)
      else:
        selectedIndex = 0

      self.moleculeSelect = PulldownList(self.frameA, texts=self.moleculeList, index=selectedIndex, sticky='ew')
      self.moleculeSelect.grid(row=self.moleculeSelectRow, column=0)
     
      #
      # Now update Deposition tab
      #
      
      molecules = []
      for chain in chains:
        if not chain.molecule in molecules:
          molecules.append(chain.molecule)
      
      numCoords = 0
      for model in models:
        numCoords += len(model.coords)

      self.updateDepositionImportLabel(molecules=len(molecules),chains=len(chains),models=len(models),coordinates=numCoords)
      
      # Add the molSystem to the entry!
      if chains and not chains[0].molSystem == self.entry.molSystem:
        self.entry.molSystem = chains[0].molSystem

    self.updateAll()

  def importShifts(self):
    
    currentChain = self.getCurrentChain()

    if not currentChain:
      showWarning('Failure','Please first read in a sequence or coordinate file and select the molecule relevant for this shift list.', parent=self)
      return
    
    elif self.currentShiftList:
      shiftListChainPair = (self.currentShiftList,currentChain)
      
      if shiftListChainPair in self.shiftListChainPairs:
        showWarning('Failure',"You already read in chemical shifts for this chain.\nPlease read in related shifts for the other chain(s), if present, or press the 'Import new set of shifts' button to read in a new set of shifts.", parent=self)
        return   
    
    dataType = 'shifts'
    
    fileTypes = self.createFileTypes([dataType])

    if self.getFileName('Import chemical shift file',fileTypes):
      
      formatNameSuggestions = self.fcWrapper.determineFormatNamesForFile(dataType,self.fileName)

      if not formatNameSuggestions:       
        showWarning('Failure','This file cannot be read by this software.\nPlease send the file with an explanation to %s.' % self.email, parent=self)
        return False
      
      #
      # Let user select if multiple options, otherwise take default
      #
            
      if len(formatNameSuggestions) == 1:
        
        formatName = formatNameSuggestions[0]

        if not showYesNo('File type detected', 'Reading as a chemical shift file in %s format. Is this correct?' % formatName, parent=self):
          showWarning('Failure','This file cannot be read by this software.\nPlease send the file with an explanation to %s.' % self.email, parent=self)
          return False
        
      else:
      
        #
        # Create a selection (hopefully user-understandable)
        #

        selectionList = []
        selectionDict = {}
        
        formatNameSuggestions.sort()

        for formatName in formatNameSuggestions:
          selectionString = "chemical shift file in %s format." % (formatName)
          selectionList.append(selectionString)
          selectionDict[selectionString] = formatName
                  
        interaction = SelectionListPopup(self, selectionList, title = 'File format selection', text = 'This is a:', selectionDict = selectionDict, dismissButton = True, modal = True)

        #
        # Check if anything was selected...
        #
        
        formatName = None
          
        if interaction.isSelected:
          formatName = interaction.selection
        else:
          showWarning('Failure','This file cannot by read without a format selection.\nIf the correct format is not available, please send the file with an explanation to %s' % self.email, parent=self)
          return False
      
      #
      # Now read the file, need to do some field updates! Also make sure to re-use shift list for other molecules...
      #
      
      (fileRead,fileInformation) = self.fcWrapper.readFile(dataType,formatName,self.fileName,addKeywords = {'measurementList': self.currentShiftList})
              
      if not fileRead:
        showWarning('Failure','This file cannot be read by this software:%s\nPlease send the file with an explanation to %s.' % (fileInformation,self.email), parent=self)
        return False
      
      (conversionInfo,conversionSuccess) = (self.fcWrapper.formatConversion.conversionInfo,self.fcWrapper.formatConversion.conversionSuccess)
      
      if not conversionSuccess:
        showWarning('Failure','This file was read by the software but contains invalid information.\nPlease send the file with an explanation to %s.' % self.email, parent=self)
        return False
      
      #
      # Set info if import worked OK
      #
      
      conversionLines = conversionInfo.split(": ")
      
      showInfo("Import chemical shift information",":\n".join(conversionLines),parent=self)
      
      self.shiftImport.setText(self.fileName) # TODO change color or something?
      self.shiftImport.configure(foreground = self.fontOkColor)
      self.shiftsLoaded = True
      self.shiftsFormatName = formatName

      shiftList = self.fcWrapper.importReturns[dataType]
      
      if not self.currentShiftList:
        self.currentShiftList = shiftList
        
      self.shiftListChainPairs.append((self.currentShiftList,currentChain))
            
      self.updateDepositionImportLabel(shiftList=shiftList)

      if not shiftList in self.entry.measurementLists:
        print shiftList
        self.entry.addMeasurementList(shiftList)
        print self.entry.measurementLists

    self.updateAll()

  
  #
  # Updaters
  #

  def selectTab(self, index):

    funcsDict = {0:(self.updateMain, ),
                 1:(self.updateCoordinates, ),
                 2:(self.updateDeposition, )}

    for func in funcsDict[index]:
      func()
      
  def updateMain(self):
  
    pass
    
  def updateCoordinates(self):
  
    pass
    
  def updateDeposition(self):
  
    pass

  def updateAll(self):

    self.selectTab(self.tabbedFrame.selected)
    
    if self.sequenceCoordinatesLoaded and self.shiftsLoaded:
      if showYesNo('Connect shifts to sequence', 'You have to check whether the chemical shift information matches the sequence. Do you want to this now?', parent=self):
        self.connectShiftsSequence()
        
    self.waiting = False


  def getCurrentChain(self):
  
    if len(self.moleculeList) == 1:
      chain = self.moleculeDict.values()[0]
    else:
      try:
        moleculeSelected = self.moleculeSelect.getText()
        chain = self.moleculeDict[moleculeSelected]
      except:
        chain = None
      
    return chain
  
  
  def connectShiftsSequence(self):
  
    if not self.linkResDone:
      
      changeResetColor = False
      
      #
      # Get chain mapping and run linkResonances
      #
      
      chain = self.getCurrentChain()
      
      forceChainMappings = self.fcWrapper.linkResonancesToSequence(chain=chain)
      self.fcWrapper.formatConversion.linkResonances(forceChainMappings=forceChainMappings,guiParent=self)        
      
      #
      # Get information about the linking process
      #
      # TODO Should only have an nmrProject (no restraint import, should be included?)
      # In any case, is always the LAST info in numResonancesLinked info
      
      numResonancesLinked = self.fcWrapper.formatConversion.numResonancesLinked
      (origLinked,origUnlinked,linked,unlinked) = (numResonancesLinked['origLinked'][-1],numResonancesLinked['origUnlinked'][-1],numResonancesLinked['linked'][-1],numResonancesLinked['unlinked'][-1])

      #
      # Track number of new resonances, and reset for new import
      #
      
      if self.fcWrapper.formatConversion.allResonancesLinked:
        status = 'All information matches (for all imports).'
        foreground = self.fontOkColor
        changeResetColor = True
      else:
        if origUnlinked - unlinked == self.fcWrapper.numNewResonances: 
          status = 'All information matches (for this import)'
          foreground = self.fontOkColor
          changeResetColor = True
        else:
          if origUnlinked != unlinked:
            status = 'Not all information matches (for this and/or another import).'
            foreground = self.fontCheckColor
            changeResetColor = True
          else:
            status = 'No information matches (for this import).'
            foreground = self.fontBadColor
          
          otherUnlinked = (origUnlinked - self.fcWrapper.numNewResonances)
          notLinked = unlinked - otherUnlinked

          status += "\nUnable to link %d out of %d imported shifts (%.2f%%)." % (notLinked,self.fcWrapper.numNewResonances,(notLinked * 100.0) / self.fcWrapper.numNewResonances)

      self.linkResCheckInfo.set("Status: %s" % status)             
      self.linkResCheckInfo.config(foreground = foreground)

      self.linkResDone = True
      
      #
      # Change the color of the reset button to indicate OK to do next one
      #
      
      if changeResetColor:
        
        self.mainButtons.buttons[0].config(foreground = self.fontOkColor)
        
        self.updateDepositionImportLabel(shiftList=None,percentageLinked=(linked * 100.0/(unlinked+linked)),connectedChain=chain)


  def updateDepositionImportLabel(self,shiftList=None,molecules=0,chains=0,models=0,coordinates=0,percentageLinked=None,connectedChain=None):
  
    if shiftList and shiftList not in self.importedShiftLists:
      self.importedShiftLists.append(shiftList)
      self.depositionImportNums[0] += 1
      
    shifts = 0
    for shiftList in self.importedShiftLists:
      shifts += len(shiftList.measurements)
    
    self.depositionImportNums[1] += shifts
    self.depositionImportNums[2] += molecules
    self.depositionImportNums[3] += chains
    self.depositionImportNums[4] += models
    self.depositionImportNums[5] += coordinates
    
    if percentageLinked != None:
      self.depositionImportNums[6] = percentageLinked
    
    if connectedChain and connectedChain not in self.connectedChains:
      self.depositionImportNums[7] += 1
      self.connectedChains.append(connectedChain)
      
    self.depositionImportLabel.destroy()
    
    finalForeground = self.fontBadColor
    if self.depositionImportNums[0] == 0 and self.depositionImportNums[2] == 0:
      # Nothing imported
      foreground = self.fontBadColor
      
    elif self.depositionImportNums[6]:
      # Linked shifts available - TODO base this on % of shifts linked?
      foreground = self.fontOkColor
      if self.depositionImportNums[5]:
        finalForeground = self.fontOkColor
      else:
        finalForeground = self.fontCheckColor
      
    else:
      # Intermediate state
      foreground = self.fontCheckColor
    
    self.depositionImportLabel = Label(self.frameD, text = self.depositionImportText % tuple(self.depositionImportNums), foreground = foreground)
    self.depositionImportLabel.grid(row=self.depositionImportLoc[0], column=self.depositionImportLoc[1], sticky='ew')

    self.eciStart.configure(foreground = finalForeground)

  def resetSequenceImport(self):

    doReset = True

    if not self.linkResDone and self.sequenceCoordinatesLoaded and self.shiftsLoaded:
      if showYesNo('Shifts not connected to sequence', 'You have not checked whether the imported chemical shift information matches the imported sequence. Do you want to this first? If not, the last imported data will be invalid.', parent=self):
        self.connectShiftsSequence()
        doReset = False
 
    if doReset:
    
      self.mainButtons.buttons[0].config(foreground = self.fontDefaultColor)

      self.sequenceCoordinatesLoaded = self.shiftsLoaded = self.linkResDone = False

      self.sequenceCoordinatesImport.setText(self.defaultSelectText) 
      self.sequenceCoordinatesImport.configure(foreground = self.fontDefaultColor)
    
      self.moleculeSelect.destroy()
      self.moleculeSelect = Label(self.frameA, text = "None available yet - import valid file first", foreground = self.fontBadColor)
      self.moleculeSelect.grid(row=self.moleculeSelectRow, column=0, sticky='ew')

      self.shiftImport.setText(self.defaultSelectText) 
      self.shiftImport.configure(foreground = self.fontDefaultColor)

      self.linkResCheckInfo.set("")

  def resetShiftImport(self):

    doReset = True

    if not self.linkResDone and self.sequenceCoordinatesLoaded and self.shiftsLoaded:
      if showYesNo('Shifts not connected to sequence', 'You have not checked whether the imported chemical shift information matches the imported sequence. Do you want to this first? If not, the last imported data will be invalid.', parent=self):
        self.connectShiftsSequence()
        doReset = False
 
    if doReset:
    
      self.mainButtons.buttons[1].config(foreground = self.fontDefaultColor)

      self.shiftsLoaded = self.linkResDone = False

      self.currentShiftList = None
      
      self.shiftImport.setText(self.defaultSelectText) 
      self.shiftImport.configure(foreground = self.fontDefaultColor)

      self.linkResCheckInfo.set("")

  def destroy(self):

    Frame.destroy(self)

  #
  # Instructions
  #
  
  def showMainInstructions(self):

    popup = getPopup(self)

    message = """Use this tab to import the chemical shifts and the coordinate and/or sequence information for your molecular chains.
    
The imported chemical shift file should contain information for only *one* molecular chain to make it easier to connect the molecule information to the chemical shift information. You therefore have to select a single chain for each shift list from the dropdown menu that will appear after you imported a coordinate or sequence file.
    
For example, when using sequence files for a dimer, import the information for the chemical shifts for each chain separately:
    
  1. Import the sequence for the first molecular chain.
  2. Import a chemical shift file with shifts only for this molecular chain.  
  3. Reset using the 'Import new sequence' button
  4. Import the sequence for the second molecular chain
  5. Import a chemical shift file with shifts only for this second chain
  
Alternatively, it is possible to read in the molecule information from a full coordinate file:

  1. Import a coordinate file with all molecular information, including coordinates.
  2. Select a molecular chain.
  3. Import a chemical shift file with shifts only for the selected molecular chain.
  4. Go back to step 2. if necessary.
  
You can also import multiple sets of chemical shifts (e.g. for the sample in different conditions). In this case, you have to import all chemical shift information that belongs together for all molecular chains, then press the 'Import new set of shifts' button.

Notes:

1. This application always creates a new CCPN project. It is not possible to import files into existing projects.
2. If your chemical shift file contains information for multiple chains, you have to edit it manually to split up the information per chain.
    """

    showHelpText(self, message, popup = popup)

  def showFormats(self):

    popup = getPopup(self)

    message = """For chemical shifts, the following formats are recognised:

*** Auremol ***

section_sequenzdefinition
_Residue_seq_code
_Atom_num_code
_Residue_label
_Atom_name
_Atom_type
_Atom_alias
_Atom_equivalent
_Atom_CSA
  1   1 MET   HN H   -      -            8.95
  1   2 MET    N N   -      -          157.00
  1   3 MET   CA C   -      -           40.00


*** Autoassign ***

AA        HN    N15    CO-1   CA-1   CB-1   HA-1          CO     CA     CB     HA

A31      8.14  121.4                                            51.3   19.6                 (GS178 115.HSQC)
D32      8.88  122.9          51.3   19.5                       55.4   39.6                 (GS271 22.HSQC)


*** CNS ***

do ( store1 = 53.13218 ) ( resid 78 and name CA )
do ( store1 = 0.7356673 ) ( resid 15 and name HD1# )
do ( store1 = 120.5381 ) ( resid 8 and name N )
do ( store1 = 121.1414 ) ( resid 78 and name N )


*** Cosmos ***

CS_VALUES 3
C_ALA 176.6
CA_ALA 51.66
CB_ALA 17.26
END


*** CSI ***

#     A       HA       CA       CO       CB       Consensus
#
1     MET       0 C      0 C      NA       0 C         0 C 
2     GLY       0 C      0 C      NA       0 C         0 C 


*** MARS ***

            H         N         CO-1      CA        CA-1   
PR_2        8.900   123.220   170.540    55.080    54.450  
PR_4        8.320   115.340   175.920      -       55.080  


*** MONTE ***

         1            102.544      8.211     45.853     54.925      0.000     18.069    180.112 
         2            103.276      8.580     45.334     54.154      0.000     35.650    175.087 
         3            103.997      7.407     45.165      0.000      0.000      0.000      0.000 


*** NMR-STAR ***

data_test

save_shifts1
   _Saveframe_category               assigned_chemical_shifts

   loop_
      _Atom_shift_assign_ID
      _Residue_seq_code
      _Residue_label
      _Atom_name
      _Atom_type
      _Chem_shift_value
      _Chem_shift_value_error
      _Chem_shift_ambiguity_code

          1     1   ASP  CA    C   52.000  0.02  1  
          2     1   ASP  HA    H    4.220  0.02  1  

   stop_

save_


*** NMRVIEW ***

  2.CG1     18.549 0
  2.CG2     18.844 0
  2.HG1#     0.800 0
  2.HG2#     0.723 0
  3.HG2      2.298 0
  3.HG1      2.298 0


*** PIPP ***

RES_ID          1
RES_TYPE        MET
SPIN_SYSTEM_ID  1
    CA        55.9920
    CB        33.1470
    HA         4.1141
    HB#        2.0492
    HG#        2.4250
END_RES_DEF


*** PISTACHIO ***

   1    1  GLY    C     C  172.621  1.000  0 
   2    1  GLY   CA     C   44.308  1.000  0 
   3    2  SER    N     N  122.241  1.000  0 


*** PRONTO ***

Spin system   HN          HA      Other:

1: Val-1                  3.76    HB: 1.945, HG1: 0.770, HG2: 0.608
2: Ile-2      8.80        4.26    HB: 1.526, HG1: 1.278, HG2: 0.728, HD: 0.918


*** SHIFTX ***

  NUM RES   HA     H       N        CA      CB       C
--- --- ------ ------ -------- ------- ------- --------
 2     T  4.4161 8.1749 111.0443 61.8324 70.3867 172.5362
 3     Y  4.9022 9.0239 120.2106 56.0493 41.4218 173.0761

 NUM RES  H    HA   HB   HB2  HB3  HD1  HD2  HD21 HD22 HD3  HE   HE1 HE2  HE21 HE22 HE3  HG   HG1  HG12 HG13 HG2  HG3  HZ
 2     T  8.17 4.42 4.24 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 1.21 0.00 0.00
 3     Y  9.02 4.90 0.00 2.22 2.20 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00


*** SPARKY ***

 Group   Atom  Nuc    Shift   SDev  Assignments

     R2     CA  13C   56.539  0.003      3
     R2     CB  13C   30.808  0.009      3


*** TALOS ***

VARS   RESID RESNAME PHI PSI DPHI DPSI DIST COUNT CLASS
FORMAT %4d %s %8.3f %8.3f %8.3f %8.3f %8.3f %2d %s

   1 Q 9999.000 9999.000    0.000    0.000    0.000  0 None
   2 N  -85.000  124.000   23.000   28.000   85.920 10 Good


***XEASY/CYANA:

   1 117.803 0.000 N       1
   2   8.208 0.002 HN      1
   3  56.508 0.055 CA      1
   4 999.000 0.000 HA      1
   5  29.451 0.004 CB      1
    
    """

    showHelpText(self, message, popup = popup)

# General functions to be moved to a general spot later.

def getTopLevelStore(memopsRoot, className):

  if not memopsRoot:
    return

  store = getattr(memopsRoot, 'current'+className) or \
          getattr(memopsRoot, 'findFirst'+className)() or \
          getattr(memopsRoot, 'new'+className)(name='eciDefault')

  return store

def getTopLevelStore2(memopsRoot, className):

  if not memopsRoot:
    return

  store = getattr(memopsRoot, 'current'+className) or \
          getattr(memopsRoot, 'findFirst'+className)() or \
          getattr(memopsRoot, 'new'+className)(namingSystem='eciDefault')

  return store

def getMethodStore(memopsRoot):

  return getTopLevelStore(memopsRoot, 'MethodStore')

