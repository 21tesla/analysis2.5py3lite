
"""
======================COPYRIGHT/LICENSE START==========================

GroupChemShiftsPopup.py: GUI window for grouping chemical shifts and checking percentage of assignment

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

from memops.universal.Io import joinPath

from memops.gui.BasePopup import BasePopup
from memops.gui.Label import Label
from memops.gui.CheckButton import CheckButton
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledFrame import ScrolledFrame
from memops.gui.Separator import Separator
from memops.gui.MessageReporter import showWarning

from ccpnmr.format.gui.FormatFilePopup import FormatFilePopup

from ccpnmr.format.general.Io import getHelpUrlDir
from ccpnmr.format.general.Constants import chemShiftClasses, newline

from ccpnmr.format.process.makeChemShiftSelections import makeChemShiftSelections

class GroupChemShiftsPopup(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'GroupChemShifts.html')

  def __init__(self, parent, project):
    
    #
    # Info for writing file...
    #
    
    self.groupText = {}
    
    #
    # Set up molSystem information
    #
    
    self.project = project
    
    self.molSysList = list(project.sortedMolSystems())
    self.molSystems = None
    
    self.molSysLabelList = []
    self.molSysLabelDict = {}
    self.molSysRelationDict = {}
    
    molSysLabel = 'All molSystems'
    
    if self.setChainLabels(molSysLabel,self.molSysList):
      self.molSysLabelList.append(molSysLabel)
      self.molSysLabelDict[molSysLabel] = self.molSysList
    
    for molSys in self.molSysList:
      molSysLabel = '%s' % molSys.code

      if self.setChainLabels(molSysLabel,[molSys]):
        self.molSysLabelList.append(molSysLabel)
        self.molSysLabelDict[molSysLabel] = [molSys]
        
    
    if not self.molSysLabelList:
      showWarning('Warning','No chemical shift lists available! Exiting...')
      return
      
    #
    # Some initializing...
    #
    
    self.chains = None
    self.shiftList = None
    self.shiftListLabel = None
    self.results = None

    # modal = true means that it won't continue unless this one returns value
    BasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + 'Group chemical shift values', modal=False, transient=True)
    
  def setChainLabels(self,molSysLabel,molSysList):
   
    #
    # Set up chain information
    #
    
    chainLabelList = []
    chainLabelDict = {}
    chainLabelShiftListDict = {}
    
    self.molSysRelationDict[molSysLabel] = [chainLabelList,chainLabelDict,chainLabelShiftListDict]

    chains = []
    
    for molSys in molSysList:
      chains.extend(list(molSys.sortedChains()))
    
    chainLabel = 'All chains'
    
    if self.setShiftListLabels(chainLabel,chains,self.molSysRelationDict[molSysLabel][2]):
      self.molSysRelationDict[molSysLabel][0].append(chainLabel)
      self.molSysRelationDict[molSysLabel][1][chainLabel] = chains
        
    for chain in chains:
      chainLabel = "'%s' (mol. '%s')" % (chain.code,chain.molecule.name)
    
      if self.setShiftListLabels(chainLabel,[chain],self.molSysRelationDict[molSysLabel][2]):
        self.molSysRelationDict[molSysLabel][0].append(chainLabel)
        self.molSysRelationDict[molSysLabel][1][chainLabel] = [chain]
    
    return self.molSysRelationDict[molSysLabel][0]
    
  def setShiftListLabels(self,chainLabel,chains,chainLabelShiftListDict):

    #
    # Set up chemical shift list information (slooooww so done at start)
    #
    
    shiftLists = []
    shiftListCount = {}
    resonanceTracker = []
    shiftListLabels = []
    shiftListLabelsDict = {}
    atomCount = 0

    for chain in chains:
      for residue in chain.residues:
        for atom in residue.atoms:
          atomCount += 1
          if atom.atomSet:
            for resonanceSet in atom.atomSet.resonanceSets:
              for resonance in resonanceSet.resonances:
                
                if resonance in resonanceTracker:
                  continue
                
                resonanceTracker.append(resonance)
                
                for shift in resonance.shifts:
                  
                  if shift.parentList not in shiftLists:
                    shiftLists.append(shift.parentList)
                    shiftListCount[shift.parentList] = 0
                    
                  shiftListCount[shift.parentList] += 1
    
    for shiftList in shiftLists:
      percentage = (shiftListCount[shift.parentList] * 100.0) / atomCount
      shiftListLabel = "%d:%s (%.1f%%)" % (shiftList.serial,shiftList.name,percentage)
      shiftListLabels.append(shiftListLabel)
      shiftListLabelsDict[shiftListLabel] = shiftList
    
    chainLabelShiftListDict[chainLabel] = [shiftListLabels,shiftListLabelsDict]
    
    return shiftListLabels
    
  def body(self, master):

    #
    # Setup header 
    #

    row = 0
    self.columnSpan = 3

    label = Label(master, text = 'Select molecular system, chains and chemical shift list to be analyzed:')
    label.grid(row=row, column=0, columnspan = self.columnSpan, sticky=Tkinter.EW)

    row += 1

    self.molSysSelect = PulldownMenu(master, entries = self.molSysLabelList, callback = self.setupChains)
    self.molSysSelect.grid(row=row, column=0, sticky=Tkinter.W)

    self.chainSelect = PulldownMenu(master, entries = self.chainLabelList, callback = self.setupShiftList)
    self.chainSelect.grid(row=row, column=1, sticky=Tkinter.W)

    self.shiftListSelect = PulldownMenu(master, entries = self.shiftListLabels, callback = self.setupPercentageFrame, do_initial_callback = False)
    self.shiftListSelect.grid(row=row, column=2, sticky=Tkinter.W)

    row += 1

    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, columnspan = self.columnSpan, sticky = Tkinter.EW)

    row += 1

    master.grid_rowconfigure(row,weight = 1)
    for i in range(self.columnSpan):
      master.grid_columnconfigure(i,weight = 1)

    self.percentageFrame = ScrolledFrame(master, height = 180, doExtraConfig = False)
    self.percentageFrame.grid(row=row, columnspan = self.columnSpan, sticky=Tkinter.NSEW)
    self.percentageFrameRow = row
    self.percentageFrameMaster = master

    row += 1

    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, columnspan = self.columnSpan, sticky = Tkinter.EW)

    row += 1

    master.grid_rowconfigure(row,weight = 1)
    self.resultsFrame = ScrolledFrame(master, height = 180, doExtraConfig = False)
    self.resultsFrame.grid(row=row, columnspan = self.columnSpan, sticky=Tkinter.NSEW)
    self.resultsFrameRow = row
    self.resultsFrameMaster = master

    row += 1

    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, columnspan = self.columnSpan, sticky = Tkinter.EW)
    
    row += 1

    texts = [ 'Recalculate', 'Write to file']
    commands = [ self.recalc, self.writeFile ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = 'Exit', help_url=self.help_url)
    buttons.grid(row=row, columnspan = self.columnSpan, sticky=Tkinter.EW)

    self.shiftListSelect.callback(0,self.shiftListLabels[0])

  def setupChains(self,listIndex,molSysLabel):
  
    molSystems = self.molSysLabelDict[molSysLabel]
    
    if self.molSystems == molSystems:
      return
      
    self.molSystems = molSystems

    self.chainLabelList = self.molSysRelationDict[molSysLabel][0]
    self.chainLabelDict = self.molSysRelationDict[molSysLabel][1]
    self.chainLabelShiftListDict = self.molSysRelationDict[molSysLabel][2]
    
    if hasattr(self,'chainSelect'):
      self.chains = None      
      self.chainSelect.callback(0,self.chainLabelList[0])

  def setupShiftList(self,listIndex,chainLabel):
    
    chains = self.chainLabelDict[chainLabel]

    if self.chains == chains:
      return
      
    self.chains = chains
    
    #
    # Reset the pulldown menu
    #
    
    if hasattr(self,'chainSelect'):
      self.chainSelect.replace(self.chainLabelList)
    
    #
    # Set up chemical shift lists
    #

    self.shiftListLabels = self.chainLabelShiftListDict[chainLabel][0]
    self.shiftListLabelsDict = self.chainLabelShiftListDict[chainLabel][1]
    
    if hasattr(self,'percentageFrame'):
      self.shiftList = None
      self.shiftListLabel = None
      self.setupPercentageFrame(0,self.shiftListLabels[0])
    
  def setupPercentageFrame(self,listIndex,shiftListLabel):
    
    if self.shiftListLabel == shiftListLabel:
      return
      
    self.shiftList = self.shiftListLabelsDict[shiftListLabel]
    self.shiftListLabel = shiftListLabel

    #
    # Reset chemical shift groups...
    #

    self.shiftGroupsList = []

    for i in range(1,10):
      shiftGroupLabel = 'Group_%s' % i
      self.shiftGroupsList.append(shiftGroupLabel)

    #
    # Reset the pulldown menu
    #
    
    if hasattr(self,'shiftListSelect'):
      if hasattr(self.shiftListSelect, 'entries') and tuple(self.shiftListLabels) != tuple(self.shiftListSelect.entries): # HACK
        # otherwise can get infinite recursion but even with above could
        # if no self.shiftListLabels because of "None" entry in PulldownMenu
        self.shiftListSelect.replace(self.shiftListLabels)
      self.shiftListSelect.setSelected(shiftListLabel)

    #
    # Reset frame...
    #
    
    self.percentageFrame.destroy()
    self.percentageFrame = ScrolledFrame(self.percentageFrameMaster, height = 180, doExtraConfig = False)
    self.percentageFrame.grid(row=self.percentageFrameRow, columnspan = self.columnSpan, sticky=Tkinter.NSEW)
    
    #
    # Recalculate results
    #
    
    self.results = makeChemShiftSelections(self.parent,self.chains,self.shiftList)
    resultKeys = self.results.keys()
    resultKeys.sort()
    
    #
    # Set up the information
    #
    
    frameRow = 0
    frame = self.percentageFrame.frame

    self.shiftGroupObjects = {}
    
    for resultKey in resultKeys:
    
      allAtoms = self.results[resultKey][0]
      shiftAtoms = self.results[resultKey][1]
      
      percentage = shiftAtoms*100.0/allAtoms
      
      label = Label(frame, text = resultKey)
      label.grid(row = frameRow, column = 0, sticky = Tkinter.W)
      
      label = Label(frame, text = "%.1f%%" % percentage)
      label.grid(row = frameRow,column = 1, sticky = Tkinter.W)
    
      label = Label(frame, text = "(%d/%d)" % (shiftAtoms,allAtoms))
      label.grid(row = frameRow,column = 2, sticky = Tkinter.W)
      
      self.shiftGroupObjects[resultKey] = PulldownMenu(frame, entries = self.shiftGroupsList)
      self.shiftGroupObjects[resultKey].grid(row = frameRow,column = 3, sticky = Tkinter.E)
    
      frameRow += 1
              
    return True
 
  def recalc(self):
    
    groups = {}
    groupsInfo = {}
    
    for resultKey in self.shiftGroupObjects:
      group = self.shiftGroupObjects[resultKey].getSelected()

      if not groups.has_key(group):
        groups[group] = []
        groupsInfo[group] = [0,0]
      
      (allAtoms,shiftAtoms) = self.results[resultKey]
      groupsInfo[group][0] += allAtoms
      groupsInfo[group][1] += shiftAtoms
     
      groups[group].append(resultKey)
    
    #
    # Reset frame...
    #
    
    self.resultsFrame.destroy()
    self.resultsFrame = ScrolledFrame(self.resultsFrameMaster, height = 180, doExtraConfig = False)
    self.resultsFrame.grid(row=self.resultsFrameRow, columnspan = self.columnSpan, sticky=Tkinter.NSEW)

    #
    # Set info in lower frame
    #
    
    frameRow = 0
    frame = self.resultsFrame.frame

    groupKeys = groups.keys()
    groupKeys.sort()
    
    self.groupText = {}
    
    for group in groupKeys:
          
      percentage = groupsInfo[group][1]*100.0/groupsInfo[group][0]

      label = Label(frame, text = group)
      label.grid(row = frameRow, column = 0, sticky = Tkinter.W)
      
      label = Label(frame, text = groups[group][0])
      label.grid(row = frameRow, column = 1, sticky = Tkinter.W)

      label = Label(frame, text = "%.1f%%" % percentage)
      label.grid(row = frameRow,column = 2, sticky = Tkinter.W)
    
      label = Label(frame, text = "(%d/%d)" % (groupsInfo[group][1],groupsInfo[group][0]))
      label.grid(row = frameRow,column = 3, sticky = Tkinter.W)
      
      self.groupText[group] = "%-10s %-20s %8.1f%% (%d/%d)%s" % (group,
                                                                 groups[group][0],
                                                                 percentage,
                                                                 groupsInfo[group][1],
                                                                 groupsInfo[group][0],
                                                                 newline)
      
      frameRow += 1
      
      for otherResultKey in groups[group][1:]:
        label = Label(frame, text = otherResultKey)
        label.grid(row = frameRow, column = 2, sticky = Tkinter.W)

        self.groupText[group] += "%-10s %-20s%s" % ('',otherResultKey,newline)

        frameRow += 1
    
    return True
  
  def writeFile(self):
  
    filePopup = FormatFilePopup(self,component = 'text')
  
    if filePopup.fileSelected:
      
      groupKeys = self.groupText.keys()
      groupKeys.sort()
      
      fout = open(filePopup.file,'w')
      fout.write("Project: %s" % self.project.name + newline)
      fout.write("Molsystems: %s" % self.molSysSelect.getSelected() + newline)
      fout.write("Chains: %s" % self.chainSelect.getSelected() + newline)
      fout.write("Shiftlist: %s" % self.shiftListSelect.getSelected() + newline * 2)
      
      for group in groupKeys:

        fout.write(self.groupText[group])
  
      print "Wrote file %s..." % filePopup.file
      
  def apply(self):
    
    return True
