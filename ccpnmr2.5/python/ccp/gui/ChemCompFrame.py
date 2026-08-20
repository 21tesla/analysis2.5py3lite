"""
======================COPYRIGHT/LICENSE START==========================

ChemCompFrame.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in license/LGPL.license
 
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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================
"""
import os, glob

import Tkinter


from memops.universal.Io import joinPath
from memops.general.Io import modifyPackageLocators

from memops.gui.ButtonList import ButtonList
from memops.gui.CheckButtons import CheckButtons
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.Separator import Separator
from memops.gui.ProgressBar import ProgressBar
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.MessageReporter import showError, showYesNo

from ccp.general.Io import getDataPath, getChemComp, getCcpCodeFromXmlFile

from ccp.general.ChemCompOverview import chemCompStdDict

import ccp.api.molecule.ChemComp

molTypeList = ('protein', 'DNA', 'RNA', 'carbohydrate', 'other')
chemCompList = ('ChemComp', 'ChemCompVar')


# Fake chemCompHead class for ScrolledTable
class ChemCompInfo:

  def __init__(self,molType,ccpCode):
  
    self.molType = molType
    self.ccpCode = ccpCode
    self.chemComp = None
   
  def setInfoFromDict(self,dictItem):

    self.code1Letter = dictItem[0]
    self.cifCode = dictItem[1]
    self.name = dictItem[2]
    # TODO this item not provided any more...
    #self.formula = dictItem[3]
  
  def setInfoFromChemComp(self,chemComp):
  
    self.chemComp = chemComp
  
    # Not fully correct, but gives an idea...
    refChemCompVar = chemComp.findFirstChemCompVar(isDefaultVar = True)
    if not refChemCompVar:
      refChemCompVar = chemComp.findFirstChemCompVar()

    formula = refChemCompVar.formula

    self.code1Letter = chemComp.code1Letter
    
    cifNamingSystem = chemComp.findFirstNamingSystem(name = 'CIF')
    if cifNamingSystem:
      self.cifCode = cifNamingSystem.findFirstChemCompSysName()
    else:
      self.cifCode = None

    self.name = chemComp.name
    self.formula = formula

class ChemCompFrame(Frame):

  chemCompClasses = {}

  def __init__(self, parent, project, path = None,
               molTypeEntries = None, chemCompEntries = None,
               selectedChemComps = None, selectLinking = None,
               *args, **kw):

    Frame.__init__(self, parent, *args, **kw)

    self.project = project
    self.molTypeEntries = molTypeEntries
    self.chemCompEntries = chemCompEntries
    self.selectLinking = selectLinking

    if (not path):
      path = getDataPath()

    self.path = path
      
    self.chemCompInfoDict = {}
    self.chemCompInfoList = []
    self.chemCompDownload = False
    self.chem_comps_shown = {}

    for entry in chemCompList:
      self.chemCompClasses[entry] = getattr(ccp.api.molecule.ChemComp, entry)

    self.grid_columnconfigure(0, weight=1)

    row = 0
    
    if (molTypeEntries is None):
      headerText = "Show residues (select molecular type(s)):"
    else:
      headerText = "Show %s residues:" % (str(molTypeEntries))
    
    #
    #
    # TODO TODO: HERE need to do some niftier stuff for displaying!
    #
    #
    
    headerTextWidget = Label(self, text = headerText)
    headerTextWidget.grid(row=row, column=0, sticky=Tkinter.W)
    row = row + 1

    if (molTypeEntries is None):
      self.mol_type_buttons = CheckButtons(self, entries=molTypeList,
                                select_callback=self.updateTables)
      self.mol_type_buttons.grid(row=row, column=0, sticky=Tkinter.EW)
      row = row + 1
    else:
      self.mol_type_buttons = None

    #
    # The chemComps to display...
    #
    
    self.showLocalText = 'Show local'
    self.showWebText = 'Show available via web'
    self.display_buttons = CheckButtons(self, entries=[self.showLocalText,self.showWebText],
                                              select_callback=self.updateTables,
                                              selected = [True,False])
    self.display_buttons.grid(row=row, column=0, sticky=Tkinter.EW)
    row = row + 1

    self.grid_rowconfigure(row, weight=2)
    headings = ('number', 'show details', 'molType', 'ccpCode', 'code1Letter', 'cifCode', 'name')
    editWidgets = 7 * [ None ]
    editGetCallbacks = [ None, self.toggleShow, None, None, None, None, None ]
    editSetCallbacks = 7 * [ None ]
    self.chem_comp_table = ScrolledMatrix(self, headingList=headings,
                                  editWidgets=editWidgets,
                                  editGetCallbacks=editGetCallbacks,
                                  editSetCallbacks=editSetCallbacks)
    self.chem_comp_table.grid(row=row, column=0, sticky=Tkinter.NSEW)

    row = row + 1
    texts = [ 'Show all in details window', 'Clear details window' ]
    commands = [ self.showAll, self.showNone ]
    buttons = ButtonList(self, texts=texts, commands=commands)
    buttons.grid(row=row, column=0, sticky=Tkinter.EW)
    
    row = row + 1

    separator = Separator(self,height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row=row, column=0, sticky=Tkinter.EW)

    row = row + 1
    
    headerTextWidget = Label(self, text = "Select the residue variant:")
    headerTextWidget.grid(row=row, column=0, sticky=Tkinter.W)

    row = row + 1
    if (chemCompEntries is None):
      self.chem_comp_buttons = CheckButtons(self, entries=chemCompList,
                                 selected=('ChemComp',),
                                 select_callback=self.updateChemCompVarTable)
      self.chem_comp_buttons.grid(row=row, column=0, sticky=Tkinter.EW)
      row = row + 1
    else:
      self.chem_comp_buttons = None

    self.grid_rowconfigure(row, weight=1)
    headings = ('number', 'molType', 'ccpCode', 'linking', 'descriptor', 'molecularMass', 'formula', 'nonStereoSmiles', 'stereoSmiles')
    self.chem_comp_var_table = ScrolledMatrix(self, headingList=headings)
    self.chem_comp_var_table.grid(row=row, column=0, sticky=Tkinter.NSEW)
    self.chem_comp_var_headings = headings[1:]
    
    if selectedChemComps:
      for chemComp in selectedChemComps:
        key = (chemComp.molType, chemComp.ccpCode)
        self.chem_comps_shown[key] = 1

    self.updateTables()

  def toggleShow(self, chemCompInfo):

    key = (chemCompInfo.molType, chemCompInfo.ccpCode)
    if (self.chem_comps_shown.has_key(key)):
      del self.chem_comps_shown[key]
    else:
      self.chem_comps_shown[key] = 1

    self.updateTables()

  def findChemComps(self):

    self.chemCompInfoDict = {}
    self.chemCompInfoList = []
    
    #
    # Get the molTypes
    #
    
    if (self.molTypeEntries is None):
      selected = self.mol_type_buttons.getSelected()
    else:
      selected = self.molTypeEntries
    
    #
    # Find out where to get the chemComp info from
    #  
    
    searchTypes = self.display_buttons.getSelected()
    
    #
    # Now get the info...
    #
    
    for molType in selected:
    
      if self.showLocalText in searchTypes:
    
        #
        # Find the files that are locally available
        #
        # Note Wim: taken out central storage, no point in adding that in (is in overview dictionary).
        # If people make their own chemComps they should add them to the data/ccp/molecule/chemComp directory.
        #
        
        ccpCodes = []
        
        localPath = self.project.findFirstRepository(name='userData').url.path
        
        for ccPath in (self.path,localPath):
  
          chemCompFileSearchString = "%s+*.xml" % molType
          chemCompFileSearchPath = os.path.join(ccPath,'ccp','molecule','ChemComp')
  
          chemCompFileNameMatches = glob.glob(joinPath(chemCompFileSearchPath,chemCompFileSearchString))
          
          #
          # Add them to the dictionaries...
          #
          
          
          for chemCompFile in chemCompFileNameMatches:
            if os.path.isfile(joinPath(chemCompFileSearchPath, chemCompFile)) and (chemCompFile[-4:] == '.xml'):
              ccpCodes.append(getCcpCodeFromXmlFile(joinPath(chemCompFileSearchPath, chemCompFile)))       
        
        for ccpCode in ccpCodes:
          
          chemCompKey = None
          
          # First try dictionary (quicker)
          if molType != 'other':
            if chemCompStdDict[molType].has_key(ccpCode):
              chemCompKey = (molType,ccpCode)
              chemCompDict = chemCompStdDict[molType]
          
          else:
            # molType == 'other'
            from ccp.general.ChemCompOverview import chemCompOverview
            if chemCompOverview[molType].has_key(ccpCode):
              chemCompKey = (molType,ccpCode)
              chemCompDict = chemCompOverview[molType]
          
          # If found, set info
          if chemCompKey:

            chemCompInfo = ChemCompInfo(molType,ccpCode)
            chemCompInfo.setInfoFromDict(chemCompDict[ccpCode])
            self.chemCompInfoList.append(chemCompKey)
            self.chemCompInfoDict[chemCompKey] = chemCompInfo
                  
          # Else try to get from local directory- TODO is this obsolete?!
          else:  
            chemCompKey = (molType,ccpCode)
            
            chemComp = self.project.findFirstChemComp(molType = molType, ccpCode = ccpCode)
            
            if not chemComp:
              chemComp = getChemComp(self.project, molType, ccpCode, 
                                     download=False)
            
            if (chemComp):
              chemCompInfo = ChemCompInfo(molType,ccpCode)
              chemCompInfo.setInfoFromChemComp(chemComp)
              self.chemCompInfoList.append(chemCompKey)
              self.chemCompInfoDict[chemCompKey] = chemCompInfo
              
      if self.showWebText in searchTypes:
       
        if molType == 'other':
          from ccp.general.ChemCompOverview import chemCompOverview
          overviewDict = chemCompOverview[molType]
        else:
          overviewDict = chemCompStdDict[molType]
         
        ccpCodes = overviewDict.keys()
        ccpCodes.sort()
       
        for ccpCode in ccpCodes:
         
          chemCompKey = (molType,ccpCode)

          if chemCompKey not in self.chemCompInfoList:
            self.chemCompDownload = True
            chemCompInfo = ChemCompInfo(molType,ccpCode)
            chemCompInfo.setInfoFromDict(overviewDict[ccpCode])
            self.chemCompInfoList.append(chemCompKey)
            self.chemCompInfoDict[chemCompKey] = chemCompInfo

  def updateTables(self, *extra):

    self.updateChemCompTable()
    self.updateChemCompVarTable()

  def updateChemCompTable(self):

    self.findChemComps()

    textMatrix = []
    n = 0
    
    chemCompInfos = []
    
    for (molType, ccpCode) in self.chemCompInfoList:
      n = n + 1
      isShown = self.chem_comps_shown.has_key((molType, ccpCode))
      chemCompInfo = self.chemCompInfoDict[(molType, ccpCode)]
      chemCompInfos.append(chemCompInfo)
      text = [n, self.booleanString(isShown), molType, ccpCode, chemCompInfo.code1Letter, chemCompInfo.cifCode, chemCompInfo.name]
      textMatrix.append(text)

    self.chem_comp_table.update(objectList=chemCompInfos, textMatrix=textMatrix)

  def updateChemCompVarTable(self, *extra):

    textMatrix = []
    allChemComps = []
    
    if (self.chemCompEntries is None):
      selected = self.chem_comp_buttons.getSelected()
    else:
      selected = self.chemCompEntries
      
    numberToLoad = self.numberChemCompsToLoad()
    if (selected and (numberToLoad > 1)):
      use_progress_bar = True
      p = ProgressBar(self, text='loading', total=numberToLoad)
    else:
      use_progress_bar = False
    try:
      n = 0
      
      for (molType, ccpCode) in self.chemCompInfoList:
        chemCompInfo = self.chemCompInfoDict[(molType,ccpCode)]
        isShown = self.chem_comps_shown.has_key((molType, ccpCode))
        has_been_loaded = False
        if (isShown):
          chemComp = chemCompInfo.chemComp 
          if not chemComp:
            chemComp = getChemComp(self.project, molType, ccpCode, 
                                   showError=showError)
            
          chemCompAndVars = [chemComp] + list(chemComp.chemCompVars)# this loads data, if needed
          for chemCompOrVar in chemCompAndVars:
            if (self.chemCompShown(chemCompOrVar)):
              n = n + 1
              text = [n]
              for heading in self.chem_comp_var_headings:
                if hasattr(chemCompOrVar,heading):
                  text.append(getattr(chemCompOrVar, heading))
                elif hasattr(chemCompOrVar,'chemComp') and hasattr(chemCompOrVar.chemComp,heading):
                  text.append(getattr(chemCompOrVar.chemComp, heading))
                else:
                  text.append("n/a")
              textMatrix.append(text)
              allChemComps.append(chemCompOrVar)
        if (use_progress_bar and has_been_loaded):
          p.increment()
    finally:
      if (use_progress_bar):
        p.destroy()

    self.chem_comp_var_table.update(objectList=allChemComps, textMatrix=textMatrix)

  def numberChemCompsToLoad(self):
    # TODO THIS NEEDS TO BE UPDATED!
    n = 0
    for (molType,ccpCode) in self.chemCompInfoList:
     
      chemCompInfo = self.chemCompInfoDict[(molType,ccpCode)]
      if chemCompInfo.chemComp:
        continue
      isShown = self.chem_comps_shown.has_key((molType, ccpCode))
      if (isShown):
        n = n + 1

    return n

  def chemCompShown(self, chemComp):

    returnStatus = False
    
    if (self.chemCompEntries is None):
      selected = self.chem_comp_buttons.getSelected()
    else:
      selected = self.chemCompEntries

    for s in selected:
      if (isinstance(chemComp, self.chemCompClasses[s])):
        if self.selectLinking and hasattr(chemComp,'linking'):
          if chemComp.linking == self.selectLinking:
            returnStatus = True
        else:
          returnStatus = True

    return returnStatus

  def showAll(self):

    loadAll = False
    if self.chemCompDownload:
      if (showYesNo('Show all variants','Are you sure you want to display all variants? This will also download the XML files for all chemComps that are not stored locally!')):
        loadAll = True
    else:
      loadAll = True

    if loadAll:
      for key in self.chemCompInfoList:
        self.chem_comps_shown[key] = 1

      self.updateTables()

      self.chemCompDownload = False
    
  def showNone(self):

    for key in self.chemCompInfoList:
      if (self.chem_comps_shown.has_key(key)):
        del self.chem_comps_shown[key]

    self.updateTables()

  def booleanString(self, boolean):
 
    if (boolean):
      return 'yes'
    else:
      return 'no'

  def getSelectedChemComp(self):

    return self.chem_comp_var_table.currentObject

if (__name__ == '__main__'):

  from memops.api.Implementation import MemopsRoot
  from memops.gui.Button import Button
  from memops.gui.MessageReporter import showInfo

  from ccp.general.Io import getDataPath

  project = MemopsRoot(name='p')

  root = Tkinter.Tk()
  root.top = root

  root.grid_rowconfigure(0, weight=1)
  root.grid_columnconfigure(0, weight=1)
  
  path = getDataPath()

  frame = ChemCompFrame(root, project, path = path)
  frame.grid(sticky=Tkinter.NSEW)

  def getSelected():

    chemComp = frame.getSelectedChemComp()
    if (chemComp):
      if hasattr(chemComp,'chemComp'):
        ccType = 'chemCompVar'
        chemComp = chemComp.chemComp
      else:
        ccType = 'chemComp'
        
      showInfo('ChemComp', 'Selected %s = (%s, %s)' % (ccType,chemComp.molType, chemComp.ccpCode), parent=root)
    else:
      showInfo('ChemComp', 'No chemComp selected', parent=root)

  button = Button(root, text='Get selected', command=getSelected)
  button.grid()

  root.mainloop()

