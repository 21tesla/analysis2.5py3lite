
"""
======================COPYRIGHT/LICENSE START==========================

AutoAssignIOCyclePopup.py: GUI window to do an export/import cycle for AutoAssign

Copyright (C) 2006 Wim Vranken (European Bioinformatics Institute)

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
import Tkinter, time

from memops.api import Implementation

from memops.universal.Io import joinPath
from memops.universal.Util import returnFloat, returnInt

from memops.gui.BasePopup import BasePopup
from memops.gui.Label import Label
from memops.gui.Entry import Entry
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.MessageReporter import showInfo, showError
from memops.gui.Separator import Separator

from ccpnmr.format.gui.ObjectEditPopup import createSelection
from ccpnmr.format.gui.FormatFilePopup import FormatFilePopup

from ccpnmr.format.general.Io import getHelpUrlDir
from ccpnmr.format.general.Constants import dateTimeStamp_kw, ioCycleTag_kw

from ccpnmr.format.converters.AutoAssignFormat import AutoAssignFormat
from ccp.format.autoAssign.chemShiftsIO import AutoAssignChemShiftFile

from ccp.general.Util import findChemAtomSysName

class AutoAssignIOCyclePopup(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'AutoAssignIOCycle.html')

  def __init__(self, parent, project):

    self.project = project
    self.guiParent = parent
      
    self.autoAssignFormat = AutoAssignFormat(project, guiParent = parent)
    
    self.format = 'autoAssign'

    self.chainStampSep = '-'

    self.dateTimeFlag = time.strftime("%Y%m%d.%H:%M") 
    
    self.getPreviousDateTimeStamps()

    self.getChainsAndChainStamps()
    
    self.importFileName = self.exportFileName = None
    
    BasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + 'AutoAssign export/import cycle.', modal=True, transient=True)

  def getPreviousDateTimeStamps(self):
    
    self.dateTimeStamps = []
    appData = self.project.findAllApplicationData(application = self.format, keyword = dateTimeStamp_kw)
    for appDatum in appData:
      self.dateTimeStamps.append(appDatum.value)

  def getChainsAndChainStamps(self):
    
    self.chains = []  
    self.chainDateTimeStamps = []
    self.chainDateTimeStampDict = {}
    
    for molSys in self.project.sortedMolSystems():
      for chain in molSys.sortedChains():
        self.chains.append(chain)
        
        appData = chain.findFirstApplicationData(application = self.format, keyword = dateTimeStamp_kw)
        if appData and appData.value in self.dateTimeStamps:
          (tlist,tdict) = createSelection([chain])
          cdtsTag = "%s%s%s" % (tlist[0],self.chainStampSep,appData.value)
          self.chainDateTimeStamps.append(cdtsTag)
          self.chainDateTimeStampDict[cdtsTag] = chain

  def body(self, master):
    
    self.geometry('600x400')

    #
    # Quick check
    #
    
    if not self.chains:
      showError("No chains","No chains available - cannot use export/import cycle.")
      self.destroy()
    
    #
    # Set it all up
    #
    
    columnspan = 2
    
    row = 0
    
    label = Label(master, text= "AutoAssign export/import cycle.")
    label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)
    
    row += 1
    
    separator = Separator(master,height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)

    row += 1
    
    label = Label(master, fg = 'red', text= "Popup to export %s data from the CCPN data model,\nrun %s, then re-import the output." % (self.format,self.format))
    label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)   
     
    #
    # Make a break...
    #

    row += 1
    
    separator = Separator(master,height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)

    #
    # Set up export file info
    #
    
    row += 1
    
    label = Label(master, text= "Export menu (using date/time label '%s')" % self.dateTimeFlag)
    label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)
    

    row += 1
    
    (chainList,self.chainDict) = createSelection(self.chains)
    
    label = Label(master, text= "Select chain to export:") 
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.chainSelect = PulldownMenu(master,entries = chainList)
    self.chainSelect.grid(row=row, column=1, sticky=Tkinter.W)

    row += 1
    
    label = Label(master, text= "Export %s project file:" % self.format) 
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.selectExportFileButton = Tkinter.Button(master, text = 'Select export file', command = (lambda : self.selectExportProjectFile()))
    self.selectExportFileButton.grid(row=row, column=1, sticky=Tkinter.W)

    row += 1

    self.exportButton = Tkinter.Button(master, text = 'Export', command = (lambda : self.doExport()))
    self.exportButton.grid(row=row, column=0, columnspan = columnspan,  sticky=Tkinter.EW)
        
    #
    # Make a break...
    #

    row += 1
    
    separator = Separator(master,height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)
    
    #
    # Set up import file info
    #
    
    row += 1
    
    label = Label(master, text= "Re-import menu")
    label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)

    row += 1
    
    #
    # Select the right chain with the right date/time flag...
    #
    
    label = Label(master, text= "Select chain with correct date/time flag:") 
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.chainDateTimeSelect = PulldownMenu(master,entries = self.chainDateTimeStamps)
    self.chainDateTimeSelect.grid(row=row, column=1, sticky=Tkinter.W)    
    # TODO UPDATE THIS WHEN EXPORT BUTTON PRESSED AND GETTING OK FROM EXPORT ITSELF!    
    # Probably also need just a message if no importDateTimeStamp available...    

    row += 1
    
    label = Label(master, text= "Import %s output file:" % self.format)
    label.grid(row=row, column=0, sticky=Tkinter.W)
    
    self.selectImportFileButton = Tkinter.Button(master, text = 'Select import file', command = (lambda : self.selectImportShiftFile()))
    self.selectImportFileButton.grid(row=row, column=1, sticky=Tkinter.W)

    row += 1

    self.importButton = Tkinter.Button(master, text = 'Import', command = (lambda : self.doImport()))
    self.importButton.grid(row=row, column=0, columnspan = columnspan,  sticky=Tkinter.EW)

    #
    # Make a break...
    #

    row += 1
    
    separator = Separator(master,height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)

    row += 1

    buttons = createDismissHelpButtonList(master, texts=[], commands=[], help_url=self.help_url)
    buttons.grid(row=row, columnspan = columnspan, column=0)

  def selectExportProjectFile(self):
    
    if self.exportFileName:
      fileName = self.exportFileName
    else:
      fileName = 'table.aat'
  
    popup = FormatFilePopup(self, file = fileName, component = 'project', format = self.format, title = 'Select name for %s project file to export.' % self.format)
    
    if popup.fileSelected:
     
      self.selectExportFileButton.config(text = popup.file)
      self.exportFileName = popup.file

  def selectImportShiftFile(self):
  
    if self.exportFileName:
      fileName = self.exportFileName
    else:
      fileName = 'autoAssign.out'
  
    popup = FormatFilePopup(self, file = fileName, component = 'shifts', format = self.format, title = 'Select name for %s output file to re-import.' % self.format)
    
    if popup.fileSelected:
     
      self.selectImportFileButton.config(text = popup.file)
      self.importFileName = popup.file
  
  def doExport(self):
    
    if self.exportFileName:
      chain = self.chainDict[self.chainSelect.getSelected()]
      returnValue = self.autoAssignFormat.writeProject(self.exportFileName, chain = chain, setTag = self.dateTimeFlag)
      if not returnValue:
        showError("No export file written","There were problems while exporting the %s project file."  % self.format)
      else:
        showInfo("File written","File written successfully")
        cdtsTag = "%s-%s" % (self.chainSelect.getSelected(),self.dateTimeFlag)
        if not cdtsTag in self.chainDateTimeSelect.entries:
          if "<None>" in self.chainDateTimeSelect.entries:
            self.chainDateTimeSelect.entries.pop(self.chainDateTimeSelect.entries.index("<None>"))
          self.chainDateTimeSelect.replace([cdtsTag] + self.chainDateTimeSelect.entries)
          self.chainDateTimeStampDict[cdtsTag] = chain
        
    else:
       showError("No export file defined","Please define a name for the %s project file to export to."  % self.format)
  
  def doImport(self):

    if self.importFileName:
      cdtsTag = self.chainDateTimeSelect.getSelected()

      if not self.chainDateTimeStampDict.has_key(cdtsTag):
        showError("No import tag defined","No chain with date/time stamp defined - cannot re-import.")
      
      else:
        chain = self.chainDateTimeStampDict[cdtsTag]
        (chainText,curDateTimeStamp) = cdtsTag.split(self.chainStampSep)
        
        #
        # Get relevant info from data model (not immediately necessary but as example)
        #
        # Note that could in principle use the normal peakNum tag, but dangerous in case
        # multiple exports were done... this is more laborious though.
        #
        
        peakLists = []

        rootPeakListTag = str((curDateTimeStamp,'isRoot'))
        rootPeakNumToPeak = {}
        
        for nmrExp in self.project.currentNmrProject.sortedExperiments():
          for ds in nmrExp.sortedDataSources():
            for peakList in ds.sortedPeakLists():
              appData = peakList.findAllApplicationData(application = self.format, keyword = dateTimeStamp_kw)
              for appDatum in appData:
                if appDatum and appDatum.value == curDateTimeStamp:
                  peakLists.append(peakList)
                  if peakList.findFirstApplicationData(application = self.format, keyword = ioCycleTag_kw, value = rootPeakListTag):
                    for peak in peakList.sortedPeaks():
                      appData = peak.findFirstApplicationData(application = self.format, keyword = ioCycleTag_kw)
                      if appData:
                        (curTag,peakNum) = eval(appData.value)
                        if curTag == curDateTimeStamp:
                          rootPeakNumToPeak[peakNum] = peak
        
        #
        # Now get the actual chemical shift info from the AutoAssign output file 
        #
        
        autoAssignChemShiftFile = AutoAssignChemShiftFile(self.importFileName)
        autoAssignChemShiftFile.read()
        
        #
        # Set the mapping between the chain residues and the seqCode in the chem shift file
        # 
        
        residues = list(chain.sortedResidues())
        seqCodeToResidue = {}
        for i in range(0,len(residues)):
          residue = residues[i]
          (seqCode,code1Letter) = autoAssignChemShiftFile.seqCodes[i]
          seqCode = returnInt(seqCode)
          seqCodeToResidue[seqCode] = residue
        
        for rawChemShift in autoAssignChemShiftFile.chemShifts:
          
          peak = None
          residue = None
          prevResidue = None
          
          seqCode = rawChemShift.seqCode
          atomName = rawChemShift.atomName
          #allValues = rawChemShift.allValues          
          #rawChemShift.spinSystemId
          
          if seqCodeToResidue.has_key(seqCode):
            residue = seqCodeToResidue[seqCode]
          else:
            # THIS SHOULD NEVER HAPPEN!
            print "  Error: no match for seqCode %d while re-importing project." % seqCode
            continue
          
          #
          # Set info on residue/atom level
          #
          
          atom = self.findMatchingAtom(residue,atomName)          
          self.autoAssignFormat.setSeqAssignTag(atom,rawChemShift.value,AppDataClass = Implementation.AppDataFloat)
          
          #
          # Set info on peak level
          #
          
          if rawChemShift.peakId:
            (peakNum,rootName) = rawChemShift.peakId.split('.') # TODO set this somewhere?
            peakNum = returnInt(peakNum)
            if rootPeakNumToPeak.has_key(peakNum):
              peak = rootPeakNumToPeak[peakNum]
              self.autoAssignFormat.setSeqAssignTag(peak,str((chain.code,residue.seqId)))

    else:
       showError("No import file defined","Please define a name for the %s shift output file to import." % self.format)  

  def findMatchingAtom(self,residue,atomName):

    atom = residue.findFirstAtom(name = atomName)
    if not atom:
      # Rough search but should be OK for bb atoms
      namingSystem = residue.chemCompVar.chemComp.findFirstNamingSystem(name = 'XPLOR')
      chemAtomSysName = findChemAtomSysName(namingSystem,{'sysName': atomName})
      atom = residue.findFirstAtom(name = chemAtomSysName.atomName)

    return atom
              
  def apply(self):
    
    if not 0:
      showError("No root spectrum","Can't continue: need a root spectrum...")
      return False
         
    return True
