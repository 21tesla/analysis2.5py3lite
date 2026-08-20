"""
======================COPYRIGHT/LICENSE START==========================

BacusPopup.py: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
import Tkinter

from memops.editor.BasePopup import BasePopup
from memops.gui.LabelFrame   import LabelFrame
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Entry        import Entry
from memops.gui.Button       import Button
from memops.gui.Label        import Label
from memops.gui.Util         import createDismissHelpButtonList
from memops.gui.FileSelect   import FileType
from memops.gui.FileSelectPopup import FileSelectPopup

from ccpnmr.clouds.ResonanceIdentification import run3dBacus, run2dBacus

from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes
from ccpnmr.analysis.core.PeakBasic       import copyPeakListNew

class BacusPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

  
    self.project    = parent.project
    self.hNoesy    = None
    self.nNoesy     = None
    self.cNoesy     = None
    self.shiftList  = None

    BasePopup.__init__(self, parent, title="BACUS control module", **kw)
  
  def body(self, guiParent):
    
    guiParent.grid_columnconfigure(0, weight=1)

    row = 0
    guiParent.grid_rowconfigure(row, weight=1)
    inputFrame = LabelFrame(guiParent, text='Inputs:')
    inputFrame.grid(row=row, column=0, sticky=Tkinter.NSEW)
 
    label = Label(inputFrame, text='NOESY:')
    label.grid(row=0, column=0, sticky=Tkinter.NW)
    
    self.hNoesyPulldown = PulldownMenu(inputFrame,callback=self.setHNoesy)
    self.hNoesyPulldown.grid(row=0, column=1, sticky=Tkinter.NW)

    label = Label(inputFrame, text='15N NOESY:')
    label.grid(row=1, column=0, sticky=Tkinter.NW)
    
    self.nNoesyPulldown = PulldownMenu(inputFrame,callback=self.setNNoesy)
    self.nNoesyPulldown.grid(row=1, column=1, sticky=Tkinter.NW)
    
    label = Label(inputFrame, text='13C NOESY:')
    label.grid(row=2, column=0, sticky=Tkinter.NW)
 
    self.cNoesyPulldown = PulldownMenu(inputFrame,callback=self.setCNoesy)
    self.cNoesyPulldown.grid(row=2, column=1, sticky=Tkinter.NW)
    
    label = Label(inputFrame, text='Shift List:')
    label.grid(row=3, column=0, sticky=Tkinter.NW)
 
    self.shiftListPulldown = PulldownMenu(inputFrame,callback=self.setShiftList)
    self.shiftListPulldown.grid(row=3, column=1, sticky=Tkinter.NW)
    
    label = Label(inputFrame, text='2d BACUS executable:')
    label.grid(row=4, column=0, sticky=Tkinter.NW)
    
    self.executableEntry = Entry(inputFrame, text='/home/tjs23/clouds/justin/SpiBacusMidge/bacus/bacus_tjs.exe')
    self.executableEntry.grid(row=4, column=1, sticky=Tkinter.NW)

    self.executableButton = Button(inputFrame, text='Choose file', command=self.chooseExecutable)
    self.executableButton.grid(row=5, column=1, sticky=Tkinter.EW)
    
    row += 1
    outputFrame = LabelFrame(guiParent, text='Output:')
    outputFrame.grid(row=row, column=0, sticky=Tkinter.NSEW)
 
    row += 1
    texts    = ['Run BACUS 2d','Run BACUS 3d']
    commands = [self.runBacus,self.runBacus3d]
    self.bottomButtons = createDismissHelpButtonList(guiParent,texts=texts,commands=commands,expands=0,help_url=None)
    self.bottomButtons.grid(row=row, column=0, sticky=Tkinter.EW)
    
    self.update()    
  
  
  def chooseExecutable(self):
  
    fileTypes = [FileType('Table', ['*.exe']), FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
                                      title = 'Choose 2d BACUS executable', dismiss_text = 'Cancel',
                                      selected_file_must_exist = False)

    fileName = fileSelectPopup.getFile() 
    if fileName:
      self.executableEntry.set(fileName)
  
  
  def update(self):
 
    # set defaults via pulldown callbacks

    names = [self.getPeakListlabel(x) for x in self.gethNoesys()]
    if names:
      self.hNoesyPulldown.setup(names,0)
    else:
      self.hNoesyPulldown.setup([],-1)

    names = [self.getPeakListlabel(x) for x in self.getNNoesys()]
    if names:
      self.nNoesyPulldown.setup(names,0)
    else:
      self.nNoesyPulldown.setup([],-1)

    names = [self.getPeakListlabel(x) for x in self.getCNoesys()]
    if names:
      self.cNoesyPulldown.setup(names,0)
    else:
      self.cNoesyPulldown.setup([],-1)
  
    names = ['Shift List %d' % x.serial for x in self.getShiftLists()]
    if names:
      self.shiftListPulldown.setup(names,0)
    else:
      self.shiftListPulldown.setup([],-1)
  
    self.updateButtons()
  
  
  def updateButtons(self):
    
    if ( self.cNoesy or self.nNoesy ) and self.shiftList:
      self.bottomButtons.buttons[1].enable()
    else:
      self.bottomButtons.buttons[1].disable()
      
      
    if self.hNoesy and self.executableEntry.get() and self.shiftList:
      self.bottomButtons.buttons[0].enable()
    else:
      self.bottomButtons.buttons[0].disable()
    
  
  def runBacus3d(self, fileRoot='CCPN1'):
  
    if ( self.cNoesy or self.nNoesy ) and self.shiftList:
      run3dBacus(fileRoot,self.shiftList,self.cNoesy,self.nNoesy)
 
  def runBacus(self, fileRoot='CCPN1'):
     
    executable = self.executableEntry.get()
    if self.hNoesy and self.shiftList and executable:
      run2dBacus(fileRoot, executable, self.hNoesy,  self.shiftList)


  def getPeakListlabel(self, peakList):
  
    return '%s:%s:%d' % (peakList.dataSource.experiment.name, peakList.dataSource.name, peakList.serial)
  
  
  def getShiftLists(self):
  
    shiftLists = []
    if self.project:
      nmrProject = self.project.currentNmrProject
      shiftLists = nmrProject.findAllMeasurementLists(className='ShiftList')
    
    return list(shiftLists)
  
  
  def setShiftList(self, index, name=None):
  
    shiftLists = self.getShiftLists()
    if shiftLists:
      self.shiftList = shiftLists[index]
    else:
      self.shiftList = None
    
  def gethNoesys(self):
  
    peakLists = []
    if self.project:
      for experiment in self.project.currentNmrProject.experiments:
        for spectrum in experiment.dataSources:
          isotopes = getSpectrumIsotopes(spectrum)
          if isotopes == ['1H', '1H']:
            for peakList in spectrum.peakLists:
              #if 'NOE' in experiment.name.upper():
              #  peakLists.insert(0, peakList)
              #else:
              peakLists.append(peakList)
    
    return peakLists
  
  def getNNoesys(self):
  
    peakLists = []
    if self.project:
      for experiment in self.project.currentNmrProject.experiments:
        for spectrum in experiment.dataSources:
          isotopes = getSpectrumIsotopes(spectrum)
          isotopes.sort()
          if isotopes == ['15N', '1H', '1H']:
            for peakList in spectrum.peakLists:
              #if 'NOE' in experiment.name.upper():
              #  peakLists.insert(0, peakList)
              #else:
              peakLists.append(peakList)
    
    return peakLists
  
  
  def getCNoesys(self):
  
    peakLists = []
    if self.project:
      for experiment in self.project.currentNmrProject.experiments:
        for spectrum in experiment.dataSources:
          isotopes = getSpectrumIsotopes(spectrum)
          isotopes.sort()
          if isotopes == ['13C', '1H', '1H']:
            for peakList in spectrum.peakLists:
              if 'NOE' in experiment.name.upper():
                peakLists.insert(0, peakList)
              else:
                peakLists.append(peakList)
    
    return peakLists
        
        
  def setHNoesy(self, index, name=None):
 
    peakLists = self.gethNoesys()
    if peakLists:
      self.hNoesy = peakLists[index]
    else:
      self.hNoesy = None
        
        
  def setNNoesy(self, index, name=None):
 
    peakLists = self.getNNoesys()
    if peakLists:
      self.nNoesy = peakLists[index]
    else:
      self.nNoesy = None

 
  def setCNoesy(self, index, name=None):
  
    peakLists = self.getCNoesys()
    if peakLists:
      self.cNoesy = peakLists[index]
    else:
      self.cNoesy = None
 

  def destroy(self):

    BasePopup.destroy(self)
