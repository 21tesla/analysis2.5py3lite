
"""
======================COPYRIGHT/LICENSE START==========================

EditPeak.py: Part of the CcpNmr Analysis program

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
from memops.general import Implementation

from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.PulldownList import PulldownList
from memops.gui.Text import Text

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import getPrimaryDataDimRef
from ccpnmr.analysis.core.PeakBasic import pickPeak, movePeak, setManualPeakIntensity
from ccpnmr.analysis.core.UnitConverter import pnt2ppm, unit_converter

class EditPeakPopup(BasePopup):
  """
  **Edit Position, Intensity & Details for a Peak**
  
  This popup window provides an means of editing peak information as an
  alternative to editing values in the main peak tables. This popup is also used
  to specify parameters for when a new peak is explicitly added to a peak list
  using a tabular display.

  The user can specify the position of the peak's dimensions in ppm, Hz or data
  point units. Also, the user can adjust the height and volume peak intensity
  values and a textual "Details" field, ,which can carry the user's comments
  about the peak.

  When editing an existing peak, no changes are made to the peak until the
  [Update] button is pressed. Likewise for a new peak the [Add Peak] button
  commits the changes. If the popup window is closed before the changes are 
  committed then the entire editing or peak addition operation is cancelled.

  """
  def __init__(self, parent, peak=None, peakList=None, *args, **kw):

    self.titleColor = '#000080'
    self.numDims    = 0
    self.peak       = peak
    
    kw['borderwidth'] = 6
    BasePopup.__init__(self, parent=parent, title='Edit Peak', **kw)

    self.registerNotify(self.deletedPeak, 'ccp.nmr.Nmr.Peak', 'delete')

    for func in ('setAnnotation','setDetails','setFigOfMerit'):
      self.registerNotify(self.updatePeak, 'ccp.nmr.Nmr.Peak', func)
    for func in ('setAnnotation','setPosition','setNumAliasing'):
      self.registerNotify(self.updatePeak, 'ccp.nmr.Nmr.PeakDim', func)
    for func in ('__init__','delete','setValue'):
      self.registerNotify(self.updatePeak, 'ccp.nmr.Nmr.PeakIntensity', func)

    self.dimensionLabels =[]
    self.dimensionEntries=[]
    self.update(self.peak, peakList)

  def body(self, guiParent):

    self.geometry("+150+150")
    
    guiParent.grid_columnconfigure(0, weight=1)
    self.master_frame = guiParent
    
    units = ('ppm','point','Hz')
    
    self.unit = 'ppm'
    
    self.specLabel = Label(guiParent, fg=self.titleColor, grid=(0,0), sticky='ew')
 
    self.peakLabel = Label(guiParent, grid=(0,1), sticky='ew')
 
    self.unit_frame = frame = Frame(guiParent, grid=(1,1), gridSpan=(1,2))

    self.unitLabel  = Label(frame, text='Current units: ', grid=(0,0))
    tipText = 'Selects which unit of measurement to display peak dimension positions with'
    self.unitSelect = PulldownList(frame, callback=self.changeUnit,
                                   texts=units, grid=(0,1), tipText=tipText)
 
    self.heightLabel = Label(guiParent, text='Height',
                             borderwidth=2, relief='groove')
    tipText = 'Sets the peak height; the value of the spectrum point intensity (albeit often interpolated)'
    self.heightEntry = FloatEntry(guiParent, borderwidth=1, tipText=tipText)
    self.volumeLabel = Label(guiParent, text='Volume',
                             borderwidth=2, relief='groove')
    tipText = 'Sets the peak volume integral; normally a summation of data point values'
    self.volumeEntry = FloatEntry(guiParent, borderwidth=1, tipText=tipText)
    self.detailLabel = Label(guiParent, text='Details',
                             borderwidth=2, relief='groove')
    tipText = 'A user-configurable textual comment for the peak, which appears an tables and occasionally on spectrum displays'
    self.detailEntry = Entry(guiParent, borderwidth=1, tipText=tipText)

    tipTexts = ['Commits the specified values to update the peak and closes the popup',]
    texts = [ 'Update' ]
    commands = [ self.commit ]
    self.buttons = UtilityButtonList(guiParent, texts=texts, commands=commands, 
                                     doClone=False, helpUrl=self.help_url,
                                     tipTexts=tipTexts)
  def open(self):
  
    self.updatePeak()
    BasePopup.open(self)
    

  def updatePeak(self, object=None):
  
    peak = None
    if object:
      if object.className == 'Peak':
        peak = object
      elif object.className == 'PeakDim':
        peak = object.peak
      elif object.className == 'PeakIntensity':
        peak = object.peak
    
    if (peak is None) or (peak is self.peak):
      self.update(peak=self.peak)
    
  def update(self, peak = None, peakList = None):

    # first destroy old labels and entries (saves grid hassles)

    for label in self.dimensionLabels:
      label.destroy()
    for entry in self.dimensionEntries:
      entry.destroy()

    # now setup required data

    if peak:
      title = 'Edit Peak'
      self.buttons.buttons[0].config(text='Update')
    else:
      title = 'Add Peak'
      self.buttons.buttons[0].config(text='Add Peak')

    self.setTitle(title)

    self.peak = peak
    self.peakList = peakList
    if not peakList:
      if peak:
        self.peakList = peak.peakList
      else:
        return

    peakList  = self.peakList
    spectrum = peakList.dataSource.name
    self.numDims  = peakList.dataSource.numDim
    self.posn     = self.numDims * [0]
    self.dataDims = peakList.dataSource.sortedDataDims()

    if self.peak:
      
      serial   = self.peak.serial
      dims     = self.peak.sortedPeakDims()
      details  = self.peak.details
      if not details:
        details = ''
      if self.peak.annotation:
        annotn = '%0.16s' % self.peak.annotation
      else:
        annotn = ''
 
      heightIntensity = self.peak.findFirstPeakIntensity(intensityType='height')
      volumeIntensity = self.peak.findFirstPeakIntensity(intensityType='volume')

      if heightIntensity:
        height = heightIntensity.value
      else:
        height = 0.0
        
      if volumeIntensity:
        volume = volumeIntensity.value
      else:
        volume = 0.0
	
      for i in range(self.numDims):
        peakDim = dims[i]
        dataDimRef = peakDim.dataDimRef
        if dataDimRef:
          self.posn[i] = peakDim.position + (peakDim.numAliasing*dataDimRef.dataDim.numPointsOrig)
        else:
          self.posn[i] = peakDim.position

      
    else:
          
      dict = peakList.__dict__.get('serialDict')
      if dict is None:
        serial = 1
      else:
        serial = dict.get('peaks',0) + 1
      
      height  = 0.0 
      volume  = 0.0
      details = ''
      annotn  = ''

    self.specLabel.set(text='Experiment: %s Spectrum: %s PeakList: %d' % (peakList.dataSource.experiment.name,spectrum,peakList.serial))
    self.peakLabel.set(text='Peak: %d' % serial)
 
    self.dimensionLabels =self.numDims*['']
    self.dimensionEntries=self.numDims*['']
    for i in range(self.numDims):
      pos = self.posn[i]
      if self.unit != 'point':
        dataDim = self.dataDims[i]
        if dataDim.className == 'FreqDataDim':
          pos = unit_converter[('point', self.unit)]( pos, getPrimaryDataDimRef(dataDim) )
      self.dimensionLabels[i] = Label(self.master_frame, text='F%d' % (i+1), borderwidth=2, relief='groove')
      tipText = 'The peak position in dimension %d, in the specified units' % (i+1)
      self.dimensionEntries[i] = FloatEntry(self.master_frame, borderwidth=1,
                                            text='%8.4f' % pos, tipText=tipText)

    self.heightEntry.set(text='%f' % height)
    self.volumeEntry.set(text='%f' % volume)
    self.detailEntry.set(text=details)

    row = 0
    self.specLabel.grid(row = row, column = 0, columnspan=2, sticky='nsew')
 
    row = row + 1
    self.peakLabel.grid(row = row, column = 0,  sticky='nsew')
    self.unit_frame.grid(row = row, column = 1, columnspan=2, sticky='nsew')

    for i in range(self.numDims):
      row = row + 1
      self.dimensionLabels[i].grid(row = row, column = 0, sticky='nsew')
      self.dimensionEntries[i].grid(row = row, column = 1, columnspan=3, sticky='e')
      
    row = row + 1
    self.heightLabel.grid(row = row, column = 0, sticky='nsew')
    self.heightEntry.grid(row = row, column = 1, columnspan=3, sticky='e')

    row = row + 1
    self.volumeLabel.grid(row = row, column = 0, sticky='nsew')
    self.volumeEntry.grid(row = row, column = 1, columnspan=3, sticky='e')

    row = row + 1
    self.detailLabel.grid(row = row, column = 0, sticky='nsew')
    self.detailEntry.grid(row = row, column = 1, columnspan=3, sticky='e')
    
    row = row + 1
    self.buttons.grid(row = row, column = 0, columnspan = 4, sticky='nsew')

  def changeUnit(self, unit):
  
    posDisp = self.numDims*[None]
    for i in range(self.numDims):
      posDisp[i] = float(self.dimensionEntries[i].get() )
      if self.unit != 'point':
        dataDim = self.dataDims[i]
        if dataDim.className == 'FreqDataDim':
          posDisp[i] = unit_converter[(self.unit,'point')](posDisp[i],getPrimaryDataDimRef(dataDim))
        
    self.unit = unit
    if self.unit != 'point':
      for i in range(self.numDims):
        dataDim = self.dataDims[i]
        if dataDim.className == 'FreqDataDim':
          posDisp[i] = unit_converter[('point',self.unit)](posDisp[i], getPrimaryDataDimRef(dataDim) )
     
    for i in range(self.numDims):
      value = posDisp[i]
      if value is None:
        self.dimensionEntries[i].set('None') 
      else:
        self.dimensionEntries[i].set('%8.4f' % posDisp[i])    

  def commit(self):
      
    posDisp = self.numDims * [0]
     
    for i in range(self.numDims):
      posDisp[i] = float(self.dimensionEntries[i].get() )
      if self.unit != 'point':
        dataDim = self.dataDims[i]
        if dataDim.className == 'FreqDataDim':
          self.posn[i] = unit_converter[(self.unit,'point')]( posDisp[i], getPrimaryDataDimRef(dataDim) )
          
      else:
        self.posn[i] = posDisp[i]
    
    if self.peak:
      movePeak(self.peak,self.posn)
    else:
      self.peak = pickPeak(self.peakList, self.posn)
   
    height = self.heightEntry.get()
    volume = self.volumeEntry.get()
    setManualPeakIntensity(self.peak, height, intensityType='height')
    setManualPeakIntensity(self.peak, volume, intensityType='volume')
 
    details    = self.detailEntry.get() or None
    
    self.peak.setDetails( details )

    self.close()
    
  def deletedPeak(self, peak):

    if self.peak is peak:
      self.close()

  def destroy(self):

    self.unregisterNotify(self.deletedPeak, 'ccp.nmr.Nmr.Peak', 'delete')

    for func in ('setAnnotation','setDetails','setFigOfMerit'):
      self.unregisterNotify(self.updatePeak, 'ccp.nmr.Nmr.Peak', func)
    for func in ('setAnnotation','setPosition','setNumAliasing'):
      self.unregisterNotify(self.updatePeak, 'ccp.nmr.Nmr.PeakDim', func)
    for func in ('__init__','delete','setValue'):
      self.unregisterNotify(self.updatePeak, 'ccp.nmr.Nmr.PeakIntensity', func)

    BasePopup.destroy(self)
