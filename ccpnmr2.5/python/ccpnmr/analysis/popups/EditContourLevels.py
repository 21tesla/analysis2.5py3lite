
"""
======================COPYRIGHT/LICENSE START==========================

EditContourLevels.py: Part of the CcpNmr Analysis program

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
from memops.gui.CheckButtons import CheckButtons
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showError
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons
from memops.gui.Separator import Separator
from memops.gui.ValueRamp import ValueRamp

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import calcContourLevels
from ccpnmr.analysis.core.Util import updateSpectrumLevelParams

notify_funcs = ('__init__', 'delete', 'setValue')

multiplier_text = 'Level multiplier'
adder_text = 'Level increment'
multiplier_changes = [ 1.2, 1.3, 1.4, 1.5, 1.6, 1.7 ]
adder_changes = [ 0.1, 0.2, 0.3, 0.5, 1.0, 2.0 ]

class EditContourLevelsPopup(BasePopup):
  """
  **Change Levels in Spectrum Contour Displays**
  
  This popup window is used to specify which contour levels (lines of constant
  intensity) are drawn for spectra within the spectrum windows of Analysis. A
  spectrum can have both positive and negative levels, and the user can set
  these on an individual basis or as a regular incremental (usually geometric)
  series. It should be noted that spectrum contour colours, which may be
  different for positive and negative levels, is set elsewhere, in the main
  Spectra_ table.

  In normal operation the user first selects the spectrum to change the contours
  for in the upper left pulldown menu. Although contour levels are specified 
  for individual spectra, if there are several spectra that share the same
  levels (for example if they form a series of some kind) then the user may
  spread the contour level information from one to many via the [Propagate
  Contours] function within the "Display Options" of the main Spectra_ table.
  The user may also affect the contour levels of multiple spectra by changing
  the "Global scale" settings. This global scale will affect *all spectra* in
  the project; the scale value multiplies all of the contour levels for all
  spectra, although in practice it is rarely changed. Using a global scale gives
  the advantage of being able to specify smaller numbers for the levels. 

  The individual contour levels for a spectrum are listed in the "positive
  levels" and "Negative levels" fields. They are multiplied by the global scale
  factor when they are used for spectrum display. The user may type in values
  for the levels directly into these field, or they can be filled using the
  "Auto contour levels" mechanism. If manual levels are specified the user
  should press [Apply manual Levels] when done to commit the changes and see the
  results. The [Apply Auto Levels] function differs in that it always applies a
  regular series of levels, according to the settings, and these will overwrite
  any manual specifications when committed; the actual levels applied are always
  displayed in the levels field.

  The setting of "Auto" levels using a series involves choosing a base level,
  which represents the first contour (closest to zero). This base level is
  usually set according to the level of noise in the spectrum. Typically it is a
  value just above most of the noise, so the noise is not seen, but it may be
  set lower in some instances to show weak signals. The base level applies to
  both positive and negative contour series, although for the negative side it
  is naturally used with negative sign. The user sets how many levels there should be in total; the subsequent levels
  start from the base level and move away from zero. Normally this involves a
  geometric series, with a constant multiplication factor applied to a value to
  generate the next in the series. This factor is set in the "Level multiplier",
  although some of the commonly used values are quickly set via the adjacent
  buttons. Under a few special circumstances it is helpful to  have a constant
  difference  between levels, in which case "Add levels" may be used instead of
  "Multiply levels". 

  **Caveats & Tips**

  Often it is useful to initially setup contours using an "Auto" series, but
  then manually remove some of the levels. For example having all the positive
  levels but only the first negative level may be helpful to show where peaks
  are truncated.

  Note this system sets the levels that would be used to make any contour files
  (an optional way of working, especially for 4D spectra etc.). However,
  once contour files are made then this system will not affect the levels
  within the files.

  .. _Spectra: EditSpectrumPopup.html

  """

  def __init__(self, parent, *args, **kw):

    BasePopup.__init__(self, parent=parent, title='Spectrum Contour Levels', **kw)

  def body(self, guiFrame):

    self.doUpdateForm = True
    self.spectrum = None

    guiFrame.grid_columnconfigure(1, weight=1)

    row = 0
 
    specFrame = LabelFrame(guiFrame, text='Spectrum', grid=(row,0))

    tipText = 'The spectrum for which you are setting the contour levels'
    self.expt_spectrum = PulldownList(specFrame, grid=(0,0), tipText=tipText,
                                      callback=self.setSpectrumDetails)
    
    tipText = 'Open the edit spectrum table to view and change other properties associated with the selected spectrum'
    button = Button(specFrame, text='Spectrum Properties',
                    command=self.editProperties, borderwidth=1,
                    grid=(1,0), sticky='ew', tipText=tipText)

    globalFrame = LabelFrame(guiFrame, text='Global scale')
    globalFrame.grid(row=row, column=1, sticky='nsew')
    globalFrame.grid_columnconfigure(0, weight=1)

    tipText = 'The value by which all contour levels in all spectra are multiplied; to give the actual contour level in terms of the stored spectrum data'
    self.global_entry = FloatEntry(globalFrame, width=10, tipText=tipText,
                                   text=self.analysisProject.globalContourScale,
                                   returnCallback=self.applyAuto, grid=(0,0),
                                   gridSpan=(1,2), sticky='ew')

    tipText = 'Divide the global scaling factor by two; moving all contour for all spectra levels closer to zero'
    button = Button(globalFrame, borderwidth=1, text='/2', sticky='ew', tipText=tipText,
                    command=lambda: self.changeGlobalScale(0.5), grid=(0,2))

    tipText = 'Multiple the global scaling factor by two; moving all contour for all spectra levels away from zero'
    button = Button(globalFrame, borderwidth=1, text='*2', sticky='ew', tipText=tipText,
                    command=lambda: self.changeGlobalScale(2.0), grid=(0,3))

    tipText = 'Set the global contour scale for all spectra to the default of 100,000'
    button = Button(globalFrame, borderwidth=1, text='10^5', tipText=tipText,
                    command=self.defaultGlobalScale, grid=(1,0))
    
    tipText = 'Click to decrease "-" or increase "+" the global contour scale by small amounts'
    frame = ValueRamp(globalFrame, callback=self.changeGlobalScale, tipText=tipText)
    frame.grid(row=1, column=1, columnspan=3, sticky='e')


    row += 1
    self.autoFrame = LabelFrame(guiFrame, text='Auto contour levels')
    self.autoFrame.grid(row=row, column=0, columnspan=5, sticky='nsew')
    self.autoFrame.grid_columnconfigure(1, weight=1)
    frow = 0
    
    label = Label(self.autoFrame, text='Base level:', grid=(frow,0), sticky='e')
    
    tipText = 'The first contour level (closest to zero) for an automated series of levels: the start of a geometric or arithmetic series defining levels'
    self.base_entry = FloatEntry(self.autoFrame, returnCallback=self.applyAuto,
                                 width=10, grid=(frow,1), sticky='ew', tipText=tipText)

    command = lambda: self.changeBaseLevel(0.5)
    tipText = 'Lower the base contour level so that it is half of the previous value; moving the series of contour levels closer to zero'
    button = Button(self.autoFrame, borderwidth=1, text='/2', tipText=tipText,
                    command=command, grid=(frow,2), sticky='ew')

    command = lambda: self.changeBaseLevel(2.0)
    tipText = 'Raise the base contour level so that it is double the previous value; moving the series of contour levels further from zero'
    button = Button(self.autoFrame, borderwidth=1, text='*2', tipText=tipText,
                    command=command, grid=(frow,3), sticky='ew')

    tipText = 'Click to decrease "-" or increase "+" the base contour level by small amounts'
    frame = ValueRamp(self.autoFrame, callback=self.changeBaseLevel, tipText=tipText)
    frame.grid(row=frow, column=4, columnspan=3, sticky='ew')

    frow += 1
    label = Label(self.autoFrame, text='Number of levels:',
                  grid=(frow,0), sticky='e')
    #self.numberEntry = IntEntry(guiFrame, text=numberLevels,
    #                            returnCallback=self.applyAuto)
    tipText = 'The number of contour levels to make in the automated series'
    self.numberEntry = IntEntry(self.autoFrame, returnCallback=self.applyAuto,
                                width=10, grid=(frow,1), sticky='ew', tipText=tipText)

   
    command = lambda w=-1: self.changeNumberLevels(w)
    tipText = 'Decrease the number of contour levels in the series by one'
    button = Button(self.autoFrame, borderwidth=1, text='-1', tipText=tipText,
                   command=command, grid=(frow,2), sticky='ew')
      
    command = lambda w=1: self.changeNumberLevels(w)
    tipText = 'Increase the number of contour levels in the series by one'
    button = Button(self.autoFrame, borderwidth=1, text='+1', tipText=tipText,
                   command=command, grid=(frow,3), sticky='ew')
 
    
    n = 4               
    for w in (5, 10, 15, 20):
      tipText = 'Set the number of contour levels in the series to %d' % w
      command = lambda w=w: self.setNumberLevels(w)
      button = Button(self.autoFrame, borderwidth=1, text='%d' % w,
                      command=command, grid=(frow,n), sticky='ew',
                      tipText=tipText)
      n += 1

    frow += 1
    self.change_label = Label(self.autoFrame, text='%s:' % multiplier_text,
                              grid=(frow,0), sticky='e')

    tipText = 'The multiplication factor (or increment if adding levels) to repetitively apply to the base level to generate the series of levels'
    self.change_entry = FloatEntry(self.autoFrame, returnCallback=self.applyAuto,
                                   width=10, grid=(frow,1), sticky='ew',
                                   tipText=tipText)

    self.change_level_buttons = []
    n = 2
    for w in multiplier_changes:
      tipText = 'Set the automated series multiplication factor or increment to %f' % w
      command = lambda w=w: self.setChange(w)
      button = Button(self.autoFrame, borderwidth=1,
                      text=str(w), command=command,
                      grid=(frow,n), sticky='ew',
                      tipText=tipText)
      self.change_level_buttons.append(button)
      n += 1

    frow += 1
    frame = Frame(self.autoFrame, grid=(frow,0), gridSpan=(1,7), sticky='ew')
   
    tipTexts = ['Toggles whether positive contour levels from the automated series will be used. Overrides any manual settings',
                'Toggles whether negative contour levels from the automated series will be used. Overrides any manual settings']
    entries = ('Positive', 'Negative')
    selected = (True, False)
    self.which_buttons = CheckButtons(frame, entries, selected=selected,
                                      select_callback=self.applyAuto,
                                      grid=(0,0), sticky='w', tipTexts=tipTexts)

    tipTexts = ['Set the contour level generation to use a geometric series, starting from the base level and using the specified factor',
                'Set the contour level generation to use an arithmetic series, starting from the base level and using the specified increment']
    entries = ('Multiply levels', 'Add levels')
    self.change_mode_buttons = RadioButtons(frame, entries, grid=(0,1), sticky='ew',
                                            select_callback=self.modeChanged,
                                            tipTexts=tipTexts)

    row += 1
    manualFrame = LabelFrame(guiFrame, text='Positive levels',
                             grid=(row,0), gridSpan=(1,2))
    manualFrame.expandGrid(None, 0)
    tipText = 'The positive contour levels that will be used for the spectrum; filled in by the automation or set/modified manually'
    self.posLevelsEntry = FloatEntry(manualFrame, isArray=True, width=60,
                                     returnCallback=self.applyManual,
                                     tipText=tipText, grid=(0,0), sticky='ew')

    row += 1
    manualFrame = LabelFrame(guiFrame, text='Negative levels',
                             grid=(row,0), gridSpan=(1,2))
    manualFrame.expandGrid(None, 0)
    tipText = 'The negative contour levels that will be used for the spectrum; filled in by the automation or set/modified manually'
    self.negLevelsEntry = FloatEntry(manualFrame, isArray=True, width=60,
                                     returnCallback=self.applyManual,
                                     tipText=tipText, grid=(0,0), sticky='ew')

    row += 1
    tipTexts = ['Set the spectrum contour levels, updating the display, to the automated series, ignoring any manual edits.',
                'Set the spectrum contour levels, updating the display, to the values displayed in the level entry fields']
    texts = ['Apply Auto Levels', 'Apply Manual Edits']
    commands = [ self.applyAuto, self.applyManual ]
    self.buttons = UtilityButtonList(guiFrame, texts=texts, commands=commands, 
                                     helpUrl=self.help_url, grid=(row,0),
                                     gridSpan=(1,2), tipTexts=tipTexts)
    guiFrame.grid_rowconfigure(row, weight=1)

    self.curateNotifiers(self.registerNotify)
    self.update()

  def update(self, spectrum=None):

    if not spectrum:
      spectrum = self.spectrum

    spectra = self.parent.getSpectra()
    if spectra:
      if spectrum not in spectra:
        spectrum = spectra[0]
      index = spectra.index(spectrum)
      names = ['%s:%s' % (x.experiment.name, x.name) for x in spectra]
    else:
      index = 0
      names = []

    self.expt_spectrum.setup(names, spectra, index)

    self.setSpectrumDetails(spectrum)

  def updateNotifier(self, *extra):

    self.update()

  def open(self):
  
    self.updateForm()
    BasePopup.open(self)
    
  def modeChanged(self, *entry):

    spectrum = self.spectrum
    if not spectrum or spectrum.isDeleted:
      return
      
    analysisSpectrum = spectrum.analysisSpectrum

    changeMode = self.change_mode_buttons.getIndex() and 'add' or 'multiply'
    if changeMode == 'add':
      text = adder_text
      changes = adder_changes
    else:
      text = multiplier_text
      changes = multiplier_changes

    self.change_label.set('%s:' % text)
    n = 0
    for button in self.change_level_buttons:
      w = changes[n]
      command = lambda w=w: self.setChange(w)
      button.config(command=command)
      button.setText(str(w))
      n = n + 1

    levelChanger = changes[2]
    self.doUpdateForm = False

    analysisSpectrum.autoLevelChanger = levelChanger
    analysisSpectrum.autoLevelMode    = changeMode
    
    self.doUpdateForm = True
    self.setContourLevels()

  def defaultGlobalScale(self):
  
    self.global_entry.set(100000)
    self.applyAuto()

  def close(self):
  
    self.applyManual()
    BasePopup.close(self)

  def destroy(self):
 
    self.curateNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  def curateNotifiers(self, notifyFunc):

    notifyFunc(self.updateContourLevels, 'ccpnmr.Analysis.AnalysisSpectrum', 'setPosLevels')
    notifyFunc(self.updateContourLevels, 'ccpnmr.Analysis.AnalysisSpectrum', 'setNegLevels')
    notifyFunc(self.updateForm, 'ccpnmr.Analysis.AnalysisSpectrum', 'setAutoLevelChanger')
    notifyFunc(self.updateForm, 'ccpnmr.Analysis.AnalysisSpectrum', 'setAutoLevelMode')
    notifyFunc(self.updateForm, 'ccpnmr.Analysis.AnalysisSpectrum', 'setAutoNumLevels')
    notifyFunc(self.updateForm, 'ccpnmr.Analysis.AnalysisSpectrum', 'setAutoBaseLevel')
    
    notifyFunc(self.updateForm, 'ccpnmr.Analysis.AnalysisProject', 'setGlobalContourScale')

    for clazz in ('Experiment', 'DataSource'):
      for func in ('__init__', 'delete', 'setName'):
        notifyFunc(self.updateNotifier, 'ccp.nmr.Nmr.%s' % clazz, func)
 
  def editProperties(self):

    self.parent.editSpectrum(self.spectrum)

  def updateForm(self, *extra):

    #print 'updateForm'
    if (not self.doUpdateForm):
      return

    self.global_entry.set(self.analysisProject.globalContourScale)
    spectrum = self.spectrum
    if spectrum and not spectrum.isDeleted:
      analysisSpectrum = spectrum.analysisSpectrum
      
      self.base_entry.set(analysisSpectrum.autoBaseLevel)
      self.numberEntry.set(analysisSpectrum.autoNumLevels)
      self.change_entry.set(analysisSpectrum.autoLevelChanger)
      
      if analysisSpectrum.autoLevelMode == 'add':
        i = 1
      else:
        i = 0  
      
      self.change_mode_buttons.setIndex(i)

  def updateContourLevels(self, analysisSpectrum):

    spectrum = self.spectrum
    if spectrum and not spectrum.isDeleted:
      analysisSpectrum = spectrum.analysisSpectrum
      
      posLevels = list(analysisSpectrum.posLevels)
      negLevels = list(analysisSpectrum.negLevels)
      
      self.posLevelsEntry.set(posLevels)
      self.negLevelsEntry.set(negLevels)
      
      self.doUpdateForm = False
      updateSpectrumLevelParams(analysisSpectrum, posLevels, negLevels)

      self.doUpdateForm = True
      self.base_entry.set(analysisSpectrum.autoBaseLevel)
      self.numberEntry.set(analysisSpectrum.autoNumLevels)
      self.change_entry.set(analysisSpectrum.autoLevelChanger)
      
      self.setWhichLevels(spectrum)

      if analysisSpectrum.autoLevelMode == 'add':
        i = 1
      else:
        i = 0  
       
      self.change_mode_buttons.setIndex(i)

  def setWhichLevels(self, spectrum):

    analysisSpectrum = spectrum.analysisSpectrum
    
    posLevels = analysisSpectrum.posLevels
    negLevels = analysisSpectrum.negLevels
    
    if posLevels:
      isSelected = True
    else:
      isSelected = False
    self.which_buttons.setIndexSelection(0, isSelected)

    if negLevels:
      isSelected = True
    else:
      isSelected = False
      
    self.which_buttons.setIndexSelection(1, isSelected)

  def setSpectrum(self, spectrum):

    if spectrum is not self.spectrum:
      self.update(spectrum)

    #if (spectrum and not spectrum.isDeleted):
    #  self.setWhichLevels(spectrum)

  speed_scale = 6.0
  speed_delay = 50 # msec

  def changeGlobalScale(self, multiplier):

    self.analysisProject.globalContourScale = multiplier * self.analysisProject.globalContourScale

  def changeBaseLevel(self, multiplier):

    spectrum = self.spectrum
    if (not spectrum or spectrum.isDeleted is True):
      return

    analysisSpectrum = spectrum.analysisSpectrum
    baseLevel = multiplier * analysisSpectrum.autoBaseLevel
    self.base_entry.set(baseLevel)
    self.doUpdateForm = False
    
    analysisSpectrum.autoBaseLevel = abs(baseLevel)
    
    self.doUpdateForm = True
    self.setContourLevels()

  def changeNumberLevels(self, change):

    spectrum = self.spectrum
    if not spectrum or spectrum.isDeleted is True:
      return

    analysisSpectrum = spectrum.analysisSpectrum
    numberLevels = analysisSpectrum.autoNumLevels + change
    self.numberEntry.set(numberLevels)
    self.doUpdateForm = False
    
    analysisSpectrum.autoNumLevels = numberLevels
    
    self.doUpdateForm = True
    self.setContourLevels()
    
  def setChange(self, levelChanger):

    spectrum = self.spectrum
    if not spectrum or spectrum.isDeleted:
      return

    self.doUpdateForm = False
    analysisSpectrum = spectrum.analysisSpectrum
    analysisSpectrum.autoLevelChanger = levelChanger
    self.doUpdateForm = True

    self.setContourLevels()

  def setNumberLevels(self, numberLevels):
    
    spectrum = self.spectrum
    if (not spectrum or spectrum.isDeleted is True):
      return

    self.doUpdateForm = False
    analysisSpectrum = spectrum.analysisSpectrum
    analysisSpectrum.autoNumLevels = numberLevels
    self.doUpdateForm = True
    
    self.setContourLevels()

  def setSpectrumDetails(self, spectrum):

    if spectrum is self.spectrum:
      return

    self.spectrum = spectrum
    if spectrum and not spectrum.isDeleted:
      analysisSpectrum = spectrum.analysisSpectrum
      
      posLevels = list(analysisSpectrum.posLevels)
      negLevels = list(analysisSpectrum.negLevels)
      
      self.posLevelsEntry.set(posLevels)
      self.negLevelsEntry.set(negLevels)
      
      self.doUpdateForm = False
      updateSpectrumLevelParams(analysisSpectrum, posLevels, negLevels)
      
      self.doUpdateForm = True
      self.base_entry.set(analysisSpectrum.autoBaseLevel)
      self.numberEntry.set(analysisSpectrum.autoNumLevels)
      self.change_entry.set(analysisSpectrum.autoLevelChanger)
      self.autoFrame.setText('Auto contour levels - %s:%s' % (spectrum.experiment.name, spectrum.name))
      self.setWhichLevels(spectrum)
      
    else:
      self.posLevelsEntry.set('')
      self.negLevelsEntry.set('')
      self.autoFrame.setText('Auto contour levels')

  def setContourLevels(self):

    spectrum = self.spectrum
    if not spectrum or spectrum.isDeleted is True:
      return

    try:
      analysisSpectrum = spectrum.analysisSpectrum
      baseLevel    = analysisSpectrum.autoBaseLevel   
      numberLevels = analysisSpectrum.autoNumLevels   
      levelChanger = analysisSpectrum.autoLevelChanger
      changeMode   = analysisSpectrum.autoLevelMode   
      
      posLevels = []
      if self.which_buttons.isIndexSelected(0):
        posLevels.extend(calcContourLevels(baseLevel, numberLevels, levelChanger, changeMode))
        
      negLevels = []
      if self.which_buttons.isIndexSelected(1):
        if changeMode == 'add':
          levelChanger = -levelChanger
        negLevels.extend(calcContourLevels(-baseLevel, numberLevels, levelChanger, changeMode))
        
      self.posLevelsEntry.set(posLevels)
      self.negLevelsEntry.set(negLevels)
      analysisSpectrum.posLevels = posLevels
      analysisSpectrum.negLevels = negLevels
      
    except Implementation.ApiError, e:
      showError('Contour levels error', e.error_msg,  parent=self)

  def applyManual(self, *extra):

    try:
      levels = self.posLevelsEntry.get()
      posLevels = [l for l in levels if l > 0]
      
      levels = self.negLevelsEntry.get()
      negLevels = [l for l in levels if l < 0]
      
      spectrum = self.spectrum
      analysisSpectrum = spectrum.analysisSpectrum
      
      analysisSpectrum.posLevels = posLevels
      analysisSpectrum.negLevels = negLevels
      
    except Implementation.ApiError, e:
      showError('Apply error', e.error_msg,  parent=self)

  def applyAuto(self, *extra):

    try:
      spectrum = self.spectrum
      scale = self.global_entry.get()
      changeMode = self.change_mode_buttons.getIndex() and 'add' or 'multiply'
      baseLevel = self.base_entry.get()
      
      if (baseLevel <= 0):
        raise Implementation.ApiError('Base level must be set to positive float')
      numberLevels = self.numberEntry.get()
      if (numberLevels < 1):
        raise Implementation.ApiError('Number of levels must be set to positive int')
        
      levelChanger = self.change_entry.get()
      if changeMode == 'add':
        if levelChanger <= 0.0:
          raise Implementation.ApiError('Level adder must be set to number > 0.0')
      else:
        if levelChanger <= 1.0:
          raise Implementation.ApiError('Level multiplier must be set to number > 1.0')
          
      self.analysisProject.globalContourScale = scale
      self.doUpdateForm = False
      
      analysisSpectrum = spectrum.analysisSpectrum
      analysisSpectrum.autoBaseLevel    = baseLevel
      analysisSpectrum.autoNumLevels    = numberLevels
      analysisSpectrum.autoLevelChanger = levelChanger
      analysisSpectrum.autoLevelMode    = changeMode
      
      self.setContourLevels()
      self.doUpdateForm = True
    except Implementation.ApiError, e:
      showError('Levels error', e.error_msg,  parent=self)
