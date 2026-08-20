"""
======================COPYRIGHT/LICENSE START==========================

EditMeasurementLists.py: Part of the CcpNmr Analysis program

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

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel, showWarning
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.core.ExperimentBasic import newShiftList
from ccpnmr.analysis.core.ChemicalShiftBasic import makeRandomCoilShiftList
from ccpnmr.analysis.core.AssignmentBasic import averageShiftValue, makeResonanceGuiName, getResonanceResidue
from ccpnmr.analysis.popups.BasePopup import BasePopup

class EditMeasurementListsPopup(BasePopup):
  """
  **The NMR Derived Measurements and Their Lists**
  
  The purpose of this popup window is to display lists of measurements within
  the current CCPN project. A "Measurement" in this regard is a value that is
  derived from NMR data and connects to one or more resonances (which may be
  assigned to atoms). The most commonly encountered kind of measurement is the
  chemical shift, which is usually derived from recording the positions of peaks
  in spectra. In this instance the chemical shift measurement of a resonance is
  made when it is assigned to the dimension of a peak. Because CCPN allows
  multiple shift lists (technically a type of measurement list) a resonance may
  have a several chemical shift measurements; useful studying different
  conditions. 

  There are several different kinds of measurement and hence measurement list that
  may be included within a CCPN project, for example J-coupling, hydrogen
  exchange rate, chemical shift anisotropy, T1 relaxation and chemical shift.
  All types will be displayed by this popup if they are available within the
  current project, although Analysis does not necessarily have tools to record
  all of these different kinds of measurement.

  The layout of this popup consists of two tabs, one for the organising lists
  that contain the measurements and another to display the individual
  measurements within a single list. With the exception of blank
  lists for chemical shifts, measurement lists are not made using this popup,
  instead they are made at the point where an analysis is performed, e.g. T1
  lists are made with data from the `Follow Intensity Changes`_
  tool. 

  The second tab, containing the "Measurements" table displays the details of
  the individual measured values within a selected list. What the values mean
  and which units they are in, if any, naturally depends on the kind of list
  being viewed. If measurements are made within CCPN software (as opposed to
  being imported) then the individual measurements usually record the spectrum
  peaks that were used in the calculation of the measured value, e.g. the peaks
  for which positions record chemical shift or intensities record T1.

  **Caveats & Tips**

  If there is no chemical shift list within a project a new one will be made
  automatically to record the shifts of any assignments.

  The shift list with which an experiment is associate may be changed via the "Shift
  List" column of the main `Experiments`_ table.

  Measurements other than chemical shift, like T1 and T2 relaxation times will only
  appear within this system if a measurement list entity is formally made. Simply
  measuring values from spectra, for example using the `Follow Intensity Changes`_
  tool

  .. _`Follow Intensity Changes`: CalcRatesPopup.html
  .. _`Experiments`: EditExperimentPopup.html
  
  """

  def __init__(self, parent, *args, **kw):

    self.guiParent   = parent
    self.experiment  = None
    self.measurementList = None
    self.measurementListB = None
    self.measurement = None
    self.waitingList = False
    self.waitingMeasurement = False
  
    BasePopup.__init__(self, parent, title="Resonance : Measurement Lists", **kw)

  def body(self, guiFrame):

    self.geometry('600x600')

    guiFrame.expandGrid(0,0)
    
    tipTexts = ['A table of all the NMR measurement lists in the project, including shift lists, T1 lists, J-coupling lists etc.',
                'A table listing all of the individual measurements within an NMR measurement list']
    options = ['Measurement Lists','Measurements Table']
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              grid=(0,0), tipTexts=tipTexts)
    self.tabbedFrame = tabbedFrame
    frameA, frameB = tabbedFrame.frames
    
    # Measurement Lists
    
    frameA.expandGrid(1,0)
    self.detailsEntry = Entry(self, text='', returnCallback=self.setDetails, width=12)
    self.nameEntry    = Entry(self, text='', returnCallback=self.setName, width=12)


    row = 0
    frame0 = Frame(frameA, grid=(row,0), gridSpan=(1,2))
    frame0.grid_columnconfigure(2, weight=1)
    
    label = Label(frame0, text='Experiment:', grid=(0,0), sticky='e')
    
    tipText = 'Selects an experiment, if required, to restrict the measurement list table display; showing only lists which were derived using the experiment'
    self.experimentPulldown = PulldownList(frame0, callback=self.setExperiment,
                                           grid=(0,1), tipText=tipText)

    row += 1
    tipTexts = ['The serial number of measurement list',
                'The type of measurement list, e.g. shift list, T1 list, J coupling list',
                'A short identifying name for the list, for graphical displays',
                'The number of measurements contained within the list',
                'The unit of measurement used for the values in the list',
                'The names of the experiments which were used to derive the measurements',
                'A user-specified textual comment for the measurement list']
    justifyList      = ['center','center','center','center','center','left']
    #colHeadings      = ['List Type','Name','Size','Unit','Experiments','Other Info','Details']
    colHeadings      = ['#', 'List Type','Name','Size','Unit','Experiments','Details']
    editWidgets      = [None,None,self.nameEntry,None,None,None,None,self.detailsEntry]
    editGetCallbacks = [None,None,self.getName,  None,None,None,None,self.getDetails  ]
    editSetCallbacks = [None,None,self.setName,  None,None,None,None,self.setDetails  ]
    
    self.listsMatrix = ScrolledMatrix(frameA, grid=(row,0), gridSpan=(1,2),
                                      justifyList=justifyList,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks,
                                      editWidgets=editWidgets,
                                      headingList=colHeadings,
                                      callback=self.selectListCell,
                                      deleteFunc=self.deleteMeasurementList,
                                      tipTexts=tipTexts)

    row += 1
    tipTexts = ['Show a table of the individual measurements within the selected measurement list',
                'Make a new, blank chemical shift list within the project',
                'Make a synthetic chemical shift list using random coil values, adjusting protein backbone values for sequence where approprate',
                'Delete the selected measurement list']
    texts    = ['Show Measurements','New Shift List','Make Random Coil Shift List','Delete']
    commands = [self.showMeasurements,self.addShiftList,
                self.makeRandomCoilShiftList,self.deleteMeasurementList]
    self.listButtons = ButtonList(frameA, texts=texts, commands=commands,
                                   grid=(row,0), gridSpan=(1,2), tipTexts=tipTexts)

    # Measurements
 
    self.measurementDetailsEntry = Entry(self, text='', returnCallback=self.setMeasurementDetails,width=12)
    self.measurementMeritEntry   = FloatEntry(self, returnCallback=self.setMeasurementMerit, width=6)

    row = 0
    frame0 = Frame(frameB, grid=(row,0))
    frame0.grid_columnconfigure(2, weight=1)
    
    label = Label(frame0, text='Measurement List:', grid=(0,0),sticky='e')
    
    tipText = 'Selects the measurement list to display measurements for'
    self.listPulldown = PulldownList(frame0, callback=self.setMeasurementList,
                                     grid=(0,1), tipText=tipText)

    row += 1
    frameB.expandGrid(row, 0)
    tipTexts = ['The serial number of the measurement within its containing list',
                'The number or assignment of the NMR resonance(s) to which the measurement applies',
                'The numeric value of the NMR measurement on the specified resonance(s), and the unit of measurement',
                'The standard deviation error in the measured value',
                'The molecular chain, if any, to which the measurement relates by virtue of atom assigned resonances',
                'The isotope type(s) of the measures resonance(s)',
                'A figure-of-merit value for the measurement indicating its quality or reliability',
                'The number of peaks in the CCPN project used to take the measurement',
                'A user-defined textual comment for the measurement']
    justifyList      = ['center','center','center','center','center',
                        'center','center','center','center','left']
    colHeadings      = ['#','Resonance','Value','SD','Chain','Isotope',
                        'Fig of\nMerit','Peaks','Details']
    editWidgets      = [None,None,None,None,None,None,
                        self.measurementMeritEntry,None,self.measurementDetailsEntry]
    editGetCallbacks = [None,None,None,None,None,None,
                        self.getMeasurementMerit ,self.showPeaks,self.getMeasurementDetails  ]
    editSetCallbacks = [None,None,None,None,None,None,
                        self.setMeasurementMerit ,None,self.setMeasurementDetails  ]

    self.measurementsMatrix = ScrolledMatrix(frameB, grid=(row,0),
                                         multiSelect=True, tipTexts=tipTexts,
                                         justifyList=justifyList,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         headingList=colHeadings,
                                         callback=self.selectMeasurementCell)

    row += 1
    tipTexts = ['Show a table containing peaks that were used to derive the selected measurements',
                'For some measurement lists (currently only shift lists) manually trigger a recalculation of values',
                'Show a table containing the resonances that relate to the selected measurements',
                'Delete the selected measurement records; cannot be done for chemical shift values still ties to peaks via assignment']
    texts    = ['Show Peaks','Recalculate','Show Resonances','Delete']
    commands = [self.showPeaks,
                self.recalculateMeasurements,
                self.showResonances,
                self.deleteMeasurements]
    self.measurementButtons = ButtonList(frameB,texts=texts, grid=(row,0),
                                         commands=commands, tipTexts=tipTexts)
   
    # Main Frame

    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url,
                                           grid=(0,0), gridSpan=(1,2), sticky='e')

    self.updateMeasurementListAfter()
    self.updateMeasurementsAfter()
    
    self.administerNotifiers(self.registerNotify)
    
  def administerNotifiers(self, notifyFunc):
    
    for func in ('__init__', 'delete','setName',
                 'setDipolarRelaxList','setHExchProtectionList',
                 'setHExchRateList','setNoeList',
                 'setJCouplingList','setRdcList',
                 'setShiftAnisotropyList','setShiftDifferenceList',
                 'setShiftList','setT1List',
                 'setT1RhoList','setT2List'):
      for clazz in ('ccp.nmr.Nmr.Experiment',):
        notifyFunc(self.updateExperiments, clazz, func)

    for func in ('__init__', 'delete','setDetails','setName','addExperiment','removeExperiment','setExperiments'):
      for clazz in ('ccp.nmr.Nmr.DipolarRelaxList', 'ccp.nmr.Nmr.HExchProtectionList',
                    'ccp.nmr.Nmr.HExchRateList', 'ccp.nmr.Nmr.NoeList',
                    'ccp.nmr.Nmr.JCouplingList', 'ccp.nmr.Nmr.RdcList',
                    'ccp.nmr.Nmr.ShiftAnisotropyList', 'ccp.nmr.Nmr.ShiftDifferenceList',
                    'ccp.nmr.Nmr.ShiftList', 'ccp.nmr.Nmr.T1List',
                    'ccp.nmr.Nmr.T1RhoList', 'ccp.nmr.Nmr.T2List'):
        notifyFunc(self.updateMeasurementListAfter, clazz, func)
 
    # Measurements
 
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.DipolarRelaxation', 'ccp.nmr.Nmr.HExchProtection',
                     'ccp.nmr.Nmr.HExchRate', 'ccp.nmr.Nmr.Noe',
                     'ccp.nmr.Nmr.JCoupling', 'ccp.nmr.Nmr.Rdc',
                     'ccp.nmr.Nmr.Shift', 'ccp.nmr.Nmr.ShiftAnisotropy',
                     'ccp.nmr.Nmr.ShiftDifference', 'ccp.nmr.Nmr.T1',
                     'ccp.nmr.Nmr.T1Rho', 'ccp.nmr.Nmr.T2'):
        notifyFunc(self.updateMeasurementAfter, clazz, func)
         
    for func in ('__init__', 'delete','setDetails','setError',
                  'setFigOfMerit','setValue','setPeaks',
                  'addPeak','removePeak'):
      for clazz in ('ccp.nmr.Nmr.DipolarRelaxation',
                     'ccp.nmr.Nmr.HExchProtection',
                     'ccp.nmr.Nmr.HExchRate',
                     'ccp.nmr.Nmr.Noe',
                     'ccp.nmr.Nmr.JCoupling',
                     'ccp.nmr.Nmr.Rdc',
                     'ccp.nmr.Nmr.Shift',
                     'ccp.nmr.Nmr.ShiftAnisotropy',
                     'ccp.nmr.Nmr.ShiftDifference',
                     'ccp.nmr.Nmr.T1',
                     'ccp.nmr.Nmr.T1Rho',
                     'ccp.nmr.Nmr.T2'):
        self.registerNotify(self.updateMeasurementAfter, clazz, func)  
         
  def open(self):
  
    self.updateMeasurementListAfter()
    self.updateMeasurementsAfter()
    BasePopup.open(self)
  
  def makeRandomCoilShiftList(self):
    
    molSystems = self.project.molSystems
    shiftList = makeRandomCoilShiftList(molSystems)
  
  def addShiftList(self):
  
    shiftList = newShiftList(self.project)
    self.setExperiment(None)
  
  def getName(self, measurementList):

    if measurementList :
      self.nameEntry.set(measurementList.name)

  def setName(self, event):

    text = self.nameEntry.get()
    if text and text != ' ':
      self.measurementList.setName( text )
  
  def getDetails(self, measurementList):

    if measurementList and measurementList.details:
      self.detailsEntry.set(measurementList.details)
  
  def setDetails(self, event):

    text = self.detailsEntry.get()
    if text and text != ' ':
      self.measurementList.setDetails( text )
  
  def setExperiment(self, experiment):

    if experiment is not self.experiment:
      self.experiment = experiment
      self.measurementList = None
      self.updateMeasurementListAfter()

  def deleteMeasurementList(self, *event):

    if self.measurementList:

      if ( len(self.measurementList.experiments) > 0) and (self.measurementList.className == 'ShiftList'):
        showWarning('Warning','Deletion of shift lists with associated experiments prohibited', parent=self)
        return
      elif len(self.measurementList.measurements) > 0:
        if showOkCancel('Confirm','List is not empty. Really delete?', parent=self):
          self.measurementList.delete()
          self.measurementList = None
        else:
          return
            
      else:
        self.measurementList.delete()
        self.measurementList = None

  def showMeasurements(self):
  
    self.updateMeasurements(self.measurementList)
    self.tabbedFrame.select(1)
    
  def selectListCell(self, object, row, col):
  
    self.measurementList = object
    if self.measurementList:
      self.listButtons.buttons[0].enable()
      self.listButtons.buttons[3].enable()
      self.setMeasurementList(self.measurementList)
    
  def updateMeasurementAfter(self, measurement):
  
    if (self.experiment is None) or (self.experiment in measurement.parentList.experiments):
      self.updateMeasurementListAfter()
   
    if measurement.parentList is self.measurementListB:
      self.updateMeasurementsAfter()


  def updateMeasurementListAfter(self, measurementList=None):
     
    self.updateListPulldown()
 
    if self.waitingList:
      return
    
    if (measurementList is None) \
     or (self.experiment is None) \
     or (self.experiment in measurementList.experiments):
      self.waitingList = True
      self.after_idle(self.updateListsTable)

  def updateExperiments(self, null=None):
  
    index = 0
    experiments = [None,] + self.nmrProject.sortedExperiments()
    names = ['<Any>',] + [e.name for e in experiments[1:]]
    experiment = self.experiment
    
    if experiments:
      if experiment not in experiments:
        experiment = experiments[0]
      
      index = experiments.index(experiment)  
    
    else:
      experiment = None
      
    if self.experiment is not experiment:
      self.experiment = experiment
      self.updateMeasurementListAfter()      
     
    self.experimentPulldown.setup(names, experiments, index)     

  def updateListsTable(self, experiment=None):

    if experiment is None:
      experiment = self.experiment
    else:
      self.experiment = experiment
    
    if self.experiment is None:
      self.measureList = None

    self.updateExperiments()

    if self.measurementList:
      self.listButtons.buttons[0].enable()
      self.listButtons.buttons[3].enable()
    else:
      self.listButtons.buttons[0].disable()
      self.listButtons.buttons[3].disable()
      

    if self.experiment == None:
      measurementLists = self.nmrProject.sortedMeasurementLists()
    else:
      measurementLists = []
      lists = self.nmrProject.sortedMeasurementLists()
      for measurementList in lists:
        if self.experiment in measurementList.experiments:
          measurementLists.append(measurementList)
    
    objectList  = []
    textMatrix  = []
    for measurementList in measurementLists:
      objectList.append(measurementList)
      names = ['%s' % (e.name) for e in measurementList.experiments]
      if len(names) > 6:
        names.insert(int(len(names)/2), '\n')
      experimentText = ' '.join(names)

      """otherInfoText  = ''
      for attribute in measurementList.metaclass.getAllAttributes():
        name = attribute.name
        otherInfoText += '%s:%s ' % (name,getattr(measurementList, name))

      for role in measurementList.metaclass.getAllRoles():
        name = role.name
        otherInfoText += '%s:%s ' % (name,getattr(measurementList, name))"""
          
      datum = [measurementList.serial,
               measurementList.className,
               measurementList.name,
               len(measurementList.measurements),
               measurementList.unit,
               experimentText,
               #otherInfoText,
               measurementList.details]
 
      textMatrix.append(datum)
    
    if not objectList:
      textMatrix.append([])
    self.listsMatrix.update(objectList=objectList, textMatrix=textMatrix)

    self.waitingList = False

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  # Measurement functions
  
  def deleteMeasurements(self):
  
    measurements = self.measurementsMatrix.currentObjects
    if measurements:
      
      className = self.measurementListB.className
      if className == 'ShiftList':
        msg = 'Really delete %s chemical shifts?' % len(measurements)
        if not showOkCancel('Confirm', msg, parent=self):
          return
          
        untouched = set()
      
        for shift in measurements:
          for peakDim in shift.peakDims:
            shiftList = peakDim.peak.peakList.dataSource.experiment.shiftList
            if shiftList is self.measurementListB:
              untouched.add(shift)
              break
          else:
            shift.delete()

        if untouched:
          if len(untouched) > 1:
            msg = 'Could not delete %d shifts because they are still assigned to peaks' % len(untouched)
         
          else:
            msg = 'Could not delete shift because it is still assigned to peaks'

          showWarning('Warning', msg , parent=self)
      
      else:
        msg = 'Really delete %d %s measurements?' % (len(measurements), className[:-4])
        if not showOkCancel('Confirm', msg, parent=self):
          return
          
        for measurement in measurements:
          measurement.delete()

  def recalculateMeasurements(self):
  
    measurements = self.measurementsMatrix.currentObjects
    if measurements and (self.measurementListB.className == 'ShiftList'):
      for shift in measurements:
        averageShiftValue(shift)

  def showPeaks(self, *event):
  
    measurements = self.measurementsMatrix.currentObjects
    if measurements:
      peaksDict = {}
      for measurement in measurements:
        for peak in measurement.peaks:
          peaksDict[peak] = 1
 
      peaks = peaksDict.keys()
      if len(peaks) > 0:
        self.guiParent.viewPeaks(peaks)

  def showResonances(self, *event):
  
    measurements = self.measurementsMatrix.currentObjects
    if measurements:
      resonances = set()
      
      for measurement in measurements:
        if hasattr(measurement, 'resonances'):
          resonances.update(measurement.resonances)
        else:
          resonances.add(measurement.resonance)
      
      if resonances:
        self.guiParent.viewSelectedResonances(resonances)
 
  def setMeasurementList(self, measurementList):
  
    if measurementList is not self.measurementListB:
      self.measurementListB = measurementList
      self.updateMeasurementsAfter()
  
  def getMeasurementMerit(self, measurement):

    if measurement:
      self.measurementMeritEntry.set(measurement.figOfMerit)
  
  def setMeasurementMerit(self, event):

    value = self.measurementMeritEntry.get()
    if value is not None:
      self.measurement.setFigOfMerit( max(0.0, min(1.0, float(value)) ) )
  
  def getMeasurementDetails(self, measurement):

    if measurement and measurement.details:
      self.measurementDetailsEntry.set(measurement.details)
  
  def setMeasurementDetails(self, event):

    text = self.measurementDetailsEntry.get()
    if text and text != ' ':
      self.measurement.setDetails( text )

  def selectMeasurementCell(self, object, row, col):
  
    self.measurement = object
    if self.measurement:
      self.measurementButtons.buttons[0].enable()
      self.measurementButtons.buttons[2].enable()
      if self.measurementListB.className == 'ShiftList':
        self.measurementButtons.buttons[1].enable()
      else:
        self.measurementButtons.buttons[1].disable()
 
  def getMeasurementListName(self, measurementList):
  
    name = measurementList.name or \
           '%d:%s' % (measurementList.serial,measurementList.className)

    return name

  def updateListPulldown(self, null=None):

    index = 0 
    measurementLists = self.nmrProject.sortedMeasurementLists()
    nameFunc = self.getMeasurementListName
    names = [nameFunc(ml) for ml in measurementLists]
    measurementList = self.measurementListB
    
    if measurementLists:
      if measurementList not in measurementLists:
        measurementList = measurementLists[0]
        
      index = measurementLists.index(measurementList)
        
    else:
      measurementList = None
    
    if self.measurementListB is not measurementList:
      self.measurementListB = measurementList
      self.updateMeasurementsAfter()

    self.listPulldown.setup(names, measurementLists, index)   



  def updateMeasurementsAfter(self, *opt):

    if self.waitingMeasurement:
      return
    else:
      self.waitingMeasurement = True
      self.after_idle(self.updateMeasurements)
 
  def updateMeasurements(self, measurementList=None):

    headingList = ['#','Resonance','Value','SD','Chain',
                   'Isotope','Fig of\nMerit','Peaks','Details'] 
      
    if measurementList is not None:
      self.measurementListB = measurementList
      self.updateListPulldown()
      
    if self.measurementListB is None:
      self.measurement = None
      measurements = []
    else:
      if self.measurementListB.unit:
        headingList[2] = 'Value\n(%s)' % (self.measurementListB.unit)
      measurements = self.measurementListB.sortedMeasurements()
      if measurements and hasattr(measurements[0], 'resonances'):
        headingList[1] = 'Resonances'
        headingList[5] = 'Isotopes'
 
    if self.measurement:
      self.measurementButtons.buttons[0].enable()
      self.measurementButtons.buttons[2].enable()
      if self.measurementListB.className == 'ShiftList':
        self.measurementButtons.buttons[1].enable()
      else:
        self.measurementButtons.buttons[1].disable()
    else:
      self.measurementButtons.buttons[0].disable()
      self.measurementButtons.buttons[1].disable()
      self.measurementButtons.buttons[2].disable()
  
    objectList  = []
    textMatrix  = []
    i = 0
    for measurement in measurements:
      i += 1
      resonanceText = ''
      isotopeText   = ''
      chainText     = ''
      if hasattr(measurement, 'resonances'):
        resonanceText = ' '.join([makeResonanceGuiName(r) for r in measurement.resonances])
        isotopeText = ' '.join([r.isotopeCode for r in measurement.resonances])
        residueDict = {}
        for resonance in measurement.resonances:
          if resonance.resonanceSet:
            residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
            residueDict[residue] = None
        
        keys = residueDict.keys()
        if len(keys) == 1:
          residue    = keys[0]
          residueNum = residue.seqCode
          chainText  = '%s:%s' % (residue.chain.molSystem.code,residue.chain.code)
      
      elif hasattr(measurement, 'resonance'):
        resonance     = measurement.resonance
        resonanceText = makeResonanceGuiName(resonance)
        isotopeText   = measurement.resonance.isotopeCode
        if resonance.resonanceSet:
          residue   = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
          chainText = '%s:%s' % (residue.chain.molSystem.code,residue.chain.code)
              
      datum = [i,
               resonanceText,
               measurement.value,
               measurement.error,
               chainText,
               isotopeText,
               measurement.figOfMerit,
               len(measurement.peaks),
               measurement.details or ' ']
 

      objectList.append(measurement)
      textMatrix.append(datum)
      
    if not objectList:
      textMatrix.append([])
                
    tipTexts = self.measurementsMatrix.tipTexts # unchanged, despite variable headings

    self.measurementsMatrix.update(headingList=headingList, objectList=objectList,
                                   textMatrix=textMatrix, tipTexts=tipTexts)

    self.waitingMeasurement = False


 
  
