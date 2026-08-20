
"""
======================COPYRIGHT/LICENSE START==========================

AcqProcParsEditPopup.py: Editor for acquisition/processing parameters.

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
import Tkinter, string, os

from memops.universal.Util import returnInt
from memops.universal.Util import returnFloat
from memops.universal.Util import isBigEndian
from memops.universal.Io import joinPath

from memops.general.Application import createAppData

from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Label import Label
from memops.gui.BasePopup import BasePopup
from memops.gui.Separator import Separator
from memops.gui.ToggleButton import ToggleButton
from memops.gui.Entry import Entry
from memops.gui.CheckButton import CheckButton
from memops.gui.Button import Button
from memops.gui.HelpPopup import showHelpUrl
from memops.gui.ObjectsTable import ObjectsTable

from ccp.general.Io import getDataSourceFileName, setDataSourceDataStore

from ccpnmr.format.gui.FormatFilePopup import FormatFilePopup
from ccpnmr.format.gui.ObjectButton import ObjectButton
from ccpnmr.format.gui.WindowFunctionSelect import WindowFunctionSelect
from ccpnmr.format.gui.ProcFilePopup import ProcFilePopup
from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

from ccpnmr.format.general.Constants import varianAcqOrdInfo
from ccpnmr.format.general.Constants import ccpNmr_kw, importGeneral_kw
from ccpnmr.format.general.Constants import IoSetupList

from ccpnmr.format.general.Util import getAzaraBlockSizes
from ccpnmr.format.general.Util import setDataDimRefs
from ccpnmr.format.general.Util import createPpmFreqDataSource
from ccpnmr.format.general.Io import getHelpUrlDir

from memops.gui.MessageReporter import showError

import memops.api.Implementation as Implementation
import ccp.api.nmr.Nmr as Nmr

class AcqProcParsEditPopup(TemporaryBasePopup):
 
  help_urls = {}
  help_urls['main'] = joinPath(getHelpUrlDir(),'AcqProcParsEditMain.html')
  
  def __init__(self, parent, experiment = None, newFidDataSources = None):
     
    #
    # General setup
    #
    
    self.objects = {}
    
    if experiment:
      self.experiment = experiment
    elif hasattr(parent,'experiment') and parent.experiment:
      self.experiment = parent.experiment
    else:
      showError("No experiment","No experiment is available... aborting edit.")

    self.applicationName = ccpNmr_kw
    self.scriptName = 'AcqProcParsEditPopup'
    
    if self.experiment.rawData:
      self.rawData = self.experiment.rawData
      self.rawDataDims = self.experiment.rawData.sortedDataDims()
    else:
      self.rawData = None
    
    self.acqOrdMenu = None

    self.guiParent = parent

    #
    # Setup expdimrefs (only FIRST ONE currently!)
    #

    self.expDimRefs = []

    for expDim in self.experiment.sortedExpDims():
      self.expDimRefs.append(expDim.findFirstExpDimRef())


    # Setup acquisition order stuff (from general.settings variable)

    self.acqOrd = varianAcqOrdInfo

    #
    # See if a format dataSource(s) already exist(s)
    #

    self.freqDataDims = {}
    self.freqFormats = {}
    
    for dataSource in self.experiment.findAllDataSources(dataType = 'processed'):
      self.freqFormats[dataSource.name] = dataSource
      self.freqDataDims[dataSource.name] = list(dataSource.sortedDataDims())
          
    #
    # Create info for data sources that have to be newly created...
    #
    
    if newFidDataSources:
      for formatName in newFidDataSources:
        if not self.freqFormats.has_key(formatName):
          dataSource = createPpmFreqDataSource(self.experiment,formatName,'processed',self.experiment.numDim)          
          self.freqFormats[formatName] = dataSource
          self.freqDataDims[formatName] = list(dataSource.sortedDataDims())
    
    #
    # Get the cases where writing processing data is possible...
    #
    
    self.writeFormats = []
    self.writeFormatAvailable = False
    
    for infoList in IoSetupList:
    
      if infoList[0] == 'processingPars':
      
        IoDefs = infoList[3]

        for formatName in IoDefs.keys():

          if IoDefs[formatName][2]:

            #setattr(self,formatName,formatName)
            self.writeFormats.append(formatName)
            self.writeFormatAvailable = True
    
    #
    # Initialize parsList dict (for creation of objectTables)
    #

    self.parsList = {}

    #
    # Create basic window
    #

    TemporaryBasePopup.__init__(self,parent = parent, title = "Project '%s': " % self.experiment.root.name + 'Edit parameters', modal = False, transient=True)

    self.popups = {}
 
  def body(self, master):
  
    if hasattr(self.parent,'top'):
      master.top = self.parent.top
    else:
      master.top = self.top
  
    #
    # Local settings for ToggleButton
    #
        
    self.toggledDict = { 'fg': 'blue', 'relief': 'sunken' }
    self.untoggledDict = { 'fg': 'red', 'relief': 'groove' }
    
    self.columnsTotal = 3

    #
    # Grid config
    #

    master.grid_rowconfigure(0,weight = 1)
    for i in range(0,self.columnsTotal):
      master.grid_columnconfigure(i,weight = 1)

    #
    # Popup window
    #

    row = 0

    label = Label(master, text= "Experiment '%s'" % self.experiment.name, fg = 'blue')
    label.grid(row=row, column=0, columnspan = self.columnsTotal, sticky=Tkinter.EW)
   
   
    #
    # Raw data spectrum location
    #
    
    if self.rawData:
    
      row += 1    

      # Spectrum location
      label = Label(master, text= "Time domain data file:")
      label.grid(row=row, column=0, sticky=Tkinter.E)

      fileName = getDataSourceFileName(self.experiment.rawData)

      self.selectFileButton = Tkinter.Button(master, text = fileName, command = (lambda : self.selectFile()))
      self.selectFileButton.grid(row=row, column=1, columnspan = 2, sticky=Tkinter.W)

      row += 1    

      #
      # Sensitivity enhanced (all but first button)
      #

      label = Label(master, text= "Sensivity enhanced dimension:")
      label.grid(row=row, column=0, sticky=Tkinter.E)

      self.sensEnhMenu = self.setSensEnhMenu(master)
      self.sensEnhMenu.grid(row=row, column=1, sticky=Tkinter.W)

      #
      # Varian acquisition order
      # 

      if self.experiment.numDim > 2 and self.experiment.rawData.dataStore.fileType == 'Varian':
        row += 1   
        # Varian acqord (if Varian data)
        label = Label(master, text= "Varian acquisition order:")
        label.grid(row=row, column=0, sticky=Tkinter.E)

        self.acqOrdMenu = self.setAcqOrdMenu(master)
        self.acqOrdMenu.grid(row=row, column=1, sticky=Tkinter.W)


      if self.writeFormatAvailable:

        row += 1    

        #
        # Swap data
        #

        label = Label(master, text= "Swap data (computer dependent):")
        label.grid(row=row, column=0, sticky=Tkinter.E)

        self.swapData = self.setSwapData(master)
        self.swapData.grid(row=row, column=1, sticky=Tkinter.W)

      row += 1    

      #
      # Help button
      #

      button = Button(master, text= "Help", command = lambda top=self.guiParent, url=self.help_urls['main']: showHelpUrl(self.guiParent, self.help_urls['main']))
      button.grid(row=row, column=0, columnspan = self.columnsTotal, sticky=Tkinter.EW)

      # Empty row
      row += 1   

      label = Label(master, text= "")
      label.grid(row=row, column=0, sticky=Tkinter.E)

    if self.rawData:

      #
      # Time domain data table
      #

      row += 1

      separator = Separator(master, height = 3)
      separator.setColor('black', bgColor = 'black')
      separator.grid(row = row, column = 0, columnspan = self.columnsTotal, sticky = Tkinter.EW)

      row += 1    

      fidParRow = row + 1
      fidParTable = 'fidPars'

      self.fidParsButton = ToggleButton(master, text= "Time domain data",
                                                toggledDict = self.toggledDict,
                                                untoggledDict = self.untoggledDict,
                                                command = lambda: self.toggleTable(master,fidParRow,fidParTable))

      self.fidParsButton.grid(row=row, column=0, columnspan = self.columnsTotal)

      row += 1    

      self.fidParsTable = self.setupObjectsTable(master,fidParTable)
    
    #
    # Reference data table
    #

    row += 1
    
    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, column = 0, columnspan = self.columnsTotal, sticky = Tkinter.EW)
    
    row += 1    

    refParRow = row + 1
    refParTable = 'refPars'

    self.refParsButton = ToggleButton(master, text= "Reference data",
                                              toggledDict = self.toggledDict,
                                              untoggledDict = self.untoggledDict,
                                              command = lambda: self.toggleTable(master,refParRow,refParTable))

    self.refParsButton.grid(row=row, column=0, columnspan = self.columnsTotal)

    row += 1

    self.refParsTable = self.setupObjectsTable(master,refParTable)
    
    #
    # Now do the processed data stuff...
    #
    
    formatNames = self.freqFormats.keys()
    formatNames.sort()
    
    for formatName in formatNames:
    
      if formatName == importGeneral_kw:
        continue

      row += 1

      separator = Separator(master, height = 3)
      separator.setColor('black', bgColor = 'black')
      separator.grid(row = row, column = 0, columnspan = self.columnsTotal, sticky = Tkinter.EW)

      row += 1
      
      row = self.createFreqDataSourceTable(row,formatName,formatNames,master)

    #
    # Exit button
    #


    row += 1
    
    separator = Separator(master, height = 3)
    separator.setColor('black', bgColor = 'black')
    separator.grid(row = row, column = 0, columnspan = self.columnsTotal, sticky = Tkinter.EW)

    row += 1    

    button = Button(master, text= "Exit/continue", command = self.ok)
    button.grid(row=row, column=0, columnspan = self.columnsTotal, sticky=Tkinter.EW)

  def apply(self):
    
    return True

  """
  
  Initialization functions

  """

  def createFreqDataSourceTable(self,row,formatName,formatNames,master):
  
    refRow = row + 1

    widgetsName = formatName + 'Widgets'

    setattr(self,widgetsName,[])

    button = ToggleButton(master, text= "%s processing" % formatName,
                                  toggledDict = self.toggledDict,
                                  untoggledDict = self.untoggledDict,
                                  command = lambda : self.toggleTable(master,refRow,formatName))

    button.grid(row=row, column=0, columnspan = self.columnsTotal)

    setattr(self,formatName + 'Button',button)

    row += 1    

    # This is where the objectsTable will appear (not created unless requested)
    
    if self.rawData and formatName in self.writeFormats:

      row += 1

      button = Tkinter.Button(master, text = "Write processing files", command = (lambda : self.procFileWrite(formatName)))
      getattr(self,widgetsName).append([button,row,0,self.columnsTotal])

    row += 1

    if len(formatNames) > 1:
    
      for otherFormat in formatNames:
        if otherFormat != formatName and otherFormat != importGeneral_kw:

          button = Tkinter.Button(master, text = "Copy data to %s" % otherFormat, command = (lambda : self.setupDataSource(otherFormat, fromFormat = formatName)))
          getattr(self,widgetsName).append([button,row,0,self.columnsTotal])
          row += 1

    return row  

  def setupDataSource(self,toFormat, fromFormat = None):

    x = 0
 
    # Just checking:
    
    dataSource = self.experiment.findFirstDataSource(name = toFormat)
    initialized = False
    
    if dataSource:
      createNew = 0
      
      appData = dataSource.findFirstApplicationData(application = self.applicationName, keyword = self.scriptName, value = True)

      if not fromFormat and appData:
        print "  Warning: %s datasource already exists..." % toFormat
        return
        
      if appData:
        initialized = True
                
    else:
      createNew = 1
    
    
    self.updateFromWidget()
    
    exp = self.experiment
    
    #
    # Application data to be set
    #
    #
    # TODO TODO TODO This is tricky - if the fid data is not complex then this has to be zero?
    # Also should really set the isComplex flag in the data model for this. The number of complex
    # points for FT should not change depending on the isComplex flag, but the output for Azara and
    # nmrPipe should!!! For now assuming imaginary data...
    # ALSO see isComplex flag setting if dataSource created..
    #
    
    freqDimApplDataItems = [['deleteImaginary',1]
                           ]

    #
    # Get reference data
    #

    if not fromFormat:
      if self.freqDataDims.has_key(importGeneral_kw):
        refFormat = importGeneral_kw
      elif self.freqDataDims.keys():
        refFormat = self.freqDataDims.keys()[0]
      else:
        print "  Error: no valid reference information for setting up %s datasource" % toFormat
        refFormat = None

    else:
      refFormat = fromFormat
    
    
    if refFormat:
      refDims =  self.freqDataDims[refFormat]
    else:
      refDims = None
    
    if createNew:
      print "  Creating %s dataSource based on %s..." % (toFormat,refFormat)

      #
      # Setup new dataSource
      #
    
      dsFreq = Nmr.DataSource(exp,name = toFormat, numDim = exp.numDim, dataType = 'processed')

    elif initialized:

      print "  Copying %s dataSource information to %s..." % (refFormat,toFormat)
      dsFreq = exp.findFirstDataSource(name = toFormat)

    else:
    
      print "  Updating %s dataSource..." % (toFormat)
      dsFreq = exp.findFirstDataSource(name = toFormat)


    #
    # Loop over freqDataDims
    #

    for i in range(0,exp.numDim):
      
      expDim = exp.sortedExpDims()[i]
      
      notSet = False
     
      if refDims:
     
        refDim = refDims[i]

        isComplex = refDim.isComplex         # Default is deleteImaginary ON! TODO: set this correctly
        valuePerPoint = refDim.valuePerPoint
        numPoints = refDim.numPoints
        numPointsOrig = refDim.numPointsOrig

        # Phase sign change between nmrPipe and azara..
        signChangeList = [['azara','nmrPipe']]
        formatList = [toFormat,fromFormat]
        formatList.sort()

        if formatList in signChangeList:
          sign = -1
        else:
          sign = 1

        for phaseKeywd in ('phase0','phase1'):
          if getattr(refDim,phaseKeywd):
            globals()[phaseKeywd] = sign * getattr(refDim,phaseKeywd)
          else:
            globals()[phaseKeywd] = 0.0

        if refDim.unit:
          unit = refDim.unit
        else:
          unit = 'Hz'
          notSet = True

      if createNew:
        freqDataDim = Nmr.FreqDataDim(dsFreq, expDim = expDim, dim = expDim.dim, isComplex = isComplex,
                                      valuePerPoint = valuePerPoint,numPoints = numPoints, unit = unit,
                                      numPointsOrig = numPointsOrig, phase0 = phase0, phase1 = phase1)

      else:
        freqDataDim = dsFreq.findFirstDataDim(expDim = expDim)
        
        if fromFormat or notSet:
          freqDataDim.valuePerPoint = valuePerPoint
          freqDataDim.numPoints = numPoints
          freqDataDim.unit = unit
          freqDataDim.numPointsOrig = numPointsOrig

          freqDataDim.phase0 = phase0
          freqDataDim.phase1 = phase1

      #
      # Set/copy application data
      #
      
      for applDataItem in freqDimApplDataItems:
        (keyword,value) = applDataItem
        
        if refDims:
          refApplData = refDim.findFirstApplicationData(keyword = keyword)

          if refApplData:
            value = refApplData.value
        
        #
        # Check if exists, create if not. Overwrite only if createNew!
        #
        
        applData = freqDataDim.findFirstApplicationData(keyword = keyword)
        
        if not applData:
          createAppData(freqDataDim,self.applicationName,keyword,value)

        elif createNew:
          applData.value = value

    #
    # Set self freqDataDims only at this stage: otherwise problems with above treatment
    #
    
    self.freqDataDims[toFormat] = dsFreq.sortedDataDims()
    
    if not fromFormat:
      
      #
      # Set default reference values if no datadimrefs exist...
      #
      
      if not self.freqDataDims[toFormat][0].dataDimRefs:
      
        setDataDimRefs(dsFreq)

    else:
      
      #
      # Copy reference values
      #
  
      for i in range(0,exp.numDim):
          
        expDim = exp.sortedExpDims()[i]
        refDataDim = refDims[i]
        dataDim = self.freqDataDims[toFormat][i]
        
        for expDimRef in expDim.expDimRefs:
          refDataDimRef = expDimRef.findFirstDataDimRef(dataDim = refDataDim)
          refRefPoint = refDataDimRef.refPoint
          refRefValue = refDataDimRef.refValue
        
          if createNew:
            dataDimRef = Nmr.DataDimRef(dataDim,refPoint = refRefPoint, refValue = refRefValue, expDimRef = expDimRef)
          else:
            dataDimRef = expDimRef.findFirstDataDimRef(dataDim = dataDim)
            dataDimRef.refPoint = refRefPoint
            dataDimRef.refValue = refRefValue
    
    #
    # Also add some special dataSource, dataDimRef attributes and application data
    # for a newly created dataSource
    #

    if toFormat == 'azara':

      dataSourceItems = {'isBigEndian': isBigEndian(),
                         'nByte': 4,
                         'numberType': "float",
                         'fileHeaderSize': 0}

      # This info will be set automatically when window updated
      #freqDataDimItems = {'blockSize': getAzaraBlockSizes(self.rawDataDims)}

      dataDimRefApplDataItems = []

    else:
      dataSourceItems = {}
      freqDataDimItems = {}
      dataDimRefApplDataItems = [['label', lambda x = value: x.expDimRef.isotopeCodes[0]]
                                ]
    #
    # DataSource attributes
    #

    for keyword in dataSourceItems.keys():
      
      if createNew or not hasattr(dsFreq,keyword):
      
        setattr(dsFreq,keyword,dataSourceItems[keyword])

    """
    #
    # freqdatadim attributes
    #

    for i in range(0,len(self.freqDataDims[toFormat])):

      freqDataDim = self.freqDataDims[toFormat][i]

      for keyword in freqDataDimItems.keys():
        #
        # If it's a list then get values depending on expDim!
        #

        if type(freqDataDimItems[keyword]) == type([]):
          value = freqDataDimItems[keyword][i]
        else:
          value = freqDataDimItems[keyword]

        setattr(freqDataDim,keyword,value)        
    """  
    #
    # DataDimRef application data - keep or overwrite?
    #

    for (keyword,function) in dataDimRefApplDataItems:

      for freqDataDim in self.freqDataDims[toFormat]:

        for dataDimRef in freqDataDim.dataDimRefs:
        
          applData = dataDimRef.findFirstApplicationData(application = toFormat, keyword = keyword)

          if createNew or not applData:
            createAppData(dataDimRef,toFormat,keyword,function(dataDimRef))
            #print toFormat, keyword

    if not createNew:
    
      self.updateFromDataModel()
    
    #
    # TODO: Should I really always remove this? Or always keep a link to the original info?      
    #
    
    if refFormat == importGeneral_kw:
      
      del self.freqDataDims[importGeneral_kw]
      self.experiment.findFirstDataSource(name = refFormat).delete()

    dsFreq.addApplicationData(Implementation.AppDataBoolean(application = self.applicationName, keyword = self.scriptName, value = True))

  def createParsLists(self,parsListKey):
    
    #
    # Setup data for objectTables...
    #
    
    dim = 0
    value = 0
      
    # Max 8k points for 2d and over
    powerTwoList = ['4','8','16','32','64','128','256','512','1024','2048','4096','8192']
    
    if self.experiment.numDim == 1:
      powerTwoList.extend(['16384','32768','65536'])
    
    allDims = range(0,self.experiment.numDim)

    #
    # ObjectTable definitions
    #
    # column[0]: row label
    # column[1]: widget type
    # column[2]: function to read data into widget
    # column[3]: function to write data from widget
    # column[4]: dimensions in which widget appears
    # column[5]: dictionary for additional widget settings
    #
    
    if parsListKey == 'fidPars':
  
      specWidthUnits = ""
    
      for rawDataDim in self.rawDataDims:
        specWidthUnits += "1/" + rawDataDim.unit + "|"
      
      specWidthUnits = specWidthUnits[:-1]
      
      self.help_urls['fidPars'] = joinPath(getHelpUrlDir(),'AcqProcParsEditFidPars.html')
      self.parsList['fidPars'] = [
                     
                     ['Complex data',
                      CheckButton,
                      lambda x = dim: self.rawDataDims[x].isComplex,
                      lambda x = dim, value = value: setattr(self.rawDataDims[x],'isComplex',returnInt(value)),
                      allDims,
                      {}],
                     ['Number of valid FID points',
                      Entry,
                      lambda x = dim: self.rawDataDims[x].numPointsValid,
                      lambda x = dim, value = value: setattr(self.rawDataDims[x],'numPointsValid',returnInt(value)),
                      allDims,
                      {}],
                     ["Spectral width \n(in %s)" % specWidthUnits,
                      Entry,
                      lambda x = dim: self.rawDataDims[x].spectralWidth,
                      lambda x = dim, value = value: setattr(self.rawDataDims[x],'valuePerPoint',returnFloat((self.rawDataDims[x].isComplex + 1) / (2 * returnFloat(value)))),
                      allDims,
                      {}],
                     ["Value per point",
                      Label,
                      lambda x = dim: str(self.rawDataDims[x].valuePerPoint) + " " + self.rawDataDims[x].unit,
                      lambda x = dim, value = value: self.void(),
                      allDims,
                      {}],
                     ['TPPI (alternate sign)',
                      CheckButton,
                      lambda x = dim: self.rawDataDims[x].alternateSign,
                      lambda x = dim, value = value: setattr(self.rawDataDims[x],'alternateSign',returnInt(value)),
                      allDims,
                      {}],
                     ['Flip axis (negate imag.)',
                      CheckButton,
                      lambda x = dim: self.rawDataDims[x].negateImaginary,
                      lambda x = dim, value = value: setattr(self.rawDataDims[x],'negateImaginary',returnInt(value)),
                      allDims,
                      {}]
                    ]    
                    
    elif parsListKey == 'refPars':
    
      self.help_urls['refPars'] = joinPath(getHelpUrlDir(),'AcqProcParsEditRefPars.html')
      self.parsList['refPars'] = [
                     
                     ["",
                      Label,
                      lambda x = dim: "Reference %s" % self.expDimRefs[x].serial,
                      lambda x = dim, value = value: self.emptyString(),
                      allDims,
                      {}],
                     ["Cycle reference",   # If changed, have to propagate to cycleExpDimRef
                      CheckButton,
                      lambda x = dim: self.void(),
                      lambda x = dim, value = value: self.void(),
                      allDims,
                      {'callback': self.cycleExpDimRef}],
                     ["Carrier frequency (in MHz)",
                      Entry,
                      lambda x = dim: self.expDimRefs[x].sf,
                      lambda x = dim, value = value: setattr(self.expDimRefs[x],'sf',returnFloat(value)),
                      allDims,
                      {}],
                     ["Base frequency (in MHz)",
                      Entry,
                      lambda x = dim: self.expDimRefs[x].baseFrequency,
                      lambda x = dim, value = value: setattr(self.expDimRefs[x],'baseFrequency',returnFloat(value)),
                      allDims,
                      {}],
                     ["Isotope",
                      PulldownMenu,
                      lambda x = dim: self.expDimRefs[x].isotopeCodes[0],
                      lambda x = dim, value = value: setattr(self.expDimRefs[x],'isotopeCodes',[value]),
                      allDims,
                      {'entries': ['1H','2H','13C','15N','31P'] }],
                     ['Folded',
                      CheckButton,
                      lambda x = dim: returnInt(self.expDimRefs[x].isFolded),
                      lambda x = dim, value = value: setattr(self.expDimRefs[x],'isFolded',returnInt(value)),
                      allDims,
                      {}]
                    ]

    # TODO: truncated FT for processed
    
    elif self.rawData and parsListKey == 'nmrPipe':
      
      # TODO Pipe: which aq2D? 
      nmrPipe = 'nmrPipe'

      self.help_urls[nmrPipe] = joinPath(getHelpUrlDir(),'AcqProcParsEditNmrPipe.html')
      self.parsList[nmrPipe] = [
                     
                     ['Number of complex points for FT',   #TODO: change for updateFromOtherWidgets if changed here!
                      PulldownMenu,
                      lambda x = dim: str(self.freqDataDims[nmrPipe][x].numPointsOrig),
                      lambda x = dim, value = value: self.changeFreqNumPoints(nmrPipe,x,value),
                      allDims,
                      {'entries': powerTwoList, 'callback': (lambda x = dim, y = dim: self.updateOwnWidget('Number of complex points for FT'))}],
                     ["Value per point", # If changed, have to propagate to changeFreqNumPoints
                      Label,
                      lambda x = dim: str(self.freqDataDims[nmrPipe][x].valuePerPoint) + " " + self.freqDataDims[nmrPipe][x].unit,
                      lambda x = dim, value = value: self.void(),
                      allDims,
                      {}],
                     ['Phase 0',
                      Entry,
                      lambda x = dim: self.freqDataDims[nmrPipe][x].phase0,
                      lambda x = dim, value = value: setattr(self.freqDataDims[nmrPipe][x],'phase0',returnFloat(value)),
                      allDims,
                      {}],
                     ['Phase 1',
                      Entry,
                      lambda x = dim: self.freqDataDims[nmrPipe][x].phase1,
                      lambda x = dim, value = value: setattr(self.freqDataDims[nmrPipe][x],'phase1',returnFloat(value)),
                      allDims,
                      {}],
                     ['Window function', # If changed, have to propagate to def windowFunction!
                      ObjectButton,
                      lambda x = dim: self.windowFunction(nmrPipe,dim = x, getValue = 1),
                      lambda x = dim, value = value: self.void(),
                      allDims,
                      {'command': lambda x = dim: self.windowFunction(nmrPipe),
                       'toggledDict': self.toggledDict,
                       'untoggledDict': self.untoggledDict}],
                     ['Delete imaginary data',
                      CheckButton,
                      lambda x = dim: self.getDeleteImaginary(nmrPipe,x),
                      lambda x = dim, value = value: setattr(self.freqDataDims[nmrPipe][x].findFirstApplicationData(application = self.applicationName, keyword = 'deleteImaginary'),'value',value),
                      allDims,
                      {}],
                     ["",
                      Label,
                      lambda x = dim: "Reference %s (%s)" % (self.expDimRefs[x].serial, self.expDimRefs[x].isotopeCodes[0]),
                      lambda x = dim, value = value: self.emptyString(),
                      allDims,
                      {}],
                     ["%s label" % nmrPipe,
                      Entry,
                      lambda x = dim: self.freqDataDims[nmrPipe][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]).findFirstApplicationData(application = nmrPipe, keyword = 'label').value,
                      lambda x = dim, value = value: setattr(self.freqDataDims[nmrPipe][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]).findFirstApplicationData(application = nmrPipe, keyword = 'label'),'value',value),
                      allDims,
                      {}],
                     ["Reference point",  # If changed, have to propagate to changeFreqNumPoints
                      Entry,
                      lambda x = dim: self.freqDataDims[nmrPipe][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]).refPoint,
                      lambda x = dim, value = value: setattr(self.freqDataDims[nmrPipe][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]),'refPoint',returnFloat(value)),
                      allDims,
                      {}],
                     ["Reference value",
                      Entry,
                      lambda x = dim: self.freqDataDims[nmrPipe][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]).refValue,
                      lambda x = dim, value = value: setattr(self.freqDataDims[nmrPipe][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]),'refValue',returnFloat(value)),
                      allDims,
                      {}]
                    ]

    elif self.rawData and parsListKey == 'azara':
      
      azara = 'azara'

      self.help_urls[azara] = joinPath(getHelpUrlDir(),'AcqProcParsEditAzara.html')
      self.parsList[azara] = [
                     
                     ['Number of complex points for FT', #TODO: change for updateOwnWidget if changed here!
                      PulldownMenu,
                      lambda x = dim: str(self.freqDataDims[azara][x].numPointsOrig),
                      lambda x = dim, value = value: self.changeFreqNumPoints(azara,x,value),
                      allDims,
                      {'entries': powerTwoList, 'callback': (lambda x = dim, y = dim: self.updateOwnWidget('Number of complex points for FT'))}],
                     ["Value per point", # If changed, have to propagate to changeFreqNumPoints
                      Label,
                      lambda x = dim: str(self.freqDataDims[azara][x].valuePerPoint) + " " + self.freqDataDims[azara][x].unit,
                      lambda x = dim, value = value: self.void(),
                      allDims,
                      {}],
                     ['Phase 0',
                      Entry,
                      lambda x = dim: self.freqDataDims[azara][x].phase0,
                      lambda x = dim, value = value: setattr(self.freqDataDims[azara][x],'phase0',returnFloat(value)),
                      allDims,
                      {}],
                     ['Phase 1',
                      Entry,
                      lambda x = dim: self.freqDataDims[azara][x].phase1,
                      lambda x = dim, value = value: setattr(self.freqDataDims[azara][x],'phase1',returnFloat(value)),
                      allDims,
                      {}],
                     ['Window function', # If changed, have to propagate to def windowFunction!
                      ObjectButton,
                      lambda x = dim: self.windowFunction(azara,dim = x, getValue = 1),
                      lambda x = dim, value = value: self.void(),
                      allDims,
                      {'command': lambda x = dim: self.windowFunction(azara),
                       'toggledDict': self.toggledDict,
                       'untoggledDict': self.untoggledDict}],
                     ['Delete imaginary data',
                      CheckButton,
                      lambda x = dim: self.getDeleteImaginary(azara,x),
                      lambda x = dim, value = value: setattr(self.freqDataDims[azara][x].findFirstApplicationData(application = self.applicationName, keyword = 'deleteImaginary'),'value',value),
                      allDims,
                      {}],
                     ["",
                      Label,
                      lambda x = dim: "Reference %s (%s)" % (self.expDimRefs[x].serial, self.expDimRefs[x].isotopeCodes[0]),
                      lambda x = dim, value = value: self.emptyString(),
                      allDims,
                      {}],
                     ["Reference point",  # If changed, have to propagate to changeFreqNumPoints
                      Entry,
                      lambda x = dim: self.freqDataDims[azara][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]).refPoint,
                      lambda x = dim, value = value: setattr(self.freqDataDims[azara][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]),'refPoint',returnFloat(value)),
                      allDims,
                      {}],
                     ["Reference value",
                      Entry,
                      lambda x = dim: self.freqDataDims[azara][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]).refValue,
                      lambda x = dim, value = value: setattr(self.freqDataDims[azara][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]),'refValue',returnFloat(value)),
                      allDims,
                      {}]
                    ]
                    
    elif parsListKey != importGeneral_kw:
      
      generic = parsListKey

      self.help_urls[generic] = joinPath(getHelpUrlDir(),'AcqProcParsEditGeneric.html')
      self.parsList[generic] = [
                     
                     ['Number of complex points',   #TODO: change for updateFromOtherWidgets if changed here!
                      Entry,
                      lambda x = dim: str(self.freqDataDims[generic][x].numPointsOrig),
                      lambda x = dim, value = value: self.changeFreqNumPoints(generic,x,value),
                      allDims,
                      {}],
                     ["Spectral width",
                      Entry,
                      lambda x = dim: self.freqDataDims[generic][x].spectralWidthOrig,
                      lambda x = dim, value = value: setattr(self.freqDataDims[generic][x],'valuePerPoint',returnFloat((2 * returnFloat(value)) / ((self.freqDataDims[generic][x].isComplex + 1) * self.freqDataDims[generic][x].numPointsOrig))),
                      allDims,
                      {}],
                     ["Value per point",
                      Label,
                      lambda x = dim: str(self.freqDataDims[generic][x].valuePerPoint) + " " + self.freqDataDims[generic][x].unit,
                      lambda x = dim, value = value: self.void(),
                      allDims,
                      {}],
                     ["",
                      Label,
                      lambda x = dim: "Reference %s (%s)" % (self.expDimRefs[x].serial, self.expDimRefs[x].isotopeCodes[0]),
                      lambda x = dim, value = value: self.emptyString(),
                      allDims,
                      {}],
                     ["Reference point",  # If changed, have to propagate to changeFreqNumPoints
                      Entry,
                      lambda x = dim: self.freqDataDims[generic][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]).refPoint,
                      lambda x = dim, value = value: setattr(self.freqDataDims[generic][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]),'refPoint',returnFloat(value)),
                      allDims,
                      {}],
                     ["Reference value",
                      Entry,
                      lambda x = dim: self.freqDataDims[generic][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]).refValue,
                      lambda x = dim, value = value: setattr(self.freqDataDims[generic][x].findFirstDataDimRef(expDimRef = self.expDimRefs[x]),'refValue',returnFloat(value)),
                      allDims,
                      {}]
                    ]
     
    # Setup table data

    parList = self.parsList[parsListKey]
    self.objects[parsListKey] = []

    for dim in range(0,self.experiment.numDim):
      
      self.objects[parsListKey].append(specDims())
            
      for row in range(0,len(parList)):
        rowName = 'row%d' % row
        # Get values from data model
        self.objects[parsListKey][-1].addRow(rowName,parList[row][2](dim))


  def getDeleteImaginary(self,format,x):
  
    deleteImaginary = self.freqDataDims[format][x].findFirstApplicationData(application = self.applicationName, keyword = 'deleteImaginary')
    
    if deleteImaginary:
      value = deleteImaginary.value
    else:
      value = True
      self.freqDataDims[format][x].addApplicationData(Implementation.AppDataBoolean(application = self.applicationName,keyword = 'deleteImaginary',value = value))
    
    return value

  def setupObjectsTable(self,master,parsListKey):
  
    self.createParsLists(parsListKey)
    
    parList = self.parsList[parsListKey]
    
    attributes = []
    attributeHeadings = {}
    widgetClasses = {}
    constructorKw = {}
    
    for row in range(0,len(parList)):
      rowName = 'row%d' % row
      attributes.append(rowName)
      attributeHeadings[rowName] = parList[row][0]
      widgetClasses[rowName] = parList[row][1]
      constructorKw[rowName] = parList[row][5]
        
    #
    # callback for updating data in data model
    #
    callback = self.parameterUpdate

    return ObjectsTable(master,objects = self.objects[parsListKey],
                               widgetClasses=widgetClasses,
                               attributes=attributes,
                               attributeHeadings=attributeHeadings,
                               constructorKw = constructorKw,
                               callback = callback,
                   	           help_url = self.help_urls[parsListKey])


  """

  General functions here

  """
  
  def void(self):
    pass

  def emptyString(self):
    return ""
      
  def selectFile(self):
  
    fileName = getDataSourceFileName(self.experiment.rawData)
    
    format = string.lower(self.experiment.rawData.dataStore.fileType)
    popup = FormatFilePopup(self, file = fileName, component = 'FID', format = format)
    
    if popup.file:
      (dataUrlPath, localPath) = os.path.split(popup.file)
      setDataSourceDataStore(self.experiment.rawData.dataStore,dataUrlPath, localPath)
      self.selectFileButton.config(text = popup.file)
      popup.destroy()

  """

  Specialized functions

  """
  
  def procFileWrite(self,format):
  
    popup = self.popups.get(format)

    if (popup):
      popup.open()
      popup.update()
      
    else:    
      self.popups[format] = ProcFilePopup(self.parent,self,format)
  
  def toggleTable(self,master,row,parName):

    buttonState = getattr(self,parName + 'Button').toggled
    
    # Create objectsTable if doesn't exist
    
    if not hasattr(self,parName + 'Table'):
      
      # Setup datadim and dict
      self.setupDataSource(parName)
      
      # Create objectTable
      setattr(self,parName + 'Table', self.setupObjectsTable(master,parName))

    table = getattr(self,parName + 'Table')
    
    if buttonState:
      table.grid(row=row, column=0, columnspan = self.columnsTotal, sticky=Tkinter.EW)
      
      # Also grid all widgets associated with this toggle
      
      if hasattr(self,parName + 'Widgets'):
      
        widgetsInfo = getattr(self,parName + 'Widgets')
        
        for (widget,row,column,columnspan) in widgetsInfo:
          widget.grid(row=row, column = column, columnspan = columnspan, sticky=Tkinter.EW)

    else:
      table.grid_forget()
      
      # Also ungrid all widgets associated with this toggle
      
      if hasattr(self,parName + 'Widgets'):
        
        for (widget,row,column,columnspan) in getattr(self,parName + 'Widgets'):
          widget.grid_forget()
      
      
  def changeFreqNumPoints(self,format,dim,value):
 
    if self.freqDataDims[format][dim].numPointsOrig != returnInt(value):

      origValPerPoint = self.freqDataDims[format][dim].valuePerPoint
      
      if self.rawData:
        origSpecWidth = self.rawDataDims[dim].spectralWidth
      else:
        origSpecWidth = self.freqDataDims[format][dim].spectralWidthOrig
    
      newValPerPoint = origSpecWidth/returnFloat(value)
      
      # Set value per point in data model
      # TODO: currently also changing number of points... have to update this
      # when ranging introduced!

      setattr(self.freqDataDims[format][dim],'valuePerPoint',newValPerPoint)
      setattr(self.freqDataDims[format][dim],'numPointsOrig',returnInt(value))
      setattr(self.freqDataDims[format][dim],'numPoints',returnInt(value))
         
      #
      # Also reset refpoints (for all dataDimRefs) and valuePerPoint in window
      #
            
      parList = self.parsList[format]
      objectTable = getattr(self,format + 'Table')
       
      for row in range(0,len(parList)):
        rowName = 'row%d' % row
        
        rowLabel = parList[row][0]
        
        #
        # Get origRefPoint value from window...
        # 

        if rowLabel == "Reference point":
         
          origRefPoint = returnFloat(objectTable.widgets[rowName][dim].get())
          #print origRefPoint

          for dataDimRef in self.freqDataDims[format][dim].dataDimRefs:
            # Also reset refPoint... this assumes that the first point is 1 (and its value is 0)
            # TODO: have to make sure that refPoint is within bounds at some stage!!
            firstPoint = 1

            refPoint = ((origRefPoint - firstPoint) * origValPerPoint)  / newValPerPoint + firstPoint

            dataDimRef.refPoint = refPoint
        
            if dataDimRef.expDimRef == self.expDimRefs[dim]:
              objectRefPoint = refPoint
              dataDimRef.refPoint = objectRefPoint

          #
          # Set back in window/object
          #
         
          setattr(self.objects[format][dim],rowName,objectRefPoint)
          objectTable.widgets[rowName][dim].set(objectRefPoint)

        
        #
        # Also quickly set valuePerPoint (otherwise not updated on screen)
        # 

        elif rowLabel == "Value per point":
          setattr(self.objects[format][dim],rowName,newValPerPoint)
          objectTable.widgets[rowName][dim].set(newValPerPoint)
          
          

  def cycleExpDimRef(self,value):

    # Find dim in which button pressed
    
    rowName = None
    parsListKey = 'refPars'
    
    for row in range(0,len(self.parsList[parsListKey])):
      if self.parsList[parsListKey][row][0] == "Cycle reference":
        rowName = 'row%d' % row
        
    if not rowName:
      print "  Error in cycling reference frequency: header name wrong"
      return

    for dim in range(0,self.experiment.numDim):
      value = self.refParsTable.widgets[rowName][dim].get()

      if value == 1:
        break

    # First update everything in case anything modified
    self.updateFromWidget()

    # Look for next expDimRef (cycle back to beginning if end of list)
    
    expDim = self.experiment.sortedExpDims()[dim]
    newExpDimRef = expDim.findFirstExpDimRef(serial = (self.expDimRefs[dim].serial + 1))
    
    if newExpDimRef:
      self.expDimRefs[dim] = newExpDimRef
    
    else:
      self.expDimRefs[dim] = expDim.findFirstExpDimRef()

    # Get values from data model, set in window objects, and update widgets
    
    self.updateFromDataModel()

  def setAcqOrdMenu(self,master):
  
    numDim = self.experiment.numDim
  
    if self.acqOrd.has_key(numDim):
      
      storageOrder = self.experiment.rawData.storageOrder

      if self.acqOrd[self.experiment.numDim][1].count(storageOrder) > 0:
        acqOrderSelIndex = self.acqOrd[self.experiment.numDim][1].index(storageOrder)
      else:
        acqOrderSelIndex = 0
      
      return PulldownMenu(master, entries = self.acqOrd[numDim][0], selected_index = acqOrderSelIndex)
    
    else:
      return None

  def setSwapData(self,master):
    
    selected = 0
    
    if self.experiment.rawData:
      swapApplData = self.experiment.rawData.findFirstApplicationData(application = 'allProcessing', keyword = 'swap')
    
      if swapApplData:
        selected = swapApplData.value
    
    return CheckButton(master, selected = selected)

  def setSensEnhMenu(self,master):
  
    numDim = self.experiment.numDim
  
    dims = ['None']
  
    for i in range(1,numDim):
      dims.append(str(i + 1))
  
    # Setup selection
    
    selectedIndex = 0
    
    method = self.getSensEnhanceMethod()
        
    if self.experiment.rawData.storageDetails == method:
    
      for i in range(1,len(self.rawDataDims)):

        rawDataDim = self.rawDataDims[i]

        appl = rawDataDim.findFirstApplicationData(application = 'allProcessing', keyword = method, value = 'on')

        if appl:
          selectedIndex = dims.index(str(i+1))

    return PulldownMenu(master, entries = dims, selected_index = selectedIndex)

  def windowFunction(self,format,dim = None, getValue = 0):
   
    if getValue == 0:
      
      #
      # Find the row, column in which button pressed
      #
      
      rowName = None
    
      for row in range(0,len(self.parsList[format])):
        if self.parsList[format][row][0] == "Window function":
          rowName = 'row%d' % row
        
      if not rowName:
        print "  Error in window function buttons: header name wrong"
        return
        
      for dim in range(0,self.experiment.numDim):
        
        objectTable = getattr(self,format + 'Table')
     
        value = objectTable.widgets[rowName][dim].toggled

        if value == 1:
          objectButton = objectTable.widgets[rowName][dim]
          break
          
      # Reset button
      
      objectButton.toggled = 0
      objectButton.callback(True)
    
    popup = WindowFunctionSelect(self,format,dim,getValue)
    
    if getValue:
      return popup.getApplDataAndReturnString()
      #print 
      
    else:
      self.updateFromWidget()
 
 
  """

  Functions for updating objecttables & window

  """
      
  def updateFromWidget(self, value = None):
    
    for parsListKey in self.parsList.keys():
      if hasattr(self, parsListKey + 'Table'):      
        getattr(self, parsListKey + 'Table').update()

  def updateOwnWidget(self, headerName):
        
    for parsListKey in self.parsList.keys():
      if hasattr(self, parsListKey + 'Table'):      

        parList = self.parsList[parsListKey]
        objectTable = getattr(self,parsListKey + 'Table')
        
        for dim in range(0,self.experiment.numDim):
          for row in range(0,len(parList)):
            rowName = 'row%d' % row
            
            if parList[row][0] == headerName:
          
              widgetValue = objectTable.widgets[rowName][dim].get()
              objectValue = getattr(self.objects[parsListKey][dim],rowName)
              #print widgetValue, objectValue
              if widgetValue != objectValue:
                setattr(self.objects[parsListKey][dim],rowName,widgetValue)
                
    self.parameterUpdate()
    
  def updateFromDataModel(self):
    #print "updateFromDataModel"
    self.updateObjectsTableObjects()
    
    for parsListKey in self.parsList.keys():
      if hasattr(self, parsListKey + 'Table'):
        getattr(self, parsListKey + 'Table').reset()
    
  def parameterUpdate(self, *args):
    
    #
    # Set sensitivity enhanced stuff
    #
    
    if hasattr(self,'sensEnhMenu'):

      method = self.getSensEnhanceMethod()
      rawData = self.experiment.rawData

      if self.sensEnhMenu.getSelected() != 'None':

        selectedRawDataDim = self.rawDataDims[returnInt(self.sensEnhMenu.getSelected()) - 1]

        if rawData.storageDetails != method:
          rawData.isNormalStorage = False     
          rawData.storageDetails = method

      else:

        selectedRawDataDim = None

        if rawData.storageDetails == method:
          # TODO: this not entirely correct... could be other stuff going on
          rawData.isNormalStorage = True
          rawData.storageDetails = None

      #
      # Make sure only set for correct dimension...
      #

      for rawDataDim in self.rawDataDims:

        appl = rawDataDim.findFirstApplicationData(application = 'allProcessing', keyword = method, value = 'on')

        if rawDataDim == selectedRawDataDim:    
          if not appl:         
            appl = Implementation.AppDataString(application = 'allProcessing', keyword = method, value = 'on')
            rawDataDim.addApplicationData(appl)


        else:
          if appl:
            appl.delete()
          
    #
    # Set acquisition order in data model
    #
    
    if hasattr(self,'acqOrdMenu') and self.acqOrdMenu:
    
      acqOrder = self.acqOrdMenu.getSelected()
      
      acqOrderSelIndex = self.acqOrd[self.experiment.numDim][0].index(acqOrder)
      self.experiment.rawData.storageOrder = self.acqOrd[self.experiment.numDim][1][acqOrderSelIndex]
    
    #
    # Set 'swap' in data model
    #
    
    if hasattr(self,'swapData'):

      swapState = str(self.swapData.get())

      swapApplData = self.experiment.rawData.findFirstApplicationData(application = 'allProcessing', keyword = 'swap')

      if swapApplData:
        swapApplData.value = swapState
      else:
        swapApplData = createAppData(self.experiment.rawData,'allProcessing','swap',swapState)

    #
    # Set table attributes in data model (using lambda's in parsList!)
    #
    
    for dim in range(0,self.experiment.numDim):
    
      for parsListKey in self.parsList.keys():
      
        parList = self.parsList[parsListKey]
      
        for row in range(0,len(parList)):
          rowName = 'row%d' % row
          # Set values in data model
          print dir(self.objects[parsListKey][dim])
          value = getattr(self.objects[parsListKey][dim],rowName)
          parList[row][3](dim,value)


    #
    # Special setting for Azara blockSize
    #
    
    if self.rawData and self.parsList.has_key('azara'):

      blockSizes = getAzaraBlockSizes(self.rawDataDims)
      format = 'azara'
      
      for i in range(0,len(self.rawDataDims)):
        expDim = self.rawDataDims[i].expDim
        
        #
        # In case different order or lower dim extraction
        #

        for j in range(0,len(self.freqDataDims[format])):
          if self.freqDataDims[format][j].expDim == expDim:
            freqDataDim = self.freqDataDims[format][j]
          
        setattr(freqDataDim,'blockSize',blockSizes[i])
  
  
    # Get values from data model and set in window objects
    #print "parameterUpdate"
    self.updateObjectsTableObjects()


  def updateObjectsTableObjects(self):
    
    for parsListKey in self.parsList.keys():

      parList = self.parsList[parsListKey]
      
      for dim in range(0,self.experiment.numDim):

        for row in range(0,len(parList)):
          rowName = 'row%d' % row
          setattr(self.objects[parsListKey][dim],rowName,parList[row][2](dim))
          
          #if parList[row][0] == "Reference point":
          #  print parsListKey, dim, parList[row][2](dim)
    
  def getSensEnhanceMethod(self):
  
    if self.experiment.rawData.dataStore.fileType == 'Varian':
      method = 'ranceKay'
    else:
      method = 'echoAntiEcho'
      
    return method
 
  def getDataDimRefs(self,format):
  
    dataDimRefs = []
    
    freqDataDimIndex = 0
    
    for expDimRef in self.expDimRefs:
    
      dataDim = expDimRef.findFirstDataDimRef(dataDim = self.freqDataDims[format][freqDataDimIndex])
      
      if dataDim:
      
        dataDimRefs.append(dataDim)
        freqDataDimIndex += 1
        
    return dataDimRefs 
 
#
# Create class with all info for objects in table
#

class specDims:
  
  def addRow(self,rowName,value):
    
    setattr(self,rowName,value)
