
"""
======================COPYRIGHT/LICENSE START==========================

ProcFilePopup.py: Popup for processing files writing.

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
import Tkinter, os

from memops.general.Application import createAppData
from memops.universal.Io import joinPath

from memops.gui.Label import Label
from memops.gui.BasePopup import BasePopup
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.MessageReporter import showError, showInfo

from ccp.general.Io import getDataSourceFileName, setDataSourceDataStore

from ccpnmr.format.gui.FormatFilePopup import FormatFilePopup
from ccpnmr.format.general.Constants import allFormatsDict
from ccpnmr.format.general.Io import getHelpUrlDir

class ProcFilePopup(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ProcFile.html')
 
  def __init__(self, parent, acqProcPars, format):
  
    self.rawData = acqProcPars.experiment.rawData
    self.freqData = acqProcPars.freqDataDims[format][0].dataSource
    self.expDimRefs = acqProcPars.expDimRefs
    self.format = format
    self.updateAcqProcPars = acqProcPars.parameterUpdate
    self.updateAcqProcParsWidget = acqProcPars.updateFromWidget
    self.getDataDimRefs = acqProcPars.getDataDimRefs
    
    self.formatObjectDict = parent.formatObjectDict
    
    self.guiParent = parent
    self.project = self.rawData.root
    
    #
    # Dict for setup
    #
    # TODO: really need two datasources for nmrPipe... one fid, one processed!?!?
    #
    
    for writeFormatName in acqProcPars.writeFormats:
      setattr(self,writeFormatName,writeFormatName)
    
    value = 0
    self.fileSetup = {}
    self.writeSetup = {}
    
    #
    #  nmrPipe definitions
    #
    
    self.fileSetup[self.nmrPipe] = [
                           ['Conversion script file\n        (outputs to...)',
                            lambda: self.rawData.findFirstApplicationData(application = self.nmrPipe,keyword = 'conversion script'),
                            lambda x = value: self.setApplData(self.rawData,self.nmrPipe,'conversion script',x),
                            ['autoconv.com'],
	        'conversion script'],
                           ['    nmrPipe FIDs\n(inputs to...)',
                            lambda: self.rawData.findFirstApplicationData(application = self.nmrPipe,keyword = 'FID'),
                            lambda x = value: self.setApplData(self.rawData,self.nmrPipe,'FID',x),
                            ['auto.fid','auto%03d.fid'],
	        'FID'],
                           ['Processing script file\n       (outputs to...)',
                            lambda: self.freqData.findFirstApplicationData(application = self.nmrPipe,keyword = 'processing script'),
                            lambda x = value: self.setApplData(self.freqData,self.nmrPipe,'processing script',x),
                            ['autoproc.com'],
	        'processing script'],
                           ['Processed data',
                            lambda: self.getDataSourceFileName(self.freqData),
                            lambda x = value: self.doDataStoreFileName(self.freqData,x),
                            ['auto.ft','auto%03d.ft'],
	        'processed data']
                          ]
      
    # TODO: can add execute nmrPipe program buttons... see dataNavigator
      
    self.writeSetup[self.nmrPipe] = [
                                  ['Write conversion file',
                                   lambda : self.writeFile(self.nmrPipe,['conversion'])],
                                  ['Write processing file',
                                   lambda : self.writeFile(self.nmrPipe,['processing'])],
                                  ['Write all',
                                   lambda : self.writeFile(self.nmrPipe,['conversion','processing'])]
                                 ]

    #
    #  azara definitions
    #
    
    self.fileSetup[self.azara] = [
                           ['Processing parameters file\n          (inputs to...)',
                            lambda: self.rawData.findFirstApplicationData(application = self.azara,keyword = 'processing parameters'),
                            lambda x = value: self.setApplData(self.rawData,self.azara,'processing parameters',x),
                            ['autoproc.par'],
	        'processing parameters'],
                           ['Processing script file\n (outputs to both...)',
                            lambda: self.freqData.findFirstApplicationData(application = self.azara,keyword = 'processing script'),
                            lambda x = value: self.setApplData(self.freqData,self.azara,'processing script',x),
                            ['autoproc.scr'],
	        'processing script'],
                           ['Processed data',
                            lambda: self.getDataSourceFileName(self.freqData),
                            lambda x = value: self.doDataStoreFileName(self.freqData,x),
                            ['spectrum.ft'],
	        'spectrum'],
                           ['Processed data parameter file\n',
                            lambda: self.freqData.findFirstApplicationData(application = self.azara,keyword = 'spectrum parameters'),
                            lambda x = value: self.setApplData(self.freqData,self.azara,'spectrum parameters',x),
                            ['spectrum.par'],
	        'spectrum parameters']
                          ]
      
    # TODO: can add execute process program buttons... see dataNavigator
      
    self.writeSetup[self.azara] = [
                                  ['Write processing parameters file',
                                   lambda : self.writeFile(self.azara,['processingPars'])],
                                  ['Write processing script file',
                                   lambda : self.writeFile(self.azara,['processingScript'])],
                                  ['Write all',
                                   lambda : self.writeFile(self.azara,['processingPars','processingScript'])]
                                 ]
    
    BasePopup.__init__(self,parent = parent, title = "Project '%s': " % self.project.name + '%s processing write' % format, modal = False, transient=False)
 
  def body(self, master):
    
    self.buttons = []

    #
    # Popup window
    #

    row = 0
    label = Label(master, text= "Experiment '%s'" % self.rawData.experiment.name, fg = 'blue')
    label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.EW)

    #
    # File name buttons
    #

    for buttonIndex in range(0,len(self.fileSetup[self.format])):

      row = row + 1
      
      (textLabel,getter,setter,defaults,popupTitle) = self.fileSetup[self.format][buttonIndex]
      getFuncOrValue = getter()
      
      if not getFuncOrValue:
        
        if len(defaults) > 1 and len(self.expDimRefs) > 2:        
          fileName = defaults[1]
        else:
          fileName = defaults[0]

      else:
        if type(getFuncOrValue) == type(""):
          fileName = getFuncOrValue
        else:
          fileName = getFuncOrValue.value

      label = Label(master, text= textLabel)
      label.grid(row=row, column=0, sticky=Tkinter.E)
        
      self.buttons.append(Tkinter.Button(master, text = fileName, command = (lambda x = buttonIndex: self.selectFile(self.format,x))))
      self.buttons[-1].grid(row=row, column=1, sticky=Tkinter.W)
 
    #
    # Buttons to write conversion and/or processing files...
    #  

    for buttonIndex in range(0,len(self.writeSetup[self.format])):

      row = row + 1
      
      (text,command) = self.writeSetup[self.format][buttonIndex]
      
      button = Tkinter.Button(master, text = text, command = command)
      button.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.EW)
 

    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 3)
  
  def setApplData(self,dataSource,format,keyword,value):
  
    applData = dataSource.findFirstApplicationData(application = format, keyword = keyword)
    
    if applData:
      applData.value = value
    else:
      createAppData(dataSource,application = format, keyword = keyword, value = value)
  
  def selectFile(self,format,buttonIndex):
   
    fileName = self.buttons[buttonIndex].__getitem__('text')
    
    popupTitle = 'Select ' + self.format + ' ' + self.fileSetup[format][buttonIndex][4] + ' file'
    
    popup = FormatFilePopup(self, file = fileName, component = 'processing', title = popupTitle, format = self.format)
    
    # TODO: some testing on file name? Is it correct (ie has % code in for fid/ft bits?)
    
    if popup.file:
      self.fileSetup[format][buttonIndex][2](popup.file)
      self.buttons[buttonIndex].config(text = popup.file)
      popup.destroy()

  def doDataStoreFileName(self,masterClass,file):
  
    (dataUrlPath, localPath) = os.path.split(file)
    setDataSourceDataStore(masterClass,dataUrlPath,localPath)

  def getDataSourceFileName(self,dataStore):
      
    return getDataSourceFileName(dataStore)
 
  def writeFile(self,format,tasks):
    
    self.updateAcqProcParsWidget()
    self.updateAcqProcPars()
  
    #
    # Get or set up format classes
    #
    
    if self.formatObjectDict[self.format]:
    
      formatObject = self.formatObjectDict[self.format]
      
    else:
    
      formatLabel = allFormatsDict[self.format]      
      formatModule = __import__('ccpnmr.format.converters.%sFormat' % formatLabel,{},{},['%sFormat' % formatLabel])
      formatClass = getattr(formatModule,'%sFormat' % formatLabel)
      formatObject = formatClass(self.project,self.guiParent,verbose = 1)

      self.formatObjectDict[self.format] = formatObject

    # Bit cumbersome... is there a better way to define writeFile, ... ?
    
    dataDimRefs = self.getDataDimRefs(format)
    
    for task in tasks:
  
      parFile = None
      addKeywds = {'dataDimRefs': dataDimRefs}
      
      if format == self.nmrPipe:
   
        if task == 'conversion':
          function = 'writeProcPars'
          self.fileName = self.buttons[0].__getitem__('text')
          outputFile = self.buttons[1].__getitem__('text')
          
          addKeywds['outputFile'] = outputFile
    
        elif task == 'processing':    
          function = 'writeProcScript'
          self.fileName = self.buttons[2].__getitem__('text')
          inputFile = self.buttons[1].__getitem__('text')
          outputFile = self.buttons[3].__getitem__('text')
      
          addKeywds['outputFile'] = outputFile
          addKeywds['inputFile'] = inputFile

      if format == self.azara:
   
        if task == 'processingPars':
          function = 'writeProcPars'
          self.fileName = self.buttons[0].__getitem__('text')
    
        elif task == 'processingScript':
          function = 'writeProcScript'
          self.fileName = self.buttons[1].__getitem__('text')
          inputFile = self.buttons[0].__getitem__('text')
          outputFile = self.buttons[2].__getitem__('text')
          parFile = self.buttons[3].__getitem__('text')

          addKeywds['outputFile'] = outputFile
          addKeywds['inputFile'] = inputFile
          addKeywds['parFile'] = parFile

      # Write for each task in list
      
      exportFunc = getattr(formatObject,function)
      
      if exportFunc(self.fileName, **addKeywds):
        showInfo("Success","Successfully exported file:\n%s." % self.fileName,self)

      else:
        showError("Error","Error while exporting file:\n%s." % self.fileName,self)

  def apply(self):
    
    return True
