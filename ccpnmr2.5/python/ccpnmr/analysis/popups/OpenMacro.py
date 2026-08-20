
"""
======================COPYRIGHT/LICENSE START==========================

OpenMacro.py: Part of the CcpNmr Analysis program

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
import re
import os.path

from memops.universal.Io import splitPath

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.FileSelect import FileSelect, FileType
from memops.gui.Label import Label
from memops.gui.Entry import Entry
from memops.gui.ScrolledMatrix import ScrolledMatrix


from ccpnmr.analysis.popups.BasePopup import BasePopup

class OpenMacroPopup(BasePopup):
  """
  **Locate Python Macro Scripts On Disk**
  
  This popup is used to locate a Python file and function definition on disk
  so that it can be loaded into Analysis as a "Macro". The top table
  is a file browser that the user can navigate to the location of the Python
  file (typically sending in ".py"). When a Python file is selected in the
  upper table the system will look through the contents of the file to find
  any macro functions that are defined within. Note that a macro function 
  is defined by the normal Python "def" keyword but must take "argSever"
  as its first input argument. 
  
  Any macro functions in the selected file will be displayed in the
  lower table. Clicking on the name of a function in this table
  selects it for opening. The user may edit the name of the macro for 
  display purposes within Analysis. Finally clicking [Load Macro]
  will add it to the CCPN project, stored in a user-specific way
  inside the current Analysis profile. Once loaded into the
  project the new macro will appear in the main macro table
  (the "Macros" tab of `User Options`_)
  
  .. _`User Options`: EditProfilesPopup.html
    
  """


  def __init__(self, parent, **kw):

    # TBD: properly
    transient = True
    BasePopup.__init__(self, parent=parent, title='Open macro', transient=transient, **kw)

  def body(self, guiParent):

    guiParent.grid_columnconfigure(1,weight=1)
    row = 0

    file_types = [  FileType('Python', ['*.py']), FileType('All', ['*']) ]
    self.file_select = FileSelect(guiParent, file_types=file_types,
                                  single_callback=self.chooseFile,
                                  double_callback=self.chooseFile)
    self.file_select.grid(row=row, column=0, columnspan=2, sticky='nsew')
    
    row = row + 1
    headingList=('Function',)
    self.scrolledMatrix = ScrolledMatrix(guiParent, initialRows=4,
                            headingList=headingList, callback=self.selectFunction)
    self.scrolledMatrix.grid(row=row, column=0, columnspan=2, sticky='nsew')
    guiParent.grid_rowconfigure(row,weight=1)

    row = row + 1
    self.moduleLabel1 = Label(guiParent, text='Module: ')
    self.moduleLabel1.grid(row=row, column=0, sticky='nw')
    self.moduleLabel2 = Label(guiParent, text=' ')
    self.moduleLabel2.grid(row=row, column=1, sticky='nw')

    row = row + 1
    self.functionLabel1 = Label(guiParent, text='Function: ')
    self.functionLabel1.grid(row=row, column=0, sticky='nw')
    self.functionLabel2 = Label(guiParent, text=' ')
    self.functionLabel2.grid(row=row, column=1, sticky='nw')

    row = row + 1
    self.nameLabel = Label(guiParent, text='Name: ')
    self.nameLabel.grid(row=row, column=0, sticky='nw')
    self.nameEntry = Entry(guiParent, text=' ', width=40)
    self.nameEntry.grid(row=row, column=1, sticky='nw')

    row = row + 1
    texts = [ 'Load Macro' ]
    commands = [ self.loadMacro ]
    buttons = UtilityButtonList(guiParent, texts=texts,
                                commands=commands, helpUrl=self.help_url)
    buttons.grid(row=row, column=0, columnspan=2, sticky='ew')
    
    self.loadButton = buttons.buttons[0]
    self.loadButton.disable()
    
    self.path     = None
    self.module   = None
    self.function = None

  def selectFunction(self, function, row, col):
  
    if function:
      self.loadButton.enable()
      self.function = function
      self.functionLabel2.set(function)
      self.nameEntry.set(function)
    
  def updateFunctions(self, file):
  
    self.loadButton.disable()
    functions = []
    fileHandle = open(file, 'r')
    line = fileHandle.readline()
    textMatrix = []
    while line :
      match = re.match('\s*def\s+(\w+)\s*\(.*argServer',line)
      if match:
        function = match.group(1)
        functions.append(function)
        textMatrix.append( [function,] )
      line = fileHandle.readline()

    self.scrolledMatrix.update(objectList=functions,textMatrix=textMatrix)

  def chooseFile(self, file):
    
    fullPath =self.file_select.getFile()
    if not os.path.isdir(fullPath) and file:
      (self.path, module) = splitPath(fullPath)
      module = re.sub('\..+?$','',module)
      self.module = module
      self.moduleLabel2.set(module)
      self.updateFunctions(file)

  def loadMacro(self):
    
    if self.module and self.function and self.path:
      name = self.nameEntry.get()
      if not name:
        name = self.function
        
      analysisProfile = self.analysisProfile
      
      m1 = analysisProfile.newMacro(name=name,path=self.path,
                                    function=self.function,
                                    module=self.module,ordering=1)
 
    self.close()


