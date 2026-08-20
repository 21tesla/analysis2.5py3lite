"""
======================COPYRIGHT/LICENSE START==========================

UpdatePopup.py: Part of the CcpNmr Update program

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
import os
import sys
import Tkinter

from memops.gui.BasePopup       import BasePopup
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.Label           import Label
from memops.gui.Entry           import Entry
from memops.gui.ButtonList      import UtilityButtonList
from memops.gui.MessageReporter import showWarning

from ccpnmr.update.UpdateAgent import UpdateAgent, UPDATE_SERVER_LOCATION, UPDATE_DIRECTORY, UPDATE_DATABASE_FILE

class UpdatePopup(BasePopup, UpdateAgent):
  """
  **Install Patches and Updates to CCPN Software**
  
  This popup window allows the user to view and automatically install the
  latest updates and patches to CcpNmr software. This is the normal mechanism
  of receiving incremental changes and fixes to bug reports between numbered
  CCPN releases. It should be noted that this system generally only provides
  fixes for the most recent version of the CcpNmr Software. If your version is
  older (but not as old as version 1.0) or if there has been a recent new
  release then the updating mechanism will automatically upgrade the whole CCPN
  distribution to the newest one. It should be noted that this mechanism will
  not automatically update some of the really old versions, e.g. from version
  1.0 to 2.0. Updates from 2.0 onwards ought to work seamlessly though.

  In normal operation the user clicks [Select All] then [Install Selected],
  whereupon the software should be restarted (not forgetting to save any data).
  Naturally, viewing and receiving updates requires that the user's computer
  has a working connection to the Internet. Although the user may toggle
  individual file updates on or off by double-clicking in the "Install?" column,
  it is generally not a good idea to choose only a subset of file updates unless
  there is a specific reason to do so. Choosing only some file updates risks
  getting problems where two or more files rely upon the most recent versions of
  one another.

  **Caveats & Tips**

  Using [Refresh List] will cause the system to re-query the CCPN update server
  if the user is awaiting an imminent fix.
  
  When an updated file is installed the old version is not removed; it is renamed
  with the "_old" suffix. 

  If for any reason this popup window cannot be opened the user may update the
  software using the "updateAll" script that sits alongside the main "analysis"
  program executable (usually in the $CCPN_HOME/bin/ directory).

  """

  def __init__(self, parent, serverLocation=UPDATE_SERVER_LOCATION,
               serverDirectory=UPDATE_DIRECTORY, dataFile=UPDATE_DATABASE_FILE,
               exitOnClose=False, helpUrl=None):

    self.exitOnClose = exitOnClose
    self.helpUrl = helpUrl
    UpdateAgent.__init__(self, serverLocation, serverDirectory, dataFile, isStandAlone=exitOnClose)
    self.fileUpdate  = None
 
    if exitOnClose:
      # quitFunc = sys.exit
      # NOTE:ED - small function to help Windows10 close properly
      quitFunc = lambda k=0: self.destroy(); os._exit(0)
    else:
      quitFunc = None

    BasePopup.__init__(self, parent=parent, title='CcpNmr Software Updates', quitFunc=quitFunc)

  def body(self, guiParent):
  
    self.geometry('600x300')
  
    guiParent.grid_columnconfigure(1, weight=1)
  

    url = ''
    if self.server:
      url = self.server.url
    label = Label(guiParent, text='Server location: %s' % url)
    label.grid(row=0, column=0, sticky='w', columnspan=2)

    label = Label(guiParent, text='Installation root: %s%s' % (self.installRoot,os.sep))
    label.grid(row=1, column=0, sticky='w', columnspan=2)
 
    editWidgets      = [ None] * 5
    editGetCallbacks = [ None, self.toggleSelected, None, None, None]
    editSetCallbacks = [ None] * 5

    guiParent.grid_rowconfigure(2, weight=1)
    headingList = ['File','Install?','Date','Relative Path','Priority','Comments']
    tipTexts = ['The name of the file for which an update is available',
                'Whether the file will be replaced when the user commits the upgrade',
                'The data that the new file version was posted on the CCPN server',
                'The path that the file will be installed at, relative to the CCPN installation',
                'The priority number of the installation; commonly not configured away from "1"',
                'Any comments added by the CCPN developers']
    self.scrolledMatrix = ScrolledMatrix(guiParent, headingList=headingList,
                                         highlightType = 0, tipTexts=tipTexts,
                                         editWidgets=editWidgets, callback=self.selectCell,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks)
                                         
    self.scrolledMatrix.grid(row=2, column=0, columnspan=2, sticky='nsew')

    if self.exitOnClose:
      #txt = 'Quit'
      # cmd = sys.exit
      # NOTE:ED - small function to help Windows10 close properly
      cmd = lambda k=0: self.destroy(); os._exit(0)
    else:
      #txt = 'Close'
      cmd = self.close

    texts = ['Refresh List','Select All','Install Selected']
    commands = [self.updateFiles, self.selectAll, self.install]
    tipTexts = ['Refresh the list of file updates by querying the CCPN update server, assuming an Internet connection is present',
                'Select all of the file updates for installation (recommended)',
                'Install newest versions of the selected files']
    self.buttonList = UtilityButtonList(guiParent, texts=texts, doClone=False, helpUrl=self.helpUrl,
                                        commands=commands, closeCmd=cmd, tipTexts=tipTexts)
    self.buttonList.grid(row=3, column=0, columnspan=2, sticky='ew')
    
    if self.server:
      for fileUpdate in self.server.fileUpdates:
        fileUpdate.isSelected = False

    #self.update()
    # need self.updateFiles, not just self.update
    # because otherwise the 2.0.5 to 2.0.6 upgrades do not work
    # (self.server not set on first pass so need call to updateFiles here)
    self.updateFiles()

  def toggleSelected(self, fileUpdate):
  
    boolean = fileUpdate.isSelected
    
    fileUpdate.isSelected = not boolean

    self.update()

  def install(self):
  
    if self.server:
      self.installUpdates()    
      self.updateFiles()

  def selectAll(self):
  
    if self.server:
      for fileUpdate in self.scrolledMatrix.objectList:
        fileUpdate.isSelected = True
        
    self.update()

  def selectCell(self, object, row, col):
  
    self.fileUpdate = object

  def updateButtons(self):
  
  
    buttons = self.buttonList.buttons
    
    if self.server:
      selected = False
      if self.server and self.server.fileUpdates:
        for fileUpdate in self.scrolledMatrix.objectList:
          if fileUpdate.isSelected:
            selected = True
            break
 
        buttons[1].enable()
      else:
        buttons[1].disable()

      if selected:
        buttons[2].enable()
        buttons[2].config(bg='#A0FFA0')
      else:
        buttons[2].disable()
        buttons[2].config(bg=buttons[1].cget('bg'))
 
    else:
      for button in buttons[:-1]:
        button.disable()
 
  def updateFiles(self):
    
    if self.server:
      if not self.server.fileUpdates:
        self.server.getFileUpdates()
        
      for fileUpdate in self.server.fileUpdates:
        fileUpdate.isSelected = False
        
    self.update()

  def update(self):
    
    self.updateButtons()
    
    self.fileUpdate = None
    textMatrix  = []
    objectList  = []
    colorMatrix = []
    
    if self.server:
      
      for fileUpdate in self.server.fileUpdates:
        
		
        if fileUpdate.getIsUpToDate():
          continue
          
        color = [None,'#A04040',None,None,None]
        yesNo = 'No'
        if fileUpdate.isSelected:
          color = [None,'#A0FFA0',None,None,None]
          yesNo = 'Yes'
        
        datum = []
        datum.append(fileUpdate.fileName)
        datum.append(yesNo)
        datum.append(fileUpdate.date)
        datum.append(fileUpdate.filePath+os.sep)
        datum.append(fileUpdate.priority)
        datum.append(fileUpdate.details)
        
        colorMatrix.append(color)
        textMatrix.append(datum)
        objectList.append(fileUpdate)
    
    
    self.scrolledMatrix.update(textMatrix=textMatrix, objectList=objectList, colorMatrix=colorMatrix)
    

if __name__ == '__main__':
  import os

  root = Tkinter.Tk()
  root.withdraw()
  top = UpdatePopup(root, exitOnClose=True)
  root.mainloop()
