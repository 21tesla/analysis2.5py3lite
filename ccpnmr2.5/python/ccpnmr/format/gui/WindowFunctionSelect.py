
"""
======================COPYRIGHT/LICENSE START==========================

WindowFunctionSelect.py: Selection for processing program window functions

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
import Tkinter
import string

from memops.universal.Util import returnFloat
from memops.universal.Io import joinPath

from memops.general.Application import createAppData

from memops.gui.BasePopup import BasePopup
from memops.gui.Label import Label
from memops.gui.Entry import Entry
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu

from ccpnmr.format.general.Io import getHelpUrlDir

class WindowFunctionSelect(BasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'WindowFunctionSelect.html')

  def __init__(self, parent, format, dim, getValue):

    #
    # General setup
    #
   
    self.maxArgs = 4
    self.nmrPipe = 'nmrPipe'
    self.azara = 'azara'
    self.format = format
    
    self.freqDataDim = parent.freqDataDims[format][dim]
    self.expDimRef =  parent.expDimRefs[dim]
    
    #
    # Window function settings (different per program)
    #
    
    self.windowFunctions = {}
    
    specWidth = self.freqDataDim.spectralWidthOrig
    
    # TODO: these for azara still have to change a bit...
    self.windowFunctions[self.azara] = [
                            ['sine','sinebell',[['Angle (deg)',0,'deg']]],
                            ['sine^2','sinebell2',[['Angle (deg)',0,'deg']]],
                            ['exponential (fixed)','decay',[['End value',0.5,'end']]],
	        ['exponential','decay_sw',[['Line broadening (Hz)',10,'lb'],['Spectral width (Hz)',specWidth,'sw']]],
	        ['gaussian (fixed)','gaussian',[['Total points fraction for max',0.25,'frac'],['End value',1,'end']]],
	        ['gaussian','gaussian_sw',[['Line broadening (Hz)',10,'lb'],['Sharpening factor',0.7,'shrp'],['Spectral width (Hz)',specWidth,'sw']]],
                            ['inverse cosine','inv_cosine',[['Frequency (Hz)',0,'freq'],['Spectral width (Hz)',specWidth,'sw']]]
	       ]

    self.windowFunctions[self.nmrPipe] = [
                            ['sine','SP',[['Start sine at (pi)',0,'off'],['End sine at (pi)',1,'end'],['Sine to power',1,'pow']]],
	        ['exponential','EM',[['Line broadening (Hz)',10,'lb']]],
	        ['lorentz-gauss','GM',[['Inverse exp width (Hz)',0,'g1'],['Gaussian width (Hz)',0,'g2'],['Center (0-1)',0,'g3']]],
	        ['gaussian','GMB',[['Exponential term',0,'lb'],['Gaussian term',0,'gb']]],
	        ['trapezoid','TM',[['Ramp length left',0,'t1'],['Ramp length right',0,'t2']]],
	        ['triangle','TRI',[['Point loc of apex',0,'loc'],['Left edge height',0,'lHi'],['Right edge height',0,'rHi']]]
	       ]
 
    #
    # Make list of names for pulldown menu
    #
    
    self.windowFuncList = []

    for (winFunc,winFuncAppl,winFuncArgs) in self.windowFunctions[format]:
      self.windowFuncList.append(winFunc)

    #
    # Get the right window function from appldata (if there)
    #
    
    self.windowFunctionKeyword = 'windowFunction'
    self.winFuncIndex = self.getWindowFunctionIndex()
    
    #
    # Set defaults depending on nucleus type if no window function set
    #
    
    if self.winFuncIndex == -1:

      if format == self.nmrPipe:
        self.winFuncIndex = 0
        if self.expDimRef.isotopeCodes[0] == '1H':
          # use SP '-off 0.3333 -pow 2.0'
          self.windowFunctions[self.nmrPipe][self.winFuncIndex][2][0][1] = 0.3333
          self.windowFunctions[self.nmrPipe][self.winFuncIndex][2][2][1] = 2.0
        else:
          # for non-H, use SP '-off 0.5 -pow 2.0'
          self.windowFunctions[self.nmrPipe][self.winFuncIndex][2][0][1] = 0.5
          self.windowFunctions[self.nmrPipe][self.winFuncIndex][2][2][1] = 2.0

      elif format == self.azara:
        self.winFuncIndex = 1
        if self.expDimRef.isotopeCodes[0] == '1H':
          # use sinebell2 60
          self.windowFunctions[self.azara][self.winFuncIndex][2][0][1] = 60
        else:
          # for non-H, use sinebell2 90
          self.windowFunctions[self.azara][self.winFuncIndex][2][0][1] = 90
      
    
    # Get actual values for arguments from appldata if they exists!
    
    self.updateWindowFunctions()
    
    # Also using this to return the current window function 'string value'
    if not getValue:
      
      # modal = true means that it won't continue unless this one returns value
      BasePopup.__init__(self, parent=parent, title='Select window function', modal=False, transient=True, )


  def body(self, master):
  
    # Master is the owner widget (not self.parent) - parent of the widget here
    
    #
    # Popup window
    #
    
    self.entries = []
    self.labels = []
    value = 0
    
    row = 0
    label = Label(master, text= "Window function:")
    label.grid(row=row, column=0, sticky=Tkinter.EW)
    
    self.menuRow = row

    self.menu = PulldownMenu(master, entries = self.windowFuncList, selected_index = self.winFuncIndex)
    self.menu.grid(row=row, column=1, sticky=Tkinter.W, ipadx = 20)

    self.menulabel = Label(master, text = "(" + self.windowFunctions[self.format][self.winFuncIndex][1] + ")")
    self.menulabel.grid(row=row, column=2, sticky=Tkinter.E)

    #
    # Setup arguments - make room for four
    #
        
    for argIndex in range(0,self.maxArgs):
      
      row = row + 1
      
      try:
        argRef = self.windowFunctions[self.format][self.winFuncIndex][2][argIndex]
      
        argName = argRef[0]
        argValue = argRef[1]
        
      except:
        argName = None
        argValue = None
                      
      self.labels.append(Label(master, text= "%s" % argName))
      self.labels[-1].grid(row=row, column=0, sticky=Tkinter.EW)

      self.entries.append(Entry(master,text = str(argValue)))
      self.entries[-1].grid(row=row, column=1, columnspan = 2, sticky=Tkinter.W)
      
      if not argName:
        self.labels[-1].grid_forget()
        self.entries[-1].grid_forget()
        
    #
    # Add callback to menu after setup
    #
   
    self.menu.callback = self.updateWindow

    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 3)
    
  def getApplDataAndReturnString(self):
  
    keyword = self.windowFunctionKeyword
    returnString = self.windowFunctions[self.format][self.winFuncIndex][0]
    winFuncName = self.windowFunctions[self.format][self.winFuncIndex][1]

    for (argName, argValue, argShort) in self.windowFunctions[self.format][self.winFuncIndex][2]:

      applData = self.freqDataDim.findFirstApplicationData(application = self.format, keyword = keyword + winFuncName + argShort)
      returnString += ' ' + argShort
    
      if not applData:
      
        returnString += ' ' + str(argValue)
                
      else:
        
        returnString += ' ' +  str(applData.value)

    return returnString

  def updateWindowFunctions(self):
  
    keyword = self.windowFunctionKeyword
    winFuncName = self.windowFunctions[self.format][self.winFuncIndex][1]
          
    #
    # Update windowfunction parameters
    #
    
    for i in range(0,len(self.windowFunctions[self.format][self.winFuncIndex][2])):
    
      (argName, argValue, argShort) = self.windowFunctions[self.format][self.winFuncIndex][2][i]
      
      applData = self.freqDataDim.findFirstApplicationData(application = self.format, keyword = keyword + winFuncName + argShort)
    
      if applData:
      
        self.windowFunctions[self.format][self.winFuncIndex][2][i][1] = str(applData.value)    

  def setApplData(self):
  
    keyword = self.windowFunctionKeyword
    winFuncName = self.windowFunctions[self.format][self.winFuncIndex][1]
    
    #
    # Set windowfunction
    #
    
    applData = self.freqDataDim.findFirstApplicationData(application = self.format, keyword = keyword)
    
    if not applData:
      
      applData = createAppData(self.freqDataDim,application = self.format, keyword = keyword, value = winFuncName)
    
    else:
      
      applData.value = winFuncName
      
    #
    # Set windowfunction parameters
    #
    
    for (argName, argValue, argShort) in self.windowFunctions[self.format][self.winFuncIndex][2]:
        
      applData = self.freqDataDim.findFirstApplicationData(application = self.format, keyword = keyword + winFuncName + argShort)
    
      if not applData:
      
        applData = createAppData(self.freqDataDim,application = self.format, keyword = keyword + winFuncName + argShort, value = str(argValue))
        
      else:
        
        applData.value = str(argValue)
    
    
  def getWindowFunctionIndex(self):

    keyword = self.windowFunctionKeyword

    applData = self.freqDataDim.findFirstApplicationData(application = self.format, keyword = keyword)
    
    if applData:
    
      origWinFuncAppl = applData.value
    
      for i in range(0,len(self.windowFunctions[self.format])):
      
        (winFunc,winFuncAppl,winFuncArgs) = self.windowFunctions[self.format][i]
        
        if winFuncAppl == origWinFuncAppl:
          return i
        
    return -1

  def updateWindow(self,*args):
  
    #
    # First do a value update of the current stuff?
    #

    #
    # Reset arguments
    #
    
    row = self.menuRow
    
    self.winFuncIndex = self.menu.getSelectedIndex()
    
    self.menulabel.set("(" + self.windowFunctions[self.format][self.winFuncIndex][1] + ")")

    for argIndex in range(0,self.maxArgs):
      
      row = row + 1
      
      try:
        argRef = self.windowFunctions[self.format][self.winFuncIndex][2][argIndex]
      
        argName = argRef[0]
        argValue = argRef[1]
                
      except:
        argName = None
        argValue = None
        
      self.labels[argIndex].set(argName)
      self.labels[argIndex].grid(row=row, column=0, sticky=Tkinter.W)

      self.entries[argIndex].set(str(argValue))
      self.entries[argIndex].grid(row=row, column=1, columnspan = 2, sticky=Tkinter.W)
      
      if not argName:
        self.entries[argIndex].grid_forget()
        self.labels[argIndex].grid_forget()

  def setWindowFunctionsValues(self):
  
    self.winFuncIndex
    
    for argIndex in range(0,len(self.windowFunctions[self.format][self.winFuncIndex][2])):
    
      argValue = returnFloat(self.entries[argIndex].get())
      
      self.windowFunctions[self.format][self.winFuncIndex][2][argIndex][1] = argValue

  def apply(self):
  
    self.winFuncIndex = self.menu.getSelectedIndex()
    
    self.setWindowFunctionsValues()
    
    self.setApplData()

    return True


if __name__ == '__main__':
  
  
  root = Tkinter.Tk()
  
  popup = WindowFunctionSelect(root)

  if popup.chemAtoms:
    print "chemAtoms %s selected" % (popup.chemAtoms)
  
  else:
    print "Nothing selected - aborting script"
    
