
"""
======================COPYRIGHT/LICENSE START==========================

ChainExportPopup.py: GUI window for setting chain names/codes to be used in exported files

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

from memops.universal.Io import joinPath

from memops.gui.Label import Label
from memops.gui.Util import createHelpButtonList
from memops.gui.MessageReporter import showError
from memops.gui.Entry import Entry

from memops.universal.Util import returnInt
from ccpnmr.format.general.Io import getHelpUrlDir

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class ChainExportPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ChainExport.html')
  
  def __init__(self, parent, ccpChainLabelDict, ccpChainSeqIDCodes, requireChainCode):

    # Constructor doesn't do much except call body
    # The parent is self.parent (parent of the popup)
   
    self.chain = None

    self.ccpChainLabelDict = ccpChainLabelDict
    self.ccpChainSeqIDCodes = ccpChainSeqIDCodes

    self.requireChainCode = requireChainCode    

    self.exportChainCode = {}
    self.exportFirstSeqCode = {}

    project = self.ccpChainLabelDict.values()[0].root

    # modal = true means that it won't continue unless this one returns value
    TemporaryBasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + 'Export chains', modal=False, transient=True)

  def body(self, master):
    
    #
    # Setup header 
    #

    row = 0

    label = Label(master, text= "Data model", fg = 'blue')
    label.grid(row=row, column=0, columnspan = 3, sticky=Tkinter.EW)

    label = Label(master, text= "Export")
    label.grid(row=row, column=3, columnspan = 2, sticky=Tkinter.EW)

    row += 1

    label = Label(master, text= "Sequence ID (code)", fg = 'blue')
    label.grid(row=row, column=1, columnspan = 2, sticky=Tkinter.EW)
    
    row += 1

    label = Label(master, text= "Ccp chain code", fg = 'blue')
    label.grid(row=row, column=0, sticky=Tkinter.W)

    label = Label(master, text= "start", fg = 'blue')
    label.grid(row=row, column=1, sticky=Tkinter.EW)

    label = Label(master, text= "end", fg = 'blue')
    label.grid(row=row, column=2, sticky=Tkinter.EW)
    
    column = 3
    columnspan = 2
    
    if self.requireChainCode:

      label = Label(master, text= self.requireChainCode)
      label.grid(row=row, column=column, sticky=Tkinter.W)
      
      column += 1
      columnspan -= 1

    label = Label(master, text= "Sequence code start")
    label.grid(row=row, column=column,columnspan = columnspan, sticky=Tkinter.W)

    #
    # Setup list of ccp chains labels
    # 
   
    ccpChainLabels = self.ccpChainLabelDict.keys()
    ccpChainLabels.sort()
                  
    #
    # Currently only linking start in ccp chain (ONLY ONE!) to start in chain read in from file
    # Assuming sequential order thereafter.
    #

    for ccpChainLabel in ccpChainLabels:

      row = row + 1
      
      chain = self.ccpChainLabelDict[ccpChainLabel]

      ccpCodeLow = "%d (%d)" % (self.ccpChainSeqIDCodes[ccpChainLabel][0][0],self.ccpChainSeqIDCodes[ccpChainLabel][1][0])
      ccpCodeHigh = "%d (%d)" % (self.ccpChainSeqIDCodes[ccpChainLabel][0][-1],self.ccpChainSeqIDCodes[ccpChainLabel][1][-1])

      label = Label(master, text= ccpChainLabel, fg = 'blue')
      label.grid(row=row, column=0, sticky=Tkinter.W)
      
      label = Label(master, text= ccpCodeLow, fg = 'blue')
      label.grid(row=row, column=1, sticky=Tkinter.W)

      label = Label(master, text= ccpCodeHigh, fg = 'blue')
      label.grid(row=row, column=2, sticky=Tkinter.W)

      column = 3
      columnspan = 2
      
      suggestedFirstSeqCode = str(self.ccpChainSeqIDCodes[ccpChainLabel][1][0])
    
      if self.requireChainCode:
      
        if self.requireChainCode == 'molecule name':
          chainCodeText = chain.molecule.name
          width = 20
        else:
          chainCodeText = chain.code
          width = 4

        self.exportChainCode[chain] = Entry(master, text = chainCodeText, width = width)
        self.exportChainCode[chain].grid(row=row, column=column, sticky=Tkinter.W)
      
        column += 1
        columnspan -= 1
        
        suggestedFirstSeqCode = '1'

      self.exportFirstSeqCode[chain] = Entry(master, text = suggestedFirstSeqCode, width = 4)
      self.exportFirstSeqCode[chain].grid(row=row, column=column,columnspan = columnspan, sticky=Tkinter.W)
      
    row = row + 1

    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 4)

  def apply(self):
  
    self.chainDict = {}
    
    exportSeqCodeRanges = {}
    
    for ccpChainLabel in self.ccpChainLabelDict.keys():
      
      chain = self.ccpChainLabelDict[ccpChainLabel]

      #
      # Get the first seq code for export (always necessary)
      #
        
      exportFirstSeqCode = self.exportFirstSeqCode[chain].get()
      
      #
      # Check if integer
      #
      
      if str(returnInt(exportFirstSeqCode)) == exportFirstSeqCode:
      
        exportFirstSeqCode = returnInt(exportFirstSeqCode)
        
      else:
        showError("Error","Sequence codes have to be integers!",self)
        return False
        
      #
      # If single chain export or chain already used, make sure doesn't fall within range other seqCodes
      #
      
      
      if self.requireChainCode:

        exportChainCode = self.exportChainCode[chain].get()
        
      else:
      
        exportChainCode = None

      #
      # Set ranges
      #
      
      seqLength = len(self.ccpChainSeqIDCodes[ccpChainLabel][0])

      if not exportSeqCodeRanges.has_key(exportChainCode):
        exportSeqCodeRanges[exportChainCode] = [[exportFirstSeqCode,exportFirstSeqCode + seqLength - 1]]
        checkChainCodes = 0
          
      else:
        
        #
        # Check chain seqcodes
        #
        
        for (lowSeqID,highSeqID) in exportSeqCodeRanges[exportChainCode]:
        
          if highSeqID >= exportFirstSeqCode >= lowSeqID:
          
            #
            # Illegal! For same chain anyway...
            #
            
            showError("Error","Overlapping sequence code numbers within the same export chain - try changing them.",self)
            return False
            
        exportSeqCodeRanges[exportChainCode].append([exportFirstSeqCode,exportFirstSeqCode + seqLength - 1])
        
      self.chainDict[chain] = [exportChainCode,exportFirstSeqCode]
        
    return True
