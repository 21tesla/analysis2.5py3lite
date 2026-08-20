
"""
======================COPYRIGHT/LICENSE START==========================

ChainLinkPopup.py: GUI window for linking chains from original format file to chains in data model

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

from memops.universal.Io import joinPath

from memops.universal.Util import returnStrings
from memops.universal.Util import returnInt

from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledListbox import ScrolledListbox
from memops.gui.MessageReporter import showError

from ccp.format.general.Util import getSeqAndInsertCode

from ccpnmr.format.general.Io import getHelpUrlDir

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class ChainLinkPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ChainLink.html')

  def __init__(self, parent, ccpChainLabelDict, ccpChainSeqIdCodes, formatChains, formatChainList, formatChainDict, defaultFormatChain):

    # Constructor doesn't do much except call body
    # The parent is self.parent (parent of the popup)
   
    self.chain = None

    self.formatChains = formatChains
    self.ccpChainLabelDict = ccpChainLabelDict
    self.ccpChainSeqIdCodes = ccpChainSeqIdCodes
    self.formatChainList = formatChainList
    self.formatChainDict = formatChainDict
    self.localFormatChainDict = {}
    
    self.defaultFormatChain = defaultFormatChain
    
    self.ccpChainSeqIdDict = {}
    self.ccpChainSeqCodeDict = {}
    
    project = self.ccpChainLabelDict.values()[0].root

    # modal = true means that it won't continue unless this one returns value
    TemporaryBasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name +'Link chains', modal=False, transient=True)

  def body(self, master):
    
    #
    # Setup header 
    #

    row = 0

    label = Label(master, text= "Data model", fg = 'blue')
    label.grid(row=row, column=0, columnspan = 2, sticky=Tkinter.EW)

    label = Label(master, text= "Information from external file")
    label.grid(row=row, column=2, columnspan = 2, sticky=Tkinter.EW)

    row += 1

    label = Label(master, text= "Ccp chain code", fg = 'blue')
    label.grid(row=row, column=0, sticky=Tkinter.W)

    label = Label(master, text= "Sequence Id (code) start", fg = 'blue')
    label.grid(row=row, column=1, sticky=Tkinter.W)

    label = Label(master, text= "Format chain code")
    label.grid(row=row, column=2, sticky=Tkinter.W)

    label = Label(master, text= "Sequence code start")
    label.grid(row=row, column=3, sticky=Tkinter.W)

    #
    # Setup list of ccp chains and selection menus for other chains
    # 
    
    self.formatChainMenu = {}

    self.ccpCodeLow = {}
    self.formatChainCodeLow = {}
    self.formatRemoveCode = {}
    self.formatRemoveLabel = {}

    ccpChainLabels = self.ccpChainLabelDict.keys()
    ccpChainLabels.sort()
    
    ccpChainSeqIdLabels = {}
    
    for ccpChainLabel in ccpChainLabels:
    
      ccpChainSeqIdLabels[ccpChainLabel] = []
      self.ccpChainSeqIdDict[ccpChainLabel] = {}
      self.ccpChainSeqCodeDict[ccpChainLabel] = {}
 
      #
      # Make a dict to list codes with Ids for ccp
      #
      
      for i in range(0,len(self.ccpChainSeqIdCodes[ccpChainLabel][0])):
        seqId = self.ccpChainSeqIdCodes[ccpChainLabel][0][i]
        (seqCode,seqInsertCode) = self.ccpChainSeqIdCodes[ccpChainLabel][1][i]
      
        label = "%s (%s)" % (seqId,seqCode + string.strip(seqInsertCode))
        ccpChainSeqIdLabels[ccpChainLabel].append(label)
        self.ccpChainSeqIdDict[ccpChainLabel][label] = seqId
        self.ccpChainSeqCodeDict[ccpChainLabel][label] = seqCode
              
    #
    # Currently only linking start in ccp chain (ONLY ONE!) to start in chain read in from file
    # Assuming sequential order thereafter.
    #

    for ccpChainLabel in ccpChainLabels:

      row = row + 1

      label = Label(master, text= ccpChainLabel, fg = 'blue')
      label.grid(row=row, column=0, sticky=Tkinter.W)

      self.ccpCodeLow[ccpChainLabel] = ScrolledListbox(master, initial_list = ccpChainSeqIdLabels[ccpChainLabel],width = 4, height = 4, xscroll = False)
      self.ccpCodeLow[ccpChainLabel].grid(row=row, column=1, sticky=Tkinter.EW)

      self.formatChainMenu[ccpChainLabel] = PulldownMenu(master, entries = self.formatChainList, selected_index = self.formatChainList.index(self.defaultFormatChain))
      self.formatChainMenu[ccpChainLabel].grid(row=row, column=2, sticky=Tkinter.EW)

      self.formatChainCodeLow[ccpChainLabel] = ScrolledListbox(master, initial_list = [], width = 4, height = 4, xscroll = False)
      self.formatChainCodeLow[ccpChainLabel].grid(row=row, column=3, sticky=Tkinter.EW)
      
      self.formatRemoveLabel[ccpChainLabel] = Label(master, text = "", fg = 'red')
      self.formatRemoveCode[ccpChainLabel] = Tkinter.Button(master, text = "Remove", command = (lambda x = ccpChainLabel, row = row: self.removeSeqCode(x,row)))

      self.formatRemoveLabel[ccpChainLabel].grid(row=row, column=4, sticky=Tkinter.N)
      self.formatRemoveLabel[ccpChainLabel].grid_forget()
      
      self.formatRemoveCode[ccpChainLabel].grid(row=row, column=4, sticky=Tkinter.S)

      self.formatChainMenu[ccpChainLabel].callback = (lambda x = -1, y = self.defaultFormatChain, z = ccpChainLabel, row = row: self.updateSeqCodes(x,y,z,row))

    row = row + 1

    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands, dismiss_text = 'Exit', help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = 5)

  def removeSeqCode(self,ccpChainLabel,row):

    selection = self.formatChainMenu[ccpChainLabel].getSelected()
      
    if selection != 'Do not link':
    
      chainDictKey = (ccpChainLabel,selection)
           
      (formatChainCode,seqCodesList) = self.localFormatChainDict[chainDictKey]
      
      selectedItems = self.formatChainCodeLow[ccpChainLabel].getSelectedItems()
      
      for selectedItem in selectedItems:
        seqCodeKey = getSeqAndInsertCode(selectedItem)
        
        try:
          seqCodeIndex = seqCodesList.index(seqCodeKey)
          seqCodesList.pop(seqCodeIndex)
        except:
          pass
        
      self.updateSeqCodes(-1,selection,ccpChainLabel,row)
       
  def updateSeqCodes(self,formatChainIndex,formatChainText,ccpChainLabel,row):

    if formatChainText == 'Do not link':	  
    
      self.formatChainCodeLow[ccpChainLabel].setItems([])
    
    else:
    
      chainDictKey = (ccpChainLabel,formatChainText)
      
      if not self.localFormatChainDict.has_key(chainDictKey):
        (formatChainCode,seqCodesList) = self.formatChainDict[formatChainText]
        self.localFormatChainDict[chainDictKey] = (formatChainCode,seqCodesList[:])
    
      (formatChainCode,seqCodesList) = self.localFormatChainDict[chainDictKey]
      
      self.formatChainCodeLow[ccpChainLabel].clear()
      
      seqCodeJumps = 0
      oldSeqCode = seqCodesList[0][0] - 1
      oldSeqInsertCode = seqCodesList[0][1]
      
      for (seqCode,seqInsertCode) in seqCodesList:
        self.formatChainCodeLow[ccpChainLabel].append(str(seqCode) + seqInsertCode)
        
        if oldSeqInsertCode == seqInsertCode and oldSeqCode + 1 != seqCode:
          seqCodeJumps = 1
        
        oldSeqInsertCode = seqInsertCode
        oldSeqCode = seqCode
        
      self.formatChainCodeLow[ccpChainLabel].setSelectedItems([str(seqCodesList[0][0]) + seqCodesList[0][1]])
      
      #
      # Try to match on sequence code...
      #
      
      foundMatch = 0
      
      (formatFirstSeqCode,formatSeqInsertCode) = seqCodesList[0]
      matchKey = (str(formatFirstSeqCode),formatSeqInsertCode)

      if matchKey in self.ccpChainSeqIdCodes[ccpChainLabel][1]:
        matchKey = str(formatFirstSeqCode) + string.strip(formatSeqInsertCode)
        for label in self.ccpChainSeqCodeDict[ccpChainLabel].keys():
          if self.ccpChainSeqCodeDict[ccpChainLabel][label] == matchKey:
            self.ccpCodeLow[ccpChainLabel].setSelectedItems([label])
            foundMatch = 1
            break
            
      if not foundMatch:
      
        #
        # Try to match on sequence id...
        # This could all be more intelligent but leave for now.
        #
        
        if str(formatFirstSeqCode) in self.ccpChainSeqIdCodes[ccpChainLabel][0]:
          for label in self.ccpChainSeqIdDict[ccpChainLabel].keys():
            if self.ccpChainSeqIdDict[ccpChainLabel][label] == str(formatFirstSeqCode):
              self.ccpCodeLow[ccpChainLabel].setSelectedItems([label])
              foundMatch = 1
              break

      #
      # Set 'remove' if seqInsertCodes present, and print warning
      #
      
      multipleSeqInsertCodes = 0
      for chainDictKey in self.localFormatChainDict.keys():
        (tempFormatChainCode,tempSeqCodesList) = self.localFormatChainDict[chainDictKey]
        if formatSeqInsertCode != tempSeqCodesList[0][1]:
          multipleSeqInsertCodes += 1
      
      #
      # Set user warning messages...
      #
          
      text = ""
      if multipleSeqInsertCodes:
        text += "Warning: insertion codes present"
      if seqCodeJumps:
        text += "Warning: jumps in seqCode - assumed sequential"
      
      self.formatRemoveLabel[ccpChainLabel].grid_forget()
      if text:
        self.formatRemoveLabel[ccpChainLabel].set(text)
        self.formatRemoveLabel[ccpChainLabel].grid(row=row, column=4, sticky=Tkinter.N)      

  def apply(self):
  
    self.chainDict = {}
    
    for ccpChainLabel in self.ccpChainLabelDict.keys():
      
      chain = self.ccpChainLabelDict[ccpChainLabel]
      selection = self.formatChainMenu[ccpChainLabel].getSelected()
 
      if selection != 'Do not link':
         
        chainDictKey = (ccpChainLabel,selection)
        (formatChainCode,seqCodesList) = self.localFormatChainDict[chainDictKey]
        
        selectedItems = self.ccpCodeLow[ccpChainLabel].getSelectedItems()
        
        if not selectedItems:
          showError("No sequence ID selected","Please select a sequence ID on the left side for chain %s!" % ccpChainLabel)
          return False
        
        ccpCodeLow = returnInt(self.ccpChainSeqIdDict[ccpChainLabel][selectedItems[0]])

        (formatChainSeqCodeLow,seqInsertCode) = getSeqAndInsertCode(self.formatChainCodeLow[ccpChainLabel].getSelectedItems()[0])

        exportSeqCodesList = seqCodesList[:]

        for (tempSeqCode,tempSeqInsertCode) in seqCodesList:
          if tempSeqCode == formatChainSeqCodeLow and tempSeqInsertCode == seqInsertCode:
            break
          else:
            exportSeqCodesList.pop(0)
        
        self.chainDict[chain] = [formatChainCode,ccpCodeLow,exportSeqCodesList]
    
    return True
