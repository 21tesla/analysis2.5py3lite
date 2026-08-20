
"""
======================COPYRIGHT/LICENSE START==========================

ProcParToExpDimPopup.py: GUI window for matching proc. par file to experiment dimensions

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

from memops.universal.Constants import dirsep
from memops.universal.Io import joinPath

from ccpnmr.format.general.Io import getHelpUrlDir

from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Label import Label
from memops.gui.BasePopup import BasePopup
from memops.gui.Util import createHelpButtonList

from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

class ProcParToExpDimPopup(TemporaryBasePopup):
 
  help_url = joinPath(getHelpUrlDir(),'ProcParToExpDim.html')

  def __init__(self, parent, experiment, fPars, linkDict = None):
  
    self.experiment = experiment
    self.fPars = fPars
    self.linkDict = linkDict
    
    TemporaryBasePopup.__init__(self,parent = parent, title = "Project '%s': " % experiment.root.name + "Experiment dim selection for procPar reading", modal = False, transient=True)
 
  def body(self, master):

    #
    # Initialize
    #
    
    #
    # Make expDimRef lists
    #
    
    expDimRefList = []
    self.expDimRefDict = {}
    expDimRefSelection = {}
    
    for expDim in self.experiment.sortedExpDims():
      
      if expDim.dim == 1:
        addInfo = " (acqu)"
      else:
        addInfo = ""
    
      for expDimRef in expDim.expDimRefs:
        
        isotopeString = string.join(expDimRef.isotopeCodes,dirsep)
        
        selectionString = "Dim %d, nucl %s %s" % (expDim.dim,isotopeString,addInfo)

        expDimRefList.append(selectionString)
        self.expDimRefDict[selectionString] = expDimRef
        expDimRefSelection[expDimRef] = len(expDimRefList) - 1
        
    #
    # Popup info
    #
    
    #
    # Header labels
    #
    
    columnspan = 4

    row = 0
    
    if self.fPars.has_key('inputFile'):
      text = "Processing file %s" % (self.fPars['inputFile'])
    else:
      text = "Processing experiment %s" % (self.experiment.name)
    
    label = Label(master, text= text)
    label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)

    row += 1

    label = Label(master, text= "NumPoints")
    label.grid(row=row, column=0, sticky=Tkinter.EW)
    
    label = Label(master, text= "SW")
    label.grid(row=row, column=1, sticky=Tkinter.EW)
    
    label = Label(master, text= "Isotope")
    label.grid(row=row, column=2, sticky=Tkinter.EW)
    
    label = Label(master, text= "ExpDimRef selection")
    label.grid(row=row, column=3, sticky=Tkinter.EW)

    #
    # Selection per dimension from the procpar file
    #
    
    self.expDimRefMenu = []
    
    for dim in range(0,len(self.fPars['nucleus'])):

      row = row + 1

      label = Label(master, text= str(self.fPars['numPoints'][dim]))
      label.grid(row=row, column=0, sticky=Tkinter.EW)

      label = Label(master, text= "%.1f" % self.fPars['spectralWidth'][dim])
      label.grid(row=row, column=1, sticky=Tkinter.EW)
      
      label = Label(master, text= self.fPars['nucleus'][dim])
      label.grid(row=row, column=2, sticky=Tkinter.EW)
      
      expDimRef = self.linkDict[dim]
      selectedIndex = expDimRefSelection[expDimRef]

      self.expDimRefMenu.append(PulldownMenu(master, entries = expDimRefList, selected_index = selectedIndex))
      self.expDimRefMenu[-1].grid(row=row, column=3, sticky=Tkinter.E, ipadx = 20)
 
 
    row = row + 1
    texts = [ 'OK' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here
    buttons = createHelpButtonList(master, texts=texts, commands=commands, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan = columnspan)   

  def apply(self):
  
    self.expDimRefs = []
    expDims = []
    
    for dim in range(0,len(self.expDimRefMenu)):
      
      expDimRef = self.expDimRefDict[self.expDimRefMenu[dim].getSelected()]
      
      if self.expDimRefs.count(expDimRef) > 0 or expDims.count(expDimRef.expDim) > 0:
        self.expDimRefs = []
        expDims = []
        return False
      
      elif expDimRef:
        self.expDimRefs.append(expDimRef)
        expDims.append(expDimRef.expDim)
    
    if len(self.expDimRefs) != self.experiment.numDim:
      self.expDimRefs = []
      return False
    
    return True
