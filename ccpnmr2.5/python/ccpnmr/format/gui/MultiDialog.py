#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

MultiDialog.py: Sets MultiDialog interaction to GUI.

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

#
# Import the non-GUI multiDialog - this is the main resource for
# functions on top of the standard ones.
#

import sys


from memops.gui.MessageReporter import showInfo

from ccp.gui.NmrExpPrototypeEditor import NmrExpPrototypePopup

from ccpnmr.format.general.Util import getNameInfo
from ccpnmr.format.general.MultiDialog import MultiDialog
from ccpnmr.format.general.userInteraction import setupMessageReporter

from ccpnmr.format.gui.SelectionListPopup import SelectionListPopup
from ccpnmr.format.gui.EntryListPopup import EntryListPopup
from ccpnmr.format.gui.MultiSelectionListPopup import MultiSelectionListPopup
from ccpnmr.format.gui.ExperimentCreatePopup import ExperimentCreatePopup
from ccpnmr.format.gui.DataDimRefSelectByPeakDimPopup import DataDimRefSelectByPeakDimPopup
from ccpnmr.format.gui.ColumnInfoSelectionPopup import ColumnInfoSelectionPopup
from ccpnmr.format.gui.PeakDimSelectPopup import PeakDimSelectPopup
from ccpnmr.format.gui.ProcParToExpDimPopup import ProcParToExpDimPopup
from ccpnmr.format.gui.LinkResonancesPopup import LinkResonancesPopup
from ccpnmr.format.gui.ChainLinkPopup import ChainLinkPopup
from ccpnmr.format.gui.SequenceLinkPopup import SequenceLinkPopup
from ccpnmr.format.gui.NamingSystemPopup import NamingSystemPopup
from ccpnmr.format.gui.AtomSelectPopup import AtomSelectPopup
from ccpnmr.format.gui.ResonanceSelectPopup import ResonanceSelectPopup
from ccpnmr.format.gui.SingleResonanceStatusPopup import SingleResonanceStatusPopup
from ccpnmr.format.gui.ResonanceGroupPopup import ResonanceGroupPopup
from ccpnmr.format.gui.ChainExportPopup import ChainExportPopup
from ccpnmr.format.gui.FormatFilePopup import FormatFilePopup
from ccpnmr.format.gui.InfoPopup import InfoPopup
from ccpnmr.format.gui.ProchiralStatusPopup import ProchiralStatusPopup
from ccpnmr.format.gui.ChemCompSelectionPopup import ChemCompSelectionPopup
from ccpnmr.format.gui.GroupChemShiftsPopup import GroupChemShiftsPopup
from ccpnmr.format.gui.ChemCompInfoPopup import ChemCompInfoPopup
from ccpnmr.format.gui.ChemCompAtomBondInfoPopup import ChemCompAtomBondInfoPopup
from ccpnmr.format.gui.AcqProcParsEditPopup import AcqProcParsEditPopup
from ccpnmr.format.gui.ObjectEditPopup import EntryEditPopup

#
#
# TODO TODO: check here whether the parent was just closed(), and if so somehow
# have to finish running Format script...
#
# TODO: Look at OpenPopup in Analysis to abstract out the popup creation
#

class GuiMultiDialog(MultiDialog):

  def __init__(self):

    self.messageReporter = setupMessageReporter(1)  # This will set it to GUI
    
  def waitWindow(self,guiParent,popup,kw = None):
    
    #
    # Wait for child to finish
    #
    
    waitStatus = False
    
    if kw and (not kw.has_key('modal') or kw['modal'] == False):
      waitStatus = True
      guiParent.wait_window(popup)
    elif not kw:
      waitStatus = True
      guiParent.wait_window(popup)
    
    #
    # If killed, then stop the current thread...
    #
    
    if waitStatus and hasattr(popup,'stopProcess') and popup.stopProcess:
      raise "Warning: Aborting project - process stopped" 
          
  def FileName(self,guiParent,*args,**kw):
  
    return FormatFilePopup(guiParent,*args,**kw)

  def SelectionListFunction(self,guiParent,*args,**kw):
  
    popup = SelectionListPopup(guiParent,*args,**kw)
    
    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup
    
  def EntryListFunction(self,guiParent,*args,**kw):
  
    popup = EntryListPopup(guiParent,*args,**kw)
    
    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup

  def MultiSelectionListFunction(self,guiParent,*args,**kw):
  
    popup = MultiSelectionListPopup(guiParent,*args,**kw)

    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup

  def ExperimentTypesListFunction(self,guiParent,project,*args,**kw):
    
    return NmrExpPrototypePopup(guiParent,project,*args,**kw)
    
  def ExperimentCreateFunction(self,guiParent,project,*args,**kw):

    popup = ExperimentCreatePopup(guiParent,project,*args,**kw)

    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup
    
  def AcqProcParsEditFunction(self,guiParent,*args,**kw):
  
    popup = AcqProcParsEditPopup(guiParent,*args,**kw)

    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup

  def EntryEditFunction(self,guiParent,*args,**kw):
  
    popup = EntryEditPopup(guiParent,*args,**kw)

    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup
    
  def DataDimRefSelectByPeakDimFunction(self,guiParent,*args,**kw):
  
    popup = DataDimRefSelectByPeakDimPopup(guiParent,*args,**kw)

    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup

  def ColumnInfoSelectionFunction(self,guiParent,*args,**kw):
 
    popup = ColumnInfoSelectionPopup(guiParent,*args,**kw)

    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup

  def ProcParToExpDimFunction(self,guiParent,*args,**kw):
  
    popup = ProcParToExpDimPopup(guiParent,*args,**kw)

    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup
    
  def ChemCompInfoFunction(self,guiParent,*args,**kw):
  
    popup = ChemCompInfoPopup(guiParent,*args,**kw)
    
    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup

  def ChemCompAtomBondInfoFunction(self,guiParent,*args,**kw):
  
    popup = ChemCompAtomBondInfoPopup(guiParent,*args,**kw)
    
    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup
    
  def PeakDimSelectFunction(self,guiParent,*args,**kw):
  
    popup = PeakDimSelectPopup(guiParent,*args,**kw)

    self.waitWindow(guiParent,popup,kw = kw)
    
    return popup
    
  def LinkResonancesSetupFunction(self,guiParent,formatNamesList,formatNamesDict,title):

    popup = LinkResonancesPopup(guiParent,formatNamesList,formatNamesDict,title = title)

    self.waitWindow(guiParent,popup)
    
    return popup

  def ChainLinkFunction(self,guiParent,formatChains):

    popup = ChainLinkPopup(guiParent,self.ccpChainLabelDict,
                                     self.ccpChainSeqIdCodes,
                                     formatChains,
                                     self.formatChainList,
                                     self.formatChainDict,
                                     self.defaultFormatChain)

    self.waitWindow(guiParent,popup)
    
    return popup

  def chainExportFunction(self,guiParent,requireChainCode):

    popup = ChainExportPopup(guiParent,self.ccpChainLabelDict,
                                       self.ccpChainSeqIdCodes,
                                       requireChainCode)

    self.waitWindow(guiParent,popup)
    
    return popup

  def SequenceLinkFunction(self,guiParent,project,createMoleculeDict):

    popup = SequenceLinkPopup(guiParent,project,createMoleculeDict)
    
    self.waitWindow(guiParent,popup)
    
    return popup

  def NamingSystemFunction(self,guiParent):

    popup = NamingSystemPopup(guiParent,self.namingSysList,self.namingSysDict)

    self.waitWindow(guiParent,popup)
    
    return popup
    
  def AtomSelectFunction(self,guiParent,resName,**kw):

    popup = AtomSelectPopup(guiParent,resName,
                                      self.chemAtomList,
                                      self.chemAtomDict,
                                      self.propagateList,
                                      self.ccpCode,
                                      self.shifts,
                                      matchInfoPopup = self.resonanceInfo,
                                      **kw)
                                      
    self.waitWindow(guiParent,popup, kw = kw)
    
    return popup

  def CoordinateAtomSelectFunction(self,guiParent,seqLabel,atomName,**kw):

    popup = AtomSelectPopup(guiParent,atomName,
                                      self.chemAtomList,
                                      self.chemAtomDict,
                                      self.propagateList,
                                      seqLabel,
                                      None,
                                      matchType = 'Atom',
                                      selectedIndex = 1,
                                      title = "Project '%s': " % self.project.name + 'Select coordinate atom match',
                                      helpUrl = 'CoordinateAtomSelect.html',
                                      **kw)

    self.waitWindow(guiParent,popup, kw = kw)
    
    return popup
    
  def ResonanceSelectFunction(self,guiParent,resName,multi,title):
   
    return ResonanceSelectPopup(guiParent,
                                self.resonanceList,
                                self.resonanceDict,
                                self.ccpCode,
                                resName,
                                self.infoText,
                                multi,
                                title = title)

  def ResonanceGroupFunction(self,guiParent,title):
   
    return ResonanceGroupPopup(guiParent,
                               self.resonanceNameList,
                               self.resonanceLineDict,
                               self.ccpCode,
                               title = title)
    
  def SingleResonanceStatusFunction(self,guiParent, **keywds):

    return SingleResonanceStatusPopup(guiParent,self.message,self.nmrRes,self.optionList,self.title, **keywds)
    
  def MeasurementSelectFunction(self,guiParent):

    popup = SelectionListPopup(guiParent,self.measurementValueList,topText = self.topMessage, text = self.question, urlFile = 'MeasurementSelect', selectionDict = self.measurementValueDict, dismissButton = False)

    self.waitWindow(guiParent,popup)
    
    return popup

  def ChemCompSelectionFunction(self,guiParent, project, **keywds):

    popup = ChemCompSelectionPopup(guiParent, project, **keywds)

    self.waitWindow(guiParent,popup, kw = keywds)
    
    return popup

  def GroupChemShiftsFunction(self,guiParent, project, **keywds):

    return GroupChemShiftsPopup(guiParent, project)

  def ProchiralStatusFunction(self,guiParent):

    popup = ProchiralStatusPopup(guiParent,self.chainList,self.chainAtomSetsDict)

    self.waitWindow(guiParent,popup)
    
    return popup

  def resonanceInfo(self,guiParent,resonanceName):
    
    infoValues = getNameInfo(resonanceName)

    infoItems = ['Original chain code:','Original sequence code:','Original atom name:']

    if infoValues[2] != None and infoValues[1] == None:
      infoValues = infoValues[0:1] + infoValues[2:3] + infoValues[4:]
      infoItems[1] = 'Original spin system code:'
    elif infoValues[3]:
      #Combine seqCode and seqInsertCode
      infoValues = infoValues[0:1] + (str(infoValues[1]) + infoValues[3],) + infoValues[4:]
    else:
      infoValues = infoValues[0:2] + infoValues[4:]
      
    infoList = []
    infoDict = {}
    
    for i in range(0,len(infoItems)):
      
      infoList.append(infoItems[i])
      infoDict[infoItems[i]] = infoValues[i]
      
    infoText = "This name reflects the information from the imported file(s)."

    popup = InfoPopup(guiParent,'Resonance name info',infoText,infoList,infoDict)

    self.waitWindow(guiParent,popup)
    
    return popup

multiDialog = GuiMultiDialog()
