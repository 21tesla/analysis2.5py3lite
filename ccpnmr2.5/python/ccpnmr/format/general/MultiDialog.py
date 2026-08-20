"""
======================COPYRIGHT/LICENSE START==========================

MultiDialog.py: Sets up generic data for user interaction and handles text-based interaction.

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

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

import string, traceback, sys

from memops.universal.Util import returnInt, returnInts
from memops.universal.Util import returnFloat
from memops.universal.Io import joinPath

from ccp.general.Util import setCurrentStore

from ccp.format.general.Constants import atomOrder
from ccp.format.general.Constants import defaultSeqInsertCode

from ccpnmr.format.general.userInteraction import setupMessageReporter
from ccpnmr.format.general.Constants import stereoText, newline, tagSep
from ccpnmr.format.general.Util import createExperiment
from ccpnmr.format.general.Util import getNameInfo
from ccpnmr.format.general.Constants import ccpNmr_kw, assign_kw

#
# Get Api and related
#

from memops.api import Implementation
from ccp.api.nmr import Nmr, NmrConstraint

#
# Get processing scripts
#

from ccpnmr.format.process.linkResonances import linkResonances

"""

NON-GUI VERSION!

Is main resource for functions (other than GUI interaction)

The ...Function defs are redefined in gui/MultiDialog.py so GUI
windows are called instead of shell interaction (if available)

"""

class MultiDialog:

  def __init__(self):

    self.messageReporter = setupMessageReporter(None)
    self.screenList = ScreenList()

  """
  
  FileName
  
  """

  def FileName(self,guiParent,*args,**kw):

    return self.FileNameFunction(guiParent,*args,**kw)
    
  def FileNameFunction(self,guiParent,*args,**kw):
    
    nonGuiClass = NonGuiClass()
    
    print kw['selectionText']
    
    nonGuiClass.file = self.screenList.getString()

    return nonGuiClass

  """
  
  SelectionList
  
  """

  def SelectionList(self,guiParent,*args,**kw):

    return self.SelectionListFunction(guiParent,*args,**kw)
    
  def SelectionListFunction(self,guiParent,*args,**kw):
    
    selectionList = args[0]    

    nonGuiClass = NonGuiClass()
    nonGuiClass.isSelected = None
    
    for argKey in ['title','text']:
    
      if kw.has_key(argKey):
        print kw[argKey]
      
    if kw.has_key('dismissText'):
      dismissText = kw['dismissText']
    else:
      dismissText = None

    selected = self.screenList.pickStringListItem(selectionList,dismissText = dismissText)
    
    if selected != None:
    
      nonGuiClass.isSelected = True
    
      if kw.has_key('selectionDict'):
      
        nonGuiClass.selection = kw['selectionDict'][selected]
      
      else:
    
        nonGuiClass.selection = selected

    return nonGuiClass

  """
  
  MultiSelectionList
  
  """

  def MultiSelectionList(self,guiParent,*args,**kw):
  
    return self.MultiSelectionListFunction(guiParent,*args,**kw)
    
  def MultiSelectionListFunction(self,guiParent,*args,**kw):
    
    selectionList = args[0]    

    nonGuiClass = NonGuiClass()
    nonGuiClass.isSelectedList = []
    
    for argKey in ['title','text']:
    
      if kw.has_key(argKey):
        print kw[argKey]
      
    if kw.has_key('dismissText'):
      dismissText = kw['dismissText']
    else:
      dismissText = None

    selectedItems = self.screenList.pickStringListMulti(selectionList)
    
    if selectedItems:
    
      for selectedItem in selectedItems:
      
        if kw.has_key('selectionDict'):
      
          nonGuiClass.isSelectedList.append(kw['selectionDict'][selectedItem])
      
        else:
    
          nonGuiClass.isSelectedList.append(selectedItem)

    return nonGuiClass

  """
  
  EntryList
  
  """

  def EntryList(self,guiParent,*args,**kw):

    return self.EntryListFunction(guiParent,*args,**kw)
    
  def EntryListFunction(self,guiParent,*args,**kw):
    
    valueList = args[0]
    valueInfo = args[1]

    nonGuiClass = NonGuiClass()
    nonGuiClass.values = []
    
    for argKey in ['title','topText']:
    
      if kw.has_key(argKey):
        print kw[argKey]
      
    if kw.has_key('dismissText'):
      dismissText = kw['dismissText']
    else:
      dismissText = None

    for i in range(0,len(valueList)):

      value = valueList[i]
      
      print "%s (currently '%s'):" % (valueInfo[i],str(value))
      value = self.screenList.getString()
      
      if kw.has_key('conversionFunc'):
        value = kw['conversionFunc'](value)
      
      nonGuiClass.values.append(value)

    return nonGuiClass

  """
  
  ExperimentTypesList
  
  """

  def ExperimentTypesList(self,guiParent,project,*args,**kw):
    
   
    if self.messageReporter.showYesNo("Edit reference experiments","You will directly edit the reference experiments. If done incorrectly this might results in problems later on. Do you want to continue?"):
           
      # Need nucleus reference info set correctly  
      setCurrentStore(project,'ChemElementStore')
      
      interaction = self.ExperimentTypesListFunction(guiParent,project,*args,**kw)
      
    else:
    
      interaction = None
    
    return interaction
    
  def ExperimentTypesListFunction(self,guiParent,project,*args,**kw):
    
    return
   
  """
  
  ExperimentCreate
  
  """

  def ExperimentCreate(self,guiParent,project,*args,**kw):
    
    #
    # Need nucleus reference information
    #
    
    setCurrentStore(project,'ChemElementStore')
    
    interaction = self.ExperimentCreateFunction(guiParent,project,*args,**kw)

    #
    # TODO: Bit of a mess here for return values: return the interaction object or
    # just a set of values?!?!?
    #
    # Only use interaction if a .expTypeName was set...
    #
    
    if hasattr(interaction,'refExperiment'):
    
      interaction.experiment = createExperiment(project,interaction.name,interaction.refExperiment)

      if hasattr(guiParent,'setMenuState'):
        guiParent.setMenuState()
        
    else:
      
      interaction = None

    return interaction
    
  def ExperimentCreateFunction(self,guiParent,*args,**kw):
    
    return
    
  """
  
  AcqProcParsEdit (only GUI)
  
  """

  def AcqProcParsEdit(self,guiParent,*args,**kw):
  
    return self.AcqProcParsEditFunction(guiParent,*args,**kw)
    
  def AcqProcParsEditFunction(self,guiParent,*args,**kw):
    
    return
    
  """
  
  EntryEdit (only GUI)
  
  """

  def EntryEdit(self,guiParent,*args,**kw):
  
    return self.EntryEditFunction(guiParent,*args,**kw)
    
  def EntryEditFunction(self,guiParent,*args,**kw):
    
    return

  """
  
  DataDimRefSelectByPeakDim (only GUI)
  
  """

  def DataDimRefSelectByPeakDim(self,guiParent,*args,**kw):
  
    return self.DataDimRefSelectByPeakDimFunction(guiParent,*args,**kw)
    
  def DataDimRefSelectByPeakDimFunction(self,guiParent,*args,**kw):
    
    return

  """
  
  ColumnInfoSelection (only GUI)
  
  """

  def ColumnInfoSelection(self,guiParent,valueList,columnContents,*args,**kw):

    columnData = []
    funcs = []
    
    for i in range(len(valueList)):
      value = valueList[i][0]
      try:
        value = float(value)
        funcs.append(float)
      except:
        funcs.append(None)
        
      columnData.append([value,value])
    
    for i in range(len(valueList)):
      func = funcs[i]
      for value in valueList[i]:
        if func:
          value = func(value)
        if columnData[i][0] > value:
          columnData[i][0] = value
        if columnData[i][1] < value:
          columnData[i][1] = value

    return self.ColumnInfoSelectionFunction(guiParent,columnData,columnContents,*args,**kw)
    
  def ColumnInfoSelectionFunction(self,guiParent,columnData,columnContents,*args,**kw):
    
    return

  """
  
  ProcParToExpDim (only GUI)
  
  """

  def ProcParToExpDim(self,guiParent,*args,**kw):
  
    kw['linkDict'] = self.ProcParToExpDimInitialize(*args,**kw)
  
    return self.ProcParToExpDimFunction(guiParent,*args,**kw)
    
  def ProcParToExpDimFunction(self,guiParent,*args,**kw):
    
    return
    
  def ProcParToExpDimInitialize(self,*args,**kw):

    linkDict = {}
    
    experiment = args[0]
    fPars = args[1]
    
    expDims = list(experiment.sortedExpDims())

    for dim in range(0,len(fPars['nucleus'])):
    
      isotope = fPars['nucleus'][dim]
      expDimFound = 0
      
      for eDim in range(0,len(expDims)):
      
        expDim = expDims[eDim]
        
        for expDimRef in expDim.expDimRefs:
          if isotope in expDimRef.isotopeCodes:
            expDimFound = 1
            break
        
        if expDimFound:
          linkDict[dim] = expDimRef
          expDims.pop(expDims.index(expDim))
          break
      
      if not expDimFound:
        linkDict[dim] = None
  
    #
    # Set up rest...
    #

    for dim in range(0,len(fPars['nucleus'])):
      
      if not linkDict[dim]:
        linkDict[dim] = expDims[0].sortedExpDimRefs()[0]
        expDims.pop(0)
        
    return linkDict

  """
  
  ChemCompInfo (only GUI)
  
  """

  def ChemCompInfo(self,guiParent,*args,**kw):
  
    return self.ChemCompInfoFunction(guiParent,*args,**kw)
    
  def ChemCompInfoFunction(self,guiParent,*args,**kw):
    
    return

  """
  
  ChemCompAtomBondInfo (only GUI)
  
  """

  def ChemCompAtomBondInfo(self,guiParent,*args,**kw):
  
    return self.ChemCompAtomBondInfoFunction(guiParent,*args,**kw)
    
  def ChemCompAtomBondInfoFunction(self,guiParent,*args,**kw):
    
    return
    
  """
  
  PeakDimSelect (only GUI)
  
  """

  def PeakDimSelect(self,guiParent,*args,**kw):
  
    return self.PeakDimSelectFunction(guiParent,*args,**kw)
    
  def PeakDimSelectFunction(self,guiParent,*args,**kw):
    
    return
    
  """
  
  LinkResonancesSetup
  
  """

  def LinkResonancesSetup(self,guiParent,project,nmrConstraintStore,*args,**kw):
  
    if self.linkResonancesFindFormat(project,nmrConstraintStore):
    
      title = "Project '%s': " % project.name + 'Link resonances setup'
      
      if (kw.has_key('fastLinking') and kw['fastLinking']) or self.messageReporter.showYesNo("Use fast linking","Click Yes if you want to link automatically with default settings."):
      
        interaction = self.DefaultLinkResonancesSetup()

      else:

        interaction = self.LinkResonancesSetupFunction(guiParent,
                                                       self.formatNamesList,
                                                       self.formatNamesDict,
                                                       title)
                                                                    
      #
      # Run linkResonances straight away (if OK pressed)
      #
      
      if interaction.status:

        keywds = {}

        for setKeyword in ['assignFormat','globalStereoAssign',
                           'useIupacMatching','useLinkResonancePopup',
                           'minimalPrompts','forceShiftMerge','makeMethylsEquivalent',
                           'useAmbiguity','verbose']:

          keywds[setKeyword] = getattr(interaction,setKeyword)

        for setKeyword in ['setSingleProchiral','setSinglePossEquiv']:

          value = getattr(interaction,setKeyword)

          if value == 'Unknown':
            keywds[setKeyword] = None
          elif value == 'Always ignore':
            keywds[setKeyword] = 1
          else:
            keywds[setKeyword] = 0

        linkRes = linkResonances(project,guiParent = guiParent,nmrConstraintStore = nmrConstraintStore,**keywds)

        status = linkRes.status

        if status:
          self.messageReporter.showInfo("Success","linkResonances ran succesfully")

        else:
          self.messageReporter.showError("Error","There were problems running linkResonances.")
        
      return interaction

    else:
    
      self.messageReporter.showError('No external format data','There are unlinked resonances present, but they do not have information from external format files (and were probably created by Analysis). LinkResonances quitting.')
      
      return None
      
  def LinkResonancesSetupFunction(self,guiParent,formatNamesList,formatNamesDict,title):
    
    nonGuiClass = self.DefaultLinkResonancesSetup()
    
    return nonGuiClass

  def DefaultLinkResonancesSetup(self):
    
    nonGuiClass = NonGuiClass()
    nonGuiClass.assignFormat = self.formatNamesDict[self.formatNamesList[0]]
    nonGuiClass.globalStereoAssign = False
    nonGuiClass.useIupacMatching = True
    nonGuiClass.useLinkResonancePopup = False
    nonGuiClass.minimalPrompts = True
    nonGuiClass.forceShiftMerge = True
    nonGuiClass.useAmbiguity = False
    nonGuiClass.verbose = True
    nonGuiClass.setSingleProchiral = True
    nonGuiClass.setSinglePossEquiv = True
    nonGuiClass.makeMethylsEquivalent = True
    
    nonGuiClass.status = 1
    
    return nonGuiClass

  def linkResonancesFindFormat(self,project,nmrConstraintStore):

    #
    # Find most likely format where resonance names come from
    #
  
    formatsFound = {}
    unlinkedResonances = 0
    
    if nmrConstraintStore:
      resonances = nmrConstraintStore.sortedFixedResonances()
    else:
      resonances = project.currentNmrProject.sortedResonances()

    for resonance in resonances:
  
       #
       # Only consider non-linked resonances
       #
     
       if not resonance.resonanceSet:
     
         unlinkedResonances += 1
     
         assignData = resonance.findAllApplicationData(keyword = assign_kw)
       
         for assignDatum in assignData:
       
           format = assignDatum.application
       
           if formatsFound.has_key(format):
             formatsFound[format] += 1
          
           else:
             formatsFound[format] = 1    
  
    if formatsFound != {}:
      
      #
      # Make dict
      #
      
      formatPercentages = {}
      
      for format in formatsFound.keys():
      
        percentage = formatsFound[format] * 100.0 / unlinkedResonances
        
        if not formatPercentages.has_key(percentage):
          formatPercentages[percentage] = []
        
        formatPercentages[percentage].append(format)
      
      #
      # Order according to %
      #
      
      percentages = formatPercentages.keys()
      percentages.sort()
      percentages.reverse()
      
      self.formatNamesList = []
      self.formatNamesDict = {}
      
      for percentage in percentages:
        
        for format in formatPercentages[percentage]:
          formatText = "%s (%.1f %%)" % (format,percentage)
        
          self.formatNamesDict[formatText] = format
          self.formatNamesList.append(formatText)
      
      return True
 
    else:
    
      #
      # If there is no application data with names it doesn't make sense to
      # link (no info available). Give error message & return.
      #
      
      return False


  """
  
  ChainLink
  
  """

  def ChainLink(self,guiParent,formatChains,molSystem,forceDefaultChainMapping = 0, doInteraction = True):
  
    self.ChainLinkInitialize(formatChains,molSystem)

    if forceDefaultChainMapping:
      
      if len(self.ccpChainLabelDict) == 1 and len(formatChains) == 1:
      
        ccpChainLabel = self.ccpChainLabelDict.keys()[0]
        formatChainCode = formatChains.keys()[0]
        
        nonGuiClass = NonGuiClass()
        nonGuiClass.chainDict = {}

        #
        # Try to match on sequence code...
        #

        foundMatch = 0
        
        (formatFirstSeqCode,formatSeqInsertCode) = formatChains[formatChainCode][0]
        
        matchKey = (str(formatFirstSeqCode),formatSeqInsertCode)

        if matchKey in self.ccpChainSeqIdCodes[ccpChainLabel][1]:
          for i in range(0,len(self.ccpChainSeqIdCodes[ccpChainLabel][1])):
            seqCodeKey = self.ccpChainSeqIdCodes[ccpChainLabel][1][i]
            if seqCodeKey == matchKey:
              ccpCodeLow = returnInt(self.ccpChainSeqIdCodes[ccpChainLabel][0][i])
              foundMatch = 1
              break

        if not foundMatch:

          #
          # Try to match on sequence id...
          # This could all be more intelligent but leave for now.
          #

          if str(formatFirstSeqCode) in self.ccpChainSeqIdCodes[ccpChainLabel][0]:
            for seqId in self.ccpChainSeqIdCodes[ccpChainLabel][0]:
              if seqId == str(formatFirstSeqCode):
                ccpCodeLow = returnInt(seqId)
                foundMatch = 1
                break
                
        if foundMatch:
          nonGuiClass.chainDict[self.ccpChainLabelDict[ccpChainLabel]] = [formatChainCode,ccpCodeLow,formatChains[formatChainCode]]
          return nonGuiClass
        
    if not doInteraction:
      # Hack to allow forceDefaultChainMapping option even if self.allowPopups is False
      nonGuiClass = NonGuiClass()
      nonGuiClass.chainDict = {}
      return nonGuiClass
  
    return self.ChainLinkFunction(guiParent,formatChains)
    
  def ChainLinkFunction(self,guiParent,formatChains):
    
    nonGuiClass = NonGuiClass()
    nonGuiClass.chainDict = {}

    ccpChainList = self.ccpChainLabelDict.keys()
    ccpChainList.sort()
    
    for formatChainText in self.formatChainList[:-1]:

      formatChainDictInfo = self.formatChainDict[formatChainText]
      formatChainCode = formatChainDictInfo[0]
      seqCodesList = formatChainDictInfo[1]
    
      print "Pick corresponding Ccp chain for '%s', %d residues, %d first seqCode:" % (formatChainCode,len(seqCodesList),seqCodesList[0][0])
    
      chosenChain = self.screenList.pickStringListItem(ccpChainList + [self.defaultFormatChain])
      
      if chosenChain != self.defaultFormatChain:
        (formatFirstSeqCode,formatSeqInsertCode) = seqCodesList[0]
        
        if formatSeqInsertCode != defaultSeqInsertCode:
          seqInsertCodeText = formatSeqInsertCode
        else:
          seqInsertCodeText = ''

        print "Sequence code %d%s for chain '%s' corresponds to which sequence Id for %s?" % (formatFirstSeqCode,seqInsertCodeText,formatChainCode,chosenChain)
        
        print "  Sequence: %s" % ('-'.join(self.ccpChainCcpCodes[chosenChain]))
        
        chosenSeqId = self.screenList.pickIntListItem(returnInts(self.ccpChainSeqIdCodes[chosenChain][0]))
    
        nonGuiClass.chainDict[self.ccpChainLabelDict[chosenChain]] = [formatChainCode,chosenSeqId,seqCodesList]

    return nonGuiClass
    
  def ChainLinkInitialize(self,formatChains,molSystem):
    
    #
    # Ccp chain list creation
    #
    
    self.ccpChainLabelDict = {}
    self.ccpChainSeqIdCodes = {}
    self.ccpChainCcpCodes = {}
                           
    for chain in molSystem.sortedChains():
      
      if len(chain.residues) == 1 and chain.findFirstResidue().chemCompVar.chemComp.ccpCode == 'Hoh':
        continue
    
      chainLabel = "'" + chain.code + "' (%d res, molecule %s)" % (len(chain.residues),chain.molecule.name)

      self.ccpChainLabelDict[chainLabel] = chain
      self.ccpChainSeqIdCodes[chainLabel] = [[],[]]
      self.ccpChainCcpCodes[chainLabel] = []
      
      for residue in chain.sortedResidues():
        self.ccpChainSeqIdCodes[chainLabel][0].append(str(residue.seqId))
        self.ccpChainSeqIdCodes[chainLabel][1].append((str(residue.seqCode),residue.seqInsertCode))
        self.ccpChainCcpCodes[chainLabel].append(residue.ccpCode)
             
    #
    # Format file chain list creation
    #
    
    self.defaultFormatChain = 'Do not link'
    self.formatChainList = []
    self.formatChainDict = {}
    
    for formatChainCode in formatChains.keys():
      
      #
      # Also keep track of seqInsertCodes and 'sequence jumps'
      #
      
      seqInsertCodes = []
      seqInsertCodeDict = {}
      seqCodeStretchDict = {}

      for (seqCode,seqInsertCode) in formatChains[formatChainCode]:
        if seqInsertCode not in seqInsertCodes:
          seqInsertCodes.append(seqInsertCode)
          seqInsertCodeDict[seqInsertCode] = []
          seqCodeStretchDict[seqInsertCode] = []
          
        seqInsertCodeDict[seqInsertCode].append((seqCode,seqInsertCode))

        if not seqCodeStretchDict[seqInsertCode] or seqCodeStretchDict[seqInsertCode][-1][-1][0] + 1 != seqCode:
          seqCodeStretchDict[seqInsertCode].append([(seqCode,seqInsertCode)])
        else:
          seqCodeStretchDict[seqInsertCode][-1].append((seqCode,seqInsertCode))
      
      if len(seqInsertCodes) > 1:
        seqInsertCodes.insert(0,None)
      
      #
      # Now create the dictionary... first for insertion codes
      #
        
      seqInsertCodes.sort()
      for seqInsertCode in seqInsertCodes:
        if len(seqInsertCodes) == 1:
          seqInsertText = ''
          seqCodesList = formatChains[formatChainCode][:]
        elif not seqInsertCode:
          seqInsertText = ' (all)'
          seqCodesList = formatChains[formatChainCode][:]
        else:
          seqInsertText = " (insertion code '%s')" % seqInsertCode
          seqCodesList = seqInsertCodeDict[seqInsertCode][:]
        
        formatChainText = "Link to code '%s'%s" % (formatChainCode,seqInsertText)
        self.formatChainList.append(formatChainText)
        self.formatChainDict[formatChainText] = (formatChainCode,seqCodesList)
        
        #
        # Add bits where numbering is not sequential
        #
      
        if seqInsertCode and len(seqCodeStretchDict[seqInsertCode]) > 1:
          for seqCodeStretch in seqCodeStretchDict[seqInsertCode]:
            startSeqCodeText = "'%d%s'" % (seqCodeStretch[0][0],seqCodeStretch[0][1])
            endSeqCodeText = "'%d%s'" % (seqCodeStretch[-1][0],seqCodeStretch[-1][1])
            formatChainText = "Link to code '%s' (range %s-%s)" % (formatChainCode,startSeqCodeText,endSeqCodeText)
            self.formatChainList.append(formatChainText)
            self.formatChainDict[formatChainText] = (formatChainCode,seqCodeStretch)
      

    self.formatChainList.sort()
    self.formatChainList.append(self.defaultFormatChain)
    
    
  """
  
  ChainExport
  
  """

  def chainExport(self,guiParent,chains,requireChainCode = None):
  
    self.chainExportInitialize(chains)
  
    return self.chainExportFunction(guiParent,requireChainCode)
    
  def chainExportFunction(self,guiParent,requireChainCode):
    
    nonGuiClass = NonGuiClass()
    nonGuiClass.chainDict = {}

    ccpChainLabelList = self.ccpChainLabelDict.keys()
    ccpChainLabelList.sort()
    
    exportChainCodes = []    
    exportSeqCodeRanges = []
    
    for ccpChainLabel in ccpChainLabelList:
    
      seqLength = len(self.ccpChainSeqIdCodes[ccpChainLabel][0])
    
      print "Export for ccpChain %s" % ccpChainLabel
      
      print "SeqId (code): start %d (%d), end %d (%d)" % (self.ccpChainSeqIdCodes[ccpChainLabel][0][0],
                                                          self.ccpChainSeqIdCodes[ccpChainLabel][1][0],
                                                          self.ccpChainSeqIdCodes[ccpChainLabel][0][-1],
                                                          self.ccpChainSeqIdCodes[ccpChainLabel][1][-1])
      #
      # Get a chain code / seg id for this ccp chain (if required)
      #
      
      exportChainCode = None
      
      if requireChainCode:
      
        while (not exportChainCode):

          print "Choose export %s corresponding to this chain" % requireChainCode
          exportChainCode = self.screenList.getString()
          
          if exportChainCode in exportChainCodes:

            print "Already selected for export - please choose another"
            exportChainCode = None

          else:

            exportChainCodes.append(exportChainCode)
        
        if exportChainCode == newline:
        
          exportChainCode = None 
      
      #
      # Get the first seq code for export (always necessary)
      #

      exportFirstSeqCode = None
      
      while (not exportFirstSeqCode):

        print "Choose sequence code corresponding to seqId 1 in export:"
        exportFirstSeqCode = self.screenList.getInteger()
        
        if exportFirstSeqCode and not requireChainCode:
          
          for (lowSeqId,highSeqId) in exportSeqCodeRanges:
          
            if highSeqId >= exportFirstSeqCode >= lowSeqId:
            
              #
              # Illegal! For same chain anyway...
              #
              
              exportFirstSeqCode = None
              break
              
          if exportFirstSeqCode:
        
            exportSeqCodeRanges.append([exportFirstSeqCode,exportFirstSeqCode + seqLength - 1])
       
      #
      # TODO: Here have to add chain linker residues option (for DIANA)
      #

      nonGuiClass.chainDict[self.ccpChainLabelDict[ccpChainLabel]] = [exportChainCode,exportFirstSeqCode]

    return nonGuiClass
    
  def chainExportInitialize(self,chains):
    
    #
    # Ccp chain list creation
    #
    
    self.ccpChainLabelDict = {}
    self.ccpChainSeqIdCodes = {}
    
    for chain in chains:
      chainLabel = "'" + chain.code + "' (%d res, molecule %s)" % (len(chain.residues),chain.molecule.name)

      self.ccpChainLabelDict[chainLabel] = chain      
      self.ccpChainSeqIdCodes[chainLabel] = [[],[]]
      
      for residue in chain.sortedResidues():
        self.ccpChainSeqIdCodes[chainLabel][0].append(residue.seqId)
        self.ccpChainSeqIdCodes[chainLabel][1].append(residue.seqCode)

  """
  
  SequenceLink (only GUI)
  
  """

  def SequenceLink(self,guiParent,project,createMoleculeDict):
  
    self.SequenceLinkInitialize(createMoleculeDict)
      
    return self.SequenceLinkFunction(guiParent,project,createMoleculeDict)
    
  def SequenceLinkFunction(self,guiParent,project,createMoleculeDict):
  
    nonGuiClass = NonGuiClass()
    nonGuiClass.createMoleculeDict = createMoleculeDict
    return nonGuiClass
    
  def SequenceLinkInitialize(self,createMoleculeDict):
    
    pass

  """
  
  NamingSystem
  
  """

  def NamingSystem(self,guiParent,resNameFound,numResonances):
  
    self.namingSystemPercentages(resNameFound,numResonances)

    return self.NamingSystemFunction(guiParent)

  def NamingSystemFunction(self,guiParent):
    
    #
    # Default select most populated (TODO: need shell selection here!?)
    #
    
    nonGuiClass = NonGuiClass()
    
    print "Choose a naming system (% match in brackets):"

    namingSysString = self.screenList.pickStringListItem(self.namingSysList)
    nonGuiClass.namingSysName = self.namingSysDict[namingSysString]

    return nonGuiClass

  def namingSystemPercentages(self,resNameFound,numResonances):
  
    if numResonances == 0:
      numResonances = 1

    #
    # Prepare list of naming systems with percentages
    #
    
    self.namingSysDict = {}
    atomPercentDict = {}
  
    for namingSystemName in resNameFound:
      atomPercent = resNameFound[namingSystemName] * 100.0 / numResonances
    
      if not atomPercentDict.has_key(atomPercent):
        atomPercentDict[atomPercent] = []
    
      atomPercentDict[atomPercent].append(namingSystemName)
  
    atomPercentList = atomPercentDict.keys()
    atomPercentList.sort()
    atomPercentList.reverse()
  
    self.namingSysList = []
  
    for atomPercent in atomPercentList:
      for namingSystemName in atomPercentDict[atomPercent]:
        listName = namingSystemName + " (%.2f%%)" % atomPercent
        self.namingSysList.append(listName)
        self.namingSysDict[listName] = namingSystemName
    
  """
  CoordinateAtomSelect
  """
  
  def CoordinateAtomSelect(self,guiParent,nmrRes,atomName,seqLabel,*args,**kw):
  
    self.coordinateAtomSelectInitialize(seqLabel)
    self.coordinateAtomSelectDefineElementAtomSets(nmrRes)
    self.coordinateAtomSelectSuggestMatch(atomName)
    
    return self.CoordinateAtomSelectFunction(guiParent,seqLabel,atomName,**kw)
    
  def CoordinateAtomSelectFunction(self,guiParent,seqLabel,atomName,**kw):
    
    nonGuiClass = NonGuiClass()
    
    if kw.has_key('headerLines'):
      for headerLine in kw['headerLines']:
        print headerLine

    print "Residue type '%s', atom '%s'." % (seqLabel,atomName)
          
    print "Pick the relevant atom match for this coordinate atom:"

    selectionText = showAllCode = 'Show all atoms'
    doNotLinkText = 'Do not link'
    
    atomSelected = 1
    
    while (selectionText == showAllCode):

      selectionText = self.screenList.pickStringListItem([showAllCode] + self.chemAtomList + [doNotLinkText])
      
      if selectionText == showAllCode:
        self.chemAtomList = self.chemAtomDict.keys()
        self.chemAtomList.sort()
        
      elif selectionText == doNotLinkText:
        atomSelected = 0
        break
    
    if atomSelected:

      nonGuiClass.chemAtomOrSets = self.chemAtomDict[selectionText]
    
    else:
    
      nonGuiClass.chemAtomOrSets = None

    print
    print "Propagate this mapping to:"
    
    propagate = self.screenList.pickStringListItem(self.propagateList)
    nonGuiClass.propagate = self.propagateList.index(propagate)

    return nonGuiClass

  def coordinateAtomSelectInitialize(self,seqLabel):
    
    #
    # Define list of mapping propagation choices
    # NOTE: if this is changed, code has to be changed in DataFormat!!
    #

    self.propagateList = ['Only this residue',
                          "All residues of %s type in this chain" % seqLabel,
                          "All residues of %s type in coordinate file" % seqLabel]

  def coordinateAtomSelectDefineElementAtomSets(self,nmrRes):
    
    #
    # Do all elements and set project
    #
    
    self.elementAtomSets = {}
    self.project = nmrRes.root

    #
    # If it has a chemAtomSet, use that as key
    # If not, then use chemAtom as key
    #
      
    for atom in nmrRes.atoms:
    
      chemAtom = atom.chemAtom
      
      elementSymbol = chemAtom.elementSymbol
      
      if not self.elementAtomSets.has_key(elementSymbol):
        self.elementAtomSets[elementSymbol] = {}
      
      self.elementAtomSets[elementSymbol][chemAtom] = ElementAtomSet(chemAtom,nmrRes)

  def coordinateAtomSelectSuggestMatch(self,atomName):
    
    #
    # TODO TODO COMBINE WITH ATOMSELECT CODE! Is same...
    #
        
    self.chemAtomList = []
    self.chemAtomDict = {}
  
    elementSymbols = self.elementAtomSets.keys()
    elementSymbols.sort()
    
    #
    # Try determine atom type and position from resName
    #
    
    atomMatchElement = None
    atomMatchPosition = None
    elementSymbolIndexLow = 999
    
    for elementSymbol in elementSymbols:
      
      elementSymbolIndex = string.find(atomName, elementSymbol)
      
      if elementSymbolIndex != -1 and elementSymbolIndex < elementSymbolIndexLow:
        elementSymbolIndexLow = elementSymbolIndex
        atomMatchElement = atomName[elementSymbolIndex]
        
        if (len(atomName) - 1) > elementSymbolIndex:
          atomMatchPosition = atomName[elementSymbolIndex + 1]
      
      for memRef in self.elementAtomSets[elementSymbol].keys():
        # Also create self.chemAtomDict for all atoms
        # This should work OK as no atom names should occur twice...
        self.chemAtomDict.update(self.elementAtomSets[elementSymbol][memRef].getAtomDict())
    
    #  
    # Try get a match from the elementAtomSet class
    #
    
    suggestElementAtomSets = []
    defaultElementAtomSets = []
    elementAtomSet = None
    
    if atomMatchElement != None:
    
      if self.elementAtomSets.has_key(atomMatchElement):

        for memRef in self.elementAtomSets[atomMatchElement].keys():
          
          elementAtomSet = self.elementAtomSets[atomMatchElement][memRef]
          
          if elementAtomSet.position == atomMatchPosition:
            suggestElementAtomSets.append(elementAtomSet)
          
          elif elementAtomSet.position == None:
            defaultElementAtomSets.append(elementAtomSet)
            
    if not suggestElementAtomSets:
      if not defaultElementAtomSets:
        # Just pick last one found...
        suggestElementAtomSets = [elementAtomSet]
      else:
        suggestElementAtomSets = defaultElementAtomSets
    
    for suggestElementAtomSet in suggestElementAtomSets:
      if suggestElementAtomSet:
        self.chemAtomList.extend(suggestElementAtomSet.getAtomDict().keys())
    
    self.chemAtomList.sort()
  

  """
  
  AtomSelect
  
  """

  def AtomSelect(self,guiParent,nmrRes,resName,resonances,*args,**kw):
  
    self.atomSelectInitialize(nmrRes,resonances)
    self.atomSelectDefineElementAtomSets(nmrRes)
    self.atomSelectSuggestMatch(resName)
    
    kw['title'] = "Project '%s': " % nmrRes.root.name + 'Select resonance atom match'
    
    return self.AtomSelectFunction(guiParent,resName,**kw)
    
  def AtomSelectFunction(self,guiParent,resName, **kw):
    
    nonGuiClass = NonGuiClass()
    
    print "Residue type '%s', resonance '%s'." % (self.ccpCode,resName)
    
    if self.shifts:
      print "Shift values: %s ppm" % (self.shifts)
      
    print "Pick the relevant atom match for this resonance:"

    selectionText = showAllCode = 'Show all atoms'
    doNotLinkText = 'Do not link'
    
    atomSelected = 1
    
    while (selectionText == showAllCode):

      selectionText = self.screenList.pickStringListItem([showAllCode] + self.chemAtomList + [doNotLinkText])
      
      if selectionText == showAllCode:
        self.chemAtomList = self.chemAtomDict.keys()
        self.chemAtomList.sort()
        
      elif selectionText == doNotLinkText:
        atomSelected = 0
        break
    
    if atomSelected:

      nonGuiClass.chemAtomOrSets = self.chemAtomDict[selectionText]
    
    else:
    
      nonGuiClass.chemAtomOrSets = None

    print
    print "Propagate this mapping to:"
    
    propagate = self.screenList.pickStringListItem(self.propagateList)
    nonGuiClass.propagate = self.propagateList.index(propagate)

    return nonGuiClass

  def atomSelectInitialize(self,nmrRes,resonances):

    self.ccpCode = nmrRes.molResidue.ccpCode
    
    self.stereoText = stereoText
    
    #
    # Get list of associated chem shift values (if any)
    #
    
    self.shifts = ""
    
    for resonance in resonances:

      if isinstance(resonance,Nmr.Resonance) and resonance.shifts:
    
        for shift in resonance.shifts:
        
          self.shifts += "%.2f," % shift.value
      
        self.shifts = self.shifts[:-1]

    #
    # Define list of mapping propagation choices
    # NOTE: if this is changed, code has to be changed in linkResonances!
    #

    self.propagateList = ['Only this residue',
                          "All residues with seqId %s in same molecule(s)" % nmrRes.seqId,
                          "All residues of %s type in this chain" % self.ccpCode,
                          "All residues of %s type in same molecule(s)" % self.ccpCode,
                          "All residues of %s type in molSystem" % self.ccpCode,
                          'All residues in chain',
                          'All residues in same molecule(s)',
                          'All residues in molSystem']

  def atomSelectDefineElementAtomSets(self,nmrRes):
    
    #
    # Do all elements...
    #
    
    self.elementAtomSets = {}

    #
    # If it has a chemAtomSet, use that as key
    # If not, then use chemAtom as key
    #
      
    for atom in nmrRes.atoms:
      chemAtom = atom.chemAtom
      
      elementSymbol = chemAtom.elementSymbol
      
      if not self.elementAtomSets.has_key(elementSymbol):
        self.elementAtomSets[elementSymbol] = {}
      
      chemAtomSet = chemAtom.chemAtomSet
            
      if chemAtomSet != None:
      
        deepChemAtomSet = chemAtomSet.chemAtomSet
        
        if deepChemAtomSet:

          if not self.elementAtomSets[elementSymbol].has_key(deepChemAtomSet):
            self.elementAtomSets[elementSymbol][deepChemAtomSet] = ElementAtomSet(chemAtom, nmrRes, chemAtomSet = chemAtomSet, deepChemAtomSet = deepChemAtomSet)
        
          else:
            self.elementAtomSets[elementSymbol][deepChemAtomSet].addChemAtomSet(chemAtom, chemAtomSet)
          
        else:
          if not self.elementAtomSets[elementSymbol].has_key(chemAtomSet):
            self.elementAtomSets[elementSymbol][chemAtomSet] = ElementAtomSet(chemAtom, nmrRes, chemAtomSet = chemAtomSet)
        
          else:
            self.elementAtomSets[elementSymbol][chemAtomSet].addChemAtom(chemAtom)

      else:
        self.elementAtomSets[elementSymbol][chemAtom] = ElementAtomSet(chemAtom, nmrRes)

  def atomSelectSuggestMatch(self,resName):
  
    self.chemAtomList = []
    self.chemAtomDict = {}
  
    elementSymbols = self.elementAtomSets.keys()
    elementSymbols.sort()
    
    #
    # Try determine atom type and position from resName
    #
    
    (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName)
    atomMatchElement = None
    atomMatchPosition = None
    elementSymbolIndexLow = 999
    
    for elementSymbol in elementSymbols:
      
      elementSymbolIndex = string.find(atomName, elementSymbol)
      
      if elementSymbolIndex != -1 and elementSymbolIndex < elementSymbolIndexLow:
        elementSymbolIndexLow = elementSymbolIndex
        atomMatchElement = atomName[elementSymbolIndex]
        
        if (len(atomName) - 1) > elementSymbolIndex:
          atomMatchPosition = atomName[elementSymbolIndex + 1]
      
      for memRef in self.elementAtomSets[elementSymbol].keys():
        # Also create self.chemAtomDict for all atoms
        # This should work OK as no atom names should occur twice...
        self.chemAtomDict.update(self.elementAtomSets[elementSymbol][memRef].getAtomDict())
    
    #  
    # Try get a match from the elementAtomSet class
    #
    
    suggestElementAtomSets = []
    defaultElementAtomSets = []
    elementAtomSet = None
    
    if atomMatchElement != None:
    
      if self.elementAtomSets.has_key(atomMatchElement):

        for memRef in self.elementAtomSets[atomMatchElement].keys():
          
          elementAtomSet = self.elementAtomSets[atomMatchElement][memRef]
          
          if elementAtomSet.position == atomMatchPosition:
            suggestElementAtomSets.append(elementAtomSet)
          
          elif elementAtomSet.position == None:
            defaultElementAtomSets.append(elementAtomSet)
            
    if not suggestElementAtomSets:
      if not defaultElementAtomSets:
        # Just pick last one found...
        suggestElementAtomSets = [elementAtomSet]
      else:
        suggestElementAtomSets = defaultElementAtomSets
    
    for suggestElementAtomSet in suggestElementAtomSets:
      if suggestElementAtomSet:
        self.chemAtomList.extend(suggestElementAtomSet.getAtomDict().keys())
    
    self.chemAtomList.sort()
    
  """
  
  ResonanceSelect
  
  """

  def ResonanceSelect(self,guiParent,existingResonances, nmrRes, resName, multi):
  
    self.resonanceSelectInitialize(existingResonances, nmrRes, multi)
    
    title = "Project '%s': " % self.project.name + 'Select resonances'
    
    while (1):
      returnClass = self.ResonanceSelectFunction(guiParent,resName,multi,title)
      
      if returnClass.selectedResonances:
        return returnClass
           
  def ResonanceSelectFunction(self,guiParent,resName,multi,title):
    
    nonGuiClass = NonGuiClass()
    
    print "Residue type '%s', resonance '%s'." % (self.ccpCode,resName)
    
    print "Pick already existing resonances to link to %s:" % self.infoText
    print
    
    if multi:
      selectedResNames = self.screenList.pickStringListMulti(self.resonanceList)
    else:
      selectedResNames = self.screenList.pickStringListItem(self.resonanceList)

    nonGuiClass.selectedResonances = []
    
    for selectedResName in selectedResNames:
      nonGuiClass.selectedResonances.append(self.resonanceDict[selectedResName])
   
    return nonGuiClass
  
  def resonanceSelectInitialize(self,existingResonances, nmrRes, multi):
  
    if multi:
      self.infoText = '(multiple choices possible)'
    else:
      self.infoText = '(single selection only)'

    self.ccpCode = nmrRes.molResidue.ccpCode
    self.project = nmrRes.root
  
    #
    # Make 'decomposed' list where resonances with multiple matches are split up
    #
    
    self.resonanceList = []
    self.resonanceDict = {}
    
    #
    # Set up a list for the names
    #
    
    matchedResonances = {}
    
    for i in range(0,len(existingResonances.resonances)):
    
      resonance = existingResonances.resonances[i]
      resNames = existingResonances.resNames[i]
      
      matchedResonances[resonance] = string.join(resNames,"/")

    #
    # Set up the list/dict for the user interaction
    #

    for resonance in matchedResonances.keys():
      listResName = matchedResonances[resonance]
      
      #
      # TODO: should there be an option to link stereospecifically to
      #       an atom as well? Or should this be sorted out by the user?
      #
      
      otherAtomsNameList = []
      for otherAtomSet in resonance.resonanceSet.atomSets:
        otherAtomsNameList.append(otherAtomSet.name)
      otherAtomsNameList.sort()
      
      listName = "Resonance '%s' (to atoms '%s', orig names '%s')" % (resonance.name,string.join(otherAtomsNameList,','),listResName)
      self.resonanceList.append(listName)
      self.resonanceDict[listName] = resonance
      
    self.resonanceList.sort()

  """
  
  ResonanceGroup
  
  """

  def ResonanceGroup(self,guiParent,resNameList, nmrRes, allowInteraction):
  
    returnClass = self.resonanceGroupInitialize(resNameList, nmrRes)
    
    if not returnClass.resNameGroups and allowInteraction:
      title="Project '%s': " % nmrRes.root.name + 'Group resonances'
      returnClass = self.ResonanceGroupFunction(guiParent,title)       
    elif returnClass.resNameGroups:
      print "  Warning: autogrouped %s and %s." % (str(returnClass.resNameGroups[0]),str(returnClass.resNameGroups[1]))
      
    return returnClass
      
  def ResonanceGroupFunction(self,guiParent,title):
    
    nonGuiClass = NonGuiClass()
    
    print title
    print "Residue type '%s'." % self.ccpCode
    
    print "Which of the following resonances should be grouped?"
    print "(remaining ones will also be combined)"
    print
    
    selectedResNames = self.screenList.pickStringListMulti(self.resonanceNameList)
 
    resNameGroups = []
    resNameGroups.append([])
    resNameGroups.append([])
    
    for i in range(0,len(self.resonanceNameList)):
      
      resonanceNameText = self.resonanceNameList[i]
      
      if resonanceNameText in selectedResNames:
        resNameGroups[0].append(self.resonanceLineDict[i])
      else:
        resNameGroups[1].append(self.resonanceLineDict[i])
    
    nonGuiClass.resNameGroups = resNameGroups
   
    return nonGuiClass

  def resonanceGroupInitialize(self,resNameList, nmrRes):
  
    self.ccpCode = nmrRes.molResidue.ccpCode
      
    self.resonanceNameList = []
    self.resonanceLineDict = {}
    
    returnClass = NonGuiClass()
    returnClass.resNameGroups = []
    
    numberString = '0123456789'
    wildcardString = '#%*+'
    
    #
    # Set up the list/dict for the user interaction
    #
  
    for i in range(0,len(resNameList)):

      resName = resNameList[i]
          
      listName = "Resonance '%s'" % resName
      self.resonanceNameList.append(listName)
      self.resonanceLineDict[i] = resName
      
    #
    # Do a quick comparison check...
    #
    
    tempGroups = []
    tempGroups.append([self.resonanceLineDict[0]])
    
    for i in range(1,len(self.resonanceLineDict)):
    
      resName = self.resonanceLineDict[i]
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName)
      
      matchFound = False
      
      for tempGroup in tempGroups:
      
        (refChainCode,refSeqCode,refSpinSystemId,refSeqInsertCode,refAtomName) = getNameInfo(tempGroup[0])
        
        if seqCode == refSeqCode and seqInsertCode == refSeqInsertCode:        
          if (refAtomName.count(atomName) or atomName.count(refAtomName)) or (atomName[-2] in numberString and atomName[:-1] == refAtomName[:-1] and atomName[-2] == refAtomName[-2]):
            if (atomName[-1] in numberString + wildcardString and \
                refAtomName[-1] in numberString + wildcardString):

              if atomName[-1] in numberString:
                tempGroup.append(resName)
              else:
                tempGroup.insert(0,resName)

              matchFound = True
              break
        
      
      if not matchFound:
        tempGroups.append([resName])
        
    if len(tempGroups) == 2:
    
      returnClass.resNameGroups = [tempGroups[0],tempGroups[1]]

    return returnClass      
      
  """
  
  SingleProchiralStatus
  
  """

  def SingleProchiralStatus(self,guiParent,refAtomSet,otherAtomName):
  
    self.singleProchiralStatusInitialize(refAtomSet,otherAtomName)
  
    return self.SingleResonanceStatusFunction(guiParent, urlFile = 'SingleProchiralStatus')
     
  def singleProchiralStatusInitialize(self,refAtomSet,otherAtomName):
  
    self.nmrRes = refAtomSet.findFirstAtom().residue
    
    self.optionList = ["%s is assigned and has the same chemical shift as %s." % (otherAtomName,refAtomSet.name),
                       "%s is not assigned, the chemical shift is unknown." % otherAtomName#,
                       #"%s not assigned, this one stereo assigned" % otherAtomName
                      ]
    
    self.message = "Atom(set) %s is prochiral, but nothing is known about the other prochiral atom(set) %s.\nPlease select its status below." % (refAtomSet.name,otherAtomName)
    self.title = "Project '%s': " % self.nmrRes.root.name + 'Prochiral status selection'

  """
  
  SingleResonanceStatus
  
  """

  def SingleResonanceStatus(self,guiParent,refAtom,otherAtomName):
  
    self.SingleResonanceStatusInitialize(refAtom,otherAtomName)
  
    return self.SingleResonanceStatusFunction(guiParent, urlFile = 'SingleResonanceStatus')
    
  def SingleResonanceStatusFunction(self,guiParent, **keywds):
    
    nonGuiClass = NonGuiClass()
    
    print
    print "Residue '%s-%d'" % (self.nmrRes.molResidue.ccpCode,self.nmrRes.seqCode)
    
    print self.message
    print
    
    selectedOption = self.screenList.pickStringListItem(self.optionList)
    
    nonGuiClass.singleResonanceStatus = self.optionList.index(selectedOption)
       
    return nonGuiClass
 
  def SingleResonanceStatusInitialize(self,refAtom,otherAtomName):
  
    self.nmrRes = refAtom.residue
    
    self.optionList = ["%s and %s have the same chemical shift (equivalent)." % (otherAtomName,refAtom.name),
                       "%s and %s have a different chemical shift (not equivalent)." % (otherAtomName,refAtom.name)
                      ]
    
    self.message = "The atoms %s and %s usually have the same chemical shift,\nbut in this case some separate information was identified.\nPlease select their status below." % (refAtom.name,otherAtomName)
    self.title = "Project '%s': " % self.nmrRes.root.name + 'Give status possible equivalent resonance'

  """
  
  MeasurementSelect
  
  """

  def MeasurementSelect(self,guiParent,resName,resonance,newValue,measurementList,measurementApiCode):
  
    selectValue = self.measurementSelectInitialize(resonance,measurementList,resName,newValue,measurementApiCode)
    
    if selectValue != None:
      returnClass = NonGuiClass()
      returnClass.selection = selectValue
    else:
      returnClass = self.MeasurementSelectFunction(guiParent)

    return returnClass

  def MeasurementSelectFunction(self,guiParent):
    
    #
    # Default select most populated (TODO: need shell selection here!?)
    #
    
    nonGuiClass = NonGuiClass()
    
    print self.topMessage
    print self.question
    
    selectedMeasurement = self.screenList.pickStringListItem(self.measurementValueList)
    
    nonGuiClass.selection = self.measurementValueDict[selectedMeasurement]

    return nonGuiClass

  def measurementSelectInitialize(self,resonance,measurementList,resName,newValue,measurementApiCode):
  
    #
    # Get the measurements
    #
    
    selectValue = None
    
    findFunction = getattr(resonance,'findFirst%s' % measurementApiCode)
    measurementValue = findFunction(parentList = measurementList).value
    measurementValueString = str(measurementValue)
    
    newValueString = str(newValue)

    mergedValue = (measurementValue + newValue) /2
    mergedValueString = str(mergedValue) + " (merged)"

    #
    # Get a list of corresponding atoms for this resonance (if exist)
    #
  
    resonanceSet = resonance.resonanceSet
    
    if resonanceSet:

      atomNames = []
      
      for atomSet in resonanceSet.atomSets:
        atomNames.append(atomSet.name)
        
      atomNameText = "(connected to atoms %s) " % string.join(atomNames,'/')

    else:

      atomNameText = ""
    
    #
    # If other resonances connected to this resonanceSet,
    # get their name and measurement for this measurement list
    #
    
    otherResonances = []
    otherMeasurementValues = []
    
    if resonanceSet:
      for otherResonance in resonanceSet.resonances:
        if otherResonance != resonance:
          findFunction = getattr(otherResonance,'findFirst%s' % measurementApiCode)
          otherMeasurement = findFunction(parentList = measurementList)
          if otherMeasurement:
            otherResonances.append(otherResonance)
            otherMeasurementValues.append(otherMeasurement.value)
    
    #
    # Check if can match automatically...
    #
    
    if otherMeasurementValues:
      if measurementValue in otherMeasurementValues and newValue not in otherMeasurementValues:
        selectValue = newValue
      elif measurementValue not in otherMeasurementValues and newValue in otherMeasurementValues:
        selectValue = measurementValue
    
    if selectValue == None:

      #
      # Setup the message list
      #
      
      measurementText = "%s value" % measurementApiCode.lower()

      self.topMessage = newline
      self.topMessage += ("Resonance '%s' %shas %s %.3f, " % (resName,atomNameText,measurementText,measurementValue))
      self.topMessage += "but %s %.3f also appears possible." % (measurementText,newValue)

      if otherMeasurementValues:

        self.topMessage += newline * 2

        for i in range(0,len(otherMeasurementValues)):
          value = otherMeasurementValues[i]
          name = otherResonances[i].name
          self.topMessage += "(other resonance %s connected to these atoms has %s %.3f)" % (name,measurementText,value) + newline

      self.question = "Which value do you want to use:"

      #
      # Setup the list of measurements (as strings), and a dict to keep track
      #

      self.measurementValueList = []
      self.measurementValueDict = {}

      self.measurementValueList.append(measurementValueString)
      self.measurementValueDict[measurementValueString] = measurementValue

      self.measurementValueList.append(newValueString)
      self.measurementValueDict[newValueString] = newValue

      self.measurementValueList.append(mergedValueString)
      self.measurementValueDict[mergedValueString] = mergedValue
    
    return selectValue

  """
  
  ChemCompSelection (only GUI)
  
  """

  def ChemCompSelection(self,guiParent,project,*args,**kw):

    return self.ChemCompSelectionFunction(guiParent,project,*args,**kw)
    
  def ChemCompSelectionFunction(self,guiParent,project,*args,**kw):

    return None

  """
  
  GroupChemShifts (only GUI)
  
  """

  def GroupChemShifts(self,guiParent,project):
  
    chainList = self.groupChemShiftsInitialize(project)
    
    if chainList and project.currentNmrProject.resonances:
      
      return self.GroupChemShiftsFunction(guiParent,project)
    
    else:
      
      text = "Cannot activate this window because:\n"
      
      if not chainList:
        text += "No chains are present (try creating a molecular system).\n"
      
      if not project.currentNmrProject.resonances:
        text += "No resonances are present (try reading in NMR data files).\n"
        
      self.messageReporter.showError("Error",text)
      
      return None
      
  def GroupChemShiftsFunction(self,guiParent,project):
   
    return None

  def groupChemShiftsInitialize(self,project):
      
    chainList = []
    
    #
    # See if any chain lists present...
    #
    
    for molSys in project.sortedMolSystems():
      for chain in molSys.sortedChains():
        chainList.append(chain)
    
    return chainList    

  """
  
  ProchiralStatus (only GUI)
  
  """

  def ProchiralStatus(self,guiParent,project):
  
    self.prochiralStatusInitialize(project)
    
    if self.chainList and self.chainAtomSetsDict:
      
      return self.ProchiralStatusFunction(guiParent)
    
    else:
      
      text = "Cannot run this window because:\n"
      
      if not self.chainList:
        text += "No chains are present (try creating a molecular system).\n"
      
      if not self.chainAtomSetsDict:
        text += "No (fixed) resonances are present and connected to atoms (try reading in NMR data files or running linkResonances).\n"
        
      self.messageReporter.showError("Error",text)
      
      return None
      
  def ProchiralStatusFunction(self,guiParent):
   
    return None

  def prochiralStatusInitialize(self,project):
      
    self.chainList = []
    self.chainAtomSetsDict = {}
    
    resParents = project.sortedNmrProjects()
    resParents.extend(project.sortedNmrConstraintStores())
    
    for resParent in resParents:
      self.chainAtomSetsDict[resParent] = {}
    
    #
    # Set up the list/dict for the user interaction
    #
    
    for molSys in project.sortedMolSystems():
      for chain in molSys.sortedChains():
        self.chainList.append(chain)
        
        for resParent in resParents:
          self.chainAtomSetsDict[resParent][chain] = []
          
        for nmrRes in chain.sortedResidues():
          atomList = list(nmrRes.atoms)
          atomSetsList = {}
          atomSetsListIndex = {}
          
          while (atomList):          
            atom = atomList[0]
            atomSetDict = {}
            
            #
            # Set atomSet for each atomSet, FixedAtomSet...
            #
            
            for resParent in resParents:
              if isinstance(resParent,NmrConstraint.NmrConstraintStore):
                atomSet = atom.findFirstFixedAtomSet(nmrConstraintStore = resParent)
              else:
                atomSet = atom.atomSet
            
              if atomSet:
                atomSetDict[resParent] = atomSet
                if not atomSetsList.has_key(resParent):
                  atomSetsList[resParent] = []
            
            #
            # If any atomSets found, continue...
            #
            
            if atomSetDict:
              #atomSet = atom.atomSet
              chemAtom = atom.chemAtom
              chemAtomSet = chemAtom.chemAtomSet
            
              if chemAtomSet:
                deepChemAtomSet = chemAtomSet.chemAtomSet
                prochiralChemAtoms = []
                
                if chemAtomSet.isProchiral:
                  prochiralChemAtoms = chemAtomSet.sortedChemAtoms()
                  
                elif deepChemAtomSet and deepChemAtomSet.isProchiral:
                  chemAtomSets = deepChemAtomSet.sortedChemAtomSets()
                  otherChemAtomSet = chemAtomSets[not chemAtomSets.index(chemAtomSet)]
                  
                  prochiralChemAtoms = otherChemAtomSet.sortedChemAtoms()
               
                if prochiralChemAtoms:
                  
                  #
                  # Do for all atomSets (Fixed and normal)
                  #
                  
                  for resParent in atomSetDict.keys():
                  
                    atomSet = atomSetDict[resParent]
                    atomPositionCode = atomSet.name[1]
                    atomSetsListIndex[resParent] = len(atomSetsList[resParent])
                  
                    for i in range(0,len(atomSetsList[resParent])):
                      if atomOrder.index(atomPositionCode) < atomOrder.index(atomSetsList[resParent][i][0].name[1]):
                        atomSetsListIndex[resParent] = i
                        break
                  
                    atomSetsList[resParent].insert(atomSetsListIndex[resParent],[atomSet])
                  
                  #
                  # Generic again... just use last one to get rid of other atoms
                  # This assumes that methyl protons are always in an atomSet!
                  #
                  
                  for atom in atomSet.atoms:
                    atomList.pop(atomList.index(atom))
                  
                  atomIndex = 0
                  
                  while (atomIndex < len(atomList)):
                    
                    atom = atomList[atomIndex]
                    
                    if atom.chemAtom in prochiralChemAtoms:
                      atomList.pop(atomIndex)
                                            
                    else:
                      atomIndex += 1

                  otherAtomList = []
                  
                  for prochiralChemAtom in prochiralChemAtoms:
                    otherAtomList.append(nmrRes.findFirstAtom(chemAtom = prochiralChemAtom))
                  
                  if otherAtomList:
                  
                    for resParent in atomSetDict.keys():
                      if isinstance(resParent,NmrConstraint.NmrConstraintStore):
                        otherAtomSet = otherAtomList[0].findFirstFixedAtomSet(nmrConstraintStore = resParent)
                        atomSetClass = NmrConstraint.FixedAtomSet
                      else:
                        otherAtomSet = otherAtomList[0].atomSet
                        atomSetClass = Nmr.AtomSet
                    
                      if not otherAtomSet:
                        if len(otherAtomList) == 1:
                          atomSetName = otherAtomList[0].name
                        else:
                          atomSetName = otherAtomList[0].name[:-1] + '*'
                        
                        otherAtomSet = atomSetClass(resParent,atoms = otherAtomList, name = atomSetName)

                      if atomSet.name < otherAtomSet.name:
                        atomSetsList[resParent][atomSetsListIndex[resParent]].append(otherAtomSet)
                      else:
                        atomSetsList[resParent][atomSetsListIndex[resParent]].insert(0,otherAtomSet)
                                      
                  continue


            atomList.pop(0)

          #
          # Set the information...
          #
          
          for resParent in atomSetsList.keys():         
            if atomSetsList[resParent]:
              self.chainAtomSetsDict[resParent][chain].extend(atomSetsList[resParent])
        
    #
    # Remove if empty
    #

    for resParent in resParents:
      if self.chainAtomSetsDict.has_key(resParent):
        for molSys in project.sortedMolSystems():
          for chain in molSys.sortedChains():
            if self.chainAtomSetsDict[resParent].has_key(chain) and not self.chainAtomSetsDict[resParent][chain]:
              del(self.chainAtomSetsDict[resParent][chain])
          if not self.chainAtomSetsDict[resParent]:
            del(self.chainAtomSetsDict[resParent])

#
# Other classes/functions
#

class ScreenList:

  def __init__(self):
  
    self.screenLen = 60

  def screenSetup(self,maxStringLen):

    self.numCols = self.screenLen / maxStringLen
    self.format = "%%-%ds" % maxStringLen
  
  def printColumnString(self,i,string):
  
    print self.format % string,
    
    if self.numCols and (i+1)/self.numCols == float(i+1) / self.numCols:
      print

  def pickChoice(self,choiceList, dismissText = None):
  
    choice = ''
    
    choiceText = 'Choose number'
    
    if dismissText:
      choiceText += ('(pick -1 to %s)' % dismissText)
      choiceList.append(-1)
      
    choiceText += ': '
    
    while (choice not in choiceList):
    
      print choiceText,
  
      choice = raw_input()

      if str(returnInt(choice)) == choice:
        choice = returnInt(choice)

    if dismissText and choice == -1:
      return None
    
    else:
      return choice
      
  def getInteger(self):
  
    choice = None
    
    while (not choice):
    
      choice = raw_input()
      
      if str(returnInt(choice)) == choice:
        choice = returnInt(choice)
      
      else:
        choice = None
        
    return choice

  def getString(self):
  
    choice = raw_input()
        
    return choice

  def pickStringListItem(self,stringList, dismissText = None):
  
    #
    # Handles strings
    #
  
    maxStringLen = -1

    for string in stringList:
      
      if maxStringLen < len(string):
        maxStringLen = len(string)
      
    self.screenSetup(maxStringLen + 4)

    print
  
    for i in range(len(stringList)):
      
      self.printColumnString(i,str(i) + ': ' + stringList[i])

    print;print
    
    choice = self.pickChoice(range(len(stringList)), dismissText = dismissText)
  
    if choice != None:
      return stringList[choice]
    else:
      return None


  def pickStringListMulti(self,stringList):
  
    #
    # Handles strings
    #
  
    maxStringLen = -1

    for string in stringList:
  
      if maxStringLen < len(string):
        maxStringLen = len(string)
      
    self.screenSetup(maxStringLen + 4)

    print
    print "(choose -1 to end multi selection)"
  
    for i in range(0,len(stringList)):
      
      self.printColumnString(i,str(i) + ': ' + stringList[i])

    print;print
    
    choiceList = []
    
    while (1):
      choice = self.pickChoice(range(-1,len(stringList)))
      
      if choice != -1:
        choiceString = stringList[choice]
        if choiceList.count(choiceString) == 0:
          print "  Added %s." % choiceString
          print
          choiceList.append(choiceString)
      else:
        break
  
    return choiceList

  def pickIntListItem(self,valueList):
  
    #
    # Handles integers
    #
  
    maxStringLen = -1

    for value in valueList:
    
      string = str(value)  

      if maxStringLen < len(string):
        maxStringLen = len(string)
      
    self.screenSetup(maxStringLen)

    print
  
    for i in range(0,len(valueList)):
      
      self.printColumnString(i,str(valueList[i]))

    print;print
    
    choice = self.pickChoice(valueList)
  
    return choice


class NonGuiClass:

  """
  
  Fake object to return info as if it comes from a GUI popup
  
  """

  def __init__(self):
  
    self.status = False
    
  def destroy(self):
    
    #
    # TODO: this doesn't work...
    #
    
    del(self)


class ElementAtomSet:

  """
  
  Data class for getting atom/group names...

  """
  
  def __init__(self,chemAtom,nmrRes,chemAtomSet = None, deepChemAtomSet = None):
    
    self.stereoText = stereoText
    
    self.nmrRes = nmrRes
    self.chemAtoms = [chemAtom]
    self.chemAtomSets = []
    self.deepChemAtomSet = deepChemAtomSet
    
    if len(chemAtom.name) > 1:
      self.position = chemAtom.name[1]
    else:
      self.position = None
    
    if chemAtomSet:
      self.chemAtomSets = [chemAtomSet]
            
  def addChemAtomSet(self,chemAtom,chemAtomSet):
  
    if not chemAtomSet in self.chemAtomSets:
      self.chemAtomSets.append(chemAtomSet)
      
    self.chemAtoms.append(chemAtom)
  
  def addChemAtom(self,chemAtom):
  
    self.chemAtoms.append(chemAtom)
  
  def getAtomDict(self):
    
    atomDict = {}

    prochiral = 0
    equivalent = 0
    addText = ""
    chemCompVar = self.nmrRes.chemCompVar
    
    if self.deepChemAtomSet:
      atomDict[self.deepChemAtomSet.name] = [self.deepChemAtomSet]
      
      if self.deepChemAtomSet.isProchiral == True:
        prochiral = 1
        
        # Check if is really prochiral for this residue (chemCompVar)
        chemAtomSets = self.deepChemAtomSet.sortedChemAtomSets()
        chemAtomSet = chemAtomSets[0]
        otherChemAtomSet = chemAtomSets[1]
            
    if prochiral == 1:
      addText = self.stereoText
      
      noStereoText = ""
      
      chemAtomSets = self.deepChemAtomSet.sortedChemAtomSets()
      for chemAtomSet in chemAtomSets:
        noStereoText += chemAtomSet.name + '/'
      
        #
        # NOT same info as in self.deepChemAtomSet: is either one or the other!
        #
      
      atomDict[noStereoText[:-1]] = chemAtomSets
    
    # This way chemAtomSets that contain only one atom (from the loc level!)
    # can be ignored
    # These should not have a deepChemAtomSet in this case so shouldn't be a problem...
    if len(self.chemAtoms) > 1:
    
      for chemAtomSet in self.chemAtomSets:
    
        atomDict[chemAtomSet.name + addText] = [chemAtomSet]
      
        if chemAtomSet.isProchiral == True:
          prochiral = 1
                 
        else:
          prochiral = 0
        
        if chemAtomSet.isEquivalent == True:
          equivalent = 1
    
    # This true only when prochiral: so will not be used for Val HG1* and the likes
    # Means that comes from one chemAtomSet only in previous loop
    if prochiral == 1:
      addText = self.stereoText
      
      noStereoText = ""

      chemAtoms = chemAtomSet.sortedChemAtoms()
      for chemAtom in chemAtoms:
        noStereoText += chemAtom.name + '/'
        
      atomDict[noStereoText[:-1]] = chemAtoms

    else:
      addText = ""
    
    # Only list single atoms if not equivalent (so don't list Ala HB1,HB2,HB3)
    if equivalent != 1:
      for chemAtom in self.chemAtoms:
        atomDict[chemAtom.name + addText] = [chemAtom]
    
    return atomDict

    
"""

NON-GUI VERSION!

"""

multiDialog = MultiDialog()
