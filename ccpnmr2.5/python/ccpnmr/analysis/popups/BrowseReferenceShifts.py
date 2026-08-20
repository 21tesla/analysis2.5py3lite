"""
======================COPYRIGHT/LICENSE START==========================

BrowseReferenceShifts.py: Part of the CcpNmr Analysis program

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
import Tkinter
import re
import os

from memops.universal.Io import joinPath, splitPath

from memops.general import Implementation


from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Color import hexRepr,hsbToRgb
from memops.gui.Label import Label
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.ScrolledDensityMatrix import ScrolledDensityMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccp.general.Constants import standardResidueCcpCodes

#from ccpnmr.analysis.macros.ArgumentServer import ArgumentServer
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.MoleculeBasic import greekSortAtomNames

AMINO_ACIDS = standardResidueCcpCodes['protein']
SOURCE_NAMES = ['RefDB','BMRB']
ATOM_TYPES = ['Hydrogen','Heavy']

class BrowseReferenceShiftsPopup(BasePopup):
  """
  **Graphs and Charts of Database Chemical Shift Data**
  
  This popup window presents the user with the distributions of known chemical
  shift values, for various kinds of atom, as they appear in the BioMagResBank
  or re-referenced RefDB databases. This information is useful when attempting
  to determine the type of a residue from its NMR spectra, for example when
  performing protein sequence assignment.

  The popup is divided into two tabs for two different kinds of graphical
  display. The first "1D graphs" tab shows the distributions of chemical shift
  value for either hydrogens or other "heavy" atoms within a given kind of
  residue. The buttons that carry atom names above the main graph can be toggled
  two switch the display for different kinds of atom on or off.

  The second "Amino Acid CA CB" tab is a special display to assist in the 
  assignment of protein sequences when using triple-resonance experiments
  like HNCA/CB. Data for all of the common amino acids is displayed, but
  only for the alpha and beta carbons.

  **Caveats & Tips**

  Two dimensional correlations, e.g. between a 1H resonance and a covalently
  bound 13C will be added to Analysis in the future.

  DNA and RNA chemical shift distributions is not present in the RefDB data, but
  are present in the BMRB data; to view the BMRB data change the pulldown menu
  at the top-right.

  Some of the atom types have a second, minor chemical shift distribution, far
  from the main peak, that is erroneous (this data is still present in the
  source databases). An example of this is for the NE atom of the amino acid
  arginine at around 115 ppm, which in this instance is probably caused by peaks
  being assigned in HSQC spectra in the normal backbone amide ppm range when the
  NE peak is really an aliased harmonic and should be a spectrum width away,
  nearer 85 ppm.

  Distributions that are notably jagged, due to lack of data for a given atom
  type should naturally not be relied upon to any great degree.

  """

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    self.sourceName = 'RefDB'
    self.chemAtomNmrRefs = {}
    
    BasePopup.__init__(self, parent=parent, title="Resonance : Reference Chemical Shifts", **kw)

  def body(self, guiFrame):
    
    self.geometry("775x500")
    guiFrame.expandGrid(0,0)
    
    options = ['1D Graphs','Amino Acid CA CB']
    tipTexts = ['Simple graphs of how chemical shift values are distributed for a given atom type',
                'Density plots of alpha a & beta carbon shift distributions for common amnio acids']
    tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(0,0), tipTexts=tipTexts)
    frameA, frameB = tabbedFrame.frames
    
    # # # # # # #  1D GRAPHS  # # # # # # # 
    
 
    row = 0
     
    MolType = self.project.metaclass.metaObjFromQualName('ccp.molecule.ChemComp.MolType')
    molTypes = MolType.enumeration
    molTypes.remove('other')
    molTypes.remove('carbohydrate')
   
    self.molType = 'protein'
    
    ccpCodes = self.getCcpCodes(self.molType) or [None,]
    self.ccpCode = 'Ala'
    
    self.atomType = ATOM_TYPES[0]
    
    self.atomNamesDict = {}
    self.chemAtomNmrRefs = self.getCcpCodeData(self.ccpCode, atomType=self.atomType)
    
    
    tipText = 'Which of the common bio-polymer types to show data for'
    self.molTypeLabel     = Label(frameA, text = 'Molecule Type:', grid=(row,0))
    self.molTypePulldown  = PulldownList(frameA, callback=self.changeMolType,
                                         texts=molTypes, grid=(row,1), tipText=tipText)

    tipText = 'Which residue code to show chemical shift distributions for'
    self.ccpCodeLabel     = Label(frameA, text = 'Residue Code:', grid=(row,2))
    self.ccpCodePulldown  = PulldownList(frameA, callback=self.changeCcpCode, texts=ccpCodes,
                                         index=ccpCodes.index(self.ccpCode), 
                                         grid=(row,3), tipText=tipText)

    tipText = 'Whether to show distributions for hydrogen atoms or other atoms'
    self.atomTypeLabel    = Label(frameA, text = 'Atom Type:', grid=(row,4))
    self.atomTypePulldown = PulldownList(frameA, callback=self.changeAtomType,
                                         texts=ATOM_TYPES, tipText=tipText,
                                         grid=(row,5))

    row += 1

    tipText = 'The selection of atom name to display distributions for'
    self.atomSelector = PartitionedSelector(frameA, self.toggleAtom,
                                            tipText=tipText, maxRowObjects=20)
    self.atomSelector.grid(row=row, column=0, columnspan=6, sticky='ew')

    row += 1

    frameA.expandGrid(row,5)
    self.scrolledGraph = ScrolledGraph(frameA, symbolSize=2, reverseX=True, width=650,
                                       height=300, title='Chemical shift distribution',
                                       xLabel='Chemical shift', yLabel='proportion',
                                       motionCallback=self.updateCrosshairs)
    self.scrolledGraph.grid(row=row, column=0, columnspan=6, sticky='nsew')
    
    # # # # # # #  PROTEIN CA CB   # # # # # # # 
    
    frameB.expandGrid(0,0)
    matrix, ppms = self.getCaCbMatrix()
    title = 'Amino Acid CA & CB Chemical Shifts'
    self.cacbMatrix = ScrolledDensityMatrix(frameB, matrix=matrix, boxSize=14,
                                            title=title,
                                            xLabels=ppms, yLabels=AMINO_ACIDS,
                                            borderColor='grey', zoom=1.0,
                                            labelAxes=True, barPlot=False,
                                            doLegend=False, grid=(0,0))
    
    sdm = self.cacbMatrix
    font = sdm.boldFont
    
    x0, y0 = (470,370)
    sdm.canvas.create_rectangle(x0+170,y0-6,x0+182,y0+6,fill='#4040A0',
                                outline=sdm.borderColor, width=1)
    sdm.canvas.create_text(x0+200, y0, text='CA', font=font)
    sdm.canvas.create_rectangle(x0+220,y0-6,x0+232,y0+6,fill='#A04040',
                                outline=sdm.borderColor, width=1)
    sdm.canvas.create_text(x0+250, y0, text='CB', font=font)
    sdm.canvas.create_text(x0, y0, text='13C PPM +/- 0.5', font=font)
    
    # # # # # # #  M A I N   # # # # # # # 
    
    tipText = 'Whether to use chemical shift data from the BMRB (unfiltered) or RefDB sets'
    label = Label(tabbedFrame.sideFrame, text='Source Database:', grid=(0,0), sticky='e')
    index = SOURCE_NAMES.index(self.sourceName)
    self.sourcePulldown = PulldownList(tabbedFrame.sideFrame, self.changeSource,
                                       texts=SOURCE_NAMES, index=index, 
                                       grid=(0,1), sticky='e', tipText=tipText)

    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame, expands=True,
                                           helpUrl=self.help_url, sticky='e',
                                           grid=(0,2))
                                               
    self.waiting = False
    self.updateAfter()

    for func in ('__init__', 'delete'):
      self.registerNotify(self.updateAfter, 'ccp.nmr.NmrReference.NmrReferenceStore', func)

  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)
  
  def changeSource(self, sourceName):
  
    if sourceName is not self.sourceName:
      self.sourceName = sourceName
      self.updateCcpCodes()
      self.updateAfter()
    
  def getAtomColor(self, i,N):
    
    h = i/float(N)
    s = 1
    v = 0.8
    [r,g,b] = hsbToRgb(h,s,v)
   
    return hexRepr(r,g,b)

  def toggleAtom(self, atomName):
  
    if self.atomNamesDict.get(atomName):
      self.atomNamesDict[atomName] = False
    else:
      self.atomNamesDict[atomName] = True

    self.updateAfter()

  def getCcpCodes(self, molType):
  
    sourceName = self.sourceName

    ccpCodes = []
    for nmrRef in self.project.findAllNmrReferenceStores(molType=molType):
      chemCompNmrRef = nmrRef.findFirstChemCompNmrRef(sourceName=sourceName)
      
      if chemCompNmrRef:
        ccpCodes.append(nmrRef.ccpCode)
    
    ccpCodes.sort()
      
    return ccpCodes
  
  def getCcpCodeData(self, ccpCode, molType=None, atomType=None):

    if not molType:
      molType = self.molType

    dataDict = {}
    project = self.project
    sourceName = self.sourceName
    
    nmrRefStore = project.findFirstNmrReferenceStore(molType=molType,ccpCode=ccpCode)
     
    chemCompNmrRef = nmrRefStore.findFirstChemCompNmrRef(sourceName=sourceName)
    if chemCompNmrRef:
    
      chemCompVarNmrRef = chemCompNmrRef.findFirstChemCompVarNmrRef(linking='any',descriptor='any')
      if chemCompVarNmrRef:
        for chemAtomNmrRef in chemCompVarNmrRef.chemAtomNmrRefs:
          atomName = chemAtomNmrRef.name
          element  = chemAtomNmrRef.findFirstChemAtom().elementSymbol
  
          if not atomType:
            dataDict[atomName] = chemAtomNmrRef
  
          elif (atomType == 'Hydrogen' and element == 'H') or \
               (atomType == 'Heavy' and element != 'H'):
            dataDict[atomName] = chemAtomNmrRef
 
    return dataDict  
 
  def changeMolType(self, molType):
  
    self.molType = molType
    self.updateCcpCodes()
  
  def changeCcpCode(self, ccpCode):
  
    self.ccpCode = ccpCode
    self.updateAtomNames()
       
  def changeAtomType(self, atomType):
  
    self.atomType = atomType
    self.updateAtomNames()
       
  def updateAfter(self, *opt):
  
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)

  def updateCcpCodes(self):
    
    ccpCodes = self.getCcpCodes(self.molType)
    if ccpCodes:
      self.ccpCodePulldown.setup(ccpCodes, ccpCodes, 0)
      self.changeCcpCode(ccpCodes[0])
    else:
      self.ccpCodePulldown.setup([], [], 0)
      self.changeCcpCode(None)
      
  
  def updateAtomNames(self):
  
    if self.ccpCode:
      self.chemAtomNmrRefs = self.getCcpCodeData(self.ccpCode, atomType=self.atomType)
      atomNames = self.chemAtomNmrRefs.keys()
      atomNames = greekSortAtomNames(atomNames, self.molType)
      N = len(atomNames)

      for atomName in atomNames:
        if self.atomNamesDict.get(atomName) is None:
          self.atomNamesDict[atomName] = True

      colors = []
      for i in range(N):
        colors.append( self.getAtomColor(i,N) )

 
    else:
      colors = []
      atomNames = []
      N = 0
 
    self.atomSelector.update(objects=atomNames,labels=atomNames,colors=colors)
        
    for i in range(N):
      if self.atomNamesDict[atomNames[i]]:
        self.atomSelector.setButtonState(i, True)
      else:
        self.atomSelector.setButtonState(i, False)

    self.updateAfter()
  
  def destroy(self):

    for func in ('__init__', 'delete'):
      self.unregisterNotify(self.updateAfter, 'ccp.nmr.NmrReference.NmrReferenceStore', func)

    BasePopup.destroy(self)
    
  def update(self):

    matrix, ppms = self.getCaCbMatrix()
    self.cacbMatrix.xLabels = ppms
    self.cacbMatrix.update(matrix)
       
    self.updateAtomNames()
    dataSets = []
    colorList = []
    atomsString = ''
    
    if self.ccpCode:
      c = 0
      atomNames = self.chemAtomNmrRefs.keys()
      N = len(atomNames)
      for atomName in greekSortAtomNames(atomNames, self.molType):
        if self.atomNamesDict[atomName]:
          atomsString += '%s ' % atomName
          chemAtomNmrRef = self.chemAtomNmrRefs[atomName]
          distribution   = chemAtomNmrRef.distribution
          refPoint       = chemAtomNmrRef.refPoint
          refValue       = chemAtomNmrRef.refValue
          valuePerPoint  = chemAtomNmrRef.valuePerPoint

          data = []
          for i in range(len(distribution)):
            x = refValue + valuePerPoint*( i-refPoint)
            y = distribution[i]
            data.append( (x,y) )
 
          dataSets.append(data)
          colorList.append(self.getAtomColor(c,N))
        c += 1
 
      self.scrolledGraph.title = '%s %s Chemical shift distribution' % (self.ccpCode,atomsString)
    else:
      dataSets = [[(0,0)],]
      self.scrolledGraph.title = 'Chemical shift distribution'
    
    self.scrolledGraph.zoom  = 1.0
    
    self.scrolledGraph.update(dataSets=dataSets, dataColors=colorList)  
    
    self.waiting = False
    
  def getCaCbMatrix(self):  
    
    ppms = range(73,14,-1)
    blankCol = [0.0]*len(ppms)
    matrix = [[0.0]*len(AMINO_ACIDS) for x in ppms]
    
    for i, ccpCode in enumerate(AMINO_ACIDS):
      atomRefDict = self.getCcpCodeData(ccpCode, molType='protein')
      valueDict = {}
      
      for atomName in ('CA','CB'):
        valueDict[atomName] = blankCol[:]
        
        chemAtomNmrRef = atomRefDict.get(atomName)
        if not chemAtomNmrRef:
          continue
          
        distribution   = chemAtomNmrRef.distribution
        refPoint       = chemAtomNmrRef.refPoint
        refValue       = chemAtomNmrRef.refValue
        valuePerPoint  = chemAtomNmrRef.valuePerPoint
        
        
        for j, ppm in enumerate(ppms):
          ppmMin = ppm-0.5
          ppmMax = ppm+0.5
 
          v = 0.0
          n = 0.0
          for k in range(len(distribution)):
            y = refValue + valuePerPoint*(k-refPoint)
 
            if ppmMin < y <= ppmMax:
              v += distribution[k]
              n += 1.0
            elif y > ppmMax:
              break
 
          if n:
            valueDict[atomName][j] = v/n
      
      for j, ppm in enumerate(ppms):
        matrix[j][i] = valueDict['CA'][j] - valueDict['CB'][j]

    return matrix, ppms

  def updateCrosshairs(self, event):

    # code below is a bit dangerous

    position = self.scrolledGraph.getPlotCoords(event)[0]

    typeLocation = []
    isHydrogen = (self.atomType == 'Hydrogen')
    for panelType in self.analysisProject.panelTypes:
      axisType = panelType.axisType
      if isHydrogen:
        if axisType.name in ('1H', '2H'):
          typeLocation.append((panelType, position))
      else:
        if axisType.name not in ('1H', '2H'):
          typeLocation.append((panelType, position))

    self.scrolledGraph.mouseEnter(event)
    self.parent.drawWindowCrosshairs(typeLocation)

  def drawCrosshairs(self, typeLocation):
    """ Draw crosshairs at specified locations.
        typeLocation = tuple of (PanelType, position)
    """

    isHydrogen = (self.atomType == ATOM_TYPES[0])

    positions = []
    for (panelType, position) in typeLocation:
      axisType = panelType.axisType
      if axisType.measurementType != 'Shift':
        continue
        
      if position is None:
        continue

      if isHydrogen:
        if axisType.name in ('1H', '2H'):
          positions.append(position)
      else:
        if axisType.name not in ('1H', '2H'):
          positions.append(position)

    self.scrolledGraph.drawVerticalLines(positions)

