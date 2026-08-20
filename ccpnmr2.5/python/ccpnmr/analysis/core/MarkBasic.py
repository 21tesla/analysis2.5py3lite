
LICENSE = """
======================COPYRIGHT/LICENSE START==========================

MarkBasic.py: Part of the CcpNmr Analysis program

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
try:
  import memops.gui.Color as Color
except ImportError:
  pass # could lead to exception when using code below, but relevant code should not be in gui in the first place

from memops.general import Implementation

from ccpnmr.analysis.core.Util        import getPeakListColor
from ccpnmr.analysis.core.WindowBasic import getSpectrumWindowView

markCounter = 0
rulerCounter = 0

def removeAllMarksAndRulers(argServer):
  """
  Macro to remove all marks and rulers in a project.

  .. describe:: Input
  
  ArgumentServer 
  
  .. describe:: Output
  
  None
  """

  project = argServer.getProject()

  removeRulers(project, removeAll=True)
  removeMarks(project,  removeAll=True)
  
  
def removeMarks(project, removeAll=False):
  """
  Remove marks in a project. If removeAll = True then remove all marks.
  Otherwise remove only marks above the user-specified limit maxMarks.
  
  .. describe:: Input

  Project, Boolean
  
  .. describe:: Output

  None
  """
 
  analysisProject = project.currentAnalysisProject

  if removeAll:
    maxMarks = 0
  else:
    maxMarks = analysisProject.maxMarks
 
  marks = analysisProject.sortedMarks()
 
  m = len(marks) - maxMarks
  for i in range(m):
    marks[i].delete()
 
 
def createPeakDimRuler(peakDim, windowPane, lineWidth=1, dashLength=3,
                       gapLength=1, color = None, remove=True):
  """
  Create an Analysis ruler at a peakDim position given a window
  
  .. describe:: Input
  
  Nmr.Peak, Analysis.SpectrumWindowPane, Int, Int, Int,
  Analysis.Color, Boolean
  
  .. describe:: Output

  Analysis.Ruler
  """
 
  global rulerCounter

  project  = peakDim.root
  analysisProject = peakDim.topObject
  spectrum = peakDim.peak.peakList.dataSource
  
  view = getSpectrumWindowView(windowPane, spectrum)
  if not view:
    return
  
  dataDimRef = peakDim.dataDimRef
  if not dataDimRef: # Sampled
    return
    
  analysisDataDim = dataDimRef.dataDim.analysisDataDim
   
  axisMapping = view.findFirstAxisMapping(analysisDataDim=analysisDataDim)
  if not axisMapping:
    return

  axisPanel = windowPane.findFirstAxisPanel(label=axisMapping.label)
  panelType = axisPanel.panelType

  if not color:
    color = Color.grey.hex
 
  ruler = analysisProject.newRuler(panelType=panelType, position=peakDim.value,
                                   lineWidth=lineWidth, dashLength=dashLength,
                                   gapLength=gapLength, color=color)
 
  if remove:
    removeRulers(project)

  return ruler


def createPeakMark(peak, lineWidth=1, dashLength=2, gapLength=2, remove=True, axisTypeDict=None):
  """
  Create a mark positioned at a given peak

  .. describe:: Input
  
  Nmr.Peak, Int, Int, Int, Boolean

  .. describe:: Output
  
  Analysis.Mark
  """

  if not axisTypeDict:
    axisTypeDict = {}

  peakList = peak.peakList
  project = peakList.root
  analysisProject = project.currentAnalysisProject
  
  analysisPeakList = peakList.analysisPeakList
  if analysisPeakList:
    color = analysisPeakList.symbolColor
  elif project.currentAnalysisProfile:
    color = project.currentAnalysisProfile.fgColor
  else:
    color = Color.grey.hex

  mark = analysisProject.newMark(lineWidth=lineWidth, dashLength=dashLength,
                                 gapLength=gapLength, color=color)

  for peakDim in peak.peakDims:
    axisType = axisTypeDict.get(peakDim)
    if not axisType:
      dataDimRef = peakDim.dataDimRef
  
      if dataDimRef:
        expDimRef = dataDimRef.expDimRef
        isotopeCodes = expDimRef.isotopeCodes
        axisType = analysisProject.findFirstAxisType(isotopeCodes=isotopeCodes)
      
        if not axisType:
          msg = 'Unknown axis type isotope codes ' + isotopeCodes
          raise Implementation.ApiError(msg)
        
    if axisType:
      markDim = mark.newMarkDim(position=peakDim.value, axisType=axisType)
 
  mark.peak = peak

  if remove:
    removeMarks(project)

  return mark


def createNonPeakMark(positions, axisTypes, lineWidth=1,
                      dashLength=2, gapLength=2, color=None, remove=True):
  """
  Create a mark at a general (i.e. non-peak) position.
  
  .. describe:: Input

  List of Floats (Analysis.Mark.positions), List of Analysis.AxisTypes,
  Int, Int, Int, Word (Analysis.Color.name)
  
  .. describe:: Output

  Analysis.Mark
  """
 
  global markCounter

  analysisProject = list(axisTypes)[0].topObject
  project = analysisProject.root

  if not color:
    analysisProfile = project.currentAnalysisProfile
    
    if analysisProfile:
      scheme = analysisProfile.marksColor
  
      if scheme:
        colors = scheme.colors
        n = len(colors)
        color = colors[markCounter % n]
        markCounter += 1
      else:
        color = analysisProfile.fgColor
    
    else:
      color = Color.red.hex      
 
  mark = analysisProject.newMark(lineWidth=lineWidth, dashLength=dashLength,
                                 gapLength=gapLength, color=color)

  for n in range(len(positions)):
    axisType = axisTypes[n]
    
    if axisType.name != 'value':
      markDim = mark.newMarkDim(position=positions[n], axisType=axisType)
 
  mark.peak = None

  if remove:
    removeMarks(project)

  return mark


def setPeakMarkColor(peakList):
  """
  Set the colour of marks positioned on peaks.
  
  .. describe:: Input
  
  Nmr.PeakList
  
  .. describe:: Output
  
  None
  """

  project = peakList.root
  analysisProject = project.currentAnalysisProject

  
  for mark in analysisProject.marks:
    if (hasattr(mark, 'peak') and mark.peak and (mark.peak.peakList == peakList)):
      analysisPeakList = peakList.analysisPeakList
      
      if analysisPeakList:
        mark.color = analysisPeakList.symbolColor
      

def setNonPeakMarkColor(project):
  """
  Set the colour of general marks (i.e. not positioned on peaks).
  
  .. describe:: Input
  
  Project, Word (Analysis.Color.name)
  
  .. describe:: Output
  
  None
  """

  global markCounter

  analysisProject = project.currentAnalysisProject
  analysisProfile = project.currentAnalysisProfile
    
  if analysisProfile:
    scheme = analysisProfile.marksColor
  
    if scheme:
      colors = scheme.colors
      n = len(colors)
      
      for i, mark in enumerate(analysisProject.marks):
        if not hasattr(mark, 'peak') or not mark.peak:
          mark.color = colors[i % n]
          
        markCounter = i
  
  
def removeRulers(project, removeAll = False):
  """
  Remove rulers in a project. If removeAll = True then remove all rulers.
  Otherwise remove only rulers above the user-specified limit maxRulers.
  
  .. describe:: Input
  
  Project, Boolean
  
  .. describe:: Output
  
  None
  """
 
  analysisProject = project.currentAnalysisProject

  if (removeAll):
    maxRulers = 0
  else:
    maxRulers = analysisProject.maxRulers

  rulers = analysisProject.sortedRulers()
 
  m = len(rulers) - maxRulers
  for i in range(m):
    rulers[i].delete()
 
 
def createRuler(position, panelType, lineWidth=1, dashLength=4,
                gapLength=2, color = None, remove=True):
  """
  Create an Analysis ruler
  
  .. describe:: Input
  
  Float (Analysis.Ruler.position), Analysis.PanelType,
  Int, Int, Int, Word (Analysis.Color.name)
  
  .. describe:: Output
  
  Analysis.Ruler
  """
 
  global rulerCounter

  project = panelType.root
  analysisProject = panelType.topObject

  if not color:
    analysisProfile = project.currentAnalysisProfile
    
    if analysisProfile:
      scheme = analysisProfile.rulersColor
  
      if scheme:
        colors = scheme.colors
        n = len(colors)
        color = colors[rulerCounter % n]
        rulerCounter += 1
      else:
        color = analysisProfile.fgColor
    
    else:
      color = Color.blue.hex 
 
  ruler = analysisProject.newRuler(panelType=panelType, position=position,
                                   lineWidth=lineWidth, dashLength=dashLength,
                                   gapLength=gapLength, color=color)
 
  if remove:
    removeRulers(project)

  return ruler


def setRulerColor(project):
  """
  Set the colour of rulers.
  
  .. describe:: Input
  
  Project, Word (Analysis.Color.name)
  
  .. describe:: Output
  
  None
  """

  global rulerCounter

  analysisProject = project.currentAnalysisProject
  analysisProfile = project.currentAnalysisProfile
    
  if analysisProfile:
    scheme = analysisProfile.marksColor
  
    if scheme:
      colors = scheme.colors
      n = len(colors)
      
      for i, ruler in enumerate(analysisProject.rulers):
        ruler.color = colors[i % n]
          
        rulerCounter = i

