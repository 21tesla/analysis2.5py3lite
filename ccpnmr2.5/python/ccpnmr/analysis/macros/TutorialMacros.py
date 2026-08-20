
"""
======================COPYRIGHT/LICENSE START==========================

TutorialMacros.py: Part of the CcpNmr Analysis program

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
#False = 0

def addMarksToPeaks(argServer, peaks=None):
  """Descrn: Adds position line markers to the selected peaks.
     Inputs: ArgumentServer, List of Nmr.Peaks
     Output: None
  """
  from ccpnmr.analysis.core.MarkBasic import createPeakMark
  
  if not peaks:
    peaks = argServer.getCurrentPeaks()

  # no peaks - nothing happens
  for peak in peaks:
    createPeakMark(peak, remove=0)

def calcAveragePeakListIntensity(argServer, peakList=None, intensityType='height'):
  """Descrn: Find the average height of peaks in a peak list.
     Inputs: ArgumentServer, Nmr.PeakList
     Output: Float
  """
  from ccpnmr.analysis.core.ConstraintBasic import getMeanPeakIntensity
  
  if not peakList:
    peakList = argServer.getPeakList()

  if not peakList:
    argServer.showWarning('No peak list selected')
    return
  
  answer = argServer.askYesNo('Use peak volumes? Height will be used otherwise.')

  if answer: # is true
    intensityType = 'volume'

  spectrum   = peakList.dataSource
  experiment = spectrum.experiment

  intensity = getMeanPeakIntensity(peakList.peaks, intensityType=intensityType)
  
  argServer.showInfo('Mean peak %s for %s %s peak list %d is %e' % (intensityType,experiment.name,spectrum.name,peakList.serial,intensity))
  
  return intensity


def openMyPopup(argServer):
  """Descrn: Opens and example popup.
     Inputs: ArgumentServer
     Output: None
  """

  peakList = argServer.getPeakList()
  popup = MyPopup(argServer.parent, peakList)


from memops.gui.BasePopup import BasePopup
from memops.gui.ButtonList import ButtonList
from memops.gui.ScrolledGraph import ScrolledGraph
from ccpnmr.analysis.core.PeakBasic import getPeakHeight, getPeakVolume

class MyPopup(BasePopup):

  def __init__(self, parent, peakList, *args, **kw):
     
    self.peakList  = peakList
    self.colours   = ['red', 'green']
    self.dataSets  = []
     
    BasePopup.__init__(self, parent=parent,  title='Test Popup', **kw)

  def body(self, guiParent):

     row = 0
     self.graph = ScrolledGraph(guiParent)
     self.graph.grid(row=row, column=0, sticky='NSEW')
     
     row += 1
     texts    = ['Draw graph','Goodbye']
     commands = [self.draw, self.destroy]
     buttons  = ButtonList(guiParent, texts=texts, commands = commands)
     buttons.grid(row=row, column=0, sticky='NSEW')

  def draw(self):
    
    self.dataSets  = self.getData()
    self.graph.update(self.dataSets, self.colours)

  def getData(self):

    peakData = [( getPeakVolume(peak) or 0.0, peak) for peak in self.peakList.peaks]
    peakData.sort()
    
    heights = []
    volumes = []
    i = 0
    for volume, peak in peakData:
      heights.append([i, getPeakHeight(peak) or 0.0])
      volumes.append([i, volume])
      i += 1
  
    return [heights, volumes]
