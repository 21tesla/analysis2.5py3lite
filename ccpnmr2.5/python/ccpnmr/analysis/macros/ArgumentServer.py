
"""
======================COPYRIGHT/LICENSE START==========================

ArgumentServer.py: Part of the CcpNmr Analysis program

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

import re

from memops.general.Implementation import ApiError

from ccp.general.ArgumentServer import ArgumentServer as GenArgumentServer
from ccpnmr.analysis.core.AssignmentBasic import makeResonanceGuiName
from ccpnmr.analysis.core.MoleculeBasic   import getResidueCode
from ccpnmr.analysis.core.Util            import getAnalysisProfile, getAnalysisProject

class ArgumentServer(GenArgumentServer):
 
  def getCurrentPeaks(self):

    peaks  = list(self.parent.currentPeaks)
    if not peaks:
      self.messageReporter.showWarning('Warning','No peaks currently selected')
      if not self.inGui:
        print 'Need to define %s.currentPeaks\n' % self.parent
    
    return peaks
  
  def getCurrentStrip(self):
        
    return self.parent.currentStrip

  def getCurrentPeak(self):
        
    return self.parent.currentPeak
  
  def getCurrentWindow(self):
  
    return self.parent.currentWindow
    
  def getCurrentWindowPane(self):
  
    return self.parent.currentWindowPane

  def getCurrentWindowPopup(self):
  
    return self.parent.currentWindowPopup

  def getCurrentSpectra(self):
  
    return list(self.parent.currentSpectra) or []

  def getCurrentCanvas(self):

    return self.parent.currentCanvas

  def getCurrentEvent(self):

    return self.parent.currentEvent

  def getCurrentPeakLists(self):
  
    return list(self.parent.currentPeakLists) or []
   
  def getCurrentPosition(self):
  
    if self.parent.currentPosition:
      return self.parent.currentPosition

    popup = self.getCurrentWindowPopup()
    canvas = self.getCurrentCanvas()
    event = self.getCurrentEvent()
  
    if popup and canvas and event:
      return popup.calcWorldCoord(canvas, event.x, event.y)[:2]

    return None

  def getCurrentSpectrum(self):
  
    peak = self.parent.currentPeak
    if peak:
      return peak.peakList.dataSource
      
  def getCurrentPeakDim(self):
  
    peakDim = self.parent.currentPeakDim
  
    if not peakDim:
  
      peak = self.getCurrentPeak()  
  
      if not peak:
        self.messageReporter.showWarning('Warning','No peak or peak dim currently selected')
        return
     
      peakDim = self.chooseObject(peak.sortedPeakDims())
  
    return peakDim

  def getAnalysisProject(self):

    return getAnalysisProject(self.getProject())

  def getAnalysisProfile(self):

    return getAnalysisProfile(self.getProject())

  def getResonanceGuiName(self, resonance, fullName=True):
    """ Analysis resonance name - overrides default version
    """
    return makeResonanceGuiName(resonance, fullName)

  def getWindow(self):

    analProject = self.getProject().currentAnalysisProject
    return self.chooseObject(analProject.sortedSpectrumWindows(),
                             key='name', objectName='Window')

  def getColorScheme(self):
  
    project = self.getProject()
    return self.chooseObject(project.currentAnalysisProfile.sortedColorSchemes())
  
  def getResidueCode(self, obj):
    """Analysis-specific version of getResidueCode - overrides generic version
    """
    return getResidueCode(obj)
