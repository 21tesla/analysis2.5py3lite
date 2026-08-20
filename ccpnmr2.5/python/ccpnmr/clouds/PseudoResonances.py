"""
======================COPYRIGHT/LICENSE START==========================

PseudoResonances.py: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
class PseudoMolSystem:

  def __init__(self):
    self.pseudoResons  = []
    self.pseudoSpinSysts = []
        
  def delete(self):
    
    for pseudoReson in self.pseudoResons:
      pseudoReson.delete()
    
    for pseudoSpinSyst in self.pseudoSpinSysts:
      pseudoSpinSyst.delete()
      
    del self

  def mergePseudoResons(self,res1,res2):

    if (res1.parent is not self) or (res2.parent is not self):
      print "PseudoResonance not in PseudoMolSystem"
      return

    for peakDim in res2.peakDims:
      res1.addPeakDim(peakDim)
      res1.recalculatePpm()
 
    ss1 = res1.pseudoSpinSyst
    ss2 = res2.pseudoSpinSyst
    if ss2:
      if ss1:
        if ss1 is not ss2:
          mergePseudoSpinSysts(ss1,ss2)
      else:
        res1.setPseudoSpinSyst(ss2)
 
      ss2.pseudoResons.remove(res2)
 
    res2.delete()
    return res1

  def mergePseudoSpinSysts(self,ss1,ss2,name):

    if (ss1.parent is not self) or (ss2.parent is not self):
      print "PseudoSpinSyst not in PseudoMolSystem"
      return

    for pseudoReson in ss2.pseudoResons:
      ss1.addPseudoReson(pseudoReson)
 
    ss2.delete()
    ss1.name = name
    return ss1

class PseudoSpinSyst:

  def __init__(self,parent,name='anonymous',pseudoResons=None):
    parent.pseudoSpinSysts.append(self)
    self.parent = parent
    self.pesudoMolSystem = parent
    self.name = name

    if pseudoResons is None:
      self.pseudoResons = []
    else:
      self.setPseudoResons(pseudoResons)
      
  def getPeaks(self):
  
    peaks = []
    for pseudoReson in self.pseudoResons:
      for peakDim in pseudoReson.peakDims:
        peak = peakDim.peak
        if peak not in peaks:
          peaks.append(peak)
          
    return peaks
  
  def getPseudoResonByName(self,name):
  
    for pseudoReson in self.pseudoResons:
      if pseudoReson.name == name:
        return pseudoReson
    
    return None
    
  def setPseudoResons(self,pseudoResons=None):
  
    self.pseudoResons = []
    if pseudoResons is not None:
      for pseudoReson in pseudoResons:
        if pseudoReson.parent is self.parent:
          pseudoReson.pseudoSpinSyst = self
          self.pseudoResons.append(pseudoReson)
     
  def addPseudoReson(self,pseudoReson):
  
    if pseudoReson.parent is self.parent:
      if pseudoReson not in self.pseudoResons:
        self.pseudoResons.append(pseudoReson)
        pseudoReson.pseudoSpinSyst = self
  
  def removePseudoReson(self,pseudoReson):
  
    if pseudoReson in self.pseudoResons:
      self.pseudoReson.remove(pseudoReson)
      pseudoReson.pseudoSpinSyst = None

  def deleteAll(self):
    for pseudoReson in self.pseudoResons:
      pseudoReson.delete()
    
    self.parent.pseudoSpinSysts.remove(self)
    del self
  
  def delete(self):
    for pseudoReson in self.pseudoResons:
      #pseudoReson.delete()
      pseudoReson.pseudoSpinSyst = None

    self.parent.pseudoSpinSysts.remove(self)
    del self

class PseudoReson:

  def __init__(self,parent,name,peakDims=None):
    parent.pseudoResons.append(self)
    self.parent = parent
    self.pesudoMolSystem = parent
    if peakDims is None:
      peakDims = []
    
    self.pseudoSpinSyst = None
    self.name     = name  
    self.peakDims = peakDims
    self.ppm      = 0.0
    if self.peakDims:
      self.recalculatePpm()

  def setPseudoSpinSyst(self,pseudoSpinSyst=None):
  
    if pseudoSpinSyst and ( pseudoSpinSyst.parent is not self.parent):
      print "Pseudo resonance has different parent to Pseudo spin system in PseudoReson.setPseudoSpinSyst"
      return

    if self.pseudoSpinSyst is pseudoSpinSyst:
      return

    if self.pseudoSpinSyst:
      self.pseudoSpinSyst.removePseudoReson(self)

    if pseudoSpinSyst:
      pseudoSpinSyst.addResonance(self)
    
    self.pseudoSpinSyst = pseudoSpinSyst
  
  def getPeaks(self):
  
    peaks = set()
    for peakDim in self.peakDims:
      peaks.add(peakDim.peak)
    
    return list(peaks)
  
  def recalculatePpm(self):
  
    n   = 0
    sum = 0.0
    for peakDim in self.peakDims:
      sum += getPeakDimPpm(peakDim)
      
    if n > 0:
      self.ppm = sum/float(n)
      
    return self.ppm
  
  def removePeakDim(self,peakDim):
  
    if peakDim in self.peakDims:
      self.peakDims.remove(peakDim)
      self.recalculatePpm
  
  def addPeakDim(self,peakDim):
    
    if peakDim not in self.peakDims:
      self.peakDims.append( peakDim )
      self.recalculatePpm
  
  def delete(self):
    if self.pseudoSpinSyst:
      self.pseudoSpinSyst.removePseudoReson(self)
      
    self.parent.pseudoResons.remove(self)

    del self
