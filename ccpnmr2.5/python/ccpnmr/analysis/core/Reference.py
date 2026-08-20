
"""
======================COPYRIGHT/LICENSE START==========================

Reference.py: Part of the CcpNmr Analysis program

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

class Reference:

  def __init__(self, npoints, sw, sf, refpt, refppm, isotopeCodes):

    self.npoints = npoints
    self.sw = sw
    self.sf = sf
    self.refpt = refpt
    self.refppm = refppm
    self.isotopeCodes = isotopeCodes

  def pnt2ppm(self, pnt):

    t = - self.npoints * self.sf / float(self.sw)
    ppm = (pnt - self.refpt)/t + self.refppm

    return ppm

  def ppm2pnt(self, ppm):

    t = - self.npoints * self.sf / float(self.sw)
    pnt = t*(ppm - self.refppm) + self.refpt

    return pnt

def getReference(dataDimRef, first = None, last = None):
  """ Return a reference object which contains the reference information
      for the dataDimRef.  Useful mainly to store old reference information.
      The parameters first and last are used if you want to restrict the
      original region.  They are specified in points (counting from 0),
      and first <= x < last.
  """

  freqDataDim = dataDimRef.parent
  npoints = freqDataDim.numPointsOrig
  sw = freqDataDim.spectralWidthOrig
  sf = dataDimRef.expDimRef.sf
  refpt = dataDimRef.refPoint
  refppm = dataDimRef.refValue
  isotopeCodes = dataDimRef.expDimRef.isotopeCodes

  if first is not None or last is not None:
    if first is None:
      first = 0
    else:
      if first < 0 or first >= npoints:
        raise Exception('first = %d, must be in range 0 to %d' % (first, npoints-1))

    if (last is None):
      last = npoints
    else:
      if last < 1 or last > npoints:
        raise Exception('last = %d, must be in range 1 to %d' % (last, npoints))
        
    if first >= last:
      raise Exception('first = %d, must be less than last = %d' % (first, last))

    n = last - first
    sw = (n * sw) / npoints
    refpt = refpt - first
    npoints = n

  return Reference(npoints, sw, sf, refpt, refppm, isotopeCodes)

def shiftDataDimRef(dataDimRef, oldReference):
  """ Shift all peaks in dim given by dataDimRef so that
      ppm value remains the same as given by oldReference
  """

  freqDataDim = dataDimRef.parent
  dim = freqDataDim.dim
  spectrum = freqDataDim.parent
  peakLists = spectrum.peakLists
  
  for peakList in peakLists:
    peaks = peakList.peaks
    
    for peak in peaks:
      peakDim = peak.findFirstPeakDim(dim=dim)
      
      if (peakDim): # should be true
        shiftPeakDim(peakDim, oldReference)

def shiftPeakDim(peakDim, oldReference):
  """ Shift peakDim so that ppm value remains the same as
      given by oldReference
  """

  dataDimRef = peakDim.dataDimRef
  if not dataDimRef:
    return
    
  """ 10 Feb 10: add in numAliasing correction """
  pnt = peakDim.position + peakDim.numAliasing*dataDimRef.dataDim.numPointsOrig
  ppm = oldReference.pnt2ppm(pnt)
  """ 30 Sep 2009: below can screw up if pnt is outside permitted range
  newReference = getReference(dataDimRef)
  pnt = newReference.ppm2pnt(ppm)
  peakDim.position = pnt
"""
  # below automatically figures out correct numAliasing
  peakDim.value = ppm
