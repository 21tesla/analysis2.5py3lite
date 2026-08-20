
"""
======================COPYRIGHT/LICENSE START==========================

UnitConverter.py: Part of the CcpNmr Analysis program

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
def ppm2pnt(ppm, dataDimRef):

  freqDataDim = dataDimRef.dataDim

  npoints = freqDataDim.numPointsOrig
  sw = freqDataDim.spectralWidthOrig
  sf = dataDimRef.expDimRef.sf
  refpt = dataDimRef.refPoint
  refppm = dataDimRef.refValue

  t = - npoints * sf / float(sw)
  pnt = t*(ppm - refppm) + refpt

  return pnt

def pnt2ppm(pnt, dataDimRef):

  freqDataDim = dataDimRef.dataDim

  npoints = freqDataDim.numPointsOrig
  sw = freqDataDim.spectralWidthOrig
  sf = dataDimRef.expDimRef.sf
  refpt = dataDimRef.refPoint
  refppm = dataDimRef.refValue

  t = - npoints * sf / float(sw)
  ppm = (pnt - refpt)/t + refppm

  return ppm

def hz2pnt(hz, dataDimRef):

  freqDataDim = dataDimRef.dataDim

  npoints = freqDataDim.numPointsOrig
  sw = freqDataDim.spectralWidthOrig
  sf = dataDimRef.expDimRef.sf
  refpt = dataDimRef.refPoint
  refppm = dataDimRef.refValue

  t = - npoints / float(sw)
  pnt = t*(hz - sf*refppm) + refpt

  return pnt

def pnt2hz(pnt, dataDimRef):

  freqDataDim = dataDimRef.dataDim

  npoints = freqDataDim.numPointsOrig
  sw = freqDataDim.spectralWidthOrig
  sf = dataDimRef.expDimRef.sf
  refpt = dataDimRef.refPoint
  refppm = dataDimRef.refValue

  t = - npoints / float(sw)
  hz = (pnt - refpt)/t + sf*refppm

  return hz

unit_converter = {}
unit_converter[('ppm', 'point')] = ppm2pnt
unit_converter[('point', 'ppm')] = pnt2ppm
unit_converter[('Hz', 'point')] = hz2pnt
unit_converter[('point', 'Hz')] = pnt2hz
