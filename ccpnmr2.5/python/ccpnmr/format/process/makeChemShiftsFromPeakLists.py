#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

makeChemShiftsFromPeakLists.py: Makes a chemical shift list from peak lists in the data model

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

import sys

from memops.universal.Util import getStandardDev

from ccpnmr.format.general.Util import getResNameText
from ccpnmr.format.general.userInteraction import setupMessageReporter
from ccpnmr.format.general.userInteraction import setupMultiDialog

import ccp.api.nmr.Nmr as Nmr

#
# Create shifts based on experiments
#

def createChemShifts(peakLists,guiParent = None, messageReporter = None, multiDialog = None, shiftList = None, useAllContribs = 1, defaultShiftError = None, shiftListName = 'None'):
  
  #
  # Set up GUI or non-GuI user interaction if not passed in
  #
    
  if not messageReporter:
    messageReporter = setupMessageReporter(guiParent)

  if not multiDialog:
    multiDialog = setupMultiDialog(guiParent)

  #
  # Initial checks 
  #

  if not peakLists:
    messageReporter.showError("Error","No peaklist(s) provided.")
    return False

  for peakList in peakLists:
    if not isinstance(peakList,Nmr.PeakList):
      messageReporter.showError("Error","Invalid peaklist provided.")
      return False
  
  #
  # Initialize some stuff
  #
  
  if not defaultShiftError:
    defaultShiftError = {}
    
  for (nucl,defValue) in (('1H',0.002),('13C',0.1),('15N',0.1)):
    if not defaultShiftError.has_key(nucl):
      defaultShiftError[nucl] = defValue

  proj = peakLists[0].root

  resonanceChemShifts = {}
  
  experiments = []
  
  #
  # Read the chemshifts and link them to the relevant resonances
  #
  
  for peakList in peakLists:
  
    if peakList.dataSource.experiment not in experiments:
      experiments.append(peakList.dataSource.experiment)
  
    for peak in peakList.sortedPeaks():
    
      for peakDim in peak.sortedPeakDims():
        
        chemShift = peakDim.value
      
        peakDimContribs = tuple(peakDim.peakDimContribs)
        if useAllContribs == 0 and len(peakDimContribs) > 1:
          
          #
          # First check if resonances for contribs connected through one resonanceSet
          # If so, then set values anyway.
          #
          
          resonanceSet = peakDimContribs[0].resonance.resonanceSet

          sameShiftGroup = 1
          
          if resonanceSet:
          
            resonanceSetResonances = resonanceSet.sortedResonances()
          
            for peakDimContrib in peakDimContribs[1:]:
            
              if peakDimContrib.resonance not in resonanceSetResonances:
                sameShiftGroup = 0
                break
            
          if not sameShiftGroup:    
          
            #
            # Ignore this peakDim - many contribs possible and not wanted (?)
            #
          
            print "  Warning: ignoring peakDim %d for peak %d... more than one resonance contribution" % (peakDim.dim,peak.serial)
            continue
          
        for peakDimContrib in peakDimContribs:
        
          resonance = peakDimContrib.resonance
          
          if not resonanceChemShifts.has_key(resonance):
            resonanceChemShifts[resonance] = []
            
          resonanceChemShifts[resonance].append(chemShift)
   
  #
  # Make a new chemical shift list if none passed in
  #
          
  if shiftList == None:
    shiftList = Nmr.ShiftList(proj.currentNmrProject, name = shiftListName, unit = 'ppm', experiments = experiments)

  #
  # Get the average (and standard dev) of the chemical shift for each resonance
  # Create the chemical shifts
  #
   
  for resonance in resonanceChemShifts.keys():
  
    chemShiftTotal = 0.0
    
    for chemShift in resonanceChemShifts[resonance]:
    
      chemShiftTotal += chemShift
    
    chemShiftAverage = chemShiftTotal / len(resonanceChemShifts[resonance])
    
    standardDev = getStandardDev(resonanceChemShifts[resonance],chemShiftTotal)
  
    if standardDev == None:
      if defaultShiftError.has_key(resonance.isotopeCode):
        shiftError = defaultShiftError[resonance.isotopeCode]
      else:
        shiftError = 0.1
    
    else:
      shiftError = standardDev * 2
    
    #
    # Set the chemical shift value
    #
    
    shift = shiftList.findFirstMeasurement(resonance = resonance)
    
    if not shift:
      Nmr.Shift(shiftList, value = chemShiftAverage, error = shiftError, resonance = resonance)
    else:
      allowedErrorMargin = shiftError / 10
      if shift.value + allowedErrorMargin < chemShiftAverage or chemShiftAverage < shift.value - allowedErrorMargin:
        resName = getResNameText(resonance)
        multiDialog.MeasurementSelect(guiParent,resName,resonance,chemShiftAverage,shiftList,'Shift')

  return shiftList
  
###################
# Main of program #
###################

if __name__ == "__main__":

  inFile = sys.argv[1]

  # Read xml file
  from memops.format.xml import XmlIO
  proj = XmlIO.loadProjectFile(inFile)

  peakLists = []
  
  for exp in proj.currentNmrProject.sortedExperiments():
    
    #
    # TODO: need selection here!
    #
    
    for ds in exp.sortedDataSources():
  
      if ds.peakLists:
        peakLists.extend(ds.sortedPeakLists())

  print peakLists  
  shiftList = createChemShifts(peakLists)
  print shiftList
  print shiftList.measurements
  
