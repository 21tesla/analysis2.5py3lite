#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

projectIO.py: I/O for nmrView project (star) file

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

import os, string

# Import general functions
from memops.universal.Io import getTopDirectory
from memops.universal.Util import returnFloat

from ccp.format.nmrView.generalIO import NmrViewGenericFile

from ccp.format.nmrStar.util import getNmrStarValue

#
# Import Jurgen Doreleijers' Star reader/writer stuff
# Has been slightly extended to encompass nmrView star
#

import ccp.format.nmrStar.bmrb.File as nmrStar

#####################
# Class definitions #
#####################

from ccp.format.nmrView.peaksIO import NmrViewPeak
from ccp.format.nmrView.peaksIO import NmrViewPeakFile

from ccp.format.nmrView.sequenceIO import NmrViewSequenceElement
from ccp.format.nmrView.sequenceIO import NmrViewSequence
from ccp.format.nmrView.sequenceIO import NmrViewSequenceFile

class NmrViewStarFile(NmrViewGenericFile):

  # Information on file level

  def initialize(self):

    self.peakFiles = []
    self.sequenceFile = None

  def read(self,verbose = 0):

    origStarFile = nmrStar.File(verbosity = 2, filename = self.name)

    #
    # Read star file in one gulp, make modifications for nmrView
    #

    if origStarFile.read(nmrView_type = 1):
      print "  Error reading nmrView star file %s" % self.name
      return

    for origSaveFrame in origStarFile.datanodes:

      saveFrameCat = None

      for tagtable in origSaveFrame.tagtables:
        try:
          tagIndex = tagtable.tagnames.index('_Saveframe_category')
          saveFrameCat = tagtable.tagvalues[tagIndex][0]
          break
        except:
          pass


      if saveFrameCat == 'peak_lists':

        for tagtable in origSaveFrame.tagtables:

          #
          # Here do tags
          #

          if tagtable.free != None:

            # Ignore tags for now use following if need them:
            #starValue = getNmrStarValue(tagtable,tagName)
            pass

          #
          # Here do tables
          #

          else:
            # set tableName for useful ones
            if tagtable.tagnames.count('_Peaklist_name') > 0:

              for i in range(0,len(tagtable.tagvalues[0])):
                # set tagname for useful ones
                peakListName = getNmrStarValue(tagtable,'_Peaklist_name',i)
                peakDimLabels = getNmrStarValue(tagtable,'_Peaklist_label',i)
                peakDimLabelList = string.split(peakDimLabels)

                peakFile = NmrViewPeakFile(self.name)
                peakFile.setSpectrumInfo(peakListName,len(peakDimLabelList))
                peakFile.dimCodes = peakDimLabelList

                self.peakFiles.append(peakFile)


      elif saveFrameCat == 'spectral_peak_list':

        #
        # Select right peakFile object
        # Warning: these should have been read in already
        #

        peakListName = string.replace(origSaveFrame.title,'spectral_peak_list_','')

        for peakFile in self.peakFiles:

          if peakFile.specNames[0] == peakListName:
            peakFile.sws = peakFile.numDims[0] * [0.0]
            peakFile.sfs = peakFile.numDims[0] * [0.0]
            numDim = len(peakFile.dimCodes)
            break

        #
        # Now get info
        #

        for tagtable in origSaveFrame.tagtables:

          #
          # Here do tags
          #

          if tagtable.free != None:

            # Ignore tags for now use following if need them:
            #starValue = getNmrStarValue(tagtable,tagName)
            pass

          #
          # Here do tables
          #

          else:

            if tagtable.tagnames.count('_Peaklist_label') > 0:

              for i in range(0,len(tagtable.tagvalues[0])):
                # set tagname for useful ones
                dimLabel = getNmrStarValue(tagtable,'_Peaklist_label',i)
                sw = getNmrStarValue(tagtable,'_Peaklist_sw',i)
                sf = getNmrStarValue(tagtable,'_Peaklist_sf',i)

                dimIndex = peakFile.dimCodes.index(dimLabel)
                peakFile.sws[dimIndex] = returnFloat(sw)
                peakFile.sfs[dimIndex] = returnFloat(sf)

            elif tagtable.tagnames.count('_Peak_list_number') > 0:

              for i in range(0,len(tagtable.tagvalues[0])):

                peaknum = getNmrStarValue(tagtable,'_Peak_list_number',i)
                volume = getNmrStarValue(tagtable,'_Intensity_volume',i)
                intensity = getNmrStarValue(tagtable,'_Intensity_height',i)
                status = getNmrStarValue(tagtable,'_Peak_status',i)

                # These not used for now
                flag = getNmrStarValue(tagtable,'_flag',i)
                comment = getNmrStarValue(tagtable,'_comment',i)
                comment = self.removeCurlyBraces(comment)

                #
                # Peakdim stuff
                #

                assign = numDim * ['']
                ppm = numDim * [0]
                width = numDim * [0]
                box = numDim * [0]
                shape = numDim * ['']

                for j in range(0,numDim):
                  dimText = '_Dim_' + str(j + 1)

                  assign[j] = getNmrStarValue(tagtable,dimText + '_label',i)
                  assign[j] = self.removeCurlyBraces(assign[j])

                  ppm[j] = getNmrStarValue(tagtable,dimText + '_chem_shift',i)
                  width[j] = getNmrStarValue(tagtable,dimText + '_line_width',i)
                  box[j] = getNmrStarValue(tagtable,dimText + '_bounds',i)

                  shape[j] = getNmrStarValue(tagtable,dimText + '_error',i)
                  shape[j] = self.removeCurlyBraces(shape[j])

                  # Ignoring other attributes! 'user' also is in curly braces

                peakFile.peaks.append(NmrViewPeak(peaknum,assign,ppm,width,box,shape,volume,intensity,status))


      elif saveFrameCat == 'monomeric_polymer':

        if not self.sequenceFile:
          sequenceFile = NmrViewSequenceFile(self.name)
          self.sequenceFile = sequenceFile

        molName = origSaveFrame.title
        self.sequenceFile.sequences.append(NmrViewSequence(molName = molName))

        for tagtable in origSaveFrame.tagtables:


          #
          # Here do tags
          #

          if tagtable.free != None:

            # Ignore tags for now use following if need them:
            #starValue = getNmrStarValue(tagtable,tagName)
            pass

          #
          # Here do tables
          #

          else:
            if tagtable.tagnames.count('_Residue_seq_code') > 0:

              for i in range(0,len(tagtable.tagvalues[0])):
                seqCode = getNmrStarValue(tagtable,'_Residue_author_seq_code',i)
                code3Letter = getNmrStarValue(tagtable,'_Residue_label',i)

                # These two ignored for now
                chainCode = getNmrStarValue(tagtable,'_Segment_definition_code',i)
                resId = getNmrStarValue(tagtable,'_Residue_seq_code',i)

                self.sequenceFile.sequences[-1].elements.append(NmrViewSequenceElement(seqCode,code3Letter))

      else:

        print "  Currently can't handle nmrView star file saveframe category %s. Ignored." % saveFrameCat


  def removeCurlyBraces(self,text):

    """
    Removes curly braces from a string
    """

    searchObj = self.patt[self.format + 'CurlyBrace'].search(text)

    if searchObj:
      return searchObj.group(1)

    else:
      return string

###################
# Main of program #
###################

if __name__ == "__main__":
                                                      
  files = ['../../reference/nmrView/Brucd44_new2.str',
           #'python/ccpnmr/format/examples/data/withass20010817.str'
           ]

  
  for inFile in files:
    
    nmrViewStarFile = NmrViewStarFile(os.path.join(getTopDirectory(),inFile))

    nmrViewStarFile.read(verbose = 1)
    
    for peakFile in nmrViewStarFile.peakFiles:
      print peakFile.specNames
      print peakFile.numDims
      print peakFile.peaks[0].ppm
      print peakFile.peaks[-1].ppm
    
    for sequence in nmrViewStarFile.sequenceFile.sequences:
      for seqEl in sequence.elements:
  
        print seqEl.seqCode, seqEl.code3Letter
