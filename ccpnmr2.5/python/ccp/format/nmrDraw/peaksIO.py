#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for nmrDraw peak list files

Copyright (C) 2006 Wim Vranken (European Bioinformatics Institute)

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
from memops.universal.Util import returnInt, returnInts
from memops.universal.Util import returnFloat

from ccp.format.nmrDraw.generalIO import NmrDrawGenericFile

#####################
# Class definitions #
#####################

class NmrDrawPeakFile(NmrDrawGenericFile):

  def initialize(self):

    self.peaks = []
    
    self.specNames = [self.name]
    
  def setDimensions(self,ndim):
  
    self.numDims = [ndim]
    
    self.refDimCodes = []
    self.refDimInfo = {}
    
    #
    # Reference data setup
    #
    
    self.refData = [['INDEX', '%5d','num',returnInt,True],
                    [['n_AXIS'],'%9.3f','points',returnFloat,True],
                    [['Dn'],'%6.3f','pointsDev',returnFloat,True],
                    [['n_PPM'],'%8.3f','ppm',returnFloat,True],
                    [['n_HZ'],'%9.3f','hz',returnFloat,True],
                    [['nW'],'%7.3f','fwhm',returnFloat,True],             # Full width half maximum (height) in pts
                    [['nW_HZ'],'%8.3f','fwhm_Hz',returnFloat,True],       # Full width half maximum (height) in Hz
                    [['n1','n3'],'%4d','limits',returnInt,True],          # Bounding region in pts
                   	['HEIGHT','%13e','intensity',returnFloat,True],       # Height by parabolic interpolation
                   	['DHEIGHT','%13e','intensityError',returnFloat,True], # Error in height due to noise (estimated)
                   	['VOL','%13e','volume',returnFloat,True],             # Volume as sum of intensities peak area
                   	['PCHI2','%.5f','chiProb',returnFloat,True],          # Prob. that peak could be due to random noise
                   	['TYPE','%d','type',returnInt,True],                  # 1 = peak, 2 = random noise, 3 = truncation artefact
                   	['ASS','%s','assignment',str,True],                   # Assigment data (which format?)
                   	['CLUSTID','%4d','clustId',returnInt,True],           # ID of cluster peak belongs to (overlap)
                   	['MEMCNT','%4d','numPeaksClust',returnInt,True],      # Total number of peaks in cluster
                    ['TROUBLE','%4d','fitSuccesful',returnInt,False]      # Info from nlinLS
                   ]
    
    self.varCodes = []
    self.formatCodes = []

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading nmrDraw peak list %s" % self.name

    fin = open(self.name, 'rU')

    refDimInfo = {}
    refDimInfo['Label'] = []
    refDimInfo['highPoint'] = []
    refDimInfo['lowPoint'] = []
    refDimInfo['highPpm'] = []
    refDimInfo['lowPpm'] = []

    refDimCodes = []

    #
    # Start reading
    #

    line = fin.readline()

    while line:

      if not self.patt['emptyline'].search(line) and not self.patt[self.format + 'Remark'].search(line):

        if self.patt[self.format + 'Dataline'].search(line):
          
          dataValues = string.split(line)
          
          if dataValues[1][-4:] == 'AXIS':

            (void,axis,axisLabel,highPoint,lowPoint,highPpm,lowPpm) = dataValues
  
            refDimCodes.append(axis[0])
  
            refDimInfo['Label'].append(axisLabel)
            refDimInfo['highPoint'].append(returnInt(highPoint))
            refDimInfo['lowPoint'].append(returnInt(lowPoint))
            refDimInfo['highPpm'].append(returnFloat(highPpm[:-3]))
            refDimInfo['lowPpm'].append(returnFloat(lowPpm[:-3]))

        elif self.patt[self.format + 'Vars'].search(line):

          #
          # All DATA lines should be read in - setup file
          #
          
          if not refDimCodes:
            print "  Warning: DATA AXIS lines missing from header - deducing AXIS values from VARS line."
            
            refDimCodes = []
            findCode = '_AXIS'
            tCols = string.split(line)
            for tCol in tCols:
              if tCol.count(findCode):
                refDimCode = tCol.replace(findCode,'')
                refDimCodes.append(refDimCode)
                dimIndex = refDimCodes.index(refDimCode)
                refDimInfo['Label'].append(refDimCode)
                refDimInfo['lowPoint'].append(1)
                refDimInfo['highPpm'].append(0)
                refDimInfo['lowPpm'].append(0)
                       
            print "  Warning: set codes %s" % refDimCodes

          n = len(refDimCodes)
          self.setDimensions(n)
          self.refDimCodes = refDimCodes

          for refDimInfoKey in refDimInfo.keys():
            self.refDimInfo[refDimInfoKey] = refDimInfo[refDimInfoKey]

          #
          # Get vars line, check if matches up with reference info
          #

          cols = string.split(line)
          colNum = 1

          for refDatum in self.refData:

            refVarCode = refDatum[0]

            if type(refVarCode) == type([]):

              for dimRefCode in self.refDimCodes:

                for refVarCodeSingle in refVarCode:

                  matchString = string.replace(refVarCodeSingle,'n',dimRefCode)
                  if refDatum[4] and matchString != cols[colNum]:
                    print "  Error: can't read nmrDraw file - headers don't match."
                    print "%6s, col %d. '%s'<->'%s'" % (refDatum[4], colNum, matchString, cols[colNum])
                    return 0
                  elif refDatum[4]:
                    self.varCodes.append(matchString)
                    colNum += 1

            else:
              if refDatum[4] and refVarCode != cols[colNum]:
                print "  Error: can't read nmrDraw file - headers don't match."
                print "%6s, col %d. '%s'<->'%s'" % (refDatum[4], colNum, refVarCode, cols[colNum])
                return 0

              elif refDatum[4]:
                self.varCodes.append(refVarCode)
                colNum += 1

        elif self.patt[self.format + 'Format'].search(line):

          self.formatCodes = string.split(line)[1:]
          self.numCols = len(self.formatCodes)

        else:

          cols = string.split(line)
          
          if len(cols) == self.numCols:

            varsDict = {}

            for varCode in self.varCodes:
              varsDict[varCode] = cols.pop(0)

            self.peaks.append(NmrDrawPeak(self))
            self.peaks[-1].setVarsDict(varsDict)

          else:

            print "  Warning: unparsed nmrDraw peak list line:\n%s" % line

      line = fin.readline()

    fin.close()

  def write(self, verbose = 0):

    if verbose == 1:
      print "Writing nmrDraw peak list %s" % self.name

    fout = open(self.name,'w')

    numDim = len(self.refDimCodes)

    #
    # Write out header
    #

    fout.write("REMARK File written by CcpNmrFormat converter.")
    fout.write(self.newline * 2)

    for dim in range(0,numDim):
      dimCode = self.refDimCodes[dim]
      fout.write("DATA  %s_AXIS %-12s 1 %5d %8.3fppm %8.3fppm" % (self.refDimCodes[dim],
                                                                  self.refDimInfo['Label'][dim],
                                                                  self.refDimInfo['lowPoint'][dim],
                                                                  self.refDimInfo['highPpm'][dim],
                                                                  self.refDimInfo['lowPpm'][dim]))
      fout.write(self.newline)

    fout.write(self.newline)

    fout.write("VARS  ")

    formatCodes = []

    for refDatum in self.refData:

      refVarCode = refDatum[0]
      refFormatCode = refDatum[1]

      if type(refVarCode) == type([]):

        for dimRefCode in self.refDimCodes:

          for refVarCodeSingle in refVarCode:

            outString = string.replace(refVarCodeSingle,'n',dimRefCode)

            fout.write(" %s" % outString)

            formatCodes.append(refFormatCode)

      else:

        fout.write(" %s" % refVarCode)
        formatCodes.append(refFormatCode)

    fout.write(self.newline)

    fout.write("FORMAT")

    formatNum = 0

    for formatCode in formatCodes:

      fout.write(" " + str(formatCode))

    fout.write(self.newline * 2)  

    #
    # Write out peaks
    #

    for peak in self.peaks:

      line = ""

      for refDatum in self.refData:

        refVarCode = refDatum[0]
        refFormatCode = refDatum[1]
        refAttrName = refDatum[2]
        refFunc = refDatum[3]

        if type(refVarCode) == type([]):

          valueList = getattr(peak,refAttrName)

          for i in range(0,numDim):

            value = valueList[i]

            if type(value) == type([]):

              for valueElement in value:

                line += (refFormatCode % refFunc(valueElement))
                line += " "

            else:

              line += (refFormatCode % refFunc(value))
              line += " "

        else:

          value = getattr(peak,refAttrName)

          line += (refFormatCode % refFunc(value))
          line += " "

      fout.write(line[:-1])
      fout.write(self.newline)

    fout.close()

class NmrDrawPeak:

  def __init__(self,parent):
  
    self.parent = parent
  
    for refDatum in parent.refData:

      refVarCode = refDatum[0]
      refFormatCode = refDatum[1]
      refAttrName = refDatum[2]
      refFunc = refDatum[3]

      if type(refVarCode) == type([]):
      
        setattr(self,refAttrName,[])
        curAttr = getattr(self,refAttrName)
      
        for dimRefCode in parent.refDimCodes:
          
          if len(refVarCode) == 1:
            
            refVarCodeSingle = refVarCode[0]
            curAttr.append(None)
            
          else:
             
            curAttr.append([])
             
            for refVarCodeSingle in refVarCode:

              curAttr[-1].append(None)
    
      else:
      
        setattr(self,refAttrName,None)

  def setVarsDict(self,varsDict):
  
    for refDatum in self.parent.refData:

      refVarCode = refDatum[0]
      refFormatCode = refDatum[1]
      refAttrName = refDatum[2]
      refFunc = refDatum[3]
      refRequired = refDatum[4]
      
      if type(refVarCode) == type([]):
      
        curAttr = getattr(self,refAttrName)
      
        for i in range(0,len(self.parent.refDimCodes)):
        
          dimRefCode = self.parent.refDimCodes[i]
          
          if len(refVarCode) == 1:
            
            refVarCodeSingle = refVarCode[0]
               
            matchString = string.replace(refVarCodeSingle,'n',dimRefCode)

            if refRequired or varsDict.has_key(matchString):     
              curAttr[i] = refFunc(varsDict[matchString])
             
          else:
             
            curAttr.append([])
             
            for j in range(0,len(refVarCode)):
              
              refVarCodeSingle = refVarCode[j]
            
              matchString = string.replace(refVarCodeSingle,'n',dimRefCode)
              if refRequired or varsDict.has_key(matchString):     
                curAttr[i][j] = refFunc(varsDict[matchString])
    
      elif refRequired or varsDict.has_key(refVarCode):
      
        setattr(self,refAttrName,refFunc(varsDict[refVarCode]))
    
    
"""

INFO from nmrDraw manual pages:

PEAK TABLE PARAMETERS
     INDEX
          the unique peak ID number.

     X_AXIS ...
          the position of the peak in points, by parabolic inter-
          polation.

     DX ...
          estimate of the error in peak position  due  to  random
          noise, in points.

     X_PPM ...
          the position of the peak in PPM, by parabolic  interpo-
          lation.

     X_HZ ...
          the position of the peak in Hz, by parabolic interpola-
          tion.

     XW ...
          the full width at half height in points,  estimated  by
          parabolic interpolation.

     XW_HZ ...
          the full width at half height in Hz, estimated by para-
          bolic interpolation.

     X1/X3 ...
          coords of the bounding region surrounding the  peak  or
          its cluster in points.

     HEIGHT
          the peak height, estimated by parabolic interpolation.

     DHEIGHT
          the error peak height due to random noise, estimated by
          parabolic interpolation.

     VOL  the peak volume, estimated as the sum of intensities in
          the neighborhood of the peak maximum.

     ASS  location to hold assignment data.

     PCHI2
          the chi-square probability for the peak (i.e. probabil-
          ity that the peak could be due to random noise alone).

     CLUSTID
          identifies the cluster of peaks that a given peak over-
          laps.  Peaks with the same CLUSTID value are considered
          to be in one overlapped group.

     MEMCNT
          the total number of peaks which are in a  given  peak's
          cluster (i.e. peaks which have the same CLUSTID value).

     TYPE An  integer   describing   the   peak   classification;
          currently  1  =  Peak, 2 = Random Noise, 3 = Truncation
          artifact.

"""
