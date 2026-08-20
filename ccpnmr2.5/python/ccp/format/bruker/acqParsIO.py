"""
======================COPYRIGHT/LICENSE START==========================

acqParsIO.py: I/O for Bruker acquisition files

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

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
# Import general functions
from memops.universal.Util import returnInt

import os

from ccp.format.bruker.generalIO import BrukerGenericFile, BrukerParFile, BrukerPar
from ccp.format.general.acqParsDefs import FreqDimReference, FidDimReference
from ccp.format.general.acqParsDefs import SpectrometerReference, GenericAcqParData

from memops.universal.Io import splitPath, joinPath

dimAttrTags = ['fidDimRefs', 'freqDimRefs', 'isComplex',
               'sf', 'baseFreq', 'nuc', 'sw', 'npts', 'numPointsValid',
               'refppm', 'refpt']

#####################
# Class definitions #
#####################

class BrukerAcqParData(BrukerGenericFile,GenericAcqParData):

  def initialize(self,numDim,parHelp):
  
    self.numDim = numDim
    self.parHelp = parHelp
    
    (self.fileDir,void) = splitPath(self.name)
    
    # TODO need to separate out pulse program!!
    self.parFiles = {}
    
    self.acqMode = None
    self.acqModeValues = {}
    
    self.pulseProgramInfo = None
    self.pulProgName = None
    
    self.fnModes = ['undef'] * (self.numDim - 1)
    self.fnModeDict = {
    
      0: ('undef',  2.0),
      1: ('qf',     1.0),
      2: ('qseq',   1.0),
      3: ('tppi',   1.0), # Correct?
      4: ('states', 2.0), 
      5: ('states-tppi',  2.0),
      6: ('echo-antiecho',2.0)
    
    }
    
    self.stdIncrements = ('0','10','11')
    
    self.software = None
    self.version = None
    
    self.dateString = None
    
    # initialise per-dim attributes
    for tag in dimAttrTags:
      setattr(self, tag, numDim*[None])

  def get(self,verbose = 0):

    ###################################
    # Read acqu(n)s and acqu(n) files #
    ###################################
    
    baseFileName = 'acqu'
    
    self.parFiles.update({baseFileName + 's': {},baseFileName:{}})

    if self.numDim >=2:
      
      for n in range(2,self.numDim + 1):

        fileName = baseFileName + str(n)
        self.parFiles.update({fileName + 's': {},fileName:{}})
        
    
    files = os.listdir(self.fileDir)

    for file in files:
    
      if self.parFiles.has_key(file):
      
        #
        # Get acqu/acqus files: most information comes from here 
        #
        
        fileName = joinPath(self.fileDir,file)
        
        self.parFiles[file] = BrukerParFile(fileName,parHelp = self.parHelp)
        self.parFiles[file].read(verbose = verbose)

      elif file == 'pulseprogram':
      
        #
        # Also read pulseprogram info if available
        #
      
        self.parFiles['pulseprogram'] = {}
        self.parFiles['pulseprogram']['increment'] = []
        self.parFiles['pulseprogram']['name'] = ""
        self.parFiles['pulseprogram']['comment'] = ""
        self.parFiles['pulseprogram']['fnMode'] = None
        self.parFiles['pulseprogram']['echoAntiEcho'] = 0  # Should be 1 or more if echo-AntiEcho, especially igrad indicative!!
        
        #
        # TODO TODO HERE have to use dictionary information, if available, to select
        # the correct parameters. If they were found, use those in code below!!!
        #
        # IS NOW PULSEPROGRAMFILE!!
        #
        
        fin = open(joinPath(self.fileDir,file))
        
        line = fin.readline()
        
        pulseLine = ""
        
        while line:

          hashSearch = self.patt[self.format + 'HashComment'].search(line)
          
          if not pulseLine and hashSearch:
            pulseLine = hashSearch.group(1).strip()
          elif not self.parFiles['pulseprogram']['comment']:
            semiColonSearch = self.patt[self.format + 'SemiColonComment'].search(line)
            if semiColonSearch:
              self.parFiles['pulseprogram']['comment'] = semiColonSearch.group(1).strip()
          
          if not self.parFiles['pulseprogram']['name'] and pulseLine and self.parFiles['pulseprogram']['comment']:
            pulseNameSearch = self.patt[self.format + 'PulseName'].search(pulseLine)
            if pulseNameSearch:
              self.parFiles['pulseprogram']['name'] = pulseNameSearch.group(1)
            else:
              self.parFiles['pulseprogram']['name'] = "generic"
          
          if not self.parFiles['pulseprogram']['fnMode']:
            fnModeSearch = self.patt[self.format + 'FnMode'].search(line)
            if fnModeSearch:
              self.parFiles['pulseprogram']['fnMode'] = fnModeSearch.group(1).strip()

          
            # TODO look for dictionary match to set parameters correctly!!
            # self.pulseProgramInfo = ...
              
          # Increment/Decrement search here - could be obsolete. Also have to add other searches
          # based on dictionary if match found!!
        
          incrSearch = self.patt[self.format + 'InOrDecrement'].search(line)
          # TODO LET 0,10 take preference over everything else?!?!
          # DO bit of after-processing in case there are non 10/0 ones!
          if incrSearch:
            incNum = incrSearch.group(2)
            self.parFiles['pulseprogram']['increment'].append(incNum)
            
          #
          # See if acquired with Echo-AntiEcho
          # TODO: is this Utrecht specific or general? Ask Daniel!?
          #
          
          echoAntiEchoSearch = self.patt[self.format + 'EchoAntiEcho'].search(line)
          
          if echoAntiEchoSearch:
            self.parFiles['pulseprogram']['echoAntiEcho'] += 1
            
          line = fin.readline()
        
        #
        # Do some after-processing of the increment/decrement stuff. Basically, if the
        # count of 'standard' increments is numdim -1, then remove all other id/dd stuff foun
        # in favour of these (see list in initialise).
        #
        # TODO: has to be dictionary-settable at one point!!
        #
        
        increments = self.parFiles['pulseprogram']['increment']
        if increments:

          if len(increments) < self.numDim - 1:
            print "  Warning: not enough increments or decrements found to determine valuePerPoint."

          elif len(increments) > self.numDim - 1:
            stdIncrCount = 0
            for stdIncr in self.stdIncrements:
              if stdIncr in increments:
                stdIncrCount += 1

            if stdIncrCount == self.numDim - 1:
              for incr in increments[:]:
                if incr not in self.stdIncrements:
                  increments.pop(increments.index(incr))
                  print "  Warning: removed increment or decrement %s in favour of standard ones." % incr

    self.process()
  
  
  def processv21(self):
    """process function for TOPSPIN versions >= 2.1
    """
    
    # Tolerance for letting SW_h override SW * sf
    swhzTolerance = 0.001
  
    mainPars = self.parFiles['acqus'].parameters
    
    # get date string
    match = self.patt[self.format + "Date"].match(self.parFiles['acqus'].infoStrings[0])
    if match:
      self.dateString = match.group()
    else:
      self.dateString = ""
    self.specRefs = specRefs = {}
    
    # default specType
    specType = None
    
    #get pulProgName:
    if mainPars.has_key('PULPROG'):
      self.pulProgName =  mainPars['PULPROG'].values[0]
    
    # convert data one dimension at a time
    for ii in range(self.numDim):
      
      # set up file names
      if ii == 0:
        acquname = 'acqu'
      else:
        acquname = 'acqu' + str(ii+1)
      acquPars =  self.parFiles[acquname].parameters
      acqusPars =  self.parFiles[acquname + 's'].parameters
      
      # set complexFactor (complex/real handling). 
      # NB assumes acq dim is alwasy complex
      if acqusPars.has_key('FnMODE') and ii != 0:
        (fnMode,complexFactor) = self.fnModeDict[acqusPars['FnMODE'].values[0]]
        self.fnModes[ii-1] = fnMode
      else:
        complexFactor = 2.0
      
      # set 
      self.isComplex[ii] = (complexFactor != 1.0)
      npts = int(acquPars['TD'].values[0])
      if npts != 1:
        # npts == 1 is default for unused/projection dimensions. Keep as ia.
        npts = int(npts/complexFactor)
      self.npts[ii] = npts
      numPointsValid = int(acqusPars['TD'].values[0])/complexFactor
      self.numPointsValid[ii] = numPointsValid
      nucleus = self.convertNucleusCode(acqusPars['NUC1'].values[0])
      self.nuc[ii] = nucleus
      self.sf[ii] = sf = float(acqusPars['SFO1'].values[0])
      self.baseFreq[ii] = baseFreq = float(acqusPars['BF1'].values[0])
      swppm = float(acqusPars['SW'].values[0])
      swhz = swppm * sf
      if acqusPars.has_key('SW_h'):
        # NB SW seems more reliable than SW_h. 
        # Use SW_h only for greater precision where the values fit.
        val = float(acqusPars['SW_h'].values[0])
        if abs(val/swhz - 1.0) < swhzTolerance:
          swhz = val
      self.sw[ii] = swhz
      if numPointsValid > npts:
        print ("  WARNING: %s dim %s has greater number of valid points than actual ones" 
               % (self.fileDir, ii))

      # set approx referencing (starting value. Formula from Bruker manual)
      self.refpt[ii] = 1.0
      self.refppm[ii] = 1.0e6 * (sf / baseFreq - 1.0) + 0.5 * swhz / baseFreq
      
      # set self.specRefs for FormatConverter purposes
      self.specRefs[nucleus] = SpectrometerReference(sf, baseFreq)
      
      #if ii != 0:
      #  from ccp.general.Constants import chemShiftRefRatios
      #  refRatio = chemShiftRefRatios.get(nucleus)
      #  if refRatio is not None:
      #    bf2 = self.baseFreq[0] * refRatio
      #    xx = 1.0e6*(sf/bf2 -1.0)
      #    print '###1', baseFreq, bf2, xx, xx + 0.5*swhz/bf2 - self.refppm[ii]
      #
      #print ('###2', ii, acquname, sf, baseFreq, swhz, swppm, self.refppm[ii], 
      #       1.0e6 * (sf / baseFreq - 1.0),  0.5 * swhz / baseFreq)
      
      # set self.fidDimRefs - for FormatCOnverter use
      timePerPoint = complexFactor / (2.0 * swhz)
      self.fidDimRefs[ii] = FidDimReference(npts, numPointsValid, 
                                            timePerPoint, nucleus)
      
      # set self.freqDimRefs - for FormatConverter use 
      # for setting up future procesasing.
      # For actual processing data read procs (etc.) file
      (nPoints,nPointsOrig,valuePerPoint) = self.getFreqDimInfo(timePerPoint,
                                                                numPointsValid)
      self.freqDimRefs[ii] = FreqDimReference(nPoints, nPointsOrig,
                                              timePerPoint, nucleus)
                  
  
  def process(self):
  
    mainPars = self.parFiles['acqus'].parameters
    
    #
    # Set parmeters that are the same for both versions
    #
    
    # get ns
    self.numScans = returnInt(mainPars['NS'].values[0])
     
        
    # get acqModCode
    
    if mainPars.has_key('AQ_mod'):
    
      acqModeCode = mainPars['AQ_mod'].values[0]
      self.acqMode = 'QSIM'
     
      if acqModeCode == 3:
        self.acqMode = 'DQD'
        
        for par in ('DECIM','DSPFVS'):
          if mainPars.has_key(par):
            self.acqModeValues[par] = mainPars[par].values[0]
        
      elif acqModeCode == 1:
        pass
        
      else:
        print "ERROR: Bruker AQ_mod %d not recognized. Using QSIM." % acqModeCode
      

    #
    # Set software and choose processing depending on version
    #
    
    if self.parFiles['acqus'].headerDict.has_key("TITLE"):
      title = self.parFiles['acqus'].headerDict["TITLE"]
      for softwareName in ("XWIN-NMR","TOPSPIN","TopSpin"):
        if softwareName in title:
          self.software = softwareName
      
      if 'TopSpin' in title: # from Bruker 3.5 or so
        import re
        patt = re.compile("TopSpin\s*(\d+\.*\d*)")
      else:
        patt = self.patt[self.format + "Version"]
      versionSearch = patt.search(title)
      
      processed = False
      
      if versionSearch:
        self.version = versionSearch.group(1)
        
        # do special processing if this is TOPSPIN >= 2.1
        if self.software in ('TOPSPIN', 'TopSpin'):
          try:
            floatVersion = float(self.version)
            if floatVersion >= 2.1:
              self.processv21()
              processed = True
          except ValueError:
            pass
    
    if not processed:
      # process with pre-21 function
      self.processpre21()
    
    
    # Reset specType if DQD data (must be DMX or higher)
    # NBNB must be done at end, where specRefs aer set
    #

    if mainPars.has_key('DECIM'):
      for nucleus in self.specRefs.keys():
        self.specRefs[nucleus].setSpecType('DMX')
  
  
  
  def processpre21(self):
  
    mainPars = self.parFiles['acqus'].parameters
    #
    # Set pulseprogram info
    #
  
    pulseProgram = None  
    if self.parFiles.has_key('pulseprogram'):
      pulseProgram = self.parFiles['pulseprogram']
      ss = pulseProgram.get('name')
      if ss:
        self.pulProgName = ss
  
    refs = {}
    
    #
    # Try get date out
    #
    
    dateInfo = self.parFiles['acqus'].infoStrings[0]
    if dateInfo.count('T '):
      dateString = dateInfo[:dateInfo.index('T ') + 1]
    else:
      dateString = dateInfo[:19].strip()
    
    self.dateString = dateString
    
    #
    # Do quick check for 3D,4D to see if dimension order in filename matches
    # the acqu2s/acqu3s, and switch their info content if necessary.
    #
    # TODO: THIS IS NOT FOOLPROOF! NUC1 for acqu2s, acqu3s not always set to 
    # nucleus used in that dimension! Is so for some setups, though (IE)
    #
    # TODO: ONLY works for 3Ds currently!
    #
    
    if (self.numDim == 3):
    
      if not hasattr(self.parFiles['acqu3s'],'parameters'):
        print "  Illegal number of dimensions!! %s is not a 3d..." % self.fileDir
        return

      if not self.parFiles['acqu3s'].parameters.has_key('NUC1'):

        #
        # No acqu3s file... basically 3D recorded as 2D (or wrong experiment put in..)
        # Assuming that L1, L2 give npts/2, getting as much info as possible from this...
        #

        print "  Warning: no acqu3s file for Bruker 3D experiment."

        pointsProduct = self.parFiles['acqu2s'].parameters['TD'].values[0]

        npts1 =  self.parFiles['acqus'].parameters['L'].values[1]
        npts2 =  self.parFiles['acqus'].parameters['L'].values[2]

        if pointsProduct == (npts1 * npts2 * 4):
        
          print "           assuming L1,L2 give number of complex points."

          # TODO: does it make sense to set this for acqu2/3?
          self.parFiles['acqu3s'].parameters['TD'] = BrukerPar('TD')
          self.parFiles['acqu3'].parameters['TD'] = BrukerPar('TD')


          self.parFiles['acqu2s'].parameters['TD'].values[0] = npts1
          self.parFiles['acqu2'].parameters['TD'].values[0] = self.parFiles['acqu'].parameters['L'].values[1]

          self.parFiles['acqu3s'].parameters['TD'].addValues(1,[npts2])
          self.parFiles['acqu3'].parameters['TD'].addValues(1,[self.parFiles['acqu'].parameters['L'].values[2]])

          #
          # Select appropriate increment per dimension
          # Bit messy... this whole script has to be cleaned up
          # TODO: very nasty to put tk stuff in here... have to:
          #        - Make whole def a class
          """
          incrementList = []

          for value in mainPars['IN'].values:
            incrementList.append(value)  

          popup = IncrementSelect(guiParent,[2,3],incrementList)

          for i in range(0,len(popup.increments)):
            increment = popup.increments[i]
            dim = i + 1

            refs['inc'][dim] = mainPars['IN'].values.index(increment)

            # TODO: nucleus!!!

          popup.destroy()
          """

        else:
        
          print "           cannot handle current case. Contact author."
        
      
    elif (self.numDim > 3):
    
      print "  ERROR: %s had more than 3 dims for acqu<n>s file order (cannot handle this)!" % (self.fileDir)

    #
    # Set increment numbers found in pulse program
    #
    
    refs['inc'] = [-1]

    if pulseProgram and self.numDim <= len(pulseProgram['increment']):
      print "  Warning: %s number of incremented delays found does not match number of dimensions!" % (self.fileDir)
    
    for i in range(0,self.numDim-1):
      if pulseProgram and len(pulseProgram['increment']) > i:
        refs['inc'].append(int(pulseProgram['increment'][i]))
      else:
        refs['inc'].append(i)
    
    #
    # If AQSEQ == 1 and experiment is a 3D, then swap the first two indirect increments 
    # to correspond to the acqu files
    #
    
    if mainPars['AQSEQ'].values[0] == 1 and self.numDim == 3:
      tempInc = refs['inc'][1]
      refs['inc'][1] = refs['inc'][2]
      refs['inc'][2] = tempInc
      print "  Warning: swapped increments for 3D based on AQSEQ setting %d!" % mainPars['AQSEQ'].values[0]
    
    #
    # TODO: Check if acqu<n>s info corresponds to reported nucleus order - reset values if so.
    #
     
    #
    # Make a list of reference values for accessing acqus file
    #
    
    self.setupBrukerRefList()
    nullValue = 'off'

    self.getSpecRefs(mainPars,nullValue)

    #
    # Now set the fid and freq info from the file
    #
    
    for i in range(0,self.numDim):
    
      #
      # Handling depends on which dim (complex points or not)
      #
    
      if i == 0:
        
        acqusPars = mainPars
        acquPars = self.parFiles['acqu'].parameters
        complexFactor = 2.0
      
        #
        # Acquisition dimension always complex
        #
        
        #
        # Spectral width in Hz per dimension and associated increment
        #swText = '(from acqus file)'
        #incText = '(calculated from SW_h, assuming complex)'
        #
        
        swh = float(acqusPars['SW_h'].values[0])        
        valuePerPoint = 1 / swh

      else:
      
        #
        # Try to get info on acquisition mode for dimension...
        #
        
        #currentAcqusFile = self.parFiles['acqu' + str(i+1) + 's']
        acqusPars = self.parFiles['acqu' + str(i+1) + 's'].parameters
        acquPars = self.parFiles['acqu' + str(i+1)].parameters

        if acqusPars.has_key('FnMODE'):
          (fnMode,complexFactor) = self.fnModeDict[acqusPars['FnMODE'].values[0]]
          self.fnModes[i-1] = fnMode
          
        else:
          complexFactor = 2.0
        #
        # Increment (based on values read in from pulseprogram!). Multiplied by 2: assuming complex.
        #incText = "(using IN%s in acqus file)" % brukerrefs['inc'][i]
        #
        
        #
        # TODO: will this cause problems if the acqu2/3s files not switched yet?
        # Keep in mind!!
        #
        
        valuePerPoint = complexFactor * float(mainPars['IN'].values[refs['inc'][i]])

      
      # set nucleus
      #nucleus = mainPars['NUC1'].values[0]
      nucleus = self.convertNucleusCode(acqusPars['NUC1'].values[0])
      self.nuc[i] = nucleus
      # Number of COMPLEX points per dimension.
      
      # set values into parameter arrays
      numPoints = int(acquPars['TD'].values[0]) / 2 
      numPointsValid = int(acqusPars['TD'].values[0]) / 2 
      self.npts[i] = numPoints
      self.numPointsValid[i] = numPointsValid
      self.sw[i] = swhz = numPoints * valuePerPoint
      self.isComplex[i] = (complexFactor != 1.0)
      self.sf[i] = sf = float(acqusPars['SFO1'].values[0])
      self.baseFreq[i] = baseFreq = float(acqusPars['BF1'].values[0])
      # set approx referencing (starting value. Formula from Bruker manual)
      self.refpt[i] = 1.0
      self.refppm[i] = 1.0e6 * (sf / baseFreq - 1.0) + 0.5 * swhz / baseFreq
      
      #
      # Now set the values for the FID
      #

      if numPointsValid > numPoints:

        print "  ERROR: %s has greater number of valid points than actual ones" % (self.fileDir)

      self.fidDimRefs[i] = FidDimReference(numPoints,numPointsValid,
                                           valuePerPoint,nucleus)
      
      #
      # And for the Freq
      #
      (numPoints,numPointsOrig,valuePerPoint) = self.getFreqDimInfo(valuePerPoint,numPointsValid)
      self.freqDimRefs[i] = FreqDimReference(numPoints,numPointsOrig,
                                             valuePerPoint,nucleus)
      
      
      
      #
      # Phase is automatically set to 0.0 each
      # TODO: Could in principle access Bruker procpars for this
      #
        

  def setupBrukerRefList(self):

    self.refList = []
    reflistTmp = ['SFO','O','NUC']
    for n in range(1,9):
      tlist = []
      for m in range(0,len(reflistTmp)):
        tlist.append(reflistTmp[m] + str(n))
      self.refList.append(tlist)
