
"""
======================COPYRIGHT/LICENSE START==========================

BrukerParams.py

Copyright (C) 2008 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

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
import os.path

from memops.general.Implementation import ApiError

from ccp.format.bruker.processingParsIO import BrukerProcessingParData
from ccp.format.bruker.acqParsIO import BrukerAcqParData
from ccp.format.bruker.acqParsHelpReader import BrukerParHelp

from ccp.format.spectra.params.ExternalParams import ExternalParams, dimParamDefs

class BrukerParams(ExternalParams):

  format = 'Bruker'
  # how many directories to keep in file path for default DataUrl split
  keepDirectories = 4

  def __init__(self, procs_file, **kw):
    
    # Accept procs directory or any file in it
    # Assume correct file is named 'procs' if there is such a file
    if os.path.isdir(procs_file):
      procs_file = os.path.join(procs_file, 'procs')
      
    else:
      dd, file = os.path.split(procs_file)
      if file != 'procs':
        ss = os.path.join(dd, 'procs')
        if os.path.isfile(ss):
          print 'Switching to input file: %s' % ss
          procs_file = ss
    
    self.procs_file = procs_file

    ExternalParams.__init__(self, **kw)

  # ExternalParams requires this to be defined
  def parseFile(self):
    
    # get processed file data
    try:
      help = BrukerParHelp('help')
      procParData = BrukerProcessingParData(self.procs_file, help.tags)
      procParData.get()
    except IOError, e:
      raise ApiError(str(e))
    
    # get acquisition file data
    procDir = procParData.fileDir
    func = os.path.dirname
    acqDir = func(func(procDir))
    files = set(os.listdir(acqDir))
    nAxes = 1
    while True:
      i = nAxes + 1
      if ('acqu%ds' % i) in files:
        nAxes = i
      else:
        break
    numDim = procParData.numDim
    if nAxes < numDim:
      nAxes = numDim
      acqParData = None
    else:
      acqParData = BrukerAcqParData(os.path.join(acqDir, 'acqus'), 
                                    nAxes, help.tags)
      acqParData.get()
      
      # check if we have a processed projection, a projection experiment or what
      if nAxes > numDim:
        extraDimHasPoints = [x for x in acqParData.npts[numDim:nAxes] if x > 1]
        if extraDimHasPoints:
          # wb104: added the npts=1 check on 11 Mar 2011
          # check if any dims have npts=1, in which case can eliminate
          onePointDims = [i for i in range(nAxes) if acqParData.npts[i]==1]
          if len(onePointDims) + numDim == nAxes:
            # can rescue this case
            # looks like don't use most of the below but what the heck, keep consistent
            for attr in ('baseFreq', 'npts', 'nuc', 'numPointsValid', 'refppm', 'refpt', 'sf', 'sw'):
              for n in reversed(onePointDims):
                try:
                  del getattr(acqParData, attr)[n]
                except:
                  # not sure should pass, perhaps should give up instead??
                  # but hopefully never in this situation
                  pass
          else:
            # something is wrong here,
            # maybe we have e.g. a 2D projection of an original 3D spectrum
            # Anyway the acquisition parameters can not be trusted to match the
            # processing dimensions.
            print (
             'WARNING, Spectrum is a computed, %sD projection or slice of a %sD Experiment'
             % (numDim, nAxes))
            acqParData = None
          nAxes = numDim
          
    self.nAxes = nAxes
    
    # set data from proc files
    aPars = procParData.aPars
    aPars['nAxes'] = max(nAxes, numDim)
    aPars['scale'] = pow(2, aPars['brukerScale'])
    self.setAttrs(aPars)
    
    # add data from acqus files
    if acqParData is not None:
    
      self.acquisitionDim = 0
    
      self.numScans = acqParData.numScans
      
      dateString = acqParData.dateString
      if dateString:
        self.date = dateString
      
      pulProgName = acqParData.pulProgName
      if pulProgName is not None:
        self.pulProgName = pulProgName
        self.pulProgType = 'bruker'
      
      # reset basefreq and sf to more accurate acqu values
      for tag in ('baseFreq', 'sf'):
        if hasattr(acqParData, tag):
          ll = getattr(acqParData, tag)
          getattr(self, tag)[:len(ll)] = ll
      
      if hasattr(acqParData, 'nuc'):
        nucs = self.nuc
        acqnucs = acqParData.nuc
        for ii, nuc in enumerate(nucs):
          #if nuc is None:
          # Added for TopSpin3.2: if nuc is 'off' also take from acqus
          if nuc in (None, 'off'):
            nucs[ii] = acqnucs[ii]
      
      # Add projection dimensions:
      # Assumes that first numDim acq dimensions are processed, and that
      # last nAxes-numDim dimensions are projection and add to last real dim
      # NB the test for projection experiment versus processed projection
      # has been done previously
      newInfoTags =  ('sf', 'baseFreq', 'nuc', 'sw', 'refpt', 'refppm')
      
      if nAxes > numDim:
        # projection experiment (e.g. Prodecomp)
        for ss in newInfoTags:
          getattr(self, ss)[numDim:nAxes] = getattr(acqParData, ss)[numDim:nAxes]  
      
      # add extra dimensions if extra displayNames present 
      # (Prodecomp exps with pos/neg peaks Ca,Cb and possibly Ha,Hb)
      posNegCaCb = False
      displayNames = self.extraData.get('displayNames')
      if displayNames:
        if len(displayNames) > nAxes:
          if 'Ca' in displayNames and 'Cb' in displayNames:
            posNegCaCb = True
            indx = [displayNames.index('Ca'), displayNames.index('Cb')]
            indx.sort()
            for ss in dimParamDefs.keys():
              ll = getattr(self, ss)
              if ll:
                ll.insert(indx[1], ll[indx[0]])
            nAxes = self.nAxes = nAxes + 1
            
            if len(displayNames) > nAxes:
              if 'Ha' in displayNames and 'Hb' in displayNames:
                indx = [displayNames.index('Ha'), displayNames.index('Hb')]
                indx.sort()
                for ss in dimParamDefs.keys():
                  ll = getattr(self, ss)
                  if ll:
                    ll.insert(indx[1], ll[indx[0]])
              nAxes = self.nAxes = nAxes + 1
      
        # reset carrier and reference from CNST for C, Cab, and Ca
        constants = acqParData.constants
        for ii, name in enumerate(displayNames):
          
          O1 = 1
          if name == 'CO':
            O1 = constants.values[21]
          elif name == 'Cab' or (posNegCaCb and name in ('Ca','Cb')):
            O1 = constants.values[23]
          elif name == 'Ca':
            O1 = constants.values[22]
            
          if O1 == 1:
            continue
          
          baseFreq = self.baseFreq[ii]
          if baseFreq:
            self.sf[ii] = baseFreq + 1.0e-6 * O1
            self.refpt[ii] = 1.0
            self.refppm[ii] = O1 + 0.5*self.sw[ii]/baseFreq

if (__name__ == '__main__'):

  import sys
  if (len(sys.argv) != 2):
    print 'Error: required argument: <procs_file>'
    sys.exit(1)

  procs_file = sys.argv[1]
  params = BrukerParams(procs_file)
