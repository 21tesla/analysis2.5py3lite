
"""
======================COPYRIGHT/LICENSE START==========================

NmrViewParams.py: A contribution to the CcpNmr Analysis program

Copyright (C) 2006 Gary Thompson (University of Leeds)

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

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk,
                       garyt@bmb.leeds.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references and acknowledge Gary S. Thompson:

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

# todo:
#   1. doesn't currently handle nmrview <spectrum>.par files      - done
#   2. there is no consistency check of file length               - done
#   3. the module doesn't check for the nmrview magic bytes       - done 
#   4. byte swapping has not been tested on a non intel machine

"""
import array
import sys
import traceback

from string import atoi
from string import atof
# from sets import Set
from os.path import splitext
from os.path import isfile
from os import stat
from copy import copy
from math import ceil
from memops.general.Implementation import ApiError

from ccp.format.spectra.params.ExternalParams import ExternalParams

class NmrViewParams(ExternalParams):

  format = 'NmrView'

  def __init__(self, file, **kw):

    self.dataFile = file
    self.head = 2048
    ExternalParams.__init__(self, **kw)

  byteOrderFlags  = {\
    ('34','18','AB','CD') : 'big',
    ('CD','AB','18','34') : 'little'
  }
  
  magicBytes = ['34','18','AB','CD']
  magicBytes.sort()
  
  def doByteSwap(self,cBuffer):
    assert cBuffer.itemsize==1
    bs0 = "%02X" % ord(cBuffer[0])
    bs1 = "%02X" % ord(cBuffer[1])
    bs2 = "%02X" % ord(cBuffer[2])
    bs3 = "%02X" % ord(cBuffer[3])
    
    
    #if __debug__:
    #  print bs0,bs1,bs2,bs3
    if sys.byteorder == self.byteOrderFlags[(bs0,bs1,bs2,bs3)]:
        return False
    else: 
        return True

  def nvParToNvNvFile(self):
    self.dataFile = splitext(self.dataFile)[0] + '.nv'
  
  def fullSanityCheck(self):
    
    headerLength = 2048
    
    rawDimLengths = []
    for i in range(0,self.ndim):
      realDimBytes = self.npts[i]
      blocks = ceil(realDimBytes/float(self.block[i]))
      rawDimBytes = blocks*self.block[i]
      rawDimLengths.append(rawDimBytes)
      #if __debug__:
      #  print i,self.npts[i],self.block[i],realDimBytes,rawDimBytes,blocks
    
    dataSize = 1
    for rawDimLength in rawDimLengths:
        dataSize *= rawDimLength
    
    expectedSize = (dataSize *4) + headerLength   
    #if __debug__:    
    #  print 'expected bytes',expectedSize
    
    
    fileSize = stat(self.dataFile).st_size
    
    if fileSize < expectedSize:
      raise ApiError('nmrview file %s\n is not the correct size (%d bytes; expected %d bytes)' % (self.dataFile,fileSize,expectedSize))
    
        
  # ExternalParams requires this to be defined
  def parseFile(self):
    
    self.nvParToNvNvFile()
    
    try:
      fp = open(self.dataFile, 'rb')
    except IOError, e:
      raise ApiError('nmrview file %s is not acessible' % self.dataFile + str(e))
    
    s = fp.read(self.head)
    if (len(s) < self.head):
      raise ApiError('nmrview file %s\n is shorter than the expected length (%d bytes) of its header' % (self.dataFile,self.head))

    fp.close()
    
    iBuffer = array.array('i')
    fBuffer = array.array('f')
    cBuffer = array.array('c')
    
    iBuffer.fromstring(s)
    fBuffer.fromstring(s)
    cBuffer.fromstring(s)
    
    self.checkMagic(cBuffer)
    
    doByteSwap = self.doByteSwap(cBuffer)
    if doByteSwap:
        iBuffer.byteswap()
        fBuffer.byteswap()
    
    ndim = self.ndim = iBuffer[6]
    self.initDims()
    
    dimBase = 256
    dimRecordSize = 128/4
    dimSizeOffset = 0
    dimBlockSizeOffest = 1
    dimSFOffest = 6
    dimSWOffest = 7
    dimRefPointOffset = 8
    dimRefValueOffset = 9
    

    for i in range(0,ndim):
      self.npts[i] = iBuffer[dimBase + (i*dimRecordSize) + dimSizeOffset]
      self.block[i] = iBuffer[dimBase + (i*dimRecordSize) + dimBlockSizeOffest]
      self.sf[i] = fBuffer[dimBase + (i*dimRecordSize) + dimSFOffest]
      self.sw[i] = fBuffer[dimBase + (i*dimRecordSize) + dimSWOffest]
      self.refpt[i]  = fBuffer[dimBase + (i*dimRecordSize) + dimRefPointOffset]+1
      self.refppm[i] = fBuffer[dimBase + (i*dimRecordSize) + dimRefValueOffset]

    self.swap = doByteSwap
    
    self.parseParFile()
    self.removeTrivialDims()

    self.fullSanityCheck()
    #if __debug__:
    #  for i in range(0,self.ndim):
    #    print i,self.npts[i],self.block[i],self.sf[i], self.sw[i],self.refpt[i], self.refppm[i]
        
    if self.ndim > 1:
      self.guessConsistentNuclei()
  
  def checkMagic(self,bytes):
    numberMagicBytes  = len(self.magicBytes)
    inMagicBytes = list(["%02X" % ord(byte) for byte in bytes[0:numberMagicBytes]])
    inMagicBytes.sort()
        
    if inMagicBytes != self.magicBytes:
      inMagicBytesString = ', '.join(inMagicBytes)
      magicBytesString = ', '.join(self.magicBytes)
      raise ApiError('nmrview file %s\n does not start with the expected magic bytes (%s; in any order) got %s instead' % (self.dataFile,magicBytesString,inMagicBytesString))
    
    
  def guessConsistentNuclei(self):
    ratios = [1.0, 0.4052, 0.2512, 0.1013]
    nuclei = ['1H', '31P', '13C', '15N']
    try:
      sfs  = list(set(self.sf))
      if len(sfs) > 1:  # at least two type of nuclei
        sfs.sort(reverse=True)
      
        #if __debug__:
        #  print 'sfs',sfs
      
        bestSets  = {}

        sfsTop =  sfs[0]
        sfs = sfs[1:]
        while len(ratios) >= self.ndim:
          topNucleus  = nuclei[0]
          ratiosTop = ratios[0]
        
        
          # remaining nuclei and spec frequencies
          ratios = ratios[1:]
          nuclei = nuclei[1:]
        
          #if __debug__:
          #  print '\n'
          #  print 'start',topNucleus, ratiosTop, sfsTop
          #  print 'remaining', nuclei, ratios,sfs
        
        
          targetRatio = sfsTop/ratiosTop
        
          #if __debug__:
          #  print 'target',targetRatio
        
          currentSet  = [(sfsTop,topNucleus)]         
          bestSum  = 0
          #if __debug__:
          #  print 'sfs',sfs
          for sf  in sfs:
        
            #if __debug__:
            #  print '\n\tsf',sf
            
            best  = int(float(sys.maxint))
            for i,pair in enumerate(zip(ratios,nuclei)):
            
              ratio =pair[0]
              nucleus = pair[1]
              error  = abs((sf/ratio) - targetRatio)
            
              #if __debug__:
              #  print '\n\t',i,nucleus,ratio,sf,sf/ratio,error
            
              if error < best:
                  best = error
                  bestIndex  = i
                
                  #if __debug__:
                  #  print '\t','update',bestIndex,best
                  
            currentSet.append((sf,nuclei[bestIndex]))
            bestSum += best
          
          bestSets[bestSum] = currentSet
        
          #if __debug__:
          #  print 'result',best,bestIndex,nuclei[bestIndex],currentSet
      
        keys = bestSets.keys()
        keys.sort()
      
        nucleiLookup = {}
        for matches in bestSets[keys[0]]:
            nucleiLookup[matches[0]]=matches[1]
      
        for i in range(0,self.ndim):
          self.nuc[i] = nucleiLookup[self.sf[i]]

      else: # len(sfs) <= 1: # just assume all 1H
        for i in range(self.ndim):
          self.nuc[i] = '1H'
      
      
    except Exception, e:
      errorSf = []
      for sf in self.sf:
        errorSf.append('%7.3f' % sf)
      print 'error: couldn\'t guess nuclei for frequencies %s' % errorSf
  
  def checkParFileNumberOfFields(self, parFile, lineIndex, fields, expected):
    if len(fields) != expected:
      tooManyFieldsError = 'line %d in nmrview .par file %s has a bad format: incorrect number of fields for %s expected %d\n (line: %s)\n'
      raise ApiError(tooManyFieldsError % (lineIndex +1 ,parFile,fields[0],expected,' '.join(fields)))
  
  def checkParFileAxisIndex(self, parFile, lineIndex, fields):
    axis = atoi(fields[1])
    if axis > self.ndim or  axis < 1:
      badAxisError = 'line %d in nmrview .par %s, has a bad axis index: %d. Permissible axis indices are 1 to %d  (line: %s)\n'
      raise ApiError(badAxisError % (lineIndex + 1, parFile, axis, self.ndim, ' '.join(fields)))
      
  
  # gst work in progress not used
  def parseParFile(self):
    try:
      parFile = splitext(self.dataFile)[0] + '.par'
      
      if isfile(parFile):
        try:
          fp = open(parFile, 'rb')
        except IOError, e:
          raise ApiError('can\'t open nmrview .par file %s: \n' % parFile + str(e))
        
        for i, line in enumerate(fp):
          fields = line.strip().split()
          if i ==0:
            if fields[0] != 'header':
              parFileFormatError = 'file %s doesn\'t appear to be an nmrview .par file the first line doesn\'t start with \'header\' %s'
              raise ApiError(parFileFormatError % (parFile,line))
          elif fields[0] == 'sw':
            axisIndex = atoi(fields[1])-1
            self.checkParFileNumberOfFields(parFile, i, fields, 3)
            self.checkParFileAxisIndex(parFile, i, fields)
            self.sw[axisIndex] = atof(fields[2])
            
          elif fields[0] == 'sf':        
            axisIndex = atoi(fields[1])-1
            self.checkParFileNumberOfFields(parFile, i, fields, 3)
            self.checkParFileAxisIndex(parFile, i, fields)
            self.sf[axisIndex] = atof(fields[2])
      
          elif fields[0] == 'ref':
            axisIndex = atoi(fields[1])-1
            self.checkParFileNumberOfFields(parFile, i, fields, 4)
            self.checkParFileAxisIndex(parFile, i, fields)
            self.refpt[axisIndex] = atof(fields[3])
            self.refppm[axisIndex] = atof(fields[2])
      
    except Exception, e:
        print e
        traceback.print_exc(file=sys.stdout)

if (__name__ == '__main__'):

  import sys
  if (len(sys.argv) != 2):
    print 'Error: correct syntax: <script> <matrix_file>'
    sys.exit(1)

  matrix_file = sys.argv[1]
  params = NmrViewParams(matrix_file)
