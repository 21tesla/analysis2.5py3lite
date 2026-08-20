"""
======================COPYRIGHT/LICENSE START==========================

processingParsIO.py: Azara processing parameter file writer

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

import os

from memops.universal.Io import getTopDirectory, normalisePath, splitPath, joinPath
from memops.universal.Constants import dirsep
from memops.universal.Util import returnInt

from ccp.format.general.Util import standardNucleusName
from ccp.format.general.formatIO import GenericProcessingParsFile
from ccp.format.azara.generalIO import AzaraGenericFile

#####################
# Class definitions #
#####################

class AzaraProcessingParsFile(AzaraGenericFile,GenericProcessingParsFile):

  def initialize(self):
    
    self.scriptName = 'processingParsIO.py'
    
    self.allPars = {
    
      'mainPars': [
    
        ('ndim','ndim',2,0),
        ('file','file',2,''),
        ('head','head',2,0),
        ('int','integer',1,False),
        ('swap','swap',1,False),
        ('big_endian','big_endian',1,None),
        ('little_endian','big_endian',1,None),
        ('varian','varian',-1,None)
    
       ],
       
      'dimPar': [
      
        ('dim','dim',2,-1)
        
       ],
    
      'dimPars': [
      
        ('npts','npts',2,None),
        ('block','block',2,None),
        ('sw','sw',2,1000.0),
        ('sf','sf',2,1000.0),
        ('refppm','refppm',2,0.0),
        ('refpt','refpt',2,1.0),
        ('nuc','nuc',2,'1H'),
        ('params','sampledValues',-1,None),
    
       ]
    }

    self.initGeneric()
                  
  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading Azara processing parameters file %s" % self.name

    try:
      fin = open(self.name, 'rU')
    except IOError, e:
      raise e

    lines = fin.readlines()

    fin.close()

    for line in lines:

      n = line.find('!')
      if (n >= 0):
        line = line[:n]

      line = line.strip()

      if (line):

        fields = line.split()
        self.parseFields(fields)

  def parseFields(self, fields):
  
    #
    # This is taken almost directly from Wayne's AzaraParams.py file
    #

    key = fields[0]
    assert key in self.allParKeys, 'key = %s' % key
    
    parType = None
    numFields = None

    for tempParType in self.allPars.keys():
      if key in self.parKeys[tempParType]:
        for parInfo in self.allPars[tempParType]:
          if key == parInfo[0]:
            (origName,analysisName,numFields,default) = parInfo
            parType = tempParType
            break
            
      if parType:
        break
      
    if not numFields:
      return

    if numFields == -1: # Basically any number of parameters!
      valueList = []
      for field in fields[1:]:
        try:
          x = int(field)
        except:
          x = float(field)
        valueList.append(x) 
      if key == 'params':
        dim = self.aPars.get('dim')
        if dim is None:
          raise Exception("'params' set before 'dim'")
        npts = self.aPars.get('npts')
        if not npts or len(npts) <= dim:
          raise Exception("'params' set before 'npts'")
        npts = npts[dim]
      if len(valueList) != npts:
        raise Exception("'params' set before 'npts'")
        
      if parType == 'mainPars':
        self.setValue(parInfo, tuple(valueList))
      elif parType == 'dimPars':
        assert self.aPars['dim'] >= 0, '"%s" must be after "dim"' % key
        dim = self.aPars['dim']
        self.setListValue(parInfo,dim,tuple(valueList))
      
      # TBD: should this be done with a setValue somehow??
      self.aPars['dimType'][dim] = self.dimTypeSampled
      
      return

    assert len(fields) == numFields, 'len(fields) = %d, should be %d' % (len(fields),numFields)

    if numFields == 1:
      if (key == 'little_endian'):  # Int is now taken care of by dictionary!
        self.setValue(parInfo,False)
      else:
        self.setValue(parInfo, True)
       
    elif numFields == 2:

      value = fields[1]

      if parType == 'mainPars':
      
        if (key == 'ndim'):
          assert self.aPars['ndim'] == 0, 'repeated "ndim"'
          self.setValue(parInfo, returnInt(value))
          assert self.aPars['ndim'] > 0, 'ndim <= 0'
          self.initDims()
        elif (key == 'file'):
          assert not self.aPars['file'], 'repeated "file"'
          path = normalisePath(value)
          if (path and (not os.path.isabs(path)) and os.path.isabs(self.name)):
            path = joinPath(splitPath(self.name)[0], path)
          self.setValue(parInfo,path)
        elif (key == 'head'):
          assert self.aPars['head'] == 0, 'repeated "head"'
          self.setValue(parInfo, 4*int(value))
          assert self.aPars['head'] >= 0, 'head < 0'
          
      elif parType == 'dimPar':
      
        if (key == 'dim'):
          assert self.aPars['ndim'] > 0, '"dim" must be after "ndim"'
          dim = self.aPars['dim']
          self.setValue(parInfo, int(value) - 1)
          assert self.aPars['dim'] != dim, 'repeated dim %d' % (dim+1)
          assert self.aPars['dim'] >= 0, 'dim = %d <= 0' % (self.aPars['dim']+1)
          assert self.aPars['dim'] < self.aPars['ndim'], 'dim = %d > ndim = %d' % (self.aPars['dim']+1, self.aPars['dim'])
          
      elif parType == 'dimPars':

        assert self.aPars['dim'] >= 0, '"%s" must be after "dim"' % key
        dim = self.aPars['dim']
        # can no longer do below because initDims() sets up default values
        #assert getattr(self, key)[dim] is None, '%s repeated for dim %d' % (key, dim+1)
        if (key  in ('npts','block')):
          v = int(value)
          assert v > 0, '%s <= 0 for dim %d' % (key, dim+1)
          self.setListValue(parInfo,dim,v)
        elif (key == 'nuc'):
          self.setListValue(parInfo,dim,standardNucleusName(value))
        else:
          self.setListValue(parInfo,dim,float(value))
  
  def writePar(self,parInfo,infoText,value,indent = "",newlines = 2):
    
    if value != None and value != False and value != '':
    
      if parInfo[2] == -1:
        newValue = ""
        for tval in value:
          newValue += "%s " % tval
        value = newValue
      elif parInfo[2] == 1:
        value = ''
        
      self.fout.write("%s%-12s%-20s !%s%s" % (indent,parInfo[0],str(value),infoText,self.newline * newlines))

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing azara processing pars file %s" % self.name

    numDim = self.fPars['numDim']
    
    self.fout = open(self.name,'w')

    #
    # Main parameters
    #

    for parInfo in self.allPars['mainPars']:
      if parInfo[1]:
        (parName,infoText) = self.parNameMapping[parInfo[1]]
        if self.fPars.has_key(parName):
          self.writePar(parInfo,infoText,self.fPars[parName])
          del (self.fPars[parName])

    #
    # Pars per dimension
    # 

    for dim in range(0,numDim):
    
      for parInfo in self.allPars['dimPar']:
        if parInfo[1]:
          (parName,infoText) = self.parNameMapping[parInfo[1]]
          if self.fPars.has_key(parName):
            self.writePar(parInfo,infoText,dim + 1)

      for parInfo in self.allPars['dimPars']:
        if parInfo[1]:
          (parName,infoText) = self.parNameMapping[parInfo[1]]
          if self.fPars.has_key(parName):
            self.writePar(parInfo,infoText,self.fPars[parName][dim],indent = '  ',newlines = 1)

      self.fout.write(self.newline)

    self.fout.write(self.newline)

    self.fout.close()

    if verbose == 1:
      print self.scriptName + " finished..."

###################
# Main of program #
###################
    
if __name__ == "__main__":

  files = [
  
      #('data/den/1/1/1/tocsy_HCCH.hhc.1/project.conversion/azara/processed/spectrum.par','local/cs.test1'),
      #('data/den/1/1/1/tocsy_HCCH.hhc.1/project.conversion/azara/autoproc.par','local/cs.test2')
      ('../../reference/ccpNmr/katCam/oxime/30.spc.par','local/cs.test3'),
      ('../../reference/ccpNmr/johnCambridge/051107/drx600/50.spc.par','local/cs.test4'),
      ('../../reference/azara/edl387_5_Sampled.spc.par','local/cs.test5'),
           
          ]

  
  for (inFile,outFile) in files:
    
    procParFile = AzaraProcessingParsFile(os.path.join(getTopDirectory(),inFile))

    procParFile.read(verbose = 1)
    
    aparsKeys = procParFile.aPars.keys()
    aparsKeys.sort()
    print 'apars'
    for aparsKey in aparsKeys:
      print "  %-20s: %s" % (aparsKey,str(procParFile.aPars[aparsKey]))
    
    fparsKeys = procParFile.fPars.keys()
    fparsKeys.sort()
    print 'fpars'
    for fparsKey in fparsKeys:
      print "  %-20s: %s" % (fparsKey,str(procParFile.fPars[fparsKey]))

    procParFile.name = outFile
 
    procParFile.write(verbose = 1)
  
  
