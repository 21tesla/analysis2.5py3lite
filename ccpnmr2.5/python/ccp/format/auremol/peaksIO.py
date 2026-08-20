"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: Peak list I/O for Auremol format files

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

from memops.universal.Util import returnInt, returnFloat

from ccp.format.auremol.generalIO import AuremolFile


class AuremolPeakFile(AuremolFile):

  def initialize(self):
  
    self.peaks = []
    self.specNames = []
    self.numDims = []
    
    self.numPeaks = 0 
    self.numSubPeaks = 0    
    
    self.spectrumInfo = {
    
      'numPoints': [],
      'swsHz':     [],
      'swsPpm':    [],
      'sfs':       [],
      'offsets':   [],
      'filters':   [],     
      'linebroadeningHz':   [],
      'acqMode':   []
      
      }
    
    self.sequence = ""
    
    
    self.headerPars = {
    
    	'EXPERIMENT':             ('specNames',None),
      'DIMENSION':              ('numDims',returnInt),
	    'MIXINGTIME_SEC':         ('mixingTime',returnFloat),
	    'AQUISITIONDELAY':        (None,None),
	    'RESOLUTION_W%d':         ('spectrumInfo.numPoints',returnInt),
	    'RESONANCEFREQUENCY_%d':  ('spectrumInfo.sfs',returnFloat),
	    'SPECTRAL_WIDTH_HZ_%d':   ('spectrumInfo.swsHz',returnFloat),
	    'SPECTRAL_WIDTH_PPM_%d':  ('spectrumInfo.swsPpm',returnFloat),
	    'OFFSET_PPM_%d':          ('spectrumInfo.offsets',returnFloat),
	    'FILTER_%d':              ('spectrumInfo.filters',None),
	    'LINEBROAD_HZ_%d':        ('spectrumInfo.linebroadeningHz',returnFloat),
	    'AQ_MODE_%d':             ('spectrumInfo.acqMode',None),
	    'SEQUENCE':               ('sequence',self.returnSequence)
      
      }

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s peak list %s" % (self.format,self.name)

    #
    # Read info and other info
    #

    headerSection = False
    peakInfo = None
    
    peakNum = 1

    fin = open(self.name, 'rU')
    line = fin.readline()

    while line:
    
      cols = line.split()
      
      if len(cols) == 0 or self.patt['hash'].search(line):
        line = fin.readline()
        continue
      
      #
      # Read information...
      #
      
      if cols[0] == 'HEADER':
        headerSection = True
        
      elif cols[0] == 'END_HEADER':
        headerSection = False

      elif headerSection:
      
        #
        # Header information
        #
        
        # Get rid of any trailing :
        parameter = cols[0].replace(':','')
        
        if parameter[-1] in '123456789':
          parameterKey = parameter[:-1] + '%d'
          dimension = int(parameter[-1])
        else:
          parameterKey = parameter
          dimension = None
        
        if self.headerPars.has_key(parameterKey):
          (mapping,returnFunc) = self.headerPars[parameterKey]
          
          if not mapping and cols[1] != '*':
            print "  Warning: not handling parameter %s, but has value '%s'!" % (parameter,cols[1])
            
          elif mapping:
            
            fileVars = mapping.split('.')
            
            if returnFunc:
              value = returnFunc(cols[1])
            else:
              value = cols[1]
              
            if len(fileVars) == 1:
            
              # Find out if this is a list
              if hasattr(self,fileVars[0]):
                attrValue = getattr(self,fileVars[0])
              else:
                attrValue = None
              
              # Set the info correctly
              if type(attrValue) == type([]):
                attrValue.append(value)
              else:
                setattr(self,fileVars[0],value)
            
            else:
              currentVar = getattr(self,fileVars[0])[fileVars[1]]             
            
              if dimension:
                currentVar.append(value)
              else:
                currentVar = value

      elif cols[0] == 'Total':
      
        self.numPeaks = int(cols[-1])
      
      elif cols[0] == 'Including':

        self.numSubPeaks = int(cols[-1])

      else:
        
        #
        # Peak level information
        #

        parameter = cols[0].replace(':','')
        
        if parameter == 'PEAKLABEL':
          
          if peakInfo:
            # Make peak here
            self.peaks.append(AuremolPeak(peakNum,peakInfo))            
            peakNum += 1
            
          peakInfo = {'s': [],'ppm': [], 'width': []}
          peakInfo['label'] = line.replace(cols[0],'').strip()
          
        elif parameter in ('COMMENT','FEATURES','COMPOUNDS'):
          peakInfo[parameter.lower()] = line.replace(cols[0],'').strip()
          
        elif parameter == 'PEAKDESCRIPTION':
        
          headerInfo = fin.readline().strip().split()
          values =     fin.readline().strip().split()

          for dim in range(1,self.numDims[0]+1):
            for (parameter,returnFunc) in (('s',returnInt),('ppm',returnFloat),('width',returnFloat)):
              dimParameter = "%s%d" % (parameter,dim)
              if dimParameter in headerInfo:
                peakInfo[parameter].append(returnFunc(values[headerInfo.index(dimParameter)]))
        
          for (parameter,returnFunc) in (('intensity',returnInt),('volume',returnInt),('prob',returnFloat)):
            if headerInfo.count(parameter):
              peakInfo[parameter] = returnFunc(values[headerInfo.index(parameter)])
                    
        elif parameter == 'SUBPEAKS':
          peakInfo['subPeaks'] = returnInt(cols[1])
          # TODO HERE HAVE TO READ NEXT LINES, BASED ON NUMBER?!?
        
        elif parameter == 'VOLERR':        
          peakInfo['volumeError'] = returnFloat(cols[1])
                      
      line = fin.readline()

    # Make peak here
    self.peaks.append(AuremolPeak(peakNum,peakInfo))
    
    fin.close()


  def returnSequence(self,seqString):
  
    seqStringElements = seqString.split('_')
  
    return seqStringElements[-1]

class AuremolPeak:

  """
{ 'subPeaks': 0,  'prob': 1.0}
  """

  def __init__(self,num,peakInfo):
    
    self.num = num
    
    for parameter in ('label','comment','volume','intensity','features','compounds','ppm','width'):
      if peakInfo.has_key(parameter):
        value = peakInfo[parameter]
      else:
        value = None
      setattr(self,parameter,value)
    
    # Assume these are assignments
    parameter = 's'
    if peakInfo.has_key(parameter) and peakInfo[parameter]:
      self.assign = peakInfo[parameter]
    else:
      self.assign = None
    
    if hasattr(self,'volume'):
      parameter = 'volumeError'
      if peakInfo.has_key(parameter) and peakInfo[parameter]:
        self.volumeError = peakInfo[parameter] / 100.0
      else:
        self.volumeError = None
    
    parameter = 'prob'
    if peakInfo.has_key(parameter) and peakInfo[parameter]:
      self.probability = peakInfo[parameter]
    else:
      self.probability = None
    

