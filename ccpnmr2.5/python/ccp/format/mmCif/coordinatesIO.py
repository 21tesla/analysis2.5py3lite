"""
======================COPYRIGHT/LICENSE START==========================

coordinatesIO.py: I/O for mmCIF coordinate file

Copyright (C) 2010 Wim Vranken (European Bioinformatics Institute)

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
from memops.universal.Util import returnInt, returnFloat

from ccp.format.mmCif.generalIO import MMCIFFile
from ccp.format.mmCif.sequenceIO import MMCIFSequenceFile

from ccp.format.general.Constants import defaultSeqInsertCode

#####################
# Class definitions #
#####################

class MMCIFCoordinateFile(MMCIFFile):

  def initialize(self):
  
    self.modelCoordinates = {}
    
    self.chains = []
    self.chainCodes = []

  def read(self, maxNum=999, ignoreResNames=None, verbose=False):
        
    #
    # Read mmCIF file and set up information...
    #

    self.readGeneric(verbose=verbose)

    self.chemCompInfo = self.mmCif.getChemCompInfo()
    
    self.coordinateInfo = self.mmCif.getCoordinateInfo()
    modelNums = self.coordinateInfo.keys()
    modelNums.sort()
    
    #
    # Run sequence code to get information on individual molecules (for CCPN)
    #
    
    sequenceFile = MMCIFSequenceFile(self.name)
    sequenceFile.read(mmCif=self.mmCif, ignoreResNames=ignoreResNames, verbose=False)

    self.residueMappings = sequenceFile.residueMappings

    #
    # Set the coordinates
    #
    
    for modelNum in modelNums:

      self.modelCoordinates[modelNum] = []
    
      for coordinateData in self.coordinateInfo[modelNum]:
      
        (serial, elementType, atomName, cifResLabel, cifChainCode, cifSeqId, entityId, pdbInsertionCode, x, y, z, occupancy, Bfactor, authSeqCode, authResLabel, authChainCode, authAtomName, atomClass) = coordinateData
          
        #
        # Handle residue names and ignore if necessary (e.g. 'HOH')
        #

        resLabel = self.convertResName(cifResLabel)

        if ignoreResNames and (resLabel in ignoreResNames or cifResLabel in ignoreResNames):
          continue
 
        if atomClass == 'HETATM':
          hetFlag = True
        else:
          hetFlag = False
          
        #
        # Set the coordinate
        #
        
        coordinate = MMCIFCoordinate(self,serial,atomName,resLabel,authChainCode,
                                     authSeqCode,x,y,z,occupancy,Bfactor,elementType,hetFlag,
                                     insertionCode=pdbInsertionCode)       
                                   
        self.modelCoordinates[modelNum].append(coordinate)
            
    #
    # Clean up empty models...
    #
    
    for modelNum in self.modelCoordinates.keys():
      if not self.modelCoordinates[modelNum]:
        print "  Warning: model %s has no coordinates - removed." % modelNum
        del(self.modelCoordinates[modelNum])

class MMCIFChain:

  def __init__(self,chainCode,refChainCode):

    self.chainId = chainCode
    self.refChainId = refChainCode

class MMCIFCoordinate:

  def __init__(self,parent,serial,atomName,resName,chainId,
                    seqCode,x,y,z,occupancy,bFactor,atomType,hetFlag, 
                    insertionCode=defaultSeqInsertCode):

    self.parent = parent

    self.serial = serial
    self.atomName = atomName
    self.origAtomName = self.atomName
    
    self.resName = resName
    
    
    # Get this ID from the sequence handling earlier... should always work
    (self.refChainId,self.seqCode) = self.parent.residueMappings[(chainId,seqCode)]

    self.chainId = chainId
    self.refSeqCode = seqCode

    # Make sure this is not None
    if not insertionCode:
      self.insertionCode = defaultSeqInsertCode

    self.x = x
    self.y = y
    self.z = z
    self.occupancy = occupancy
    self.bFactor = bFactor
    self.atomType = atomType
    self.hetFlag = hetFlag

    #
    # Track chain codes...
    #
    
    if not self.refChainId in self.parent.chainCodes:
      self.parent.chainCodes.append(self.refChainId)
      self.parent.chains.append(MMCIFChain(chainId,self.refChainId))
