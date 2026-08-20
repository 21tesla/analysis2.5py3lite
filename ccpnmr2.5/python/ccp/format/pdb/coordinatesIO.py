"""
======================COPYRIGHT/LICENSE START==========================

coordinatesIO.py: I/O for PDB coordinate file

Copyright (C) 2005-2008 Wim Vranken (European Bioinformatics Institute)

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

# Import header reader
from ccp.format.pdb.generalIO import PdbFile

from ccp.format.general.Constants import defaultSeqInsertCode, bondTypes, bioPolymerCodes
from ccp.format.general.Constants import defaultAltLoc

#####################
# Class definitions #
#####################

class PdbCoordinateFile(PdbFile):

  coordFixDict = {
  
    'Asp': [['HD1'],[{'HD1':'HD2','OD1':'OD2','OD2':'OD1'}]],
    'Arg': [['HH2'],[{'HH2':'HH1','1HH1':'1HH2','HH11':'HH21','2HH1':'2HH2','HH12':'HH22','NH1':'NH2','NH2':'NH1'}]],  
    'Glu': [['HE1'],[{'HE1':'HE2','OE1':'OE2','OE2':'OE1'}]],

  }

  def doProblemAtomFix(self,modelNum,fixProblemAtom):
    
    refCoord = self.modelCoordinates[modelNum][-1]
    resName = refCoord.resName
    
    print "  Warning: fixing problems associated with %s atom %s..." % (resName,fixProblemAtom)

    coordFixIndex = self.coordFixDict[resName][0].index(fixProblemAtom)
    curFixDict = self.coordFixDict[resName][1][coordFixIndex]
    
    chainId = refCoord.chainId
    seqCode = refCoord.seqCode
    
    for i in range(len(self.modelCoordinates[modelNum])-1,-1,-1):
      tempPdbCoord = self.modelCoordinates[modelNum][i]
      if tempPdbCoord.seqCode == seqCode and tempPdbCoord.chainId == chainId:
        if tempPdbCoord.atomName in curFixDict:
          #print "Fixed %s to %s for %s,%d!" % (tempPdbCoord.atomName,curFixDict[tempPdbCoord.atomName],tempPdbCoord.resName,tempPdbCoord.seqCode)
          tempPdbCoord.atomName = curFixDict[tempPdbCoord.atomName]
      else:
        break

  def read(self,maxNum = 999,readHeader = True, ignoreResNames = None, verbose = 0):

    # Read in header info - TODO always do this?!? Not necessary if
    # already molsystem available...

    if readHeader:
      self.readHeader(verbose = verbose)

      if not self.chains and not self.hetGroups:
        return
      
    # TODO: create pdb file with no other info if molsys available?
    # have to import PdbFile class!

    # This reader assumes files downloaded from PDB (fixed column format)
    atomRead = 1
    referenceRead = 1
    modelNum = ''
    modelCount = 0
    
    fixProblemAtom = None
    oldSeqCode = None
    chainTerminated = False
    
    self.serialToCoord = {}

    #
    # Read coordinates
    #

    fin = open(self.name, 'rU')
    line = fin.readline()

    while line:

      if self.patt[self.format + 'NewModel'].search(line): # Not necessarily there!! Is empty if so
        modelNum = returnInt(self.patt[self.format + 'NewModel'].search(line).group(1))
        modelCount += 1
        
        if modelCount > 1:
          referenceRead = 0

        if modelCount > maxNum:
          atomRead = 0
        else:
          self.modelCoordinates[modelNum] = []
          chainTerminated = False

      elif atomRead == 1 and self.patt[self.format + 'AllAtom'].search(line):
        
        # Case one model only
        if modelNum == '':
          modelNum = 1
          self.modelCoordinates[modelNum] = []

        serial = line[6:11]
        atomName = line[12:16]
        altLoc = line[16:17]
        resName = line[17:20]
        chainId = line[21:22]
        seqCode = line[22:26]
        insertionCode = line[26:27]
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        occupancy = line[54:60]
        tempFactor = line[60:66]
        segId = line[72:76]
        atomType = line[76:78]
        charge = line[78:80]

        if line[0:6] == 'HETATM':
          hetFlag = 1
        else:
          hetFlag = 0
          
        #
        # Ignore residue names (e.g. 'HOH')
        #

        if ignoreResNames and resName.strip() in ignoreResNames:
          line = fin.readline()
          continue
          
        #
        # Fix problems (see coordFixDict beginning of class)
        #
        
        if self.coordFixDict.has_key(resName):
          if atomName.strip() in self.coordFixDict[resName][0]:
            fixProblemAtom = atomName.strip()
            fixModelNum = modelNum
        
        #
        # Here clean up when seqCode changes... could be a problem but is much faster
        #
        
        if oldSeqCode != seqCode:
          oldSeqCode = seqCode
          if fixProblemAtom:
            self.doProblemAtomFix(fixModelNum,fixProblemAtom)
            fixProblemAtom = None
      
        coordinate = PdbCoordinate(self,serial,atomName,resName,chainId,seqCode,
                                   x,y,z,occupancy,tempFactor,atomType,hetFlag,
                                   insertionCode=insertionCode, altLoc=altLoc,
                                   chainTerminated=chainTerminated)       
        self.modelCoordinates[modelNum].append(coordinate)
        
        if referenceRead:
          self.serialToCoord[coordinate.serial] = coordinate

        # TODO: add in occupancy, tempFactor, ... fields if needed!

      elif self.patt[self.format + 'Bonds'].search(line):
        
        serial = line[6:11]
        bondedSerials = {}
        bondedSerials['covalent'] = [line[11:16],line[16:21],line[21:26],line[26:31]]
        bondedSerials['hydrogen'] = [line[31:36],line[36:41],line[46:51],line[51:56]]
        bondedSerials['salt'] = [line[41:46],line[56:61]]

        centralAtom = self.searchCoordinate(serial,verbose = 0)
        
        if centralAtom:
          for bondType in bondTypes:
            if bondedSerials.has_key(bondType):
              for bondSerial in bondedSerials[bondType]:
                if not self.patt['emptyline'].search(bondSerial):
                  bondedCoordinate = self.searchCoordinate(bondSerial, verbose = 0)
                  centralAtom.setBond(bondType,bondedCoordinate)
                  
      elif self.patt[self.format + 'Termination'].search(line):
        
        chainTerminated = True
      
      line = fin.readline()

    fin.close()    
    
    #
    # Fix remaining problems if any
    #
    
    if fixProblemAtom:
      self.doProblemAtomFix(fixModelNum,fixProblemAtom)
    
    #
    # Clean up empty models...
    #
    
    for modelNum in self.modelCoordinates.keys():
      if not self.modelCoordinates[modelNum]:
        print "  Warning: model %s has no coordinates - removed." % modelNum
        del(self.modelCoordinates[modelNum])
    
    #
    # Set SS bonds and links between het groups/other residues
    #
    
    for (fileAttr,bondType,defaultAtom) in (('ssBonds','disulfide','SG'),('links','link',None)):
    
      bonds = getattr(self,fileAttr)

      if bonds:
      
        atoms = []
        atomCoords = []

        # Combine info on ss bonds in separate variable (so only have to loop
        # over coordinate once)

        for bond in bonds:
          if defaultAtom:
            ((resName1,chainId1,seqCode1,iCode1),(resName2,chainId2,seqCode2,iCode2)) = bond
            atomName1 = atomName2 = defaultAtom
          else:
            ((atomName1,resName1,chainId1,seqCode1,iCode1),(atomName2,resName2,chainId2,seqCode2,iCode2)) = bond

          # Ignoring iCode for now
          atoms.append([(atomName1,resName1,chainId1,seqCode1,iCode1),(atomName2,resName2,chainId2,seqCode2,iCode2)])
          atomCoords.append(['',''])
          
        atomRange = range(len(atoms))

        # Loop over coordinates, check if there's a match
        # Put coordinate memory reference in atomCoords
        for modelNum in self.modelCoordinates:
          for pdbCoordinate in self.modelCoordinates[modelNum]:

            checkData = (pdbCoordinate.atomName,pdbCoordinate.resName,pdbCoordinate.chainId,pdbCoordinate.seqCode,pdbCoordinate.insertionCode)

            for i in atomRange:
              atom = atoms[i]
              for j in range(0,2):
                atomEl = atom[j]
                if atomEl == checkData:
                  atomCoords[i][j] = pdbCoordinate

        # Set the bonds between the coordinates
        for i in range(0,len(atomCoords)):
          atomCoord = atomCoords[i]
          atomInfo = atoms[i]
          
          badAtomCoord = False
          for j in range(0,2):
            if not atomCoord[j]:
              print "  Error: no coordinate for linkBond atom %s.%s.%s.%s%s." % (atomInfo[j][0],atomInfo[j][1],atomInfo[j][2],atomInfo[j][3],atomInfo[j][4])
              badAtomCoord = True
              
          if badAtomCoord:
            continue

          atomCoord[0].setBond(bondType,atomCoord[1])
          atomCoord[1].setBond(bondType,atomCoord[0])

    # Test routine to print out identified chains/hetgroups
    #self.printChainsAndHetGroups(pdb)

  def write(self,addPdbHeader = False,verbose = 0):
    
    #
    # Set up reference info
    #
    
    self.bioPolymerCodeList = []
    for molType in bioPolymerCodes.keys():
      for codeList in bioPolymerCodes[molType]:
        self.bioPolymerCodeList.extend(codeList)
    
    #
    # Start write
    #

    if verbose == 1:
      print "Writing PDB coordinate file %s" % self.name

    fout = open(self.name,'w')

    #
    # TODO: CURRENTLY NO HEADER!
    #
    
    if addPdbHeader:
      fout.write("""CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000      \n""")
    
    #
    # Now start writing coordinates
    # 
    
    multipleModels = 0
    
    if len(self.modelCoordinates) > 1:
      multipleModels = 1
    
    modelNums = self.modelCoordinates.keys()
    modelNums.sort()
    
    oldCoord = None
      
    for modelNum in modelNums:
      
      oldChainId = self.modelCoordinates[modelNum][0].chainId

      if multipleModels:
      
        fout.write("MODEL %8s" % modelNum + self.newline)
        
      for coord in self.modelCoordinates[modelNum]:
      
        if oldChainId != coord.chainId or (oldCoord and oldCoord.hetFlag != coord.hetFlag):

          #
          # TODO: is this correct?
          #

          oldChainId = coord.chainId
          self.writeTerRecord(fout,oldCoord)
      
        atomText = 'ATOM'
      
        if coord.hetFlag:
          atomText = 'HETATM'
        
        atomFormat = " %-3s"
       
        if self.patt['onlydigit'].search(coord.atomName[0]) or len(coord.atomName) == 4:
          
          atomFormat = "%-4s"
      
        lineFormat = "%-6s%5d " + atomFormat + "%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s"
      
        fout.write(lineFormat %
        
                    (atomText,coord.serial,coord.atomName,'',coord.resName,coord.chainId,
                     coord.seqCode,coord.insertionCode,coord.x,coord.y,coord.z,coord.occupancy,
                     coord.bFactor,'',coord.atomType,''))
                     
        fout.write(self.newline)

        oldCoord = coord

      self.writeTerRecord(fout,oldCoord)
      
      if multipleModels:
      
        fout.write("ENDMDL" + self.newline)

      #
      # TODO: Add bonds!! CONECT stuff... which level though?
      #


    fout.close()

  def writeTerRecord(self,fout,coord):
    
    # Only write this for polymer chains - this is hacky but should work, generally
    
    if coord.resName in self.bioPolymerCodeList:
  
      fout.write("TER    9999      %-3s %1s%4d%1s" % (coord.resName,coord.chainId,coord.seqCode,coord.insertionCode) + self.newline)

  def printChainsAndHetGroups(self):

    print self.code

    indent = "    "

    for chain in self.chains:
      print "  Chain: " + chain.chainId
      print indent + "MolId,name: %s,%s" % (chain.molId,chain.molName)
      print indent + "MolSynonyms: %s" % (chain.synonyms)
      print indent + "ecNums, fragment: %s,%s" % (chain.ecNums,chain.fragment)
      print indent + "engineered, mutation: %s,%s" % (chain.engineered,chain.mutation)
      print indent + "biologicalUnit, details: %s,%s" % (chain.biologicalUnit,chain.details)
      print indent + "numRes: %s" % (chain.numRes)
      #print chain.resNames

    for hetGroup in self.hetGroups:
      print "  Hetgroup: " + hetGroup.Id
      print indent + "chainId,seqNum: %s,%s" % (hetGroup.chainId,hetGroup.seqNum)
      print indent + "iCode, numHetAtoms: %s,%s" % (hetGroup.iCode,hetGroup.numHetAtoms)
      print indent + "description, chemicalName: %s,%s" % (hetGroup.description,hetGroup.chemicalName)
      print indent + "synonyms, formula: %s,%s" % (hetGroup.synonyms,hetGroup.formula)

  def searchCoordinate(self,serial,verbose = 1):

    serial = returnInt(serial)
    coordinate = None
    
    if self.serialToCoord.has_key(serial):
      coordinate = self.serialToCoord[serial]

    if not coordinate and verbose:
      print "Error: Serial %s not found" % serial
  
    return coordinate

class PdbCoordinate:

  def __init__(self,parent,serial,atomName,resName,chainId,seqCode,x,y,z,
               occupancy,tempFactor,atomType,hetFlag, 
               insertionCode=defaultSeqInsertCode, altLoc=defaultAltLoc,
               chainTerminated=False):

    self.parent = parent

    self.serial = returnInt(serial)
    self.atomName = atomName.strip()
    self.origAtomName = self.atomName
    
    self.resName = parent.convertResName(resName)
    
    self.chainId = chainId                    # DO NOT strip the chainId - has to be space if empty!
    self.refChainId = chainId
    self.seqCode = returnInt(seqCode)
    self.insertionCode = insertionCode
    self.altLoc = altLoc
    self.x = returnFloat(x)
    self.y = returnFloat(y)
    self.z = returnFloat(z)
    self.occupancy = returnFloat(occupancy,default = 1.0)
    self.bFactor = returnFloat(tempFactor)
    self.atomType = atomType
    self.hetFlag = hetFlag
    
    self.chainTerminated =chainTerminated # To detect HET atoms that are hanging off end of file with same chain code but lower sequence codes...
    
    self.bonds = {}
    
  def setBond(self, bondType, bondedCoordinate):
  
    if bondType not in bondTypes:
      print "Error: Unrecognized bond type %s." % bondType
      return
  
    if not bondedCoordinate:
      print "Error: %s bond cannot be created (coordinate for one partner not found)" % bondType
      return
      
    for tBondType in self.bonds.keys():
      if bondedCoordinate in self.bonds[tBondType]:
        if tBondType == 'covalent': # less informative than link or disulfide!
          self.bonds[tBondType].pop(self.bonds[tBondType].index(bondedCoordinate))
      
    if not self.bonds.has_key(bondType):
      self.bonds[bondType] = []
  
    self.bonds[bondType].append(bondedCoordinate)
