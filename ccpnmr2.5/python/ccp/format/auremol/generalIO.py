
"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information for Auremol format files

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
from memops.universal.Util import returnInt, returnFloat
from memops.universal.Io import getTopDirectory

from ccp.format.general.formatIO import FormatFile
from ccp.format.general.Constants import defaultMolCode

#####################
# Class definitions #
#####################

class AuremolFile(FormatFile):

  def setGeneric(self):
    
    self.format = 'auremol'
    self.defaultMolCode = defaultMolCode


class AuremolGenericFile(AuremolFile):
  
  def initialize(self):
  
    pass
  
  def setSubSection(self,name,value = None):
  
    self.subSectionType = name
  
  def setReadSection(self,sectionSearch):
  
    readSection = sectionSearch.group(1).strip()
    
    if readSection[:18] == 'compounddefinition':

      cols = readSection.split()
      readSection = cols[0]
      compoundIndex = returnInt(cols[1]) - 1
      self.compound = self.compounds[compoundIndex]
        
    return readSection
  
  def setSection(self,readSection,cols,subSectionSearch, value = None):
  
    pass

  def returnColumnValue(self,value):
  
    if value == '-':
      value = None
      
    return value
    
  def returnFloatValue(self,value):
  
    value = self.returnColumnValue(value)
    
    if value:
      value = returnFloat(value)
    
    return value
    
  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s molecular information file %s" % (self.format,self.name)

    #
    # Read info and other info
    #

    readSection = None

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
      
      sectionSearch = self.patt[self.format + 'Section'].search(line)
      subSectionSearch = self.patt[self.format + 'SubSection'].search(line)

      if len(cols) == 2 and subSectionSearch:
      
        self.setSubSection(subSectionSearch.group(1),value = cols[1])
        
      elif cols[0] == 'end_section':
      
        readSection = None

      elif sectionSearch:
        
        readSection = self.setReadSection(sectionSearch)
        self.colDefs = []
      
      elif readSection:
      
        if subSectionSearch:
          self.setSubSection(subSectionSearch.group(1))

        elif self.subSectionType and cols[0] == 'END_' + self.subSectionType:
          self.colDefs = []
          self.subSectionType = None

        elif self.patt['underscore'].search(line):
      
          self.colDefs.append(cols[0])
          
        elif cols[0][:4] == 'END_':
        
          print "Warning: double END code %s" % cols[0]
        
        else:
        
          self.setSection(readSection,cols,subSectionSearch)
            
      line = fin.readline()

    fin.close()

class AuremolMetaFile(AuremolGenericFile):
  
  def initialize(self):
  
    self.compounds = []
    self.parameters = {}
    
    self.parameterInfo = {
    
      '_Variable_type':       ('type',None),
      '_Variable_value':      ('value',self.returnColumnValue),
      '_Variable_value_unit': ('unit',self.returnColumnValue)}
    
    self.compoundInfo = {
    
      '_Mol_label':                 ('label',None),
      '_Compound_file':             ('compFile',None),
      '_Concentration_value':       ('concentration',self.returnFloatValue),
      '_Concentration_value_units': ('concentrationUnit',self.returnColumnValue),
      '_Isotopic_labeling':         ('labelling',self.returnColumnValue)}

    self.shiftInfo = {
    
      '_Residue_seq_code':          ('seqCode',returnInt),
      '_Atom_num_code':             ('atomNum',returnInt),
      '_Atom_alias':                ('atomSetName',self.returnColumnValue),
      '_Chem_shift_value':          ('value',self.returnFloatValue),
      '_Chem_shift_value_error':    ('error',self.returnFloatValue),
      '_Chem_shift_ambiguity_code': ('ambCode',returnInt),
      '_Atom_class':                ('class',returnInt),
      '_Linewidth':                 ('linewidth',self.returnFloatValue)}

  def setSection(self,readSection,cols,subSectionSearch,value = None):

    if readSection == 'probedefinition':

      if self.subSectionType == 'COMPOUNDS':

        self.compounds.append(AuremolCompoundInfo(self,cols,self.compoundInfo))

      elif self.subSectionType == 'PARAMETER':

        self.parameters[cols[0]] = AuremolParameter(self,cols,self.parameterInfo)

    elif readSection == 'compounddefinition':

      if self.subSectionType == 'CLASSDEF':
      
        pass
      
      elif self.subSectionType == 'SHIFTS':
      
        self.compound.shifts.append(AuremolShiftInfo(self,cols,self.shiftInfo))
      
      elif self.subSectionType == 'J_COUPL':
      
        pass
        
      else:
        
        pass

class AuremolCompFile(AuremolGenericFile):
  
  def initialize(self):
  
    self.compound = None
    self.atoms = []
    self.atomDict = {}
    self.couplings = []
    
    self.atomInfo = {
    
      '_Residue_seq_code': ('seqCode',returnInt),
      '_Atom_num_code':    ('atomNum',returnInt),
      '_Residue_label':    ('resLabel',None),
      '_Atom_name':        ('atomName',None),
      '_Atom_type':        ('atomType',None),
      '_Atom_alias':       ('atomSetName',self.returnColumnValue),
      '_Atom_equivalent':  ('equivalentAtoms',self.returnEquivalentAtoms),
      '_Atom_CSA':         ('csa',self.returnColumnValue)}
    
    self.couplingInfo = {
    
      '_Coupled_atom1':    ('atom1Key',None),
      '_Coupled_atom2':    ('atom2Key',None),
      '_Coupling_value':   ('value',self.returnColumnValue),
      '_Angle_name':       ('name',self.returnColumnValue)}

  def returnEquivalentAtoms(self,value):
  
    values = self.returnColumnValue(value)
    
    if values:
      values = values.split('/')
      
    return values

  def setSubSection(self,name,value = None):
  
    self.subSectionType = name
  
    if name == 'COMPOUND':
      self.compound = value
  
  def setSection(self,readSection,cols,subSectionSearch):

    if readSection == 'sequenzdefinition':

      self.atoms.append(AuremolAtomInfo(self,cols,self.atomInfo))
      self.atomDict[self.atoms[-1].atomKey] = self.atoms[-1]

    elif readSection == 'bond':

      refAtom = self.atomDict[cols[0]]
      refAtom.addBonds(cols[1:],self.atomDict)

    elif readSection == 'j-coupling':
    
      #
      # Ignoring define for now...
      #

      if self.subSectionType == 'DEFINE':

        pass

      elif self.subSectionType == 'DATA':

        self.couplings.append(AuremolCouplingInfo(self,cols,self.couplingInfo))

class AuremolGenericInfo:

  def __init__(self,parent,cols,parentInfo):
  
    self.parent = parent
  
    for i in range(0,len(parent.colDefs)):
    
      colTag = parent.colDefs[i]
      value = cols[i]
      
      if parentInfo.has_key(colTag):
      
        func = parentInfo[colTag][1]

        if func:
          value = func(value)

        setattr(self,parentInfo[colTag][0],value)
    
    self.setSpecificInfo()

class AuremolAtomInfo(AuremolGenericInfo):

  def setSpecificInfo(self):
  
    self.bondedAtoms = []
    self.atomKey = "%d,%d" % (self.seqCode,self.atomNum)
    
  def addBonds(self,atomKeys,atomDict):
  
    for atomKey in atomKeys:
      atomKey = self.parent.returnColumnValue(atomKey)
      if atomKey:
        self.bondedAtoms.append(atomDict[atomKey])

class AuremolCouplingInfo(AuremolGenericInfo):

  def setSpecificInfo(self):
    
    self.atoms = [self.parent.atomDict[self.atom1Key],self.parent.atomDict[self.atom2Key]]

class AuremolCompoundInfo(AuremolGenericInfo):

  def setSpecificInfo(self):
  
    self.shifts = []
    
class AuremolParameter(AuremolGenericInfo):

  def setSpecificInfo(self):
  
    for parameterType in ['pH','Temperature','Pressure']:
      
      if hasattr(self,parameterType):
        setattr(self,parameterType,returnFloat(getattr(self,parameterType)))
        
class AuremolShiftInfo(AuremolGenericInfo):

  def setSpecificInfo(self):
  
    self.atomKey = "%d,%d" % (self.seqCode,self.atomNum)
        
###################
# Main of program #
###################

if __name__ == "__main__":  

  file = os.path.join(getTopDirectory(),'../../reference/auremol/csp_korr2_iupac.RIGHTCOMPFILE.meta')
  auremolFile = AuremolMetaFile(file)
  
  auremolFile.read()

  for compound in auremolFile.compounds:
  
    print compound.label, compound.compFile
    
    for cs in compound.shifts:
      print "%s: %.3f| " % (cs.atomKey, cs.value),
    print
  
  for par in auremolFile.parameters.keys():
    print par, auremolFile.parameters[par].value, auremolFile.parameters[par].unit

  """
  file = os.path.join(getTopDirectory(),'../../reference/auremol/csp.comp')
  auremolFile = AuremolCompFile(file)
  
  auremolFile.read()

  for atomInfo in auremolFile.atoms:
  
    print atomInfo.atomKey, atomInfo.atomName, atomInfo.bondedAtoms
  """
