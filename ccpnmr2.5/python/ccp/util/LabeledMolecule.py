try:
  from memops.gui.MessageReporter import showWarning
except:
  from memops.universal.MessageReporter import showWarning

# TBD
#
# Still need equivalents of getSchemeAtomFractions() and
# getSchemeAtomPairFractions() that work with labelledMixtures
# Presumable this will use the sample that is associated with
# a given experiment (or peak)
#

# # # # # # # # Setup functions # # # # # # # # # 


def getLabeledMolecule(molecule):
  """Descrn: Make, or find a top-level object to store molecule labelling
             patterns and mixtures fro the input molecule.
     Inputs: Molecule.Molecule
     Output: LabeledMolecule.LabeledMolecule
  """      
  
  project = molecule.memopsRoot
  
  topObj = project.findFirstLabeledMolecule(name=molecule.name)
  
  if not topObj:
    topObj = project.newLabeledMolecule(name=molecule.name)
  
  return topObj

def getDefaultMolPattern(topObj):
  """Descrn: Make a MolLabel object to store a single molecule labelling
             pattern.
     Inputs: LabeledMolecule.LabeledMolecule
     Output: LabeledMolecule.MolLabel
  """      
  
  molLabel = topObj.findFirstMolLabel()
  
  if not molLabel:
    molLabel = topObj.newMolLabel()

  return molLabel
  
    
def clearResidueScheme(molLabel, molResidue):
  """Descrn: Remove all labelling information for a residue in a specified
             molLabel (molecule labelling pattern).
     Inputs: LabeledMolecule.MolLabel, Molecule.molResidue
     Output: None
  """

  resLabel = molLabel.findFirstResLabel(resId=molResidue.serial)

  if resLabel:
    resLabel.delete()
                               
    
def setResidueScheme(molLabel, molResidue, scheme):
  """Descrn: Link a given input residue to a standard set of itotopomers
             specified in a given labelling scheme. The residue labelling
             is only set within a given molLabel instance and so may 
             correspond to only one component of a mixture.
     Inputs: LabeledMolecule.MolLabel, Molecule.molResidue,
             ChemCompLabel.LabelingScheme
     Output: LabeledMolecule.ResLabel
  """

  chemComp = molResidue.chemComp
  chemCompLabel = scheme.findFirstChemCompLabel(chemComp=chemComp)

  if (not chemCompLabel) or not chemCompLabel.isotopomers:
    data = (chemComp.molType, chemComp.ccpCode, scheme.name)
    msg = 'No isotopomer information found for %s %s in labelling scheme %s' % data
    showWarning('Failure', msg)
    return

  resLabel = molLabel.findFirstResLabel(resId=molResidue.serial)

  if resLabel:
    for oldFraction in resLabel.resLabelFractions:
      oldFraction.delete()
    
  else:  
    resLabel = molLabel.newResLabel(resId=molResidue.serial)

  # Weight 1.0 is he default in any case
  # Isotopomer serial zero means standard set  

  resLabel.newResLabelFraction(schemeName=scheme.name,
                               weight=1.0, isotopomerSerial=0)
                               
  return resLabel


# # # # # # # # QUERY FUNCTIONS # # # # # # # # #

def getExperimentAtomFractions(experiment, atom):
  """Descrn: Get the isotope proportions for an atom
             in a given experiment (label mixture)
     Inputs: Nmr.Experiment, MolSystem.Atom
     Output: Dict of Word:Float - IsotopeCode:fraction
  """
  
  labelledMixtures = experiment.labeledMixtures
  molResidue = atom.residue.molResidue
  molecule = molResidue.molecule
      
  fractionDict = {}
  
  for mixture in labelledMixtures:
    if mixture.labeledMolecule.molecule is molecule:
      resId = molResidue.serial
      atomName = atom.name
      fractionDict = singleAtomFractions(mixture, resId, atom.name)
      break

  return fractionDict

def getSchemeAtomFractions(scheme, atom):
  """Descrn: Get the isotope proportions for a given molecular
             system atom within a labelling scheme.
     Inputs: ChemCompLabel.LabelingScheme, MolSystem.Atom
     Output: Dict of Word:Float - IsotopeCode:fraction
  """

  fractionDict = {}

  residue = atom.residue
  ccpCode = residue.ccpCode
  molType = residue.molResidue.molType
  
  atomName = atom.name
  subType  = atom.chemAtom.subType
  
  chemCompLabel = scheme.findFirstChemCompLabel(ccpCode=ccpCode,
                                                molType=molType)
  
  if chemCompLabel:
    isotopomers  = chemCompLabel.isotopomers
    fractionDict = getIsotopomerSingleAtomFractions(isotopomers,
                                                    atomName,
                                                    subType)
  else:
    warnChempCompLabelFailure(ccpCode,molType,scheme)

  return fractionDict


# To work out atom-level pair incorportaions etc...

def getExperimentAtomPairFractions(experiment, atom1, atom2):
  """Descrn: Get the combined isotope proportions for a given pair of molecular
             system atoms for a given experiment (label mixture).
     Inputs: Nmr.Experiment, MolSystem.Atom, MolSystem.Atom
     Output: Dict of Word:Float - IsotopeCode:fraction
  """
  residue1 = atom1.residue
  residue2 = atom2.residue
  molResidue1 = residue1.molResidue
  molResidue2 = residue2.molResidue
  molecule1 = molResidue1.molecule
  molecule2 = molResidue2.molecule
  
  resIds = (molResidue1.serial, molResidue2.serial)
  atomNames = (atom1.name, atom2.name)
  labelledMixtures = experiment.labeledMixtures
      
  fractionDict = {}
  
  if molecule1 is molecule2:
    for mixture in labelledMixtures:
      if mixture.labeledMolecule.molecule is molecule1:
        fractionDict = atomPairFractions(mixture, resIds, atomNames)
        break
  
  else:
    fractionDict1 = getExperimentAtomFractions(experiment, atom1)
    fractionDict2 = getExperimentAtomFractions(experiment, atom2)
    
    for isotope1 in fractionDict1.keys():
      fraction1 = fractionDict1[isotope1]
      
      for isotope2 in fractionDict2.keys():
        fraction2 = fractionDict2[isotope2]
        key = (isotope1, isotope2)
        fractionDict[key] = fraction1 * fraction2

  return fractionDict

def getSchemeAtomPairFractions(scheme, atom1, atom2):
  """Descrn: Get the combined isotope proportions for a given pair of molecular
             system atoms within a labelling scheme.
     Inputs: ChemCompLabel.LabelingScheme, MolSystem.Atom, MolSystem.Atom
     Output: 
  """

  fractionDict = {}

  residue1 = atom1.residue
  ccpCode1 = residue1.ccpCode
  molType1 = residue1.molResidue.molType
  
  atomName1 = atom1.name
  subType1  = atom1.chemAtom.subType

  residue2 = atom2.residue
  ccpCode2 = residue2.ccpCode
  molType2 = residue2.molResidue.molType
  
  atomName2 = atom2.name
  subType2  = atom2.chemAtom.subType
  
  if residue1 is residue2:
    # Abundances are correlated by isotopomer site incorporation
    chemCompLabel = scheme.findFirstChemCompLabel(ccpCode=ccpCode1,
                                                  molType=molType1)
 
    if chemCompLabel:
      atomNames    = (atomName1, atomName2)
      subTypes     = (subType1, subType2)
      isotopomers  = chemCompLabel.isotopomers
      fractionDict = getIsotopomerAtomPairFractions(isotopomers, atomNames, subTypes)
    else:
      warnChempCompLabelFailure(ccpCode1,molType1,scheme)
 
  else:
    # Abundances are independent, and presumably equilibrium, for the two sites
    chemCompLabel1 = scheme.findFirstChemCompLabel(ccpCode=ccpCode1,
                                                   molType=molType1)
    chemCompLabel2 = scheme.findFirstChemCompLabel(ccpCode=ccpCode2,
                                                   molType=molType2)
 
    if chemCompLabel1 and chemCompLabel2:
      isotopomers1  = chemCompLabel1.isotopomers
      fractionDict1 = getIsotopomerSingleAtomFractions(isotopomers1,
                                                       atomName1,
                                                       subType1)
      isotopomers2  = chemCompLabel2.isotopomers
      fractionDict2 = getIsotopomerSingleAtomFractions(isotopomers2,
                                                       atomName2,
                                                       subType2)
 
      for isotope1 in fractionDict1.keys():
        fraction1 = fractionDict1[isotope1]
        
        for isotope2 in fractionDict2.keys():
          fraction2 = fractionDict2[isotope2]
          
          key = (isotope1, isotope2)
          
          fractionDict[key] = fraction1 * fraction2
 
    elif not chemCompLabel1:
      warnChempCompLabelFailure(ccpCode1,molType1,scheme)
 
    elif not chemCompLabel2:
      warnChempCompLabelFailure(ccpCode2,molType2,scheme)
     

  return fractionDict
  

def warnChempCompLabelFailure(ccpCode,molType,scheme):
  """Descrn: Display a general warning for not being able to find a chemComplLabel
     Inputs: Word, Word, ChemCompLabel.LabelingScheme, MolSystem.Atom
     Output: None
  """

  data = (ccpCode,molType,scheme.name)
  msg  = 'Could not find chemp comp labelling information for %s(%s) in scheme %s'
  showWarning('Failure',msg % data)
  

def singleAtomFrac(labeledMixture, resId, atName, isotopeCode):
  """get the isotope fraction for a given isotope at a given position 
  
  labeledMixture:  LabeledMixture object
  resid, atName: residue serial and atom name for atom
  isotopeCode: isotope code you want fraction for
  Returns isotopeCode:fraction dictionary with fractions normalised to 1.0
  """
  return singleAtomFractions(labeledMixture, resId, atName).get(isotopeCode)


def atomPairFrac(labeledMixture, resIds, atNames, isotopeCodes):
  """get the isotope pair : fraction dictionary for a labeledMixture
  
  labeledMixture:  LabeledMixture object
  resIds: length-two tuple of residue serials
  atNames: length-two tuple of atom names
  isotopeCodes: length-two tuple of isotopeCodes
  
  Returns (isotopeCode1, isotopeCode2):fraction dictionary 
  with fractions normalised to 1.0
  """
  
  if len(isotopeCodes) != 2:
    raise("Error: length of isotopeCodes %s must be 2" % isotopeCodes)
  
  
  return atomPairFractions(labeledMixture, resIds, atNames).get(
          tuple(isotopeCodes)
         )


def singleAtomFractions(labeledMixture, resId, atName):
  """get the isotopeCode:fraction dictionary for a labeledMixture
  
  labeledMixture:  LabeledMixture object
  resid, atName, residue serial and atom name for atom
  Returns isotopeCode:fraction dictionary with fractions normalised to 1.0
  """
  
  # set up
  molResidue = labeledMixture.labeledMolecule.molecule.findFirstMolResidue(
               serial=resId)
  if molResidue is None:
    return None
  
  chemAtom = molResidue.chemCompVar.findFirstChemAtom(name=atName)
  if chemAtom is None:
    return None
  elementName = chemAtom.elementSymbol
  subType = chemAtom.subType
  
  # get MolLabelFractions
  useLabel = labeledMixture.averageComposition
  if not useLabel:
    molLabelFractions = list(labeledMixture.molLabelFractions)
    if len(molLabelFractions) == 1:
      useLabel = molLabelFractions[0]
  
  if useLabel:
    # average composition
    resLabel = useLabel.molLabel.findFirstResLabel(resId=resId)
    result = molLabelFractionsDict(resLabel, atName, subType, elementName)
    
  else:
    # average over multiple molLabels
    molWeightSum = sum([x.weight for x in molLabelFractions])
    result = {}
  
    for mlf in molLabelFractions:
      resLabel = mlf.molLabel.findFirstResLabel(resId=resId)
      partResult = molLabelFractionsDict(resLabel, atName, subType, elementName)
      for key, val in partResult.items():
        result[key] = result.get(key, 0.0) + val * mlf.weight / molWeightSum
    
  #
  return result
  

def atomPairFractions(labeledMixture, resIds, atNames,):
  """get the isotope pair : fraction dictionary for a labeledMixture
  
  labeledMixture:  LabeledMixture object
  resIds: length-two tuple of residue serials
  atNames: length-two tuple of atom names
  
  Returns (isotopeCode1, isotopeCode2):fraction dictionary 
  with fractions normalised to 1.0
  """
  
  result = {}

  if len(resIds) != 2:
    raise("Error: length of resIds %s must be 2" % resIds)
  if len(atNames) != 2:
    raise("Error: length of atNames %s must be 2"
          % atNames)
  
  # calculate starting parameters
  elementNames = []
  subTypes = []
  for ii in (0,1):
    molResidue = labeledMixture.labeledMolecule.molecule.findFirstMolResidue(
                 serial=resIds[ii])
    if molResidue is None:
      return None
 
    chemAtom = molResidue.chemCompVar.findFirstChemAtom(name=atNames[ii])
    if chemAtom is None:
      return None
    elementName = chemAtom.elementSymbol
    elementNames.append(elementName)
    subTypes.append(chemAtom.subType)
  
  # get MolLabelFractions
  avLabel = labeledMixture.averageComposition
  if avLabel:
    molLabelFractions = [avLabel]
  else:
    molLabelFractions = labeledMixture.molLabelFractions
  molWeightSum = sum([x.weight for x in molLabelFractions])
  
  # calculate result
  for mlf in molLabelFractions:
    molFactor = mlf.weight / molWeightSum
    molLabel = mlf.molLabel
    
    uncorrelatedAtoms = True
    if resIds[0] == resIds[1]:
      oneResLabel = molLabel.findFirstResLabel(resId=resIds[0])
      for ii in (0,1):
        if (oneResLabel.findAllAtomLabels(atomName=atNames[ii]) or
            oneResLabel.findAllAtomLabels(elementName=elementNames[ii])):
          uncorrelatedAtoms = False

    if uncorrelatedAtoms:
      # isotope frequencies are uncorrelated at the residue level
      dds = []
      for ii in (0,1):
        resLabel = molLabel.findFirstResLabel(resId=resIds[ii])
        dds.append(molLabelFractionsDict(resLabel, atNames[ii], 
                                     subTypes[ii], elementNames[ii]))
      
      for iso0 in dds[0]:
        for iso1 in dds[1]:
          key = (iso0, iso1)
          contrib = dds[0][iso0] * dds[1][iso1] * molFactor
          result[key] = result.get(key, 0.0) + contrib
  
    else:
      # isotope frequencies are correlated at the residue level
      # Only happens if both are from the same residue and neither has
      # any AtomLabels. Loop over ResLabelFractions only
      resLabelFractions = oneResLabel.resLabelFractions
      rlfWeightSum = sum([x.weight for x in resLabelFractions])
      for rlf in resLabelFractions:
        partResult = getIsotopomerAtomPairFractions(rlf.iotopomers, atNames, subTypes)
        for key, val in partResult.items():
          contrib = val * rlf.weight * molFactor / rlfWeightSum
          result[key] = result.get(key, 0.0) + contrib
  
  #
  return result


def molLabelFractionsDict(resLabel, atName, subType, elementName):
  """get the isotopeCode:fraction dictionary for a single molLabel
  
  molLabel: molLabel object
  resid residue serial 
  atName, subType, elementName: atom name subType and element name for atom
  
  Returns isotopeCode:fraction dictionary with fractions normalised to 1.0
  """
  # set up
  result = {}
  
  # get atomLabels, if any
  atomLabels = resLabel.findAllAtomLabels(atomName=atName)
  if not atomLabels:
    atomLabels = resLabel.findAllAtomLabels(elementName=elementName)
  
  # calculate fractions for AtomLabels
  if atomLabels:
    atWeightSum = sum([x.weight for x in atomLabels])
    for atomLabel in atomLabels:
      isotopeCode = '%s%s' % (atomLabel.massNumber, elementName)
      result[isotopeCode] = (result.get(isotopeCode, 0.0) + 
                             atomLabel.weight / atWeightSum)
      
    
  else:
    # calculate fractions for ResLabelFractions
    resLabelFractions = resLabel.resLabelFractions
    rlfWeightSum = sum([x.weight for x in resLabelFractions])
      
    for rlf in resLabelFractions:
      isotopomers = rlf.isotopomers
      isoFactor = rlf.weight  / rlfWeightSum
      
      fractionDict = getIsotopomerSingleAtomFractions(isotopomers, atName, subType)
      
      for isotopeCode in fractionDict.keys():
        contrib = fractionDict[isotopeCode]
        result[isotopeCode] = result.get(isotopeCode,0.0) + (isoFactor * contrib)
 
  #
  return result
  
 
def getIsotopomerAtomPairFractions(isotopomers, atomNames, subTypes):
  """Descrn: Get the combined isotope proportions for a given pair of named
             atoms within a set of isotopomers. Fractions normalised to 1.0
     Inputs: List of ChemCompLabel.Isotopomers, 2-Tuple of Words (ChemAtom.name), 2-Tuple of Words (ChemAtom.subType)
     Output: Dict of Tuple:Float - (IsotopeCode,IsotopeCode):fraction
  """

  fractionDict = {}

  isoWeightSum = sum([x.weight for x in isotopomers])
  
  atLabels   = [None, None]
  sumWeights = [None, None]
  
  for isotopomer in isotopomers:
    for i in (0,1):
      atLabels[i] = isotopomer.findAllAtomLabels(name=atomNames[i], 
                                                 subType=subTypes[i])
      sumWeights[i] = sum([x.weight for x in atLabels[i]])
    
    # Done this way to guard against the divisor becoming zero
    factor  = isotopomer.weight / isoWeightSum 
    divisor = sumWeights[0] * sumWeights[1]
    
    for atl0 in atLabels[0]:
      for atl1 in atLabels[1]:
      
        if atl0 is atl1:
          contrib = atl0.weight * factor / sumWeights[0]  
        else:
          contrib = atl0.weight * atl1.weight * factor / divisor
        
        key = (atl0.isotopeCode, atl1.isotopeCode)
        fractionDict[key] = fractionDict.get(key, 0.0) + contrib

  return fractionDict


def getIsotopomerSingleAtomFractions(isotopomers, atomName, subType=1):
  """Descrn: Get the isotope proportions for a names atom in over a set
             of isotopomers. Fractions normalised to 1.0
     Inputs: List of ChemCompLabel.Isotopomers, Word (ChemAtom.name), Word (ChemAtom.subType)
     Output: Dict of Word:Float - IsotopeCode:fraction
  """

  fractionDict = {}
  isoWeightSum =  sum([x.weight for x in isotopomers])
  
  for isotopomer in isotopomers:
    atomLabels  = isotopomer.findAllAtomLabels(name=atomName,subType=subType)
    atWeightSum = sum([x.weight for x in atomLabels])
    atFactor    = isotopomer.weight  / isoWeightSum
    
    for atomLabel in atomLabels:
      isotopeCode = atomLabel.isotopeCode
      contrib     = atFactor * atomLabel.weight / atWeightSum
      fractionDict[isotopeCode] = fractionDict.get(isotopeCode,0.0) + contrib

  return fractionDict



def testMacro(argServer):

  project = argServer.getProject()
  scheme  = project.findFirstLabelingScheme(name='13Glycerol')
  
  molSystem = project.findFirstMolSystem()

  for chain in molSystem.chains:
    chainCode = chain.code
    
    for residue in chain.residues:
      resId = '%d%s' % (residue.seqCode,residue.ccpCode)
      
      for atom in residue.atoms:
        if atom.chemAtom.elementSymbol != 'C':
          continue
      
        atomName = atom.name

        print '%s %s %s:' % (chainCode, resId, atomName)

        fractions = getSchemeAtomFractions(scheme, atom)
        
        #for isotope in fractions.keys():
        #  print '  %s %.3f' % (isotope,fractions[isotope])

        isotope = '13C'
        print '  %s %.3f' % (isotope,fractions[isotope])

  for chain in molSystem.chains:
    chainCode = chain.code
    
    for residue in chain.residues:
      resId = '%d%s' % (residue.seqCode,residue.ccpCode)

      print '%s %s:' % (chainCode, resId)
      
      atoms = []
      for atom in residue.atoms:
        if atom.chemAtom.elementSymbol != 'C':
          continue
      
        atoms.append((atom.name, atom))

      residue2 = chain.findFirstResidue(seqId = residue.seqId-1)
      if residue2:
        atom = residue2.findFirstAtom(name='C')

        if atom:
          atoms.append(('C_i-1',atom))
      
      for name, atom in atoms:
        for name2, atom2 in atoms:
          fracs  = getSchemeAtomPairFractions(scheme, atom, atom2)
      
          print '  %s %s' % (name, name2)
          
          #for key in fracs:
          key = ('13C','13C')
          
          iso1, iso2 = key
          print '    %s %s %.3f' % (iso1,iso2,fracs[key])
