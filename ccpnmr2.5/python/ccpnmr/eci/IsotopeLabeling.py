import re

from memops.gui.MessageReporter import showError

def getStarIsotopeLabeling(refComp):

  bmrbLabelDict = {'NatAbun':      'natural abundance',
                   'uni_15N':      'U-95% 15N',
                   'uni_15N13C':   'U-95% 13C; U-95% 15N',
                   'uni_15N13C2H': 'U-95% 13C; U-95% 15N; U-95% 2H'}

  molecule = None

  if hasattr(refComp, 'molecule'):
    molecule = refComp.molecule

  else:
    return

  labelMix = refComp.labeledMixture
  if not labelMix:
    return

  # wb104 9 Jul 2013: do not just go up and down from labeledMixture to find a random molLabel
  # instead use averageComposition, if it exists
  molLabel = labelMix.averageComposition
  if not molLabel:
    labelMol = labelMix.parent
    molLabel = labelMol.findFirstMolLabel()

  firstBmrbLabel = ''
  ccpCodes = []
  resLabFracFlag = False
  atomLabFlag = False

  for resLabel in molLabel.sortedResLabels():
    bmrbLabel = ''

    resLabFrac = resLabel.findFirstResLabelFraction()

    if resLabFrac:
      if atomLabFlag:
        print '  Warning: mixture of labelling types in molecule %s' % molecule.name
        break

      resLabFracFlag = True

      bmrbLabel = bmrbLabelDict[resLabFrac.schemeName]
      #print 'BMRB: [%s]' % bmrbLabel

    else:
      if resLabFracFlag:
        print '  Warning: mixture of labelling types in molecule %s' % molecule.name
        break

      atomLabFlag = True


      atomLabel13 = resLabel.findFirstAtomLabel(massNumber=13)

      if atomLabel13:
        if bmrbLabel:
          bmrbLabel += '; '

        if atomLabel13.className == 'UniformAtomLabel' and atomLabel13.elementName == 'C':

          if atomLabel13.weight == 0.00:
            bmrbLabel += 'U-%s%s' % (atomLabel13.massNumber, atomLabel13.elementName)
          else:
            percent = atomLabel13.weight*100
            bmrbLabel += 'U-%d%% %s%s' % (percent, atomLabel13.massNumber, atomLabel13.elementName)

        else:

          if atomLabel13.weight == 0.00:
            bmrbLabel += '%s%s' % (atomLabel13.massNumber, atomLabel13.atomName)
          else:
            percent = atomLabel13.weight*100
            bmrbLabel += '%d%% %s%s' % (percent, atomLabel13.massNumber, atomLabel13.atomName)


      atomLabel15 = resLabel.findFirstAtomLabel(massNumber=15)

      if atomLabel15:
        if bmrbLabel:
          bmrbLabel += '; '

        if atomLabel15.className == 'UniformAtomLabel' and atomLabel15.elementName == 'N':

          if atomLabel15.weight == 0.00:
            bmrbLabel += 'U-%s%s' % (atomLabel15.massNumber, atomLabel15.elementName)
          else:
            percent = atomLabel15.weight*100
            bmrbLabel += 'U-%d%% %s%s' % (percent, atomLabel15.massNumber, atomLabel15.elementName)

        else:

          if atomLabel15.weight == 0.00:
            bmrbLabel += '%s%s' % (atomLabel15.massNumber, atomLabel15.atomName)
          else:
            percent = atomLabel15.weight*100
            bmrbLabel += '%d%% %s%s' % (percent, atomLabel15.massNumber, atomLabel15.atomName)


      atomLabel2 = resLabel.findFirstAtomLabel(massNumber=2)

      if atomLabel2:
        if bmrbLabel:
          bmrbLabel += '; '

        if atomLabel2.className == 'UniformAtomLabel' and atomLabel2.elementName == 'H':

          if atomLabel2.weight == 0.00:
            bmrbLabel += 'U-%s%s' % (atomLabel2.massNumber, atomLabel2.elementName)
          else:
            percent = atomLabel2.weight*100
            bmrbLabel += 'U-%d%% %s%s' % (percent, atomLabel2.massNumber, atomLabel2.elementName)

        else:

          if atomLabel2.weight == 0.00:
            bmrbLabel += '%s%s' % (atomLabel2.massNumber, atomLabel2.atomName)
          else:
            percent = atomLabel2.weight*100
            bmrbLabel += '%d%% %s%s' % (percent, atomLabel2.massNumber, atomLabel2.atomName)


    if bmrbLabel:

      if firstBmrbLabel:

        if firstBmrbLabel != bmrbLabel:
          print '  Warning: multiple labels for this molecule %s' % molecule.name
          break

      else:
        firstBmrbLabel = bmrbLabel

      if molecule and len(molecule.molResidues) > 1:
        ccpCode = molecule.findFirstMolResidue(serial=resLabel.resId).ccpCode

        if molecule.molType in ('DNA', 'RNA', 'DNA/RNA'):

          nucTlcDict = {'A': 'Ade',
                        'C': 'Cyt',
                        'G': 'Gua',
                        'T': 'Thy',
                        'U': 'Ura'}

          ccpCode = nucTlcDict[ccpCode]

        ccpCodes.append(ccpCode)

  if firstBmrbLabel and firstBmrbLabel != 'natural abundance':
    firstBmrbLabel = '[' + firstBmrbLabel + ']'

    if ccpCodes:
      firstCcpCode = ccpCodes[0]
      resFlag = True

      for ccpCode in ccpCodes:
        if ccpCode != firstCcpCode:
          resFlag = False
          break

      if resFlag:
        firstBmrbLabel += '-' + firstCcpCode

  return firstBmrbLabel

def parseBmrbLabelName(bmrbLabelName):

  specificResidue = None

  resLabPatt = re.compile("\]\-([A-Za-z\d]+)$")

  searchObj = resLabPatt.search(bmrbLabelName)

  if searchObj:
    specificResidue = searchObj.group(1)

  bmrbLabelName = resLabPatt.sub(']', bmrbLabelName)

  #print 'SPEC: [%s] [%s]' % (bmrbLabelName, specificResidue)

  if bmrbLabelName.count(';'):
    labels = bmrbLabelName.split(';')

  else:
    labels = [bmrbLabelName]

  labelInfo = []

  for label in labels:
    label = label.strip('[] \t')

    #print 'LABEL: [%s]' % label

    if label == 'natural abundance':
      return None # pass ???

    else:
      labelPatt = re.compile("^(U\-)?(\d+%)? ?(\d+)([HCN])([A-Z]*)$")

      searchObj = labelPatt.search(label)

      if searchObj:
        uniform  = searchObj.group(1)
        percent  = searchObj.group(2)
        mass     = searchObj.group(3)
        atomType = searchObj.group(4)
        specificAtom = searchObj.group(5)

        labelInfo.append( (uniform, percent, mass, atomType, specificResidue, specificAtom) )

  return labelInfo

def makeUniformLabels(resLabel, atomType, mass, percent):

  isotopeMasses = {'C': '13',
                   'N': '15',
                   'H': '2'}

  isotopeNuclei = {'C': 'carbon',
                   'N': 'nitrogen',
                   'H': 'hydrogen'}

  if mass != isotopeMasses[atomType]:
    print '  Warning: labelled %s isotope does not have mass of %s' % (isotopeNuclei[atomType], isotopeMasses[atomType])
    return

  if not percent:
    weight = 0.00
  else:
    weight = float(percent[:-1])/100.0

  if weight < 0 and weight > 1:
    print '  Warning: degree of isotope labelling not in the correct range (0-100%)'
    return


  uniformAtomLabel2 = None

  uniformAtomLabel1 = resLabel.newUniformAtomLabel(elementName=atomType, massNumber=int(mass), weight=weight)

  if weight != 1.00:
    uniformAtomLabel2 = resLabel.newUniformAtomLabel(elementName=atomType, massNumber=int(mass)-1, weight=1.00-weight)

  return (uniformAtomLabel1, uniformAtomLabel2)

def makeSpecAtomLabels(resLabel, atomType, mass, percent, specificAtom):

  isotopeMasses = {'C': '13',
                   'N': '15',
                   'H': '2'}

  isotopeNuclei = {'C': 'carbon',
                   'N': 'nitrogen',
                   'H': 'hydrogen'}

  if mass != isotopeMasses[atomType]:
    print '  Warning: labelled %s isotope does not have mass of %s' % (isotopeNuclei[atomType], isotopeMasses[atomType])
    return

  if not percent:
    weight = 0.00
  else:
    weight = float(percent[:-1])/100.0

  if weight < 0 and weight > 1:
    print '  Warning: degree of isotope labelling not in the correct range (0-100%)'
    return

  specAtomLabel2 = None

  specAtomLabel1 = resLabel.newSingleAtomLabel(atomName=atomType + specificAtom, massNumber=int(mass), weight=weight)
  if weight != 1.00:
    specAtomLabel2 = resLabel.newSingleAtomLabel(atomName=atomType + specificAtom, massNumber=int(mass)-1, weight=1.00-weight)

  #print 'SPECS: [%s] [%s]' % (specAtomLabel1, specAtomLabel2)

  return (specAtomLabel1, specAtomLabel2)

def doLabels(labelInfo, bmrbLabelName, molLabel, molecule=None, molRes=None):

  schemeNameDict = {'natural abundance':                'NatAbun',
                    '[U-95% 15N]':                      'uni_15N',
                    '[U-95% 13C; U-95% 15N]':           'uni_15N13C',
                    '[U-95% 13C; U-95% 15N; U-95% 2H]': 'uni_15N13C2H'}

  if molRes is None:
    molResId = 1
  else:
    molResId = molRes.serial

  resLabel = molLabel.findFirstResLabel(resId=molResId)
  if not resLabel:
    resLabel = molLabel.newResLabel(resId=molResId)

  #for key in resLabel.__dict__['resLabelFractions'].keys():
  #  del(resLabel.__dict__['resLabelFractions'][key])
  #for key in resLabel.__dict__['atomLabels'].keys():
  #  del(resLabel.__dict__['atomLabels'][key])

  for rlf in resLabel.resLabelFractions:
    rlf.delete()
  for al in resLabel.atomLabels:
    al.delete()

  if bmrbLabelName in schemeNameDict:
    #resLabelFrac = resLabel.findFirstResLabelFraction(schemeName=schemeNameDict[bmrbLabelName])

    #if not resLabelFrac:
    #  pass
    resLabelFrac = resLabel.newResLabelFraction(schemeName=schemeNameDict[bmrbLabelName])

  else:
    uniformLabels = None

    if labelInfo:
      for label in labelInfo:

        (uniform, percent, mass, atomType, specificResidue, specificAtom) = label

        ccpCode = ''

        if molecule and molRes:

          if molecule.molType in ('DNA', 'RNA', 'DNA/RNA'):

            nucTlcDict = {'A': 'Ade',
                          'C': 'Cyt',
                          'G': 'Gua',
                          'T': 'Thy',
                          'U': 'Ura'}

            ccpCode = nucTlcDict[molRes.ccpCode]

          elif molecule.molType == 'protein':

            ccpCode = molRes.ccpCode

          #elif sugar  # TODO: ???

        #print "DATA: [%s] [%s] [%s] [%s] [%s] [%s] [%s] [%s] [%s]" % (uniform, percent, mass, atomType, specificResidue, molRes, molRes.ccpCode, ccpCode, specificAtom)

        if specificResidue is not None and specificResidue != ccpCode:
          continue

        if uniform == 'U-':
          #print 'Making uniform labels'
          uniformLabels = makeUniformLabels(resLabel, atomType, mass, percent)

          if uniformLabels is None:
            continue

        elif specificAtom:
          #print 'Making specific labels'
          specAtomLabels = makeSpecAtomLabels(resLabel, atomType, mass, percent, specificAtom)

        else:
          continue

    #print 'RES2: [%s] [%s] [%s]' % (resLabel, uniformLabels[0], uniformLabels[1])

def makeLabelObjects(mr, refComp, bmrbLabelName):

  # TBD: wb104 15 Aug 2012: before the refComp.className did not have to be
  # 'MolComponent' but that then causes a problem with creating a LabeledMolecule,
  # which has to have a Molecule with that name.  Look at again.
  if not refComp.className == 'MolComponent':
    showError('Not a MolComponent', 'refComponent is not a MolComponent')
    return

  molecule = refComp.molecule
  if not molecule:
    showError('No Molecule', 'refComponent has no Molecule')
    return

  name = molecule.name
    
  labMol = mr.findFirstLabeledMolecule(name=name)
  if not labMol:
    labMol = mr.newLabeledMolecule(name=name)

  # wb104 9 Jul 2013: do not just use first MolLabel found
  """
  molLabel = labMol.findFirstMolLabel()
  if not molLabel:
    molLabel = labMol.newMolLabel()
"""

  # wb104 9 Jul 2013: check if labeledMixture set and if not always create a new one
  labMix = refComp.labeledMixture
  if not labMix:
    labMix = labMol.newLabeledMixture()

  # wb104 9 Jul 2013: set up link from LabeledMixture to new MolLabel if needed
  molLabel = labMix.averageComposition
  if not molLabel:
    molLabel = labMix.averageComposition = labMol.newMolLabel()
  """
  # previous way of doing it
  labMix = labMol.findFirstLabeledMixture()
  if not labMix:
    labMix = labMol.newLabeledMixture()
"""

  # wb104 9 Jul 2013: old fashioned way to do this (and makes grepping worse)
  #refComp.setLabeledMixture(labMix)
  refComp.labeledMixture = labMix

  molLabelFrac = labMix.findFirstMolLabelFraction(molLabel=molLabel)
  if not molLabelFrac:
    molLabelFrac = labMix.newMolLabelFraction(molLabel=molLabel)

  labelInfo = parseBmrbLabelName(bmrbLabelName)

  molecule = None

  if refComp.className == 'MolComponent':
    molecule = refComp.molecule

    if molecule and not molecule.isFinalised:
      molecule.isFinalised = True

  if molecule:
    for molRes in molecule.sortedMolResidues():
      doLabels(labelInfo, bmrbLabelName, molLabel, molecule, molRes)

  else:
    doLabels(labelInfo, bmrbLabelName, molLabel)

if __name__ == '__main__':

  pass
