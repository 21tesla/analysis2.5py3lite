"""Mapping of Resonances and ResonancGroups in version 2 to new 4-string assignment style

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================

__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - $Date: 2014-10-23 14:41:29 +0100 (Thu, 23 Oct 2014) $"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Simon Skinner, Geerten Vuister"
__license__ = ("CCPN license. See www.ccpn.ac.uk/industry/licencefees"
               "or ccpncore.memops.Credits.CcpnLicense for license text")
__reference__ = ("For publications, please use reference from www.ccpn.ac.uk/industry/licencefees"
                 " or ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification:
#=========================================================================================
__author__ = "$Author: rhfogh $"
__date__ = "$Date: 2014-10-23 14:41:29 +0100 (Thu, 23 Oct 2014) $"
__version__ = "$Revision: 7875 $"

#=========================================================================================
# Start of code
#=========================================================================================

import os
import operator

from ccpncore.lib import MoleculeQuery
from ccpncore.lib.chemComp.Util import getChemAtomSetFromAtoms


def mapAllAssignments(topObject, assignmentMap=None, molSystem=None, chainMap=None):
  """V2: make/extend resonance:assignmentList map for all assignments in NmrProject/NmrConstraintStore
  Either molSystem or chainMap must be passed in"""

  if assignmentMap is None:
    assignmentMap = {}

  assert ((molSystem is None) != (chainMap is None)), "Pass in molSystem or chainMap, but not both."

  # Map assigned resonances
  mapAssignedResonances(topObject, assignmentMap, chainMap=chainMap, molSystem=molSystem)

  # For NmrProject map partly assigned resonances
  if topObject.className == 'NmrProject':
    _mapUnAssignedResonances(topObject, assignmentMap, molSystem=molSystem, chainMap=chainMap)
  else:
    # map unassigned resonances
    _mapUnAssignedFixedResonances(topObject, assignmentMap)

  _testAssignmentMap(assignmentMap)

  #
  return assignmentMap

def _testAssignmentMap(assignmentMap):
  dd = {}
  for res,ass in assignmentMap.items():
    if res.className == 'Resonance':
      # Fixed resonance duplications do not matter
      # NB this test must be standard - these problems likely arise from incorrect assignments
      ass = tuple(ass)
      ll = dd.get(ass,[])
      ll.append(res)
      dd[ass] = ll

  for ass,ll in sorted(dd.items()):
    if len(ll) > 1:
      print ('### DUPLICATE RESONANCE %s %s' % (ass, ll))


def mapResonanceGroups(nmrProject, molSystem=None, chainMap=None, defaultChainCode=None):
  """Map ResonanceGroups to three-string assignments"""

  assert ((molSystem is  None) != (chainMap is None)),  "Pass in molSystem or chainMap, not both"

  result = {}

  assignedGroups = []

  # Handle properly assigned ResonanceGroups
  for resonanceGroup in nmrProject.sortedResonanceGroups():

    # Find residue
    residue = resonanceGroup.residue
    if residue is None:
      ll = []
      for residueProb in resonanceGroup.residueProbs:
        if residueProb.weight:
          ll.append(residueProb)
      if len(ll) == 1:
        # In principle this should never happen, but it does not hurt to be careful
        residue = ll[0].possibility

    if residue:

      # Remap chain and residue if necessary,

      chain = residue.chain
      if chainMap:
        # If there is a chainMap the chain MUST be in it.
        chain = chainMap[chain]
        # get residue in new chain
        residue = chain.findFirstResidue(seqCode=residue.seqCode,
                                         seqInsertCode=residue.seqInsertCode)

      elif residue.topObject is not molSystem:
          raise ValueError("Cannot generate consistent assignment names from mixed MolSystems - 1")

      # set residue assignment strings
      chainCode = chain.code
      sequenceCode = str(residue.seqCode)+ (residue.seqInsertCode or '').strip()
      residueType = residue.molResidue.chemComp.code3Letter

      result[resonanceGroup] = (chainCode, sequenceCode, residueType)
      assignedGroups.append(resonanceGroup)


  # Now add dependent ResonanceGroups
  for resonanceGroup in assignedGroups:
    assignment = result[resonanceGroup]

    # Add sequential stretches either way
    for direction in (+1,-1):
      stretch = _findSpinSystemStretch(resonanceGroup, delta=direction)
      offset = 0
      for rg in stretch:
        if rg in result:
          break
        else:
          # Identify by seq offset:
          chainCode = assignment[0]
          offset += direction
          sequenceCode = '%s%+s' % (assignment[1], offset)

          chemComp = rg.root.findFirstChemComp(molType=rg.molType, ccpCode=rg.ccpCode)
          if chemComp:
            residueType = chemComp.code3Letter
          else:
            residueType = None

          result[rg] = (chainCode, sequenceCode, residueType)


  # Now look for stretches of ResonanceGroups
  newChainCounter = 0
  for resonanceGroup in nmrProject.sortedResonanceGroups():

    if resonanceGroup not in result:

      if not  _findSpinSystemStretch(resonanceGroup, delta=+1):
        # Start only at +1 end of stretches

        stretch =  _findSpinSystemStretch(resonanceGroup, delta=-1)

        if not stretch:
          # Do nothing if no residue stretch found
          pass

        elif len(stretch) > 2:
          # proper assigned stretch - treat at pseudochain
          stretch.reverse()
          stretch.append(resonanceGroup)
          newChainCounter += 1
          chainCode = '@%s' % newChainCounter
          sequenceCode = 0
          for rg in enumerate(stretch):
            sequenceCode += 1
            chemComp = rg.root.findFirstChemComp(molType=rg.molType, ccpCode=rg.ccpCode)
            if chemComp:
              residueType = chemComp.code3Letter
            else:
              residueType = None

            result[rg] = (chainCode, str(sequenceCode), residueType)

        elif stretch[0] not in result:
          # 2-3 residue stretch, not assigned.
          # Set as a residue with +/- markers on neighbours
          if len(stretch) == 1 or stretch[-1] in result:
            # Two-residue stretch.
            rg = resonanceGroup
            rgprev = stretch[0]
            rgnext = None
          else:
            # assert len(stretch) == 2
            rg = stretch[0]
            rgprev = stretch[-1]
            rgnext = resonanceGroup

          # set codes for rg
          chainCode = defaultChainCode
          sequenceCode = '@%s' % rg.serial
          chemComp = rg.root.findFirstChemComp(molType=rg.molType, ccpCode=rg.ccpCode)
          if chemComp:
            residueType = chemComp.code3Letter
          else:
            residueType = None
          result[rg] = (chainCode, sequenceCode, residueType)


          # Now set for rgprev
          chemComp = rgprev.root.findFirstChemComp(molType=rgprev.molType,
                                                   ccpCode=rgprev.ccpCode)
          if chemComp:
            residueType = chemComp.code3Letter
          else:
            residueType = None
          result[rgprev] = (chainCode, '%s-1' % sequenceCode, residueType)

          if rgnext is not None:
            # Now set for rgnext
            chemComp = rgnext.root.findFirstChemComp(molType=rgnext.molType,
                                                     ccpCode=rgnext.ccpCode)
            if chemComp:
              residueType = chemComp.code3Letter
            else:
              residueType = None
            result[rgnext] = (chainCode, '%s+1' % sequenceCode, residueType)


  # Finally deal with ResonanceGroups neither assigned nor in sequential stretches
  for resonanceGroup in nmrProject.sortedResonanceGroups():
    if resonanceGroup not in result:

      # Is it identity-linked to one there:
      rg0 = _findIdentityResonanceGroup(resonanceGroup)
      assignment = result.get(rg0)
      if (assignment is not None and rg0 is not None and
          not (assignment[1][-2] in '+-' and assignment[1][-1].isdigit())):
        # Identity-linked to existing residue - set accordingly
        result[resonanceGroup] = (assignment[0], assignment[1] + '+0', assignment[2])

      else:
        # Not linked to anything - set on its own merits
        chemComp = resonanceGroup.root.findFirstChemComp(molType=resonanceGroup.molType,
                                                         ccpCode=resonanceGroup.ccpCode)
        if chemComp:
          residueType = chemComp.code3Letter
        else:
          residueType = None
        sequenceCode = '@%s' % resonanceGroup.serial
        result[resonanceGroup] = (defaultChainCode, sequenceCode, residueType)


def mapAssignedResonances(topObject, assignmentMap, molSystem=None, chainMap=None):
  """Make/extend {resonance:assignmentTuple} map in V2 for either Resonances or fixedResonances
  chainMap remaps chains to new ones with different chainCodes (for V2-V3 upgrade)."""

  if topObject.className == 'NmrProject':
    resonanceSets = topObject.sortedResonanceSets()
  else:
    resonanceSets = topObject.sortedFixedResonanceSets()

  for resonanceSet in resonanceSets:

    # set up loop-level parameters
    resonances = list(resonanceSet.resonances)
    resonances.sort(key=operator.attrgetter('name'))
    atomSets = list(resonanceSet.atomSets)
    atomSets.sort(key=operator.attrgetter('name'))
    allAtoms = [x for y in atomSets for x in y.atoms]
    chemAtomSet = getChemAtomSetFromAtoms(allAtoms)
    if len(set(x.residue for x in allAtoms)) == 1:
      residue = allAtoms[0].residue
    else:
      residue = None

    if residue:

      # Remap chain and residue if necessary
      chain = residue.chain
      if chainMap:
        # If there is a chainMap the chain MUST be in it.
        chain = chainMap[chain]
        # get residue in new chain
        residue = chain.findFirstResidue(seqCode=residue.seqCode,
                                         seqInsertCode=residue.seqInsertCode)

      elif chain.molSystem is not molSystem:
        raise ValueError("Cannot generate consistent assignment names from mixed MolSystems - 2")

      # set residue assignment strings and some more variables
      chainCode = chain.code
      sequenceCode = str(residue.seqCode) + (residue.seqInsertCode or '').strip()

      chemComp = residue.molResidue.chemComp
      chemCompVar = residue.chemCompVar or chemComp.findFirstChemCompVar(isDefaultVar=True)
      residueType = chemComp.code3Letter

      # Now for the atom name
      if chemAtomSet and len(atomSets) == 2 and len(resonances) <= 2:
        # prochiral pair. Use _getAmbigProchiralLabel for priority, and set 'X' for 'a', 'Y' for 'b'

        # get non-stereo names
        atomSetNames = [x.name for x in atomSets]
        starpos = chemAtomSet.name.find('*')
        newNames = []
        for ii, newChar in enumerate('XY'):
          chars = list(atomSetNames[ii])
          chars[starpos] = newChar
          newNames.append(''.join(chars))

        # select new name to use
        if topObject.className == 'NmrProject':
          # Real resonance - use normal procedure

          # First one resonance. NB we test prochiral label on one resonance only
          # - if data are inconsistent both might give 'a' and we must avoid a name clash.
          if _getAmbigProchiralLabel(resonances[0]) != 'a':
            # First resonance does not match first name - reverse name order
            # NB we test on one resonance only. If data are inconsistent and both give 'a'
            # we still want the names to be different
            newNames.reverse()

          for ii,resonance in enumerate(resonances):
            assignmentMap[resonance] = (chainCode, sequenceCode, residueType,
                                        newNames[ii].replace('*','#'))


        else:
          # NmrConstraintStore - these are fixed resonances

          if len(resonances) == 2:
            # resonances are sorted by name, as are newNames. Match in order.
            # NB this being in an NmrConstraintStore, we have to assume that names  are
            # consistent, so e.g. HG1* is bound to CG1 and not CG2
            for ii in range(2):
              assignmentMap[resonances[ii]] = (chainCode, sequenceCode, residueType,
                                               newNames[ii].replace('*','#'))

          else:
            # Only one resonance. Must choose which.
            # Use various heuristics, as we can not assume assignments in NmrProject still match.

            # assert len(resonances) == 1
            resonance = resonances[0]
            for ii in range(2):
              if (resonance.name in(atomSetNames[ii], newNames[ii]) or
                  resonance.name.upper() in(atomSetNames[ii], newNames[ii])):
                # name matches one of the atomSets - use matching name
                indx = ii
                break

            else:
              # Name did not match either possibility, so these are non-standard. Try anyway
              if 'b' in resonance.name or '3' in resonance.name:
                # Heuristic - these names are likely to be second in sorting
                # e.g. 'HBb' or 'HB3'
                indx = 1
              elif 'a' in resonance.name or '1' in resonance.name:
                # Heuristic = this name is likely to be first in sorting (e.g. 'HBa'
                # NB both XY1, XY1*, XY11, XY21 (Asn, Gln,Arg) XY12 (Ile) should sort first
                # Arg HH12 should sort second, but this cannot be helped
                indx = 0
              elif '2' in resonance.name:
                # This will not sort second for e,g, HB2/HB3, but that cannot be helped.
                # The only one that matters is
                # where e.g. HD1 must bind to CD1 and HD2 to CD2, and it gets those right
                # Anyway the real atom names are caught before this
                indx = 1
              else:
                realResonance = resonance.resonance
                if realResonance:
                  # We can not be sure that assignments have not changed, but
                  # better have it match the resonance than not
                  if _getAmbigProchiralLabel(realResonance) == 'a':
                    indx = 0
                  else:
                    indx = 1
                else:
                  # Stuff it. We just do not know and pick the first one.
                  # Anyway the cases where is makes a difference are cared for above.
                  indx = 0

            resonanceName = newNames[indx]
            assignmentMap[resonance] = (chainCode, sequenceCode, residueType,
                                        resonanceName.replace('*','#'))

      elif len(resonances) == 1:
        # Single resonance, not assigned to prochiral

        resonance = resonances[0]
        if len(atomSets) == 1:
          # simple one-to-one stereospecific assignment
          if chemAtomSet:
            # assignment to atomSet
            resonanceName = chemAtomSet.name.replace('*', '#')
          else:
            #asssignment to single atom
            resonanceName = allAtoms[0].name

        else:
          # multiple atomSets
          # NB we do it this way because it must work in pure V2 as well as the intermediate model
          nuc = allAtoms[0].elementSymbol
          residueChemAtoms = [chemCompVar.findFirstChemAtom(name=x.name) for x in residue.atoms]
          residueChemAtoms = [x for x in residueChemAtoms if x is not None]
          if len(allAtoms) == len(residueChemAtoms):
            # All single atoms in residue
            resonanceName = '*'

          elif (all((x.elementSymbol == nuc) for x in allAtoms) and
              len([x for x in residueChemAtoms if x.elementSymbol == nuc]) == len(allAtoms)):
            # All atoms of a given nucleus
            resonanceName = nuc + '*'

          else:
            # random multiple atom selection
            resonanceName = '/'.join(sorted((str(x.name) for x in atomSets)))

        assignmentMap[resonance] = (chainCode, sequenceCode, residueType, resonanceName)

      else:
        # multiple resonances not matching chemAtomSet
        atomsName = '/'.join(sorted((str(x.name) for x in atomSets)))
        for resonance in resonances:
          # NB this name can not be in use already, so we do not need to check
          resonanceName = '%s@%s' % (atomsName, resonance.serial)

          assignmentMap[resonance] = (chainCode, sequenceCode, residueType,
                                      resonanceName.replace('*','#'))

    else:
      # assigned to multiple residues - cannot be helped
      # Same naming style for single and multiple resonances
      partNames = []
      residueTypes = set()
      for atomSet in atomSets:
        rr = atomSet.findFirstAtom().residue
        residueTypes.add(rr.molResidue.chemComp.code3Letter)
        partNames.append('%s%s-%s' % (rr.seqCode, (rr.seqInsertCode or '').strip(), atomSet.name))
      if len(residueTypes) == 1:
        residueType = residueTypes.pop()
      else:
        residueType = None
      chainCode = None
      sequenceCode = None

      ss = '/'.join(sorted(partNames))
      for resonance in resonances:
        resonanceName = '%s@%s' % (ss, resonance.serial)

        assignmentMap[resonance] = (chainCode, sequenceCode, residueType,
                                    resonanceName.replace('*','#'))

###################################################################################
#
#  Functions for V2/upgrade only:
#
###################################################################################

def regularisedResonanceName(resonance):
  """V2: Get resonance name, starting with element type and adding @serial to impossible names"""

  # NB names like '*', 'C*' etc. will not make it through regularisation
  # but it is too dangerous to arrive at those from name strings only
  # For that you need to have them properly assigned.

  resonanceName = resonance.name or ''

  if resonance.className == 'Resonance':
    # Exclude fixedResonances
    assignNames = list(set(resonance.assignNames))
    if len(assignNames) == 1:
      # One assignName - use for assignment
      name = assignNames[0]
      return name.replace('*', '#')


    elif len(assignNames) == 2:
      # 2 assignNames - if they match prochiral return nonstereo variant
      assignNames.sort()
      prefix = os.path.commonprefix(assignNames)
      lenPrefix = len(prefix)
      if assignNames[0][lenPrefix].isdigit() and assignNames[1][lenPrefix].isdigit():
        if assignNames[0][lenPrefix+1:] == assignNames[1][lenPrefix+1:]:

          # Heuristics, try to get X or Y right
          # NB 'X' will not always match the lowest sorting - it fails for
          # ordinary HB2/HB3 methylenes and must be caught later
          # But the crucial isopropyl and aromatic cases should work
          newChar = None
          for char in reversed(resonanceName):
            if char.isdigit():
              break
          else:
            char = None
          if assignNames[0][lenPrefix] == char:
            newChar = 'X'
          elif assignNames[1][lenPrefix] == char:
            newChar = 'Y'
          elif 'a' in resonanceName and not 'b' in resonanceName:
            newChar = 'X'
          elif 'b' in resonanceName and not 'a' in resonanceName:
            newChar = 'Y'
          if newChar is not None:
            ll = list(assignNames[0])
            ll[lenPrefix] = newChar
            return ''.join(ll).replace('*', '#')

  # If we are still here, assignNames did not help. Use resonanceName only

  # get elementCode
  isotope = resonance.isotope
  if isotope is None:
    elementCode = ''
  else:
    elementCode = isotope.chemElement.symbol.upper()

  upperName = resonanceName.upper()

  result = None
  if elementCode and upperName.startswith(elementCode):
    # name is OK
    if not resonanceName.startswith(elementCode):
      # except for casing - fix the casing
      result = elementCode + resonanceName[len(elementCode):]

      if 'X' in result or 'Y' in result:
        # Necessary to avoid potential clashes with XY names set above
        # No actual assigned names contain 'X' or 'Y' anyway
        result = '%s@%s' % (result, resonance.serial)
      else:
        # Name might be proper assignment name. Change to new wildcard convention
        result.replace('*', '#')


  else:
    # Set unique default name
    result = '%s%s@%s' % (elementCode, resonanceName, resonance.serial)

  #
  return result




def _mapUnAssignedFixedResonances(nmrConstraintStore, assignmentMap):
  """Map unassigned resonances for NmrConstraintStore"""

  separator1 = '@'
  separator2 = '@@'   # To distinguish fixedResonance serial from resonance serial

  for resonance in nmrConstraintStore.sortedFixedResonances():
    if not resonance.resonanceSet:
      # unassigned - treat it

      name =  regularisedResonanceName(resonance)

      # Add resonance serial to name, if it is not in the name already.
      separator = separator1
      # use resonanceSerial if available
      serial = resonance.resonanceSerial
      if not serial:
        # use FixedResonance.serial instead, and use '@@' to distinguish
        serial = resonance.serial
        separator = separator2

      ss = '%s%s' % (separator1, serial)
      if ss not in name:
        name = '%s%s%s' % (name, separator, serial)

      #
      assignmentMap[resonance] = (None, None, None, name)


def _mapUnAssignedResonances(nmrProject, assignmentMap, molSystem, chainMap=None):
  """Map part- or un-assigned resonances for NmrProject"""

  if molSystem is None:
    # Set molSystem from chainMap
    for key,val in chainMap:
      molSystem = val.molSystem
      break

  if molSystem and len(molSystem.chains) == 1:
    defaultChainCode = molSystem.findFirstChain().code
  else:
    defaultChainCode = None

  resonanceGroupMap =  mapResonanceGroups(nmrProject, molSystem=molSystem,
                                           defaultChainCode=defaultChainCode)

  resonanceGroupMap[None] = (None, None, None)

  for resonance in nmrProject.sortedResonances():
    if resonance not in assignmentMap:
      groupAssignment = resonanceGroupMap[resonance.resonanceGroup]
      name =  regularisedResonanceName(resonance)
      assignmentMap[resonance] = groupAssignment + (name,)


def _findSpinSystemStretch(spinSystem, direction=1):
  """
  V2: Find unambiguous stretch of spin systems sequentially connected to the input one
  in direction +1/-1.

  .. describe:: Input

  Nmr.ResonanceGroup, Int
  .. describe:: Output

  Nmr.ResonanceGroup
  """

  stretch = [spinSystem]

  while True:

    ll = []
    for link in stretch[-1].findAllResonanceGroupProbs(linkType='sequential',
                                                       sequenceOffset=direction,
                                                       isSelected=True):
      ll.append(link.possibility)

    for link in stretch[-1].findAllFromResonanceGroups(linkType='sequential',
                                                       sequenceOffset=-direction,
                                                       isSelected=True):
      ll.append(link.parent)

    if len(ll) == 1:
      stretch.append(ll[0])
    else:
      break
  #
  return stretch[1:]


def  _findIdentityResonanceGroup(resonanceGroup):
  """V2: find unique resonanceGroup linked as identical to teh input"""
  ll = resonanceGroup.findAllResonanceGroupProbs(linkType='identity',
                                                  isSelected=True)
  identicals = list(x.possibility for x in ll)
  ll2 = resonanceGroup.findAllFromResonanceGroups(linkType='identity',
                                                  isSelected=True)
  identicals.extend(x.parent for x in ll2)
  if len(identicals) == 1:
    return identicals[0]
  else:
    return None

def _getAmbigProchiralLabel(resonance):
  """
  V2: Deterimine if an ambigous prochiral resonance (non-stereospecifically assigned)
  Has an "a" label or a "b" label. "a" is reserved for the upfield proton and any
  other nulceus bound to it.

  .. describe:: Input

  Nmr.Resonance

  .. describe:: Output

  Character
  """

  letter = ''
  if hasattr(resonance, 'onebond'):
    del resonance.onebond

  resonanceSet = resonance.resonanceSet

  if resonanceSet:
    if resonance.isotopeCode == '1H':
      data = []
      for resonance2 in resonanceSet.sortedResonances():
        if resonance2.shifts:
          data.append( ('%f%d' % (resonance2.findFirstShift().value,resonance2.serial),resonance2) )
        else:
          data.append( (resonance2.serial,resonance2) )

      data.sort()
      resonances = [x[1] for x in data]
      i = resonances.index(resonance)
      letter = chr(ord('a')+i)

    else:
      resonance2 = _getOnebondResonance(resonance, isotopeCode='1H')

      if resonance2 and resonance2.resonanceSet and (len(resonance2.resonanceSet.atomSets) > 1):
        letter = _getAmbigProchiralLabel(resonance2)
        resonance2.onebond = resonance

      elif (len(resonanceSet.resonances) > 1) and (len(resonanceSet.atomSets) > 1):
        for resonance2 in resonanceSet.resonances:
          if resonance2 is not resonance:
            resonance3 = _getOnebondResonance(resonance2)
            if resonance3 and resonance3.resonanceSet and (len(resonance3.resonanceSet.atomSets) > 1):
              letter = 'b'
            break

      if not letter:
        data = []
        for resonance2 in resonanceSet.resonances:
          if resonance2.shifts:
            data.append( (resonance2.findFirstShift().value,resonance2) )
          else:
            data.append( (resonance2.serial,resonance2) )

        data.sort()
        resonances = [x[1] for x in data]
        i = resonances.index(resonance)
        letter = chr(ord('a')+i)

  #keyword = 'ambigProchiralLabel'
  #app     = 'Analysis'
  #appData = resonance.findFirstApplicationData(application=app, keyword=keyword)
  #
  #if appData and (appData.value != letter):
  #  appData.delete()
  #  appData = None
  #
  #if not appData:
  #  AppDataString(resonance,application=app,keyword=keyword, value=letter)

  return letter


def _getOnebondResonance(resonance, isotopeCode=None):
  """
  V2: Find any resonance that may have a single bond connection to the input resonance
  Option to specify the isotope type

  .. describe:: Input

  Nmr.Resonance, Nmr.Resonance.isotopeCode

  .. describe:: Output

  Nmr.Resonance
  """

  resonances = _getBoundResonances(resonance)
  if resonances:
    if isotopeCode:
      for resonance1 in resonances:
        if resonance1.isotopeCode == isotopeCode:
          return resonance1

    else:
      return resonances[0]

  resonance2 = None

  for contrib in resonance.peakDimContribs:
    peakDim      = contrib.peakDim
    expDimRef1   = peakDim.dataDimRef.expDimRef
    expTransfers = expDimRef1.expTransfers

    for expTransfer in expTransfers:
      if expTransfer.transferType in ('onebond','CP'):
        expDimRef2 = None

        for expDimRef in expTransfer.expDimRefs:
          if expDimRef is not expDimRef1:
            expDimRef2 = expDimRef
            break

        if expDimRef2:
          if (not isotopeCode) or (isotopeCode in expDimRef2.isotopeCodes):
            for peakDim2 in peakDim.peak.peakDims:
              if peakDim2.dataDimRef and (peakDim2.dataDimRef.expDimRef is expDimRef2):
                for contrib2 in peakDim2.peakDimContribs:
                  if (not isotopeCode) or (contrib2.resonance.isotopeCode == isotopeCode):
                    resonance2 = contrib2.resonance

                break

    if resonance2:
      break

  return resonance2



def _getBoundResonances(resonance, recalculate=False, contribs=None, recursiveCall=False):
  """
  V2: Find all resonances that have a single bond connection to the input resonance
  Option to recalculate given assignment status (e.g. if something changes)
  Option to specify peakDimContribs to search

  .. describe:: Input

  Nmr.Resonance, Boolean, List of Nmr.PeakDimContribs

  .. describe:: Output

  List of Nmr.Resonances
  """

  if (not recalculate) and resonance.covalentlyBound:
    return list(resonance.covalentlyBound)

  resonances = set() # Linked by bound atoms irrespective of spectra
  pairResonances = set() # prochiral or other pairs that can not be determined imemdiately
  resonanceSet   = resonance.resonanceSet

  funnyResonances = set()

  if resonanceSet:
    #residue  = resonanceSet.findFirstAtomSet().findFirstAtom().residue
    atomSets = resonanceSet.atomSets

    for atomSet in atomSets:
      #for atom in atomSet.atoms:
      atom = atomSet.findFirstAtom()

      for atom2 in MoleculeQuery.getBoundAtoms(atom):
        atomSet2 = atom2.atomSet

        if atomSet2 and atomSet2.resonanceSets:

          usePaired = False
          if len(atomSets) > 1:
            chemAtomSet = atom2.chemAtom.chemAtomSet
            if chemAtomSet:
              usePaired = (chemAtomSet.isProchiral or
                           (chemAtomSet.chemAtomSet and chemAtomSet.chemAtomSet.isProchiral))

          for resonanceSet2 in atomSet2.resonanceSets:
            for resonance2 in resonanceSet2.resonances:
              if resonance2 is resonance: # should not happen
                if resonance not in funnyResonances:
                  resonance.root._logger.warning( 'in getBoundResonances():'
                                                  'resonance %d tried to be linked to itself'
                                                  % resonance.serial)
                  funnyResonances.add(resonance)
              elif usePaired:
                pairResonances.add(resonance2)
              else:
                resonances.add(resonance2)

  if not contribs:
    contribs = resonance.peakDimContribs

  expResonances = set()
  foundBothPaired = False
  for contrib in contribs:
    peakDim      = contrib.peakDim
    expDimRef1   = peakDim.dataDimRef.expDimRef
    expTransfers = expDimRef1.expTransfers

    for expTransfer in expTransfers:
      if expTransfer.transferType == 'onebond':
        expDimRef2 = None

        for expDimRef in expTransfer.expDimRefs:
          if expDimRef is not expDimRef1:
            expDimRef2 = expDimRef
            break

        if expDimRef2:
          for peakDim2 in peakDim.peak.peakDims:
            if peakDim2.dataDimRef and (peakDim2.dataDimRef.expDimRef is expDimRef2):
              expBound = set()

              for contrib2 in peakDim2.peakDimContribs:
                if (not contrib.peakContribs) and (not contrib2.peakContribs):
                  resonance2 = contrib2.resonance

                  if resonance is not resonance2:
                    expBound.add(resonance2)

                else:
                  for peakContrib in contrib.peakContribs:
                    if peakContrib in contrib2.peakContribs:
                      resonance2 = contrib2.resonance

                      if resonance is not resonance2:
                        expBound.add(resonance2)

                      break

              if len(expBound) > 1:
                # Ambiguity
                for bound in expBound:
                  # Leave the covalently bound one
                  if bound in resonances:
                    break

                else:
                  aSet = set(x for x in expBound if x in resonance.covalentlyBound)
                  if aSet and aSet != pairResonances:
                    # Resonances found. Previously linked.
                    # Not the pairResonances. Use them
                    expResonances.update(aSet)

                  else:
                    # check presence of prochiral pairs
                    ll = [x for x in pairResonances if x in expBound]
                    if len(pairResonances) == 2 and len(ll) == 2:
                      foundBothPaired= True
                    elif ll:
                      # found some prochiral pair resonances - use them
                      expResonances.update(ll)
              else:
                expResonances.update(expBound)

  if foundBothPaired and not [x for x in expResonances if x in pairResonances]:
    # particular special case.
    # Resonance is bound to both prochiral altrnatives but always as a pair.

    if recursiveCall:
      # This was called from elsewhere. We could resolve nothing, so send back to caller
      pass

    else:
      # call for sister resonances and see
      resons = resonanceSet.sortedResonances()
      newResonances = set()
      if len(resons)> 1:
        # there are sister resonances
        resons.remove(resonance)
        for reson in resons:
          boundResons = _getBoundResonances(reson, recalculate=True, contribs=contribs, recursiveCall=True)
          ll = [x for x in pairResonances if x not in boundResons]
          if not ll:
            # One sister was bound to both. Incorrect data. Bind to both here too
            newResonances.update(pairResonances)
            break
          elif len(ll) < len(pairResonances):
            # Some resonances were taken. Use the free ones.
            newResonances.update(ll)

      if newResonances:
        expResonances.update(newResonances)
      else:
        # No data anywhere to resolve which is which. Match on serials
        pairResonList = list(sorted(pairResonances, key=operator.attrgetter('serial')))
        rr = pairResonList[resonanceSet.sortedResonances().index(resonance)]
        expResonances.add(rr)


  resonances.update(expResonances)

  #if doWarning and (resonance.isotopeCode == '1H') and (len(resonances) > 1):
  #  pass

  if resonances:
    resonance.setCovalentlyBound(resonances)
  else:
    resonance.setCovalentlyBound([])

  return list(resonances)

