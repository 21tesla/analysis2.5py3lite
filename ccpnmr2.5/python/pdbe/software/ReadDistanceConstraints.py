import os
import sys

from memops.general.Io import loadProject
from pdbe.software.violationStatistics import ViolationStatistics

class ReadDistanceConstraints:

  #os.environ['PYTHONPATH'] = '.:' + os.getcwd() + '/ccpn/python'
  #print 'ENVIRONMENT: [%s]' % os.environ.get('PYTHONPATH', '')

  def __init__(self, mr):

    self.mr = mr
    self.entry = getEntry(mr)

    self.strucGens = self.entry.sortedStructureGenerations()

    self.distConstLists = self.getDistanceConstraintLists(self.strucGens)

    self.assign_kw = 'assign'
    #self.format = 'nmrStar'


  def getViolationData(self):

    nmrProject = self.mr.currentNmrProject
    if not nmrProject:
      nmrProject = self.mr.findFirstNmrProject()

    strucGen = self.entry.findFirstStructureGeneration()

    violationStatistics = ViolationStatistics(nmrProject, strucGen=strucGen)
    violData = violationStatistics.getSummaryText()

    return violData


  def getDistanceConstraintData(self, violData):

    distConstData = []

    constType = 0 # Unknown constraint type

    for dcl in self.distConstLists:

      violDataDcl = None
      if dcl.serial in violData:
        violDataDcl = violData[dcl.serial]

      #print violDataDcl

      constType = self.getDistanceConstraintType(dcl)

      distConsts = dcl.sortedConstraints()

      ambi = None # Unknown ambiguity

      # Distance length (short < 2.5 A, medium < 4.0 A, long < 6.0 A)

      distType = 0 # Distance type not set

      for dc in distConsts:
        violDataDclDc = None

        if violDataDcl and dc.serial in violDataDcl:
          violDataDclDc = violDataDcl[dc.serial]

        items = dc.sortedItems()

        if dc.targetValue < 2.5:
          distType = 1 # short
        elif dc.targetValue < 4.0:
          distType = 2 # medium
        elif dc.targetValue < 6.0:
          distType = 3 # long

        if len(items) > 1:
          ambi = True
        else:
          ambi = False

        """
        violFlag = False
        violVal = 0.0
        violCalc = 0.0
        violError = 0.0
        violFrac = 0.0

        if dc.violations:
          violFlag = True
          viol = dc.findFirstViolation()

          violVal   = viol.violation
          violCalc  = viol.calcValue
          violError = viol.calcValueError
          violFrac  = viol.fractionViolated
        """

        for item in items:
          resonances = list(item.orderedResonances)
          if not resonances:
            resonances = item.sortedResonances()

          firstRes  = resonances[0]
          secondRes = resonances[1]

          authNames1 = [appData.value for appData in firstRes.findAllApplicationData(application='ccpNmr',keyword=self.assign_kw)]
          authNames2 = [appData.value for appData in secondRes.findAllApplicationData(application='ccpNmr',keyword=self.assign_kw)]

          firstAtoms  = self.getAtoms(firstRes)
          secondAtoms = self.getAtoms(secondRes)

          #if dc.serial == 9:
            #print 'ITEM: [%s]' % item
            #print 'RES: [%s] --- [%s]' % (firstRes, secondRes)
            #print 'APP: [%s] *** [%s]' % (firstRes.applicationData,
            #                              secondRes.applicationData)
            #print 'ATOMS: [%s] +++ [%s]' % (firstAtoms, secondAtoms)

          # Contact range (intra, sequential, medium range (i->i+4 and less), long range)

          rangeType = 0 # Contact type not set

          for atom1 in firstAtoms:
            res1 = atom1.residue
            chain1 = res1.chain
            for atom2 in secondAtoms:
              res2 = atom2.residue
              chain2 = res2.chain
              if chain1 != chain2:
                rangeType = 5 # inter range
              elif res1.seqCode == res2.seqCode:
                rangeType = 1 # intra
              elif abs(res1.seqCode - res2.seqCode) == 1:
                rangeType = 2 # sequential
              elif abs(res1.seqCode - res2.seqCode) <= 4:
                rangeType = 3 # medium range (i->i+4 and less)
              elif abs(res1.seqCode - res2.seqCode) > 4:
                rangeType = 4 # long range

              if violDataDclDc:
                violEns, violMean, violSd, violMax, count, count_0, count_1, count_3, count_5 = violDataDclDc

              else:
                violEns, violMean, violSd, violMax, count, count_0, count_1, count_3, count_5 = [0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0]

              rogScore = 'green'

              try:
                if float(violMax) >= 0.3 or float(violSd) >= 0.15:
                  rogScore = 'orange'

                if float(violMax) >= 0.5 or float(violSd) >= 0.3:
                  rogScore = 'red'

              except:
                pass

              targetValue = 0.0

              if type(dc.targetValue) == type(float() ):
                targetValue = dc.targetValue

              upperLimit = 0.0

              if type(dc.upperLimit) == type(float() ):
                upperLimit = dc.upperLimit

              lowerLimit = 0.0

              if type(dc.lowerLimit) == type(float() ):
                lowerLimit = dc.lowerLimit

              dataTuple = (dcl.serial, dcl.className,
                           dc.serial, constType, ambi,
                           distType, rangeType,
                           '%.3f' % targetValue,
                           '%.3f' % upperLimit,
                           '%.3f' % lowerLimit,
                           #authNames1,
                           atom1.name, res1.seqCode, res1.ccpCode, chain1.code,
                           #authNames2,
                           atom2.name, res2.seqCode, res2.ccpCode, chain2.code,
                           violEns, violMean, violSd, violMax, rogScore,
                           count, count_0, count_1, count_3, count_5)

              if dataTuple not in distConstData:
                distConstData.append(dataTuple)

    return distConstData


  def getAtoms(self, resonance):

    atoms = []

    resSet = resonance.resonanceSet
    if resSet:
      atomSets = resSet.sortedAtomSets()

      for atomSet in atomSets:
        for atom in atomSet.sortedAtoms():
          if atom not in atoms:
            atoms.append(atom)

    return atoms


  def getDistanceConstraintLists(self, structureGenerations):

    return self.getConstraintLists(structureGenerations,
                                   ['DistanceConstraintList',
                                    'HBondConstraintList'])


  def getConstraintLists(self, structureGenerations, classNames):

    constraintLists = []

    if structureGenerations:
      for strucGen in structureGenerations:
        for className in classNames:
          if hasattr(strucGen, 'nmrConstraintStore') and strucGen.nmrConstraintStore:
            constraintLists.extend(list(strucGen.nmrConstraintStore.findAllConstraintLists(className=className) ) )

    return constraintLists


  def getDistanceConstraintType(self, constraintList):

    constraintType = 1 # 'NOE' or 'general distance'?

    if constraintList.className == 'HBondConstraintList':
      constraintType = 2 # H-Bond

    else:
      for experiment in constraintList.sortedExperiments():
        if experiment.refExperiment:
          if experiment.refExperiment.name.count('ROESY'):
            constraintType = 3 # 'ROE'

    return constraintType


  def getAllConstraintItems(self, constraintList):

    constraintItems = []

    for constraint in constraintList.sortedConstraints():
      for constraintItem in constraint.sortedItems():
        constraintItems.append(constraintItem)

    return constraintItems


  def getGeneralConstraintLogic(self, constraint):

    constraintLogic = None

    if len(constraint.items) > 1:
      constraintLogic = 'OR'

    return constraintLogic


def getEntry(memRoot):

  if not memRoot.currentNmrEntryStore and memRoot.findFirstNmrEntryStore():
    memRoot.currentNmrEntryStore = memRoot.findFirstNmrEntryStore()

  if memRoot.currentNmrEntryStore is not None:
    return memRoot.currentNmrEntryStore.findFirstEntry()


if __name__ == '__main__':

  ccpnProjectName, = sys.argv[1:]

  print ccpnProjectName

  mr = loadProject(ccpnProjectName)

  readDistConst = ReadDistanceConstraints(mr)

  violData = readDistConst.getViolationData()

  distConstData = readDistConst.getDistanceConstraintData(violData)

  for line in distConstData:
    print line
