import sys
import smtplib

from email.MIMEText import MIMEText

def getShiftList(entry):

  shiftList = entry.findFirstMeasurementList(className='ShiftList')

  return shiftList


def getXEasyObject(project):

  from ccpnmr.format.converters.XEasyFormat import XEasyFormat

  xeasyObj = XEasyFormat(project, None, allowPopups=False)

  return xeasyObj


def writeShiftList(xeasyObj, xeasyFileName, shiftList, entry):

  xeasyObj.writeShifts(xeasyFileName,
                       measurementList=shiftList,
                       chains=ccpnEntry.molSystem.sortedChains(),
                       forceDefaultChainMapping=True,
                       minimalPrompts=True,
                       verbose=False)


def getPeakList(entry):

  peakList = entry.findFirstPeakList()

  return peakList


def writePeakList(xeasyObj, xeasyPeakFileName, peakList):

  xeasyObj.writePeaks(xeasyPeakFileName,
                      peakLists=[peakList],
                      minimalPrompts=True,
                      verbose=False)


def getStructureGeneration(entry):

  strucGen = None

  if entry.structureGenerations:
    strucGen = entry.findFirstStructureGeneration()

  return strucGen


def getStructures(strucGen):

  strucEns = strucGen.structureEnsemble

  structures = []

  if strucEns:
    structures = strucEns.sortedModels()

  return structures


def getPdbObject(project):

  from ccpnmr.format.converters.PdbFormat import PdbFormat

  pdbObj = PdbFormat(project, None, allowPopups=False)

  return pdbObj


def writeStructures(pdbObj, pdbFileName, structures):

  pdbObj.writeCoordinates(pdbFileName,
                          structures=structures,
                          resetMapping=False,
                          minimalPrompts=True,
                          verbose=True,
                          addPdbHeader=True,
                          version='3.1')


def getDistConstrainList(strucGen):

  nmrConstraintStore = strucGen.nmrConstraintStore

  distConstraintList = None

  if nmrConstraintStore:
    distConstraintList = nmrConstraintStore.findFirstConstraintList(className='DistanceConstraintList')

  return distConstraintList


def getCnsObject(project):

  from ccpnmr.format.converters.CnsFormat import CnsFormat

  cnsObj = CnsFormat(project, None, allowPopups=False)

  return cnsObj


def writeRestraints(cnsObj, cnsFileName, distConstraintList):

  cnsObj.writeDistanceConstraints(cnsFileName,
                                  constraintList=distConstraintList,
                                  minimalPrompts=True,
                                  verbose=False)


def writeNmrStarFile(entry, nmrStarFileName, nmrStarVersion = '3.0'):

  from pdbe.nmrStar.IO.NmrStarExport import NmrStarExport

  nmrStarExport = NmrStarExport(entry, nmrStarVersion=nmrStarVersion)

  nmrStarExport.createFile(nmrStarFileName)
  nmrStarExport.writeFile()


def checkCurrentTopObjects(project):

  if not project.currentNmrEntryStore:
    if project.findFirstNmrEntryStore():
      project.currentNmrEntryStore = project.findFirstNmrEntryStore()

  if not project.currentNmrProject:
    if project.findFirstNmrProject():
      project.currentNmrProject = project.findFirstNmrProject()

  if not project.currentStructureEnsemble:
    if project.findFirstStructureEnsemble():
      project.currentStructureEnsemble = project.findFirstStructureEnsemble()

  if not project.currentNmrConstraintStore:
    if project.findFirstNmrConstraintStore():
      project.currentNmrConstraintStore = project.findFirstNmrConstraintStore()


def mailErrorAndExit(errValue, errString):

  msg = MIMEText(errString)

  me = 'penkett@ebi.ac.uk'
  you = 'penkett@ebi.ac.uk'

  msg['Subject'] = 'Error in exporting CCPN project to NmrStar file'
  msg['From'] = me
  msg['To'] = you

  server = smtplib.SMTP()
  server.connect()

  #server.sendmail(me, you.split(', '), msg.as_string() )
  server.sendmail(me, [you], msg.as_string() )

  server.close()

  sys.exit(errValue)


def mailUploadInfo(logFile):

  fp = open(logFile, 'r')
  msg = MIMEText(fp.read() )
  fp.close()

  me = 'penkett@ebi.ac.uk'
  you = 'penkett@ebi.ac.uk'

  msg['Subject'] = 'CCPN project exported to NmrStar file'
  msg['From'] = me
  msg['To'] = you

  server = smtplib.SMTP()
  server.connect()

  #server.sendmail(me, you.split(', '), msg.as_string() )
  server.sendmail(me, [you], msg.as_string() )

  server.close()


if __name__ == '__main__':

  from memops.general.Io import loadProject

  project = None

  errValue = 0
  errString = ''

  if len(sys.argv[1:]) != 6:
    errString = 'There are not six arguments, finishing.\n'
    errValue = -2

    mailErrorAndExit(errValue, errString)

  (ccpnProjectDir, nmrStarFileName, xeasyFileName, pdbFileName, cnsFileName, logFile) = sys.argv[1:]

  print 'Trying %s' % (ccpnProjectDir)

  sys.stdout = open(logFile, 'w')

  print 'Args: %s\n' % sys.argv[1:]

  project = loadProject(ccpnProjectDir)

  if not project:
    errString = 'Cannot load this CCPN project.\n'
    errValue = -3

  checkCurrentTopObjects(project)

  if not project.currentNmrEntryStore:
    errString = """No CCPN nmrEntryStore/Entry defined - \n
  use ECI (http://www.ebi.ac.uk/pdbe/docs/pdbe_nmr_deposition/eci.html)\n
  to create one within your project, then upload again.\n"""
    errValue = -4

    mailErrorAndExit(errValue, errString)

  ccpnEntries = None

  if project.currentNmrEntryStore:
    ccpnEntries = project.currentNmrEntryStore.sortedEntries()

  else:
    errString = 'Cannot read CCPN Entries.\n'
    errValue = -5

    mailErrorAndExit(errValue, errString)

  ccpnEntry = None

  if not ccpnEntries:
    errString = """No CCPN Entry defined - \n
  use ECI (http://www.ebi.ac.uk/pdbe/docs/pdbe_nmr_deposition/eci.html)\n
  to create one within your project, then upload again.\n"""
    errValue = -6

    mailErrorAndExit(errValue, errString)

  elif len(ccpnEntries) == 1:
    ccpnEntry = ccpnEntries[0]

  else:
    ccpnEntry = ccpnEntries[len(ccpnEntries)-1]
    #ccpnEntry = ccpnEntries[0]

    #for tempCcpnEntry in ccpnEntries:
      # TODO: set this tag in Fc/Analysis!
      #if tempCcpnEntry.findFirstApplicationData(application = self.appName, keyword = 'deposit', value = True):
        #ccpnEntry = tempCcpnEntry
        #break

    if not ccpnEntry:
      print 'Please select the relevant Entry for this deposition.\n'
      ccpnEntry = ccpnEntries[0]

    if not ccpnEntry:
      errString = 'Cannot find any CCPN Entries.\n'
      errValue = -7

      mailErrorAndExit(errValue, errString)

  if ccpnEntry:
    print '\nDoing Ccpn2NmrStar for entry: %s' % ccpnEntry.name

    #sys.__stdout__.write('NAME: [%s]\n' % ccpnEntry.name)

    shiftList = getShiftList(ccpnEntry)

    if shiftList:
      xeasyObj = getXEasyObject(project)

      writeShiftList(xeasyObj, xeasyFileName, shiftList, ccpnEntry)

      #peakList = getPeakList(entry)

      #if peakList
      #  writePeakList(xeasyObj, xeasyPeakFileName, peakList)

    else:
      print 'No chemical shifts.\n'

    strucGen = getStructureGeneration(ccpnEntry)

    if strucGen:
      structures = getStructures(strucGen)

      if structures:
        pdbObj = getPdbObject(project)

        writeStructures(pdbObj, pdbFileName, structures)

      else:
        print 'No structures.\n'

      distConstraintList = getDistConstrainList(strucGen)

      if distConstraintList:
        cnsObj = getCnsObject(project)

        writeRestraints(cnsObj, cnsFileName, distConstraintList)

      else:
        print 'No distance restraints.\n'

    else:
      print 'Cannot find a structureGeneration.\n'

    writeNmrStarFile(ccpnEntry, nmrStarFileName, nmrStarVersion = '3.0')

    project.saveModified()

    mailUploadInfo(logFile)

    sys.stdout = sys.__stdout__

  else:
    sys.stdout = sys.__stdout__

  sys.exit(errValue)
