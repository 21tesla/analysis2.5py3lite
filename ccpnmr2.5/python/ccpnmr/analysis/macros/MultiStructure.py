from ccpnmr.integrator.core import Io as intIo

prelimProtocolName = 'MultipleStructures'


def runMultiProtocol(argServer, procNames=('CYANA_SS1','CYANA_PEAKLIST', 'ROSETTA_MULTI',
                                           'UNIO_CANDID','ASDP_XPLOR',),
                     prelimProtocolName=prelimProtocolName):
  """ Run multiple protocols in one run
  """
  intIo.setupMultiInteractive(argServer, procNames, prelimProtocolName,
                                executeProc=True)


def setupMultiProtocol(argServer, procNames=('CYANA_SS1','CYANA_PEAKLIST', 'ROSETTA_MULTI',
                                             'UNIO_CANDID','ASDP_XPLOR',),
                       prelimProtocolName=prelimProtocolName,
                       executeProc=False):
  """ Run multiple, connected protocols in one run
  """

  intIo.setupMultiInteractive(argServer, procNames, prelimProtocolName)



def runCandid(argServer):

  intIo.runSingleInteractive(argServer, protocolName='UNIO_CANDID_XPLOR_TEST',
                             prelimProtocolName=prelimProtocolName)


# def runCyanaPeaklist(argServer):
#
#   intIo.runSingleInteractive(argServer, protocolName='CYANA_PEAKLIST',
#                              prelimProtocolName=prelimProtocolName)

def setupCyanaCalculation(argServer):

  intIo.setupSingleInteractive(argServer, protocolName='CYANA_SS4',
                             )
def setupCyanaCalculationDialogue(argServer):

  intIo.setupCyana2CcpnDialogue(argServer, protocolName='CYANA_SS4',
                             )
def runCyana2CcpnDialogue(argServer):

  calculationData = intIo.runCyana2CcpnDialogue(argServer, protocolName='CYANA_SS4',
                             )
  return calculationData

def runPreviousCalculation(argServer):

  intIo.runPreviousCalculation(argServer, protocolName='CYANA_SS4')

def setupPreviousCalculation(argServer):

  intIo.setupPreviousCalculation(argServer, protocolName='CYANA_SS4')

def importDataFromCyana(argServer, calculationData=None):
  if calculationData == None:
    dataSources = intIo.importDataFromCyana(argServer)
  else:
    dataSources = intIo.importDataFromCyana(calculationData)
  return dataSources

def runCyanaPeaklist(argServer):

  intIo.runCyana2Ccpn(argServer, protocolName='CYANA_SS4',
                             )


def runAsdpXplor(argServer):

  intIo.runSingleInteractive(argServer, protocolName='ASDP_XPLOR',
                             prelimProtocolName=prelimProtocolName)

def setupCandid(argServer):

  project = argServer.getProject()
  run = None
  #runIndex = None
  runIndex = 55
  if runIndex:
    for store in reversed(project.sortedNmrCalcStores()):
      run = store.findFirstRun(serial=runIndex)
      if run:
        break

  intIo.setupSingleInteractive(argServer, protocolName='UNIO_CANDID',
                               prelimProtocolName=prelimProtocolName,
                               masterRun=run)


def setupCyanaPeaklist(argServer):

  project = argServer.getProject()
  run = None
  #runIndex = None
  runIndex = 88
  if runIndex:
    for store in reversed(project.sortedNmrCalcStores()):
      run = store.findFirstRun(serial=runIndex)
      if run:
        break

  intIo.setupSingleInteractive(argServer, protocolName='CYANA_PEAKLIST',
                               prelimProtocolName=prelimProtocolName,
                               masterRun=run)



def setupAsdpXplor(argServer):

  project = argServer.getProject()
  run = None
  #runIndex = None
  runIndex = 88
  if runIndex:
    for store in reversed(project.sortedNmrCalcStores()):
      run = store.findFirstRun(serial=runIndex)
      if run:
        break

  intIo.setupSingleInteractive(argServer, protocolName='ASDP_XPLOR',
                               prelimProtocolName=prelimProtocolName,
                               masterRun=run)


def setupRosetta(argServer):

  intIo.setupSingleInteractive(argServer, protocolName='ROSETTA_MULTI',
                               prelimProtocolName=prelimProtocolName)

