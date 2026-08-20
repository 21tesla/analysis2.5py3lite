
"""
======================COPYRIGHT/LICENSE START==========================

Io.py: code for CCPN data model and code generation framework

Copyright (C) 2011  (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

A copy of this license can be found in ../../../license/LGPL.license

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

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""

import operator, os, json, sys, subprocess, traceback, time, shutil, tarfile
from memops.universal import Io as uniIo
from memops.general import Io as genIo
from ccpnmr.integrator.core import Util as intUtil
from ccpnmr.integrator.core import ParameterEditor
from cyana2ccpn.cyana2ccpn import importFromCyana

propFileName = 'Properties.ccpn.json'
propIndent = 2


integratorDir = 'integrator'
dirNameSeparator = '.'

timeFormat = '%Y-%b-%d-%X'


setupLogFile = 'ccpnSetup.log'

def makeJsonDict(containerObj, direction='input'):
  """Make dictionary of simple type data, selecting those with a code value.
  Values are grouped in dictionaries, listed by Data object type. Each
  dictionary has a serial key, with the Data object serial.  Non-linked
  parameters are in the top dictionary. Will treat repeated codes as a list of
  parameters. No way to specify lists with a single member. NB Should use an
  ordered dict, but this is not yet available.

  Input:
    containerObj: NmrCalc.Run or NmrCalc.ParameterGroup
    direction: ('input'|'output'|None)

  Output: JSON-ready dictionary

  """

  print 'Prepring JSON dict'

  result = {}
  result['RunParameter'] = paramDict = {}

  # CCPN-specific entries:
  if containerObj.className == 'Run':
    topLevel = True
    result['CCPN.NmrCalcStore.name'] = containerObj.nmrCalcStore.name
    result['CCPN.NmrProject.name'] = containerObj.nmrCalcStore.nmrProjectName
    result['CCPN.nmrCalcId'] = intUtil.getNmrCalcIdentifier(containerObj)
    result['CCPN.Run.wmsProtocolName'] = containerObj.wmsProtocolName
    masterRun = containerObj.masterRun
    if masterRun is not None:
      result['CCPN.masterCalcId'] = intUtil.getNmrCalcIdentifier(masterRun)

  else:
    topLevel = False

  # Unconnected parameters
  if direction is None:
    params = list(containerObj.findAllRunParameters(data=None))
  else:
    params = list(containerObj.findAllRunParameters(ioRole=direction, data=None))
  params.sort(key=operator.attrgetter('serial'))
  for param in params:
    code = param.code
    if code is not None:
      # Has code. Put in RunParameters dictionary
      intUtil.addParamToDict(code, intUtil.getParameterValue(param), paramDict)
    elif topLevel:
      # No code. Put in top level.
      name = param.name
      intUtil.addParamToDict(name, intUtil.getParameterValue(param), result)

  # Data objects
  # paramGroups = []
  for datum in containerObj.sortedData():
    if datum.parameterGroup is None or not topLevel:

      # get className and make list to put data
      className = datum.className
      ll = result.get(className)
      if ll is None:
        ll = result[className] = []

      if (className == 'ParameterGroup'
          and datum.parameterGroup in (None, containerObj)):
        # ParameterGoup object that belongs here.
        ll.append(makeJsonDict(datum, direction=direction))

      else:
        # Normal Data object
        paramDict= {'serial':datum.serial}
        ll.append(paramDict)

        # put parameter values in paramDict
        if direction is None:
          params = list(datum.runParameters)
        else:
          params = list(datum.findAllRunParameters(ioRole=direction))
        params.sort(key=operator.attrgetter('serial'))
        for param in params:
          code = param.code
          if code is not None:
            # Has code. Put in RunParameters dictionary
            intUtil.addParamToDict(code, intUtil.getParameterValue(param),
                                   paramDict)
  #
  return result

def makeParameterDict(nmrCalcRun, direction='input'):
  """Make code:value dictionary with simple type data only,
     and list of parameter codes
  Input:
    run: NmrCalc.Run
    direction: ('input'|'output'|None)

  Will treat repeated codes as a list of parameters.
  No way to specify lists with a single member.
  NB Should use an ordered dict, but this is not yet available.

  NBNB OBSOLETE

  """

  raise Exception("CCPN Obsolete function called")

  paramDict = {}

  # CCPN-specific entries:
  paramDict['CCPN.NmrCalcStore.name'] = nmrCalcRun.nmrCalcStore.name
  paramDict['CCPN.NmrCalcIdentifier'] = intUtil.getNmrCalcIdentifier(nmrCalcRun)

  codes = ['CCPN.NmrCalcStore.name', 'CCPN.NmrCalcIdentifier']

  if direction:
    ll = list(nmrCalcRun.findAllRunParameters(ioRole=direction))
    ll.sort(key=operator.attrgetter('serial'))
  else:
    ll = nmrCalcRun.sortedRunParameters()

  for param in ll:
    code=param.code

    if code:
      codes.append(code)
      value = intUtil.getParameterValue(param)
      xx = paramDict.get(code)
      if xx:
        # list of parameters
        if isinstance(xx, list):
          xx.append(value)
        else:
          paramDict[code] = [xx, value]
      else:
        paramDict[code] = value

  #
  return paramDict, codes

def nmrCalcDir(nmrCalcRun, dirSuffix='in', topDir=None, compact=False):
  """ get target directory for nmrCalc.Run input or output

  Input:  nmrCalcRun: NmrCalc.Run
          dirSuffix:  ('in'|'out'|None)
          topDir:     Directory to contain data
                      (defaults to projectDir/integratorDir)
  """

  if nmrCalcRun is None:
    raise Exception('cannot handle null nmrCalcRun')

  if topDir is None:
    # get default location
    project = nmrCalcRun.root
    topDir = project.getPackageLocator().findFirstRepository().url.getDataLocation()
    topDir = uniIo.joinPath(topDir,integratorDir)

  # make run-specific directory
  if compact:
    dataDirName = intUtil.objectFileName(nmrCalcRun, compact=True)
  else:
    dataDirName = intUtil.getNmrCalcIdentifier(nmrCalcRun)


  if dirSuffix:
    dataDirName = dirNameSeparator.join((dataDirName,dirSuffix))
  targetDir = uniIo.joinPath(topDir, dataDirName)

  #
  return targetDir

def getNmrCalcRunFromId(project, nmrCalcRunId):
  """ Get NmrCalc.Run corresponding to nmrCalcRunId, if any, from project.
  """
  nmrCalcRun = None
  tt = nmrCalcRunId.split('+')
  if len(tt) == 2:
    guid, ss = tt
    runSerial = int(ss)

    nmrCalcStore = project.findFirstNmrCalcStore(guid=guid)
    if nmrCalcStore is not None:
      nmrCalcRun = nmrCalcStore.findFirstRun(serial=runSerial)
  #
  return nmrCalcRun

def getNmrCalcRun(projectDir, nmrCalcRunId, pluginName=None):
  """ Open project and get NmrCalc.run from project dir and NmrCalcId string
  If pluginName is set, replace run with daughter run of correct protocol,
  if any
  """

  # find project (first it comes across)
  for dirpath, dirnames, filenames in os.walk(projectDir):
    if dirpath.endswith('/memops/Implementation'):
      projDir = (os.path.dirname(os.path.dirname(dirpath)))
      project = genIo.loadProject(projDir)
      break
  else:
    return

  # get NmrCalc.Run
  nmrCalcRun = getNmrCalcRunFromId(project, nmrCalcRunId)

  if nmrCalcRun is not None and pluginName is not None:
    derivedRuns = nmrCalcRun.sortedDerivedRuns()
    if derivedRuns:

      # Check if this run fits the protocol name
      tag = pluginName.split('.')[-1].upper()
      #if tag not in nmrCalcRun.wmsProtocolName.upper():

      # see if there is a derived run that does fit.
      for run in derivedRuns:
        if tag in run.wmsProtocolName.upper():
          newRun = run
          break
      else:
        newRun = None

      if newRun is not None:
        # we have a run that fits. Use it.
        # NB the output is NECESSARY.
        # It is used to transfer information to WMS Java
        print "CCPN_new_calcId = '%s'" % intUtil.getNmrCalcIdentifier(newRun)

        return newRun

  #
  return nmrCalcRun


def setStdParameters(protocol, nmrCalcRun, interface=None):
  """ Set NmrCalc parameters with standard types and values
  """

  for paramObj in protocol.sortedProtocolParameters():
    interfaceObj = (interface and
                interface.findFirstInterfaceParameter(protocolParameter=paramObj))
    setStdParameter(paramObj, nmrCalcRun, interfaceObj)

def setStdParameter(paramObj, nmrCalcRun, interfaceObj=None):
  """ Set NmrCalc parameter with standard types and values
  """

  # NBNB TBD to be expanded. This covers a subset of cases
  defaultStrings = None
  if interfaceObj:
    defaultStrings = interfaceObj.defaultStrings
  if not defaultStrings:
    defaultStrings = paramObj.defaultStrings

  if len(defaultStrings) == 1 and paramObj.container is None:
    default = defaultStrings[0]
    paramType = paramObj.paramType
    converter = intUtil.paramTypeConverters.get(paramType)
    if default is not None and converter is not None:
      param = nmrCalcRun.newRunParameter(name=paramObj.name, code=paramObj.code)
      setattr(param, paramType, converter(default))


def loadProtocol(memopsRoot, jsonFile, interfaceName=None, overwrite=False):
  """ Load a json Protocol definition into a WmsProtocol

  NB Assumes all default strings are single values

  NB copies hicard (but not locard) from ProtocolParameter to InterfaceParameter
  if not set in the latter.

  NB modified RHF 3/10/2013 for Integrator requirements
  """

  jsonObject = json.load(open(jsonFile))

  protocoldd = jsonObject['protocol']
  name = protocoldd['name']

  print 'Loading new protocol: %s' % name

  protocol = memopsRoot.findFirstWmsProtocol(name=name)
  if protocol is not None:
    if overwrite:
      protocol.delete()
    else:
      raise Exception('protocol named %s already exists!' % name)
  details = protocoldd.get('details')
  if details:
    # Cut down to allowed maximum length
    details = details[:254]
  protocol = memopsRoot.newWmsProtocol(name=name, details=details)

  # load parameters
  params = protocoldd['protocolParameters']
  for pardd in params:
    dd = {}
    for tag in ('name', 'code', 'locard', 'hicard',):
      ss = pardd.get(tag)
      if ss is not None:
        dd[tag] = ss

    ss = pardd['paramType']
    dd['paramType'] = intUtil.paramTypeMap.get(ss,ss)

    ss = pardd.get('value')
    if ss is not None:
      dd['defaultStrings'] = (str(ss),)

    newPar = protocol.newProtocolParameter(**dd)

#    #name = pardd.get('name')
    #if name is not None:
    #  parObjs[name] = newPar

  # set parameter relations:
  for pardd in params:
    relatedParameter = pardd.get('relatedParameter')
    if relatedParameter is not None:
      newPar = protocol.findFirstProtocolParameter(name=pardd['name'])
      container = protocol.findFirstProtocolParameter(name=relatedParameter)
      newPar.container = container

#  # set interfaceParameters
  interfaces = protocoldd['protocolInterfaces']
  for interfacedd in interfaces:
    if interfaceName is not None and interfacedd['name'] != interfaceName:
      continue

    name = interfacedd['name']
    title = interfacedd.get('title') or name
    interface = protocol.newProtocolInterface(name=name, title=title)
    for tag in ('info', 'details'):
      val = interfacedd.get(tag)
      if val is not None:
        setattr(interface, tag, val)

    #
    for tabdd in interfacedd['tabs']:
      print tabdd
      #if tabdd.get('io') != 'out':
      interfaceGroup = tabdd.get('name')

      ll = tabdd.get('interfaceParameters', ())
      for pardd in ll:

        dd = {'interfaceGroup':interfaceGroup}
        for tag in ('isOrdered', 'isEditable', 'locard', 'hicard',
                    'row', 'col', 'rowspan', 'colspan'):
          val = pardd.get(tag)
          if val is not None:
            dd[tag] = val
        param = protocol.findFirstProtocolParameter(name=pardd['name'])
        dd['protocolParameter'] = param

        if not param:
          raise Exception("Parameter named %s not found" % pardd['name'])

        if 'hicard' not in pardd:
          dd['hicard'] = param.hicard

        ss = pardd.get('value')
        if ss is not None:
          dd['defaultStrings'] = (str(ss),)
        #
        xx = interface.newInterfaceParameter(**dd)

      ll = tabdd.get('interfaceLabels', ())
      for pardd in ll:

        dd = {'interfaceGroup':interfaceGroup}
        for tag in ('label', 'row', 'col', 'rowspan', 'colspan'):
          val = pardd.get(tag)
          if val is not None:
            dd[tag] = val
        #
        xx = interface.newInterfaceLabel(**dd)
  #
  return protocol


def interfaceTabs(protocolInterface):
  """ Return list of interfaceGroup strings in protocol, in order of parameter
  creation. Assumes that
  1) interface group names are unique.
  2) Parameters are created in interface group order.
  3) groups without parameters (if any) are appended at the end of the list
  """

  result = []

  for xx in (protocolInterface.sortedInterfaceParameters() +
             protocolInterface.sortedInterfaceLabels()):
    ss = xx.interfaceGroup
    if ss not in result:
      result.append(ss)
  #
  return result

def parseResidueExpr(molSystem, selector):
  """Convert residue Expression (e.g. A:;B:-4,6,9,11-13,45-;C:)
     to list of residues from molSystem
  """

  result = []

  data = {}

  chainExprs = selector.split(';')
  for chainExpr in chainExprs:
    code, resExpr = (x.strip() for x in chainExpr.split(':',1))

    chain = molSystem.findFirstChain(code=code)
    if chain is None:
      raise Exception("%s contains no chain with code %s" % (molSystem, code))

    if resExpr:

      seqIds = parseNumberExpr(resExpr,
                               endat=max(x.seqId for x in chain.residues))
      result.extend(chain.findFirstResidue(seqId=ii) for ii in seqIds)
    else:
      result.extend(chain.sortedResidues())
  #
  return result

def parseNumberExpr(selector, startat=1, endat=None):
  """Parse numer selection string of the form
  '-5, 9, 11, 13-17, 19, 25-'
  """

  result = []

  strs = [x.strip() for x in selector.split(',')]
  if strs[0].startswith('-'):
    strs[0] = str(startat)+strs[0]
  if strs[-1].endswith('-'):
    if endat is None:
      raise Exception(
            "max mumber must be specified to allow selector ending with '-'")
    strs[-1] = strs[-1] + str(endat)


  for ss in strs:
    ll = ss.split('-',1)
    if len(ll) == 1:
      result.append(int(ll[0]))
    else:
      result.extend(range(int(ll[0]),int(ll[-1])+1))
  #
  return result

def initRunInteractive(argServer, protocolName=None, prelimProtocolName=None,
                        masterRun=None):
  """ Interactive initialisation of NmrCalc Run starting from protocol names
  """
  if masterRun is not None:
    # copy template run -
    # e.g. because we have entered data already for multistructure run
    nmrCalcRun = intUtil.makeDerivedRun(masterRun)

  elif prelimProtocolName is None:
    nmrCalcRun = None

  else:
    # First set up using preliminary protocol
    nmrCalcRun = setupRunInteractive(argServer, protocolName=prelimProtocolName)

  if protocolName is None:
    protocolName = argServer.askString("Program/protocol name", '?')

  nmrCalcRun = setupRunInteractive(argServer, protocolName=protocolName,
                                   nmrCalcRun=nmrCalcRun)
  #
  return nmrCalcRun

def initRunDialogue(argServer, protocolName=None, prelimProtocolName=None,
                        masterRun=None):
  """ Interactive initialisation of NmrCalc Run starting from protocol names
  """
  if masterRun is not None:
    # copy template run -
    # e.g. because we have entered data already for multistructure run
    nmrCalcRun = intUtil.makeDerivedRun(masterRun)

  elif prelimProtocolName is None:
    nmrCalcRun = None

  else:
    # First set up using preliminary protocol
    nmrCalcRun = setupRunDialogue(argServer, protocolName=prelimProtocolName)

  if protocolName is None:

    protocolName = argServer.askString("Program/protocol name", '?')

  nmrCalcRun = setupRunDialogue(argServer, protocolName=protocolName,
                                   nmrCalcRun=nmrCalcRun)
  #
  return nmrCalcRun

def getWmsInteractive(argServer, protocolName=None):
  """ Interactive get WmsProtocol from protocol name
  """
  from ccp.general.Io import getDataPath

  project = argServer.getProject()

  entryProtocol = project.findFirstWmsProtocol(name=protocolName)

  if entryProtocol is None:
    dataPath = getDataPath()
    jsonFile = os.path.join(dataPath,'ccpnmr', 'integrator', protocolName+'.json')
    # print 'Enter name of %s protocol definition file' % protocolName
    # jsonFile = argServer.getFile()
    entryProtocol = loadProtocol(project, jsonFile)
  #
  return entryProtocol

def initNmrCalcRun(entryProtocol, nmrCalcRun=None, nmrCalcStore=None,
                   interface=None):
  """initialise NmrCalc.Run from entryProtocol, creating new if necessary
  """

  project = entryProtocol.root

  if nmrCalcRun is None:
    if not nmrCalcStore:
      nmrCalcStore = intUtil.getNmrCalcStore(project)

    nmrCalcRun = nmrCalcStore.newRun(details=entryProtocol.name + " Run",
                                   wmsProtocolName=entryProtocol.name)
  else:
    # reset protocol name, as we are overriding
    nmrCalcRun.wmsProtocolName = entryProtocol.name

  setStdParameters(entryProtocol, nmrCalcRun, interface=interface)
  #
  return nmrCalcRun


def setupRunDialogue(argServer, protocolName=None, nmrCalcRun=None):
  """ Interactive set up (if necessary create) nmrCalcRun
  according to named protocol
  If nmrCalcRun is passed in, new parameters are added to it.
  """

  project = argServer.getProject()

  entryProtocol = getWmsInteractive(argServer, protocolName=protocolName)

  nmrCalcRun = initNmrCalcRun(entryProtocol, nmrCalcRun)

  interface = (entryProtocol.findFirstProtocolInterface(name=entryProtocol.name) or
               entryProtocol.findFirstProtocolInterface())
  # map to handle protocol/interfaceparameter connections
  paramMap = {}
  for ip in interface.sortedInterfaceParameters():
    pp = ip.protocolParameter
    if pp in paramMap:
      raise Exception("Parameter %s in interface twice" % pp)
    paramMap[pp] = ip

  popup = ParameterEditor.ParameterEditor(argServer.parent, nmrCalcRun, interface, modal=True)
  nmrCalcRun = popup.run # ParameterEditor sets to None if Cancel button hit; otherwise it is above nmrCalcRun if Ok button hit
  popup.destroy()

  return nmrCalcRun
  
  
def setupRunInteractive(argServer, protocolName=None, nmrCalcRun=None):
  """ Interactive set up (if necessary create) nmrCalcRun
  according to named protocol
  If nmrCalcRun is passed in, new parameters are added to it.
  """

  project = argServer.getProject()

  entryProtocol = getWmsInteractive(argServer, protocolName=protocolName)

  nmrCalcRun = initNmrCalcRun(entryProtocol, nmrCalcRun)

  interface = (entryProtocol.findFirstProtocolInterface(name=entryProtocol.name) or
               entryProtocol.findFirstProtocolInterface())
  # map to handle protocol/interfaceparameter connections
  paramMap = {}
  for ip in interface.sortedInterfaceParameters():
    pp = ip.protocolParameter
    if pp in paramMap:
      raise Exception("Parameter %s in interface twice" % pp)
    paramMap[pp] = ip

  # Get input:
  nmrCalcObj = None

  for ip in interface.sortedInterfaceParameters():
    print ip, ip.protocolParameter.name
    if ip.isEditable:

      pp = ip.protocolParameter

      if pp.container is None:

        # get contained parameters
        ll = [ip]
        count = 1
        for pp2 in pp.sortedContent():
          ip2 = paramMap.get(pp2)
          if ip2 is not None:
            ll.append(ip2)

        # set values
        while ip.hicard < 0 or count <= ip.hicard:
          if count > 1 and count > ip.locard:

            xx = argServer.askYesNo("Set another %s?" % pp.name)
            if not xx:
              break
          count += 1
          dataObj = None
          for ip2 in ll:
            nmrCalcObj = makeNmrCalcData(argServer, nmrCalcRun,
                                         ip2.protocolParameter,
                                         dataObj=dataObj)
            if nmrCalcObj is None:
              break
            elif nmrCalcObj.className == 'RunParameter':
              if dataObj is not None:
                nmrCalcObj.data = dataObj
            else:
              dataObj = nmrCalcObj

          if nmrCalcObj is None:
            break

  #
  return nmrCalcRun

def makeNmrCalcData(argServer, run, protocolParameter, ioRole='input',
                    dataObj=None):
  """ Make NmrCalc data or parameter for single ProtocolParameter
  """


  paramType = protocolParameter.paramType
  name = protocolParameter.name
  parDD = {'name':name, 'code':protocolParameter.code, 'ioRole':ioRole}

  if paramType in ('intValue', 'floatValue', 'textValue', 'booleanValue', ):

    defaults = protocolParameter.defaultStrings
    if defaults:
      defValue = intUtil.paramTypeConverters[paramType](defaults[0])
    else:
      defValue = None

    if paramType == 'intValue':
      value = argServer.askInteger("%s?" % name, defValue)

    elif paramType == 'floatValue':
      value = argServer.askFloat("%s?" % name, defValue)

    elif paramType == 'textValue':
      value = argServer.askString("%s?" % name, defValue)

    elif paramType == 'booleanValue':
      value = argServer.askYesNo("%s Boolean value?" % name)

    else:
      value=None
      # NBNB should not get here

    if value is None:
      return

    parDD[paramType] = value
    result = None
    # check if alredy set, and overwrite:
    if protocolParameter.hicard == 1:
      # NBNB there could also be problemd with hicard != 1,
      # but we ignore that for now
      result = run.findFirstRunParameter(name=name, data=dataObj)

    if result is None:
      result = run.newRunParameter(**parDD)

    else:
      # NB will cause trouble if the parameter had data set for wrong paramType,
      # but that should never happen
      for tag, val in parDD.items():
        setattr(result, tag, val)

  elif paramType == 'molResidues':
    molSystem = argServer.getMolSystem()
    if molSystem is None:
      return

    print 'Residue specification style: A:;B:-4,6,9,11-13,45-;C:'
    ss = ':;'.join(x.code for x in molSystem.sortedChains()) + ':'
    selector = argServer.askString("Residue selection", ss)
    if selector == ss:
      result = run.newMolResidueData(chains=molSystem.sortedChains(), **parDD)

    elif selector is None:
      return

    else:
      residues = parseResidueExpr(molSystem, selector)
      result = run.newMolResidueData(residues=residues, **parDD)

  elif paramType == 'measurementList':
    measurementList = argServer.getMeasurementList()
    if measurementList is None:
      return
    result = run.newMeasurementListData(measurementList=measurementList, **parDD)

  elif paramType == 'peakList':
    peakList = argServer.getPeakList()
    if peakList is None:
      return
    result = run.newPeakListData(peakList=peakList, **parDD)

  elif paramType == 'molSysChains':
    molSystem = argServer.getMolSystem()
    if molSystem is None:
      return

    print 'Chain code style: A,B,D'
    ss = ','.join(x.code for x in molSystem.sortedChains())
    selector = argServer.askString("Chain codes", ss)
    if selector == ss:
      result = run.newMolSystemData(chains=molSystem.sortedChains(), **parDD)

    elif selector is None:
      return

    else:
      chainCodes = [x.strip() for x in selector.split(',')]
      result = run.newMolSystemData(molSystemCode=molSystem.code,
                                    chainCodes=chainCodes, **parDD)

  elif paramType == 'structureModels':
    structureEnsemble = argServer.getStructure()
    if peakList is structureEnsemble:
      return
    print "Model number style: '-4,6,9,11-13,45-'"
    maxModels = max(x.serial for x in structureEnsemble.models)
    ss = '1-%d' % maxModels
    selector = argServer.askString("Model selection", ss)
    if selector == ss:
      result = run.newStructureEnsembleData(
                      models=structureEnsemble.sortedModels(), **parDD)
    elif selector is None:
      return

    else:
      modelNums = parseNumberExpr(selector, startat=1, endat=maxModels)
      result = run.newStructureEnsembleData(structureEnsemble=structureEnsemble,
                                            modelSerials=sorted(modelNums),
                                            **parDD)

  elif paramType == 'constraintLists':
    constraintStore = argServer.getConstraintSet()
    if constraintStore:
      print "ConstraintList number style: '-4,6,9,11-13,45-'"
      ss = ','.join(str(x.serial) for x in constraintStore.sortedConstraintLists())
      selector = argServer.askString("Constraint list selection", ss)
      if selector == ss:
        result = run.newConstraintStoreData(
                        constraintLists=constraintStore.sortedConstraintLists(),
                        **parDD)
      else:
        maxClists = max(x.serial for x in constraintStore.constraintLists)
        clistSerials = parseNumberExpr(selector, startat=1, endat=maxClists)
        result = run.newConstraintStoreData(nmrConstraintStore=constraintStore,
                                            constraintListSerials=clistSerials,
                                            **parDD)
    else:
      return
  else:
    raise Exception("Unknown parameter type: %s for %s" %
                    (paramType, protocolParameter))

  #
  return result


def writeDataFiles(nmrCalcRun, targetDir):
  """ Write input data files for Program run
    Input:
      nmrCalcRun: NmrCalc.Run
      targetDir: destination directory.
  """

  print 'Writing data files to', targetDir

  # initialise shared resonanceToAtoms mapping
  # resonanceToAtoms = None
  # fixResonanceToAtoms = None

  # get I/O handler
  dataIoHandler = intUtil.DataIoHandler(nmrCalcRun.root, nmrCalcRun)

  # set local keywords
  myKeywds = {}
  obj = nmrCalcRun.findFirstRunParameter(name='atomNamingSystem', data=None)
  if obj is not None:
    myKeywds['forceNamingSystemName'] = obj.textValue

  for tag in ('useXeasyDimCodes','skipMultiAssignments'):
    obj = nmrCalcRun.findFirstRunParameter(name=tag, data=None)
    if obj is not None:
      myKeywds[tag] = obj.booleanValue

  # set variables for resonance mapping
  presetInput = {}
  presetInput['namingSystemName'] = myKeywds.get('forceNamingSystemName')
  seqDataObj = nmrCalcRun.findFirstData(name='definedResidues')
  if seqDataObj:
    presetInput['molSystem'] = seqDataObj.molSystem

  # write sequence file
  xx = nmrCalcRun.findFirstRunParameter(name='fileFormatSequence')

  if xx is not None:
    fileFormat = xx.textValue
    fileName = nmrCalcRun.findFirstRunParameter(name='fileNameSequence').textValue
    filePath = uniIo.joinPath(targetDir, fileName)
    seqDataObj = nmrCalcRun.findFirstData(name='definedResidues')
    chain = seqDataObj.chain

    handler = dataIoHandler.getIoHandler(fileFormat, dataType='seq')
    keywds = intUtil.getFileFormatData(fileFormat, dataType='seq').get(
                'IOkeywords', {}).copy()
    keywds.update(myKeywds)
    handler.writeSequence(filePath, chains=(chain,), **keywds)

    print '### seq', fileName, fileFormat


  # write shiftLists
  for obj in nmrCalcRun.sortedData():
    if obj.name == 'shiftList' and obj.className == 'MeasurementListData':
      # ShiftList data

      xx = obj.findFirstRunParameter(name='fileFormat')
      if xx is not None:
        fileFormat = xx.textValue

        #for fileFormat in ('XEASY', 'NMRVIEW','SPARKY', 'NMRSTAR', ):

        #fileName = fileFormat + obj.findFirstRunParameter(name='fileName').textValue
        fileName = obj.findFirstRunParameter(name='fileName').textValue
        filePath = uniIo.joinPath(targetDir, fileName)

        handler = dataIoHandler.getIoHandler(fileFormat, dataType='shift')
        keywds = intUtil.getFileFormatData(fileFormat, dataType='shift').get(
                   'IOkeywords',{}).copy()
        keywds.update(myKeywds)

        if keywds.get('presetResonanceMapping'):
        # get global resonanceToAtoms map

          if 'resonanceToAtoms' not in presetInput:
            # get resonanceToAtoms mapping for presetting
            namingSystemName = keywds.get('forceNamingSystemName', 'IUPAC')
            # get nmrConstraintStore
            constraintData = nmrCalcRun.findFirstData(name='assignment',
                                           className='ConstraintStoreData')
            if constraintData:
              nmrConstraintStore = constraintData.nmrConstraintStore
            else:
              constraintStores = set(x.nmrConstraintStore
                                     for x in nmrCalcRun.findAllData(className='ConstraintStoreData'))
              if len(constraintStores) == 1:
                nmrConstraintStore = constraintStores.pop()
              else:
                raise Exception("No nmrConstraint not identified in nmrCalcRun")

            dd = {}
            for tag in ('individualAtomsIfNoSet', 'individualAtoms',
                        'ignoreChemCompsWithNoSysNames'):
              dd[tag] = keywds[tag]

            resonanceToAtoms = intUtil.getResonanceAtomMap(namingSystemName, nmrConstraintStore, shiftListObj=obj,
                                                           **dd)
            presetInput['resonanceToAtoms'] = resonanceToAtoms
            # NB resonanceToAtoms has BOTH resonances and fixedResonances as keys


          for tag,val in presetInput.items():
            setattr(handler, tag, val)

        handler.writeShifts(filePath,
                            measurementList=obj.measurementList,
                            **keywds)

        print '### Done shift', fileName, fileFormat, 'presetResonanceMapping' in keywds


  # Write data sets
  for dataObj in nmrCalcRun.sortedData():
    print dataObj.name
    if dataObj.name == 'noesyPeakList':

      # write actual peak list
      xx = dataObj.findFirstRunParameter(name='fileFormat')
      if xx is not None:

        fileFormat = xx.textValue
        print fileFormat
        fileName = dataObj.findFirstRunParameter(name='fileName').textValue
        filePath = uniIo.joinPath(targetDir, fileName)
        peakList = dataObj.peakList
        # dataSource = peakList.dataSource
        dataDimRefs = [x.findFirstDataDimRef()
                       for x in peakList.dataSource.sortedDataDims()]

        handler = dataIoHandler.getIoHandler(fileFormat, dataType='peak')
        keywds = intUtil.getFileFormatData(fileFormat, dataType='peak').get(
                  'IOkeywords',{}).copy()
        keywds.update(myKeywds)

        if keywds.get('presetResonanceMapping'):
          for tag,val in presetInput.items():
            setattr(handler, tag, val)

        handler.writePeaks(filePath, peakLists=[peakList],
                           dataDimRefs=dataDimRefs, **keywds)

        print '### Done peak', fileName, fileFormat

    elif dataObj.name == 'constraintLists':

      # write actual constraint list
      xx = nmrCalcRun.findFirstRunParameter(name='restraintFormat')
      print 'restraintFileObject', xx
      if xx is not None:

        fileFormat = xx.textValue
        print 'restraintFileFormat', fileFormat
        fileNames = [x.textValue for x in dataObj.sortedRunParameters()
                     if x.name == 'fileName' and x.ioRole == 'input']
        for ii,constraintList in enumerate(dataObj.constraintLists):
          filePath = uniIo.joinPath(targetDir, fileNames[ii])
          constraintTag = constraintList.className[:-14]
          constraintTag = constraintTag[0].lower() + constraintTag[1:]
          handler = dataIoHandler.getIoHandler(fileFormat, dataType=constraintTag)
          keywds = intUtil.getFileFormatData(fileFormat, dataType=constraintTag).get(
                    'IOkeywords',{}).copy()
          keywds.update(myKeywds)

          if keywds.get('presetResonanceMapping'):
            for tag,val in presetInput.items():
              setattr(handler, tag, val)

          print (' '.join(('### writing', constraintTag, fileName, str(constraintList))))
          handler.writeConstraints(filePath, constraintList=constraintList, **keywds)



          print '### Done constraintList', fileName, fileFormat


def setupCalculation(nmrCalcRun, targetDir, protocolName=None,
                     pluginModule=None, doGeneralAdapt=None):
  """ Set up NMR calculation
  """

  print 'Writing input to:'
  print targetDir

  if protocolName is None:
    protocolName = nmrCalcRun.wmsProtocolName

  if pluginModule is None:
    pluginModule = intUtil.getIntegratorPlugin(protocolName)
  print ('Starting setupCalculation; protocol %s; plugin %s' #
         % (protocolName, pluginModule and pluginModule.__name__))

  # set env variable for passing
  os.environ['CCPN_CURRENT_RUN'] = targetDir

  if doGeneralAdapt is True:
    intUtil.adaptNmrCalcRun(nmrCalcRun)

    pluginModule.Util.adaptNmrCalcRun(nmrCalcRun, protocolName)

  pluginModule.write.write(nmrCalcRun, targetDir)


def setupMultiInteractive(argServer, protocolNames, prelimProtocolName,
                          executeProc=False):
  """ Set up or run multiple, connected protocols in one run
  """

  # Get master data. shared between protocols
  masterRun = initRunInteractive(argServer,
                                       protocolName=prelimProtocolName)

  executeScripts = []
  for protocolName in protocolNames:
    executeScript = setupSingleInteractive(argServer, protocolName=protocolName,
                                         masterRun=masterRun)
    executeScripts.append(executeScript)

  if executeProc:
    for ii,executeScript in enumerate(executeScripts):
      if executeScript:
        pid = subprocess.Popen(['python', executeScript]).pid
        print 'CCPN executing %s process: %s' % (protocolNames[ii],pid)
      else:
        print 'CCPN not executing %s - no script generated' % protocolName


def runSingleInteractive(argServer, protocolName, prelimProtocolName=None,
                      masterRun=None):
  """ Run single protocol for one run
  """
  executeScript = setupSingleDialogue(argServer, protocolName,
                                       prelimProtocolName,
                                       masterRun)
  if executeScript:
    process = subprocess.Popen(['python', executeScript])
    pid = process.pid
    print 'CCPN executing %s process: %s' % (protocolName,pid)
    process.communicate()

    while process.returncode is not None:
       print process.returncode
       continue
    else:
      targetDir = nmrCalcDir(nmrCalcRun)
      project = argServer.getProject()
      convert(project,targetDir)
  else:
     print 'CCPN not executing %s - no script generated' % protocolName


def runCyana2Ccpn(argServer, protocolName, prelimProtocolName=None,
                           masterRun=None):
    nmrCalcRun = initRunInteractive(argServer, protocolName=protocolName,
                                        prelimProtocolName=prelimProtocolName,
                                        masterRun=masterRun)

    executeScript=prepareSingleRun(nmrCalcRun, protocolName)

    if executeScript:
      process = subprocess.call(['python', executeScript])
      #print 'CCPN executing %s process: %s' % (protocolName,pid)
      # process.poll()
      #
      # if process.returncode is not None:
      targetDir = nmrCalcDir(nmrCalcRun)
      importFromCyana(nmrCalcRun,targetDir)
    else:
     print 'CCPN not executing %s - no script generated' % protocolName


def setupCyana2CcpnDialogue(argServer, protocolName, prelimProtocolName=None,
                           masterRun=None):
    nmrCalcRun = setupSingleDialogue(argServer, protocolName=protocolName,
                                        prelimProtocolName=prelimProtocolName,
                                        masterRun=masterRun)

    executeScript = prepareSingleRun(nmrCalcRun, protocolName, doGeneralAdapt=True)

def runCyana2CcpnDialogue(argServer, protocolName, prelimProtocolName=None,
                           masterRun=None):
    nmrCalcRun = setupSingleDialogue(argServer, protocolName=protocolName,
                                        prelimProtocolName=prelimProtocolName,
                                        masterRun=masterRun)

    executeScript = prepareSingleRun(nmrCalcRun, protocolName, doGeneralAdapt=True)

    xx = argServer.askYesNo("Run Calculation")
    if xx:
      if executeScript:
        process = subprocess.call([sys.executable, executeScript])
        #print 'CCPN executing %s process: %s' % (protocolName,pid)
        # process.poll()
        #
        # if process.returncode is not None:
        # yy = argServer.showInfo("Calculation complete")
        # if yy:
        targetDir = nmrCalcDir(nmrCalcRun)
        #   importFromCyana(nmrCalcRun,targetDir)
        calculationData=(nmrCalcRun, targetDir)

      else:
       print 'CCPN not executing %s - no script generated' % protocolName
    else:
      pass
    print "calculationData",calculationData
    return calculationData

def setupPreviousCalculation(argServer, protocolName, prelimProtocolName=None,
                           masterRun=None):

    currentProject=argServer.getProject()
    currentNmrCalcStore = currentProject.findFirstNmrCalcStore()
    masterRun = currentNmrCalcStore.sortedRuns()[-1]

    nmrCalcRun = intUtil.makeDerivedRun(masterRun)

    prepareSingleRun(nmrCalcRun, protocolName, doGeneralAdapt=False)

def runPreviousCalculation(argServer, protocolName, prelimProtocolName=None,
                           masterRun=None):

    currentProject=argServer.getProject()
    currentNmrCalcStore = currentProject.findFirstNmrCalcStore()
    masterRun = currentNmrCalcStore.sortedRuns()[-1]

    nmrCalcRun = intUtil.makeDerivedRun(masterRun)

    executeScript=prepareSingleRun(nmrCalcRun, protocolName, doGeneralAdapt=False)

    if executeScript:
      process = subprocess.call(['python', executeScript])
      #print 'CCPN executing %s process: %s' % (protocolName,pid)
      # process.poll()
      #
      # if process.returncode is not None:
      targetDir = nmrCalcDir(nmrCalcRun)
      importFromCyana(nmrCalcRun,targetDir)
    else:
      print 'CCPN not executing %s - no script generated' % protocolName


def importDataFromCyana(argServer):

    targetFile=argServer.getFile()
    jsonObj=json.load(open(targetFile))
    currentProject=argServer.getProject()
    nmrCalcId=jsonObj['CCPN.nmrCalcId']
    nmrCalcRun = getNmrCalcRunFromId(currentProject, nmrCalcId)
    targetDir = nmrCalcDir(nmrCalcRun)
    dataSources = importFromCyana(nmrCalcRun,targetDir)
    return dataSources


def setupSingleInteractive(argServer, protocolName, prelimProtocolName=None,
                           masterRun=None):
  """ Set up single protocol for one run
  """

  # get data
  nmrCalcRun = initRunInteractive(argServer, protocolName=protocolName,
                                        prelimProtocolName=prelimProtocolName,
                                        masterRun=masterRun)

  return prepareSingleRun(nmrCalcRun, protocolName, doGeneralAdapt=True)


def setupSingleDialogue(argServer, protocolName, prelimProtocolName=None,
                           masterRun=None):
  """ Set up single protocol for one run
  """

  # get data

  nmrCalcRun = initRunDialogue(argServer, protocolName=protocolName,
                                        prelimProtocolName=prelimProtocolName,
                                        masterRun=masterRun)

  return nmrCalcRun

def prepareSingleRun(nmrCalcRun, protocolName, doGeneralAdapt=None):
  """ Set up single protocol for one run
  """

  # Set up calculation and write data files
  pluginModule = intUtil.getIntegratorPlugin(protocolName)

  # Create target directory. NB must not exist before this point,
  targetDir = nmrCalcDir(nmrCalcRun)
  os.makedirs(targetDir)
  outfp = open(uniIo.joinPath(targetDir, setupLogFile),'w')
  print 'Redirecting output to %s' % setupLogFile

  try:
    sys.stdout = outfp

    if doGeneralAdapt is True:
      setupCalculation(nmrCalcRun, targetDir=targetDir, protocolName=protocolName,
                       pluginModule=pluginModule, doGeneralAdapt=True)
    else:
      setupCalculation(nmrCalcRun, targetDir=targetDir, protocolName=protocolName,
                       pluginModule=pluginModule, doGeneralAdapt=False)

  except:

    traceback.print_exc(file=outfp)
    raise

  finally:
    sys.stdout = sys.__stdout__
    outfp.close()
    print 'Done'

  nmrCalcRun.root.saveModified()
  print "CCPN: Set up %s run in %s" % (protocolName,targetDir)

  # Prepare for local execution
  if hasattr(pluginModule.write, 'prepareLocalExecution'):
    commandList, scriptTargetDir, commandList2 = pluginModule.write.prepareLocalExecution(nmrCalcRun,
                                                                            targetDir)

    return writeExecuteScript(commandList, scriptTargetDir, commandList2,
                              protocolName=nmrCalcRun.wmsProtocolName)

  else:
    print ("WARNING, function 'prepareLocalExecution' not found in %s"
           % pluginModule.__path__)


def writeExecuteScript(commandList, targetDir, commandList2, scriptName='ccpnExecute.py', logFile='ccpnCalc.out',
                       protocolName='unspecified'):
  """ Write script to execute commands, and return full script path.
  """

  scriptFile = uniIo.joinPath(targetDir, scriptName)


  fp = open(scriptFile,'w')
  try:
    fp.write("""\"\"\" CCPN automatic execution script
Protocol: %s
\"\"\"
import os, subprocess
from memops.universal import Io as uniIo
if __name__ == '__main__':
  targetDir = "%s"
  logFile = "%s"
  procargs = %s
  procargs2 = %s

  # set env variable for passing
  os.environ['CCPN_CURRENT_RUN'] = targetDir

  outfp = open(uniIo.joinPath(targetDir, logFile),'w')
  print 'Redirecting output to:', logFile
  cyanatableFile = open(uniIo.joinPath(targetDir, 'cyanatable.txt'), 'w')

  # Execute commands
  curDir = os.getcwd()
  try:
    os.chdir(targetDir)

    print 'Calling calculation program:'
    print '  ' + ' '.join(procargs)
    print '  in ', targetDir
    subprocess.call(procargs, cwd=targetDir, stdout=outfp, stderr=outfp)
    subprocess.call(procargs2, cwd=targetDir, stdout=cyanatableFile, stderr=outfp)


  finally:

    outfp.close()
    os.chdir(targetDir)
    print os.getcwd()
    print 'Done calculation'


    """ % (protocolName, targetDir, logFile, commandList, commandList2))


  finally:
    fp.close()

  #
  return scriptFile




def doPrepareStdWmsRun(nmrCalcRun, pluginModule, targetDir=None):
  """ Prepare for Wms run, starting from command line args (sys.args)
  """

  print '### doPrepareStdWmsRun', nmrCalcRun, pluginModule, targetDir

  if targetDir is None:
    targetDir = nmrCalcDir(nmrCalcRun)

  # Create target directory. NB must not exist before this point,
  os.makedirs(targetDir)

  outfname = uniIo.joinPath(targetDir, setupLogFile)
  outfp = open(outfname,'w')
  print 'Redirecting output to %s' % outfname

  try:
    sys.stdout = outfp

    setupCalculation(nmrCalcRun, targetDir=targetDir,
                     pluginModule=pluginModule)

    nmrCalcRun.root.saveModified()

    #pluginModule.write.write(nmrCalcRun,targetDir)

  except:
    # NB Doing this instead of redirecting stderr,
    # as the latter does not seem to work when called from WMS server.
    traceback.print_exc(file=outfp)

  finally:
    sys.stdout = sys.__stdout__
    outfp.close()
    print 'Done'

    return targetDir



def prepareStdWmsRun(pluginName, projectDir, nmrCalcRunId, targetDir=None):
  """ Prepare for Wms run, starting from command line args (sys.args)
  """

  print '### prepareStdWmsRun', pluginName, projectDir, nmrCalcRunId

  nmrCalcRun = getNmrCalcRun(projectDir, nmrCalcRunId, pluginName)

  print '### got nmrCalcRun', nmrCalcRun

  if nmrCalcRun is None:
    print "No NmrCalcRun found. Aborting"

  else:
    pluginModule = __import__(pluginName, globals(), locals(),
                              ['read', 'write', 'Util'])

    targetDir = doPrepareStdWmsRun(nmrCalcRun, pluginModule,
                                   targetDir=targetDir)
  #
  return targetDir


def getDataSummary(nmrCalcRun):
  """Get string with summary of NmrCalcRun data:
  """

  text = []
  text.append('Run: %s %s %s %s' % (nmrCalcRun.wmsProtocolName,
                                    nmrCalcRun.softwareName,
                                    nmrCalcRun.status, nmrCalcRun))
  text.append('    RunParameters :')
  for pp in nmrCalcRun.sortedRunParameters():
    text.append('%s %s %s %s' %
                (pp.name, pp.code, pp.ioRole, intUtil.getParameterValue(pp)))
    if pp.data:
      text[-1] += ' dataSerial %s' % pp.data.serial

  text.append('    Data :')
  for obj in nmrCalcRun.sortedData():
    text.append('%s %s %s %s %s' %
                (obj.name, obj.code, obj.ioRole, obj.serial, obj.className))

  return '\n'.join(text)

def mergeParallelRuns(calcId, projectFiles, targetDir=None):
  """ Merge output of multiple NmrCalcRuns into a single project.
  Assumes input and result data are in projdir/ccpnmr/integrator/
  as will be the case with WMS output
  calcId is nmrCalcId of masterRun
  For each projectFile merges in data from most recent daughter run that has
  a '.out' file, using integrator.plugins.xyz.read
  failing that usesmost recent daughter run
  Project files can be either 'tgz' zipped projects, or project directories
  NB All but one projectFIle must have a match with a '.out' directory.
  """

  calcStoreGuid, masterRunId = calcId.split('+',1)
  masterCalcId = calcId

  # Make temporary directory
  if targetDir is None:
    targetDir = os.getcwd()
  targetDir = os.path.join(targetDir, "CCPN-" + time.strftime(timeFormat))
  os.makedirs(targetDir)
  tmpDir = os.path.join(targetDir,'tmp')
  os.makedirs(tmpDir)

  # copy or extract input project files
  mergeData = [None]*len(projectFiles)
  for ii,path in enumerate(projectFiles):

    # copy or extract input project files
    target = os.path.join(tmpDir, str(ii))
    if os.path.isdir(path):
      shutil.copytree(path, target)
    elif os.path.isfile(path) and (path.endswith('.tgz') or
                                   path.endswith('.tar.gz')):
      tarObj = tarfile.open(path)
      os.makedirs(target)
      tarObj.extractall(path)
    else:
      raise Exception("Project path is neither directory nor .tgz file : " + path)

    # get extracted project directory
    ll = genIo.findCcpnProjectDirs(target)
    if len(ll) == 1:
      projDir = ll[0]
    else:
      raise Exception(" %s should contain 1 CCPN project, %s found"
                      % (target, len(ll)))

    # find active run datadir

    # First make dictionary to sort in runId order
    integratorDir = os.path.join(projDir,'integrator')
    subdirs = [x for x in os.listdir(integratorDir)
               if x.startswith(calcStoreGuid)]
    subdirDict = {}
    for ss in subdirs:
      if ss.endswith('.in'):
        runId = int(ss.split('+')[1][:-3])
        subdirDict[runId] = ss

    # Now get highest runId with an output directory
    backupData = None
    for runId,inDir in sorted(subdirDict.items()):
      propFile = os.path.join(integratorDir, inDir,  propFileName)
      if os.path.isfile(propfile):
        propData = json.load(open(propFile))
        if masterCalcId == propData.get('CCPN.masterCalcId'):
          wmsProtocol = propData.get('CCPN.Run.wmsProtocolName')
          outDir = inDir[:-3]+'.out'
          if outDir in subdirs:
            mergeData[ii] = (runId, wmsProtocol, outDir, projDir)
          else:
            # no '.out' file. Likely an ARIA run that does nto produce them.
            backupData = (runId, wmsProtocol, None, projDir)

    # Nothing found use the backupData instead (likely an ARIA run)
    if mergeData[ii] is None:
      if backupData is None:
        raise Exception(" %s no output project found" % (target,))
      else:
        mergeData[ii] = backupData

  # Get the index of the project to use as base:
  baseIndex = None
  for ii, tt in enumerate(mergeData):
    if tt[2] is None:
      if baseIndex is None:
        baseIndex = ii
      else:
        raise Exception("Two directories lacked mergeable data: %s, %s"
                        % (mergeData[baseIndex][3],mergeData[ii][3]))
  if baseIndex is None:
    baseIndex = 0

  # Get base project
  tt = mergeData[baseIndex]
  project = genIo.loadProject(tt[3])
  nmrCalcStore = project.findFirstNmrCalcStore(guid=calcStoreGuid)
  calcRun = nmrCalcStore.findFirstRun(serial=tt[0])
  result = tt[3]
  target = os.path.join(result, 'integrator')
  if not calcRun.wmsProtocolName:
    # Necessary as ARIA projects lose their protocol names
    # in conversion to/from oler data model version
    calcRun.wmsProtocolName = tt[1]

  for ii, tt in enumerate(mergeData):
    if ii != baseIndex:
      calcRun = nmrCalcStore.findFirstRun(serial=tt[0])
      pluginModule = intUtil.getIntegratorPlugin(tt[1])
      pluginModule.read.read(calcRun, tt[2])
      if not calcRun.wmsProtocolName:
        # Necessary as ARIA projects lose their protocol names
        # in conversion to/from older data model version
        calcRun.wmsProtocolName = tt[1]
      #
      shutil.copytree(tt[2], target)
  #
  project.saveModified()

  print 'Finished Data Merge. Results are in ', result
  return result


