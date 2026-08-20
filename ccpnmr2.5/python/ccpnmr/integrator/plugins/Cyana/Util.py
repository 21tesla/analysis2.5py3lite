
"""
======================COPYRIGHT/LICENSE START==========================

Util.py: code for CCPN data model and code generation framework

Copyright (C) 2012  (CCPN Project)

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
import random

from ccpnmr.integrator.core import Io as intIo
from ccpnmr.integrator.core import Util as intUtil


defaultConfiguration = {
  "protocol": {
    "protocolParameters": [
      {"name":"fileNameSequence","paramType":"String","value":"protein.seq","code":"fileNameSequence"},
      {"name":"fileFormatSequence","paramType":"String","value":"CYANA"},
      {"name":"peakFormat","paramType":"String","value":"XEASY"},
      {"name":"shiftFormat","paramType":"String","value":"CYANA"},
      {"name":"restraintFormat","paramType":"String","value":"CYANA"},
      {"name":"atomNamingSystem","paramType":"String","value":"CYANA2.1"},
      {"name":"useXeasyDimCodes","paramType":"Boolean","value":True}
    ]
  }

}

def adaptNmrCalcRun(nmrCalcRun, protocolName):
  """ Modify nmrCalcRun from generic MultiStructure form to protocol-specific
    Input:
      nmrCalcRun: NmrCalc.Run
      String protocolName
  """

  if protocolName.startswith('CYANA'):
    cyanaNmrCalcRun(nmrCalcRun)

  else:
    raise Exception("Protocol %s not recognized" % protocolName)

def cyanaNmrCalcRun(nmrCalcRun):
  """ Modify nmrCalcRun from generic MultiStructure form to CANDID-specific
    Input:
      nmrCalcRun: NmrCalc.Run
  """

  print ('### Cyana adaptNmrCalcRun')

  ###Reset calcMode values
  valueMap = {
  'none':'default',
  'assignedPeaks':'pks',
  'unassignedPeaks':'upks',
  }
  runParam = nmrCalcRun.findFirstRunParameter(name='calcMode')
  val = valueMap.get(runParam.textValue)
  print val
  if val is None:
    raise ValueError("Cyana calculation mode %s must be one of %s"
                     % (runParam.textValue, tuple(valueMap.keys())))
  else:
    runParam.textValue = val


  # Sequence. Set RMSDresidues
  seqDataObj = nmrCalcRun.findFirstData(name='definedResidues')
  #seqIds = seqDataObj.residueSeqIds # will not work - seqIds is empty when takingwhole chain
  seqIds = [x.seqId for x in seqDataObj.sortedResidues()]
  ss = intUtil.integerListExpression(seqIds,fieldSep='..')
  nmrCalcRun.newRunParameter(name='definedResidueString', code='rmsdrange',
                             textValue=ss)

  # first and last residue numbers, for one-chain web version
  nmrCalcRun.newRunParameter(name='firstSeqId', code='firstSeqId',
                             intValue= min(seqIds))
  nmrCalcRun.newRunParameter(name='lastSeqId', code='lastSeqId',
                             intValue= max(seqIds))

  # Spectrum tolerance
  defaultTolPoints = 1.0
  runParam = nmrCalcRun.findFirstRunParameter(name='defaultTolPoints')
  if runParam:
    val = runParam.floatValue
    if val:
      defaultTolPoints = val

  minTol = 0.02
  runParam = nmrCalcRun.findFirstRunParameter(name='minTolPpm')
  if runParam:
    val = runParam.floatValue
    if val:
      minTol = val

  # random number seed
  nmrCalcRun.newRunParameter(name='ranseed', code='seed',
                             intValue= random.randint(1,10000000))

  #set 'peaks' and 'tolerance' parameters

  # First get data
  peakFiles = []
  dataSources = []
  peakListObjs = [x for x in nmrCalcRun.sortedData()
                  if x.className == 'PeakListData']


  if peakListObjs:

    for peakListData in peakListObjs:
      runParameter = peakListData.findFirstRunParameter(name='fileName')
      peakFiles.append(runParameter.textValue)

      dataSources.append(peakListData.peakList.dataSource)

    # set 'peaks' parameter
    nmrCalcRun.newRunParameter(name='peakFileNames', code='peaks',
                               textValue=','.join(peakFiles))

  customTolerances = nmrCalcRun.findFirstRunParameter(name='assignmentTolerances')
  if customTolerances is None:
    tolValues = intUtil.getAmalgamatedTolerances(dataSources, defaultTolPoints,
                                                   minTol,
                                                   tagOrder=('hx2', 'hx1',
                                                             'x1', 'x2'))
    nmrCalcRun.newRunParameter(name='toleranceString', code='tolerance',
                               textValue=','.join(str(x) for x in tolValues))
  else:
    tolValues = str(customTolerances.textValue)
    print('tolvalues',tolValues)

    nmrCalcRun.newRunParameter(name='toleranceString', code='tolerance',
                               textValue=tolValues)

  # Check number of shift lists
  shiftListObjs = [x for x in nmrCalcRun.sortedData()
                   if x.className == 'MeasurementListData'
                   and x.measurementList.className == 'ShiftList']
  if len(shiftListObjs) == 1:

    shiftListObj = shiftListObjs[0]
    fileNameObj = shiftListObj.findFirstRunParameter(name='fileName')
    nmrCalcRun.newRunParameter(name='shiftFileName', code='prot',
                                textValue=fileNameObj.textValue,
                                data=shiftListObj)
  else:
    # raise Exception("Required 1 shiftList, %s found" % len(shiftListObjs))
    pass
