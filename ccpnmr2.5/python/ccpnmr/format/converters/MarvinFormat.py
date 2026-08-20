"""
======================COPYRIGHT/LICENSE START==========================

MarvinFormat.py: Contains functions specific to Marvin conversions.

Copyright (C) 2007 Wim Vranken    (European Bioinformatics Institute)
Copyright (C) 2007 Gary Thompson  (Astbury Centre, University of Leeds)


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

- contact Gary Thompson (garyt@bmb.leeds.ac.uk)
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

import sys

import string, copy, textwrap

from ccp.general.Util import findChemAtomSysName

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords
from ccpnmr.format.general.Util import getNameInfo

from ccpnmr.format.general.Constants import peakNum_kw
from ccpnmr.format.general.Constants import volume_kw, height_kw
from ccpnmr.format.general.Constants import tagSep

from ccpnmr.format.general.Util import getResName, getResNameText

#
# Additional IOkeywords definitions
#

IOkeywords = copy.deepcopy(IOkeywords)
# note bug Wim had default and mandatory wrong way round? in probability cutoff

IOkeywords['readPeaks']['probabilityCutoff']          = (0.0,  False, 'Can set the probability cutoff to import peaks/assignments.')
IOkeywords['readPeaks']['preservePeaks']              = (True, False, 'Preserve peaks in the original list that don\'t appear in the imported list.')
IOkeywords['readPeaks']['preserveAssignments']        = (True, False, 'Preserve assignments already in the peak list.')
IOkeywords['readPeaks']['removeDuplicateAssignments'] = (True, False, 'Remove duplicate assignments present in peaks (unimplemented)')
IOkeywords['readPeaks']['defaultChainCode']           = (None, False, 'Default CCPN chain code')
IOkeywords['readPeaks']['defaultMolSystem']           = (None, False, 'Default CCPN  molecular system name')
IOkeywords['readPeaks']['createNew']                  = (False,False, 'Create a new peak if a matching peak can\'t be found')
IOkeywords['readPeaks']['axisMapping']                = ('',   False, 'Change the axis marvin axes are mapped to, 132 would map axes marvin-> analysis as 1->1,3->2,2->3 for a peak list with 3 axes, default no mapping')
IOkeywords['readPeaks']['useAppData']                 = (False,False, 'Use applicationData keyword (peakNum) on the analysis peak list to match marvin peak ids')
IOkeywords['readPeaks']['appDataKey']                 = (',',  False, 'Application,keyword pair used to define which appdata to lookup to find marvin peak ids; keyword is typically peakNum')

class MarvinFormat(DataFormat):

  def setFormat(self):

    self.format = 'marvin'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):

    self.getPeaks = self.getPeaksGeneric

  #
  # Deviations from generic import stuff
  #


  #
  # Functions different to default functions in DataFormat
  #

  def thisPeakValid(self):

    #print 'assignList', self.rawPeak.assignList
    #print 'probabilities', self.rawPeak.probabilities
    #print 'labels', self.rawPeak.label

    # filter out invalid data with the current assignment probability cutoff
    validList = [prob >= self.probabilityCutoff for prob in self.rawPeak.probabilities]
    #print 'validList', validList

    print 'thisPeakValid, rawPeak.num', self.rawPeak.num
    #print 'thisPeakValid, rawPeak.assignList', self.rawPeak.assignList
    self.rawPeak.assignList  = self.filterByValid(self.rawPeak.assignList, validList)
    self.rawPeak.label  = self.filterByValid(self.rawPeak.label, validList)
    self.rawPeak.probabilities  = self.filterByValid(self.rawPeak.probabilities, validList)
    #print 'thisPeakValid, rawPeak.assignList after', self.rawPeak.assignList
    #print 'assignList', self.rawPeak.assignList
    #print 'probabilities', self.rawPeak.probabilities
    #print 'labels', self.rawPeak.label

    isValid = False

    for probability in self.rawPeak.probabilities:
      if probability > self.probabilityCutoff:
        isValid = True

    return isValid

  def setPeakIntensity(self):

    # No info?

    pass

  def filterByValid(self,elementArray,truthArray):
      if len(truthArray) != len(elementArray):
          print elementArray
          print truthArray
          raise self.FormatConverterError('internal error the length of the truth array and element array dont agree!')

      result  = []
      for elem,truth in zip (elementArray,truthArray):
          if truth:
              result.append(elem)
      return result

  def setPeakExtras(self):

    #
    # Make peakContribs for this one...
    #

    if self.peak==None:
        return


    numContribs = len(self.rawPeak.assignList)

    for i in range(numContribs):
      self.peakContribs.append(self.peak.newPeakContrib(weight = self.rawPeak.probabilities[i]))

  def findFirstChemAtomWildCard(self,chemComp,atomName):
      result = None
      #print atomName
      atomNameLength = len(atomName)
      while atomName[-1] == '*':
          #print 'pre',atomName
          atomName = atomName[:-1]
          #print 'post',atomName
      #print atomName

      for chemAtom in chemComp.chemAtoms:
          if chemAtom.name.startswith(atomName) and len(chemAtom.name) == atomNameLength:
              result=chemAtom
              break

      return result

  def convertLightToHeavy(self,assign):
      (chain,seqId,codedAtomName)=assign.split('.')
      (typeCode,resType,atomName)=codedAtomName.split(':')
      #print ''
      #print 'assign:',assign
      targetChemComp =  self.project.findFirstChemComp(ccpCode = resType)

      #print 'target chemc comp:', targetChemComp
      #for chemAtom in targetChemComp.chemAtoms:
      #    print chemAtom

      resultAtoms=None
      resultAssign=None
      if targetChemComp != None:
          targetAtom = None
          if typeCode == 'HEAVY':
              #print 'looking for heavy for %s %s %s %s' % (chain,resType,seqId,atomName)

              targetAtom = self.findFirstChemAtomWildCard(targetChemComp, atomName)
              #targetChemComp.findFirstChemAtom(name=atomName)
              #print '\tfound target atom: ', targetAtom
              #print 'targetAtom:', targetAtom,targetAtom.chemBonds,targetAtom.chemBonds[0].chemAtoms

              targetBondAtoms = None
              if targetAtom != None:
                  targetBondAtoms = list(targetAtom.findFirstChemBond().chemAtoms)
                  #print 'bond atoms:', targetBondAtoms,targetAtom
                  if targetAtom in targetBondAtoms:
                      targetBondAtoms.remove(targetAtom)
              #print 'heavy bonds atom:', targetBondAtoms,targetBondAtoms.__class__.__name__
                  resultAssign = '.'.join((chain,seqId,targetBondAtoms[0].name),)
                  #print '\tNote: converted %s to %s'  % (assign,resultAssign)

          else:
              resultAssign =  '.'.join((chain,seqId,atomName),)

      #print assign, resultAssign
      return resultAssign
      #print targetChemComp

  def convertAtomNameByType(self,assign):
      (chain,seqId,codedAtomName)=assign.split('.')
      (typeCode,resType,atomName)=codedAtomName.split(':')
      
      newAssign = None
      
      # Warning: this is slightly dangerous because molType not defined - could in principle find wrong chemComp (Wim 080207)
      targetChemComp =  self.project.findFirstChemComp(ccpCode = resType)
      
      if targetChemComp:

        namingSystem = targetChemComp.findFirstNamingSystem(name = self.namingSystemName)
        
        if namingSystem:

          searchDict = {'sysName': atomName}
          chemAtomSysName = findChemAtomSysName(namingSystem,searchDict)


          newAtomName=chemAtomSysName.atomName
          #print chemAtomSysName,atomName,newAtomName

          newCodedAtomName=':'.join((typeCode,resType,newAtomName),)
          newAssign =  '.'.join((chain,seqId,newCodedAtomName),)

      #if assign != newAssign:
          #print '\tNote: changed %s %s %s %s to %s' % (chain,resType,seqId,atomName,newAtomName)
      return newAssign

  def getPeakResNames(self):

    self.resNames = []

    #print
    #print self.rawPeak.assignList
    for assign in self.rawPeak.assignList:

      # note is there a way to abort here and keep going...
      renamedAssign = self.convertAtomNameByType(assign[self.rawPeakDimIndex])
      if renamedAssign ==  None:
          msg = "Error: couldn't convert assignment %s (%s) to iupac nomenclature!\n stopping..."
          msg = msg % (`self.rawPeakDimIndex`,`assign[self.rawPeakDimIndex]`)
          raise self.FormatConverterError(msg)

      #print 'renamed assign: ', renamedAssign
      heavySwappedAssign = self.convertLightToHeavy(renamedAssign)
      if heavySwappedAssign == None:
          msg = "Error: couldn't convert assignment %s (%s) to it's  heavy parent!\n stopping..."
          msg = msg % (`self.rawPeakDimIndex`,`assign[self.rawPeakDimIndex]`)
          raise self.FormatConverterError(msg)

      self.resNames.append(heavySwappedAssign)



  def setPeakDim(self):

    if self.peak == None:
        return

    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef

    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]

  def setFitMethodSerial(self):

    pass

  def getPresetChainMapping(self,chainList):

    return self.getSingleChainFormatPresetChainMapping(chainList)


  def setPeak(self):

    """

    Creates a peak object inside the data model or uses it to copy information
    onto. Also creates/checks application data with original peak number.

    """
    #print
    #print self.rawPeak.num
    newPeak = 1

    #print self.existingPeaks

    if self.overwrite:
      #print
      #print self.rawPeak.num
      target_peak = None
      if self.useAppData:


          if self.origPeakNumberDict.has_key(self.rawPeak.num):
              target_peak = self.origPeakNumberDict[self.rawPeak.num]


      else:
          target_peak = self.peakList.findFirstPeak(serial=self.rawPeak.num)

      if target_peak!=None:




        self.peak = target_peak
        #print self.existingPeaks.index(target_peak), target_peak.serial
        #print self.existingPeaks.pop(self.existingPeaks.index(target_peak)).serial
        ###print len(self.existingPeaks)
        self.existingPeaks.pop(self.existingPeaks.index(target_peak))
        ###print len(self.existingPeaks)

        #
        # Get rid of intensities and peakcontribs
        # TODO: note that peakDims are kept... could in principle create problems with dataDimRefs
        #

        #for peakIntensity in target_peak.peakIntensities:
        #  peakIntensity.delete()
        if not self.preserveAssignments:
            if len(target_peak.sortedPeakContribs()) > 0:
                #peakName = "%s:%d[%d]" % (target_peak.peakList.dataSource.name,target_peak.peakList.serial, target_peak.serial)

                print 'Warning: replacing assignments for peak %s' % peakName
                for peakContrib in target_peak.sortedPeakContribs():
                    peakContrib.delete()
        ###print self.rawPeak

        newPeak = 0

    if newPeak:

      print 'Warning: ignored  marvin peak assignmnet %s no corresponding analysis peak found' % self.rawPeak.num

      #note I have guarded all the other peak reading routines so they ignore peaks of type None
      # it think this shlould work but I may need some more work elsewhere
      self.peak = None

###      #self.peakList.newPeak( details = "Original number %d" % self.rawPeak.num)
###
###      #
###      # Also set original number in applicationdata so can be reused for writing
###      #
###
###      self.peak.addApplicationData(Implementation.AppDataInt(application = self.format, keyword = peakNum_kw, value = self.rawPeak.num))
  ### can be removed when removeDuplicateAssignments is implimented
  def setPeaks(self):
    if self.removeDuplicateAssignments:
      print 'Warning: removal of duplicate assignments is selected but is not implimented!'

    DataFormat.setPeaks(self)

  def deleteRemainingPeaks(self):

    """

    If overwrite mode is on, delete remaining existing peaks that were
    not in new list.

    """

    #print self.preservePeaks
    if not self.preservePeaks:
        DataFormat.deleteRemainingPeaks(self)

  def findResonances(self,resName):
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName)

      #print '\tsearch for resonance', (chainCode,seqCode,spinSystemId,seqInsertCode,atomName)
      resonance = None

      targetMolSystem =  None
      if len(self.project.molSystems) == 1 and self.defaultMolSystem == None:
          targetMolSystem = self.project.molSystems[0]
      elif self.defaultMolSystem != None:
          targetMolSystem = self.project.findFirstMolSystem(name = self.defaultMolSystem)


      targetChain =  None
      if targetMolSystem != None:
          if len(targetMolSystem.chains) == 1 and self.defaultChainCode == None:
              targetChain = targetMolSystem.chains[0]
          elif self.defaultChainCode != None:
              targetChain = targetMolSystem.findFirstChain(code = self.defaultChainCode)


      targetResidue =  None
      if seqCode != None and targetChain != None:
          targetResidue = targetChain.findFirstResidue(seqCode = seqCode)

      targetAtoms = []
      targetAtomSets = []

      #print atomName
      atomNameLength = len(atomName)
      while atomName[-1] == '*':
          #print 'pre',atomName
          atomName = atomName[:-1]
          #print 'post',atomName
      #print atomName

      if targetResidue != None:
          for atom in targetResidue.atoms:
              if atom.name.startswith(atomName) and len(atom.name) == atomNameLength:
                  targetAtoms.append(atom)
                  targetAtomSets.append(atom.atomSet)


      targetResonances = set()

      for atom in targetAtoms:
          for resonanceSet in atom.atomSet.resonanceSets:
              #for resonance in resonanceSet.resonances:
              targetResonances.update(resonanceSet.resonances)
      #print ''

      if len(targetResonances) > 1:
          selectName = '.'.join((chainCode,`seqCode`,atomName),)
          resonanceIds = [`resonance.serial` for resonance in targetResonances]
          print 'Warning: found more than one resonance for peak (%s) %s' % (selectName,self.rawPeak.num)
          print '         resonances are: ', ','.join(resonanceIds)


      #if targetResonances == None:
      #    print '**** Oh bugger!!'
      #    sys.exit()
      return list(targetResonances)






  def setPeakDimResonance(self):

    """

    Set the connection to (existing) resonances.

    Not assuming that the peak list is for one chain only... but if no chainCode
    is given by the external reader then assuming that it is all 'one chain'.
    Ambiguity then has to be handled by linkResonances, not by this script.

    Note that for XEasy this whole thing is different - have to use atomSerial or
    nothing at all

    """

    # this indicates that the raw peak was not found in the peak list or vice versa
    # could be tackled ata higher level
    if self.peak == None:
        return
    if not hasattr(self,'peakDim'):
        return

    self.getPeakResNames()

    for resNameIndex in range(len(self.resNames)):


      self.resName = self.resNames[resNameIndex]

      #
      # Ignore if None
      #

      if self.resName == None:

        continue

      #
      # Link peakDim to resonance
      #

      if not self.resonanceNames.has_key(self.resName):
        resonances=None
        #
        # find the resonance if we have enough information
        #
        resonances =  self.findResonances(self.resName)
        #
        # Create resonance if doesn't exist.
        #
        if resonances == None:
            resonance = self.createResonance(self.resName)
            resonances = [resonance]

        #
        # Keep track of it...
        #


        self.resonanceNames[self.resName] = resonances

      else:

        #
        # Resonance(s) exists, get info from self.resonanceNames
        #
        resonances = self.resonanceNames[self.resName]


      #
      # Link resonance(s) to peakDim via peakDimContrib
      # Note that will be linked to peakContrib automatically in this case:
      # no distinction is made between the contributions (by default)
      #

      keywds = {}
      #print self.peakContribs
      #print self.rawPeak.num
      if self.peakContribs:
        keywds['peakContribs'] = [self.peakContribs[resNameIndex]]
      #print resonances
      for resonance in resonances:
        peakDimContrib = self.peakDim.findFirstPeakDimContrib(resonance = resonance)

        if not peakDimContrib:
          self.peakDim.newPeakDimContrib(resonance = resonance, **keywds)
        elif keywds:
          peakContribs = peakDimContrib.peakContribs
          for peakContrib in keywds['peakContribs']:
            if not peakContrib in peakContribs:
              peakDimContrib.addPeakContrib(peakContrib)

  APPDATA_APPLICATION =  0
  APPDATA_KEYWORD =  1
  def getAppDataKey(self):

      result = None
      if self.appDataKey != ',':
          appDataApplicationKeyPair = self.appDataKey.split(',')
          if len(appDataApplicationKeyPair) > 2 or len(appDataApplicationKeyPair) ==0:
              msg = """Error: the applicationData application,keyword pair provided in appDataKey\n
                       must be a pair of values: of the form  'application,keyword\n
                       got: '%s'"""
              msg = textwrap.dedent(msg)
              raise self.FormatConverterError(msg % `self.appDataKey`)
          if len(appDataApplicationKeyPair) == 1:
              result=(appDataApplicationKeyPair,'peakNum')
          else:
              result=appDataApplicationKeyPair

      return result


  def setPeaks(self):
  
      #
      # Get naming system name (Wim 080207)
      #
      
      self.namingSystemName = self.getFormatNamingSystemName()
      
      # 
   
      if self.useAppData == True:

          appDataKeyPair = self.getAppDataKey()
          if  appDataKeyPair != None:
              #print 'using app data %s' % appDataKeyPair
              self.origPeakNumberDict = {}
              for peak in self.peakList.sortedPeaks():
                  application = appDataKeyPair[self.APPDATA_APPLICATION]
                  keyword = appDataKeyPair[self.APPDATA_KEYWORD]
                  applData = peak.findFirstApplicationData(application = application, keyword = keyword)
                  #if applData != None:
                  #    print peak.serial,applData.value
                  #else:
                  #    print peak.serial,applData
                  if applData:
                      self.origPeakNumberDict[applData.value] = peak
                  else:

                      peakSerial= peak.serial
                      peakPeakList = peak.peakList.serial
                      peakDataSource=peak.peakList.dataSource.name
                      peakExperiment=peak.peakList.dataSource.experiment.name
                      peakId = '%s:%s[%s].%d' % (peakExperiment,peakDataSource,peakPeakList,peakSerial)
                      msg  = 'Warning: using application data %s but peak %s doesn\'t have integer appdata'
                      print  msg % (self.appDataKey,peakId)
      DataFormat.setPeaks(self)
