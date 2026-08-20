"""
======================COPYRIGHT/LICENSE START==========================

projectIO.py: I/O for Sparky project files

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

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

# Import general functions
from memops.universal.Util import returnFloat, returnFloats
from memops.universal.Util import returnInt, returnInts
from ccp.format.sparky.generalIO import SparkyGenericFile
from memops.universal.Io import splitPath, joinPath

#####################
# Class definitions #
#####################

from ccp.format.sparky.peaksIO import SparkyPeak
from ccp.format.sparky.peaksIO import SparkyPeakFile

from ccp.format.sparky.sequenceIO import SparkySequenceFile
from ccp.format.sparky.sequenceIO import SparkySequence
from ccp.format.sparky.sequenceIO import SparkySequenceElement

from ccp.format.sparky.chemShiftsIO import SparkyChemShiftFile
from ccp.format.sparky.chemShiftsIO import SparkyChemShift
     
class SparkyProjectFile(SparkyGenericFile):
  
  #
  # Information on file level
  #
  
  def initialize(self):

    self.molNames = []
    self.moleculeConditions = {}
    
    self.chemShiftLocation = {}
    self.chemShiftFiles = []
    
    self.sequenceLocation = {}
    self.sequenceFile = None
    
    self.saveFiles = []
    self.peaksLocation = {}
    self.peaksCondition = {}

    # TODO: split this up in project and save file reference?
    self.refData = [
    
           ['savefiles',[]],
           ['molecule',[['attached data',[]],
                        ['condition',[['resonances',[]]]],
                       ]
           ],
           ['user',[]],
           ['spectrum',[['attached data',[]],
                        ['view',[['params',[]]]],
                        ['ornament',[]]
                       ]
           ]

                   ]


  def read(self,fileType = 'project',saveFileKey = None,readSaveFiles = True,verbose = 0):

    if verbose == 1:
      print "Reading sparky %s file %s" % (fileType,self.name)

    fin = open(self.name, 'rU')

    #
    # Start reading, look for header line
    #

    location = []
    level = -1

    locationSubs = self.setLocationSubs(level,location)

    line = fin.readline()

    if not line.count("<sparky %s file>" % fileType):
      print "  File %s is not a %s file - ignored." % (self.name,fileType)
      return

    line = fin.readline()

    while line:
    
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        line = fin.readline()
        continue

      bracketSearch = self.patt[self.format + 'SharpBracketBetween'].search(line)

      if bracketSearch:

        bracketName = bracketSearch.group(1)
        bracketSplit = bracketName.split()

        if bracketSplit[0] == 'version':
          pass

        elif bracketSplit[0] == 'end':
          bracketName = ' '.join(bracketSplit[1:])
          
          #
          # Resets if no data was found...
          #

          if location == ['molecule','attached data'] and not self.sequenceFile.sequences[0].elements:

            self.sequenceFile = None
            del(self.sequenceLocation[self.molNames[-1]])
            
          #
          # Make sure navigation is correct
          #

          if level == -1:
            pass

          elif bracketName != location[level]:
            print "  Bracketing mismatch for %s in sparky project. Exiting." % bracketName
            return

          else:
            level -= 1
            location.pop(-1)
            locationSubs = self.setLocationSubs(level,location)

        else:

          if bracketName not in locationSubs:

            print "  Unrecognized subbracket name %s in %s!" % (bracketName,str(location))

            if level != -1:
              location.append(bracketName)
              level += 1

          else:
            location.append(bracketName)
            level += 1

            locationSubs = self.setLocationSubs(level,location)

            if location == ['spectrum']:

              peakFile = SparkyPeakFile(self.name)
              self.peaksLocation[saveFileKey] = []
              self.peaksLocation[saveFileKey].append(peakFile)    
              tempMolName = None
              
            elif location == ['molecule','condition']:
            
              self.moleculeConditions[self.molNames[-1]] = []

            elif location == ['molecule','attached data']:
            
              if not self.sequenceFile:
                self.sequenceFile = SparkySequenceFile(self.name)
              
              self.sequenceFile.sequences.append(SparkySequence())
              self.sequenceLocation[self.molNames[-1]] = self.sequenceFile.sequences[-1]

            elif location == ['molecule','condition','resonances']:

              chemShiftFile = SparkyChemShiftFile(self.name)
              self.chemShiftFiles.append(chemShiftFile)
              self.chemShiftLocation[self.molNames[-1]][self.moleculeConditions[self.molNames[-1]][-1]] = chemShiftFile
            
              
      else:

        if location == ['savefiles']:

          self.saveFiles.append(cols[0])

        elif location == ['molecule']:

          if cols[0] == 'name':
            molName = ' '.join(cols[1:])
            self.molNames.append(molName)
            self.chemShiftLocation[self.molNames[-1]] = {}

        elif location == ['molecule','attached data']:
        
          if cols[0] == 'sequence':
            code1LetterSequence = cols[1]
            
            for i in range(0,len(code1LetterSequence)):
              code1Letter = code1LetterSequence[i]
              self.sequenceLocation[self.molNames[-1]].elements.append(SparkySequenceElement(i+1,code1Letter = code1Letter))

          elif cols[0] == 'first_residue_number':
            firstSeqNum = returnInt(cols[1])
            for i in range(0,len(self.sequenceLocation[self.molNames[-1]].elements)):
              seqEl = self.sequenceLocation[self.molNames[-1]].elements[i]
              seqEl.seqCode = firstSeqNum + i

        elif location == ['molecule','condition']:

          if cols[0] == 'name':
            self.moleculeConditions[self.molNames[-1]].append(' '.join(cols[1:]))

        elif location == ['molecule','condition','resonances']:
          
          atomInfo = cols[0].split(self.resonanceSep)
          
          # Sometimes split up badly...
          if len(atomInfo) != 4 and (cols[1][0] == self.resonanceSep or cols[0][-1] == self.resonanceSep):
            atomInfo = (cols[0]+cols[1]).split(self.resonanceSep)
            if len(atomInfo) != 4:
              atomInfo = []
          
          if len(atomInfo) != 4:
            print "  Error: bad '%s' atom information for chemical shift. Ignored." % cols[0]
          else:
            (void1,residue,atom,void2) = atomInfo
  
            values = {}
            values['residue'] = [residue]
            values['atom'] = [atom]
            values['chemShift'] = [cols[1]]
            values['atomType'] = [cols[2]]
  
            curShift = SparkyChemShift(chemShiftFile,values)
            chemShiftFile.chemShifts.append(curShift)

        elif location == ['spectrum']:

          if cols[0] == 'name':
            peakFile.specNames = [cols[1]]  # TODO: or can have space?
          elif cols[0] == 'pathname':
            peakFile.specFiles = [cols[1]]
          elif cols[0] == 'dimension':
            peakFile.numDims = [returnInt(cols[1])]
          elif cols[0] == 'points':
            peakFile.specNumPoints = [returnInts(cols[1:])]
          elif cols[0] == 'molecule':
            tempMolName = ' '.join(cols[1:])
            if tempMolName not in self.molNames:
              self.molNames.append(tempMolName)
          elif cols[0] == 'condition':
            tempCondName = ' '.join(cols[1:])

            if tempMolName == None:
              tempMolName = ''
              self.molNames.append(tempMolName)
            if tempCondName not in self.moleculeConditions[tempMolName]:
              self.moleculeConditions[tempMolName].append(tempMolName)
            self.peaksCondition[peakFile] = [tempMolName,tempCondName]

        elif location == ['spectrum','ornament']:

          if cols[0] == 'type':
            curType = cols[1]

            if curType == 'peak':
              curPeak = SparkyPeak(peakFile,{})
              peakFile.peaks.append(curPeak)

          elif curType == 'peak':
            if cols[0] == 'id':
              curPeak.setValue('num',[cols[1]])
            elif cols[0] == 'pos':
              curPeak.setValue('freqDimsPpm',cols[1:])
            elif cols[0] == 'height':
              curPeak.setValue('intensity',[cols[1]])
              curPeak.setValue('intensityFit',[cols[2]])
            elif cols[0] == 'linewidth':
              curPeak.setValue('freqDimsWidth',cols[1:])
            elif cols[0] == 'integral':
              curPeak.setValue('volume',[' '.join(cols[1:3])])
            elif cols[0] == 'rs':
              assignment = '-'.join(cols[1:])
              assignment = assignment.replace(self.resonanceSep,'')
              curPeak.setValue('assignment',[assignment])
            elif cols[0] == 'fr':
              pass
              # TODO: what does this mean?

      line = fin.readline()

    if fileType == 'project' and readSaveFiles:
    
      (projectDir,fileName) = splitPath(self.name)

      for saveFile in self.saveFiles:
        self.name = joinPath(projectDir,saveFile)
        self.read(fileType = 'save', saveFileKey = saveFile, verbose = verbose)

  def write(self, verbose = 0):
  
    #
    # HOW ARE SEQUENCE FILES SAVED?!?! With molecule name?
    #
    
    if verbose == 1:
      print "Writing Sparky project file %s." % self.name

    fout = open(self.name,'w')

    fout.write('<sparky project file>' + self.newline)
    fout.write('<version 3.104>' + self.newline)
    fout.write('<savefiles>' + self.newline)
    
    for saveFile in self.saveFiles:
    
      fout.write(saveFile + self.newline)
      
    fout.write('<end savefiles>' + self.newline)
    
    #
    # Make sure that there's a chemShiftFile for each molName/condition combination!
    #
        
    for molName in self.molNames:
    
      fout.write('<molecule>' + self.newline)
      fout.write('name %s' % molName + self.newline)

      if self.sequenceLocation.has_key(molName):
      
        if self.sequenceLocation[molName].elements:
        
          fout.write('<attached data>' + self.newline)

          code1LetterSequence = ''        
          for seqEl in self.sequenceLocation[molName].elements:
            code1LetterSequence += seqEl.code1Letter
          fout.write('sequence %s' % code1LetterSequence + self.newline)

          fout.write('first_residue_number %d' % self.sequenceLocation[molName].elements[0].seqCode + self.newline)
          fout.write('<end attached data>' + self.newline)
      
      for conditionName in self.moleculeConditions[molName]:
      
        fout.write('<condition>' + self.newline)
        fout.write('name %s' % conditionName + self.newline)
        
        #
        # Check if chemical shift info available
        #
      
        if self.chemShiftLocation.has_key(molName):
      
          if self.chemShiftLocation[molName].has_key(conditionName):
            
            chemShiftFile = self.chemShiftLocation[molName][conditionName]
            
            fout.write('<resonances>' + self.newline)
            
            for chemShift in chemShiftFile.chemShifts:
            
              #
              # TODO: problem here: for case |R23,24|HA| will create both R23 and R24 HA
              # resonance: are then written out separately. This is NOT correct but how
              # else can this be handled on a chemical shift level?!?!
              #
            
              resonanceCode = self.getExportAssignCode(chemShift.resLabel,chemShift.seqCode,chemShift.atomName)
              
              fout.write("%-12s %-13.3f %s" % (resonanceCode,chemShift.value,chemShift.atomType))
              fout.write(self.newline)
              
            fout.write('<end resonances>' + self.newline)
          
        fout.write('<end condition>' + self.newline)
    
      fout.write('<end molecule>' + self.newline)
    
    fout.close()
    
    #
    # Then write saveFiles
    #
    
    for saveFile in self.saveFiles:
       
      if verbose == 1:
        print "Writing associated save file %s." % saveFile

      fout = open(saveFile,'w')

      fout.write('<sparky save file>' + self.newline)
      fout.write('<version 3.23>' + self.newline)
      
      for i in range(0,len(self.peaksLocation[saveFile])):
      
        peakFile = self.peaksLocation[saveFile][i]
      
        fout.write('<spectrum>' + self.newline)
        
        if self.peaksCondition.has_key(peakFile):
          fout.write('molecule %s' % self.peaksCondition[peakFile][0] + self.newline)
          fout.write('condition %s' % self.peaksCondition[peakFile][1] + self.newline)
        
        fout.write('name %s' % peakFile.specNames[0] + self.newline)
        fout.write('pathname %s' % peakFile.specFiles[0] + self.newline)
        fout.write('dimension %d' % peakFile.numDims[0] + self.newline)
        
        fout.write('points')
        for dim in range(0,peakFile.numDims[0]):
          fout.write(' %d' % peakFile.specNumPoints[0][dim])
        fout.write(self.newline)
          
        fout.write('assignFormat %a1-%a2' + self.newline)
        
        fout.write('<ornament>' + self.newline)
        
        for peak in peakFile.peaks:
        
          fout.write('type peak' + self.newline)
          fout.write('color 0 2 0 white' + self.newline)
          fout.write('flags 0 0 0 0' + self.newline)
          fout.write('id %d' % peak.num + self.newline)

          fout.write('pos')
          for ppm in peak.ppm:
            fout.write(' %.3f' % ppm)
          fout.write(self.newline)
          
          if peak.intensity:
            fout.write('height %.3f %.3f' % (peak.intensity,peak.intensityFit) + self.newline)
          
          if peak.widths:

            # TODO: here is probably *not* always fitted!!
            fout.write('linewidth')
            for lw in peak.widths:
              fout.write(' %.3f' % lw)
            fout.write(' fit')
            fout.write(self.newline)
          
          volumeMethod = peak.getVolumeMethodCode()          
          if volumeMethod and peak.volume != None:         
            fout.write('integral %e %s' % (peak.volume,volumeMethod) + self.newline)
          
            
          assignText = "rs"
            
          for assign in peak.assign:
            assignText += ' %s' % self.getAssignCode(assign)
            
          fout.write(assignText)
              
          fout.write(self.newline)

        fout.write('<end ornament>' + self.newline)
        fout.write('<end spectrum>' + self.newline)

    #
    # TODO TODO: sequence 
    #          
      
  def setLocationSubs(self,level,location):

    if len(location) - 1 != level :
      print "  Error: level/location mismatch"
      return None

    if level == -1:
      locationSubs = []
      for (bracketName,bracketSubs) in self.refData:
        locationSubs.append(bracketName)
      return locationSubs 

    else:
      return self.findLocationSubs(location,self.refData)

  def findLocationSubs(self,curLocation,curRefData):

    curLocationName = curLocation[0]
    curBracketSubs = None

    for (bracketName,bracketSubs) in curRefData:
      if curLocationName == bracketName:
        curBracketSubs = bracketSubs
        break

    if curBracketSubs == None:
      print "  Error: could not find subbrackets for %s" % curLocationName
      return

    else:
      if len(curLocation) == 1:
        locationSubs = []
        for (bracketName,bracketSubs) in bracketSubs:
          locationSubs.append(bracketName)
        return locationSubs 
      else:
        return self.findLocationSubs(curLocation[1:],curBracketSubs)
