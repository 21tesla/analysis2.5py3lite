"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information for XPLOR/CNS files

Copyright (C) 2005-2009 Wim Vranken (European Bioinformatics Institute)

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

import os

# Import general functions
from memops.universal.Util import returnInt

from ccp.format.general.formatIO import FormatFile

from ccp.format.general.Constants import defaultMolCode
from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class CnsGenericFile(FormatFile):

  def setGeneric(self):
    
    self.format = 'cns'
    self.defaultMolCode = defaultMolCode

  def interpreteAssignment(self,assiLine,constraint,numMembers):

    #
    # Start analyzing restraints - work with brackets!
    #

    innerElementPatt = self.patt[self.format + 'InnerElementPatt']
    
    innerElements = []
    outerElements = []
    
    depth = 0
    
    for search in innerElementPatt.finditer(assiLine):
    
      if not innerElements:
        lastOuterElement = assiLine[:search.start()]
      else:
        lastOuterElement = assiLine[prevSearch.end():search.start()]
      
      outerElements.append(lastOuterElement)
      
      outerElementLen = len(lastOuterElement)
      for oei in range(outerElementLen):
        if lastOuterElement[oei] == '(':
          depth += 1
        elif lastOuterElement[oei] == ')':
          depth -= 1
        
        if oei < outerElementLen-1 and lastOuterElement[oei:oei+2].lower() == 'or':
          innerElements.append(("OR",depth))    
      
      # If an OR in the inner element, split this up as well!
      innerElementString = search.group(1).upper()
      depth += 1
      if innerElementString.count('OR'):
        innerElementList = []
        for innerElement in innerElementString.split('OR'):
          innerElementList.append(innerElement)
          innerElementList.append("OR")
        innerElementList.pop(-1)
      else:
        innerElementList = [innerElementString]
        
      for innerElement in innerElementList:     
        innerElements.append((innerElement,depth))
      depth -= 1
      
      prevSearch = search
    
    outerElements.append(assiLine[search.end():])
    
    #print innerElements
    elementList = []
    itemIndex = 0
    
    #print innerElements
    
    for i in range(len(innerElements)):
      innerElement = innerElements[i]
      
      if i == 0:
        elementList.append([])
        for i in range(numMembers):
          elementList[-1].append([])

        #print i, innerElement, itemIndex
        elementList[0][0].append(innerElement[0])
        continue
      
      if innerElement[0] != 'OR':
        if innerElements[i-1][0] == 'OR':   
          if innerElements[i-1][1] == 0:  # TODO this might not be fully correct. Do OR's need to be on ASSI level or not?
             # Back on highest level, need new element
            elementList.append([])
            for i in range(numMembers):
              elementList[-1].append([])
            itemIndex = 0
          else:
            # If within this itemIndex, just add there
            pass
                  
        else:
          itemIndex += 1
      
        #print i, innerElement, itemIndex
        elementList[-1][itemIndex].append(innerElement[0])
    
    #print
    #print elementList
    #print
    
    #print assiLine
    
    for element in elementList:
      
      # Count how many single-member elements.
      numElements = 1      
      for i in range(numMembers):
        numElements *= len(element[i])
      
      # Set up singleElements.
      singleElements = []
      for j in range(len(element[0])):
        divider = numElements/len(element[0])
        for k in range(divider):
          singleElements.append([element[0][j]])
      
      # Now add the rest, scattered order.
      for i in range(1,numMembers):
        curNumMembers = len(element[i])
        divider = numElements/curNumMembers
        index = 0
        for k in range(divider):
          for j in range(curNumMembers):
            singleElements[index].append(element[i][j])
            index += 1
        
      for singleElement in singleElements:      
        self.setConstraintItem(constraint,singleElement)
        #print singleElement

    #print
           
    """

    while (openBracket):
    
      print "OPEN", elementNum, elementInnerNum, assiLine
      print constraintElements
      print

      #
      # Inside an element
      #

      assiLine = assiLine[openBracket.end():]

      #
      # If there's another opening bracket, go in...
      #

      openInnerBracket = self.patt['bracketOpen'].match(assiLine)

      while (openInnerBracket):
      
        print 'OPEN INNER', elementNum, elementInnerNum, assiLine
        print

        assiLine = assiLine[openInnerBracket.end():]

        closeBracket = self.patt['bracketClose'].search(assiLine)
        element = assiLine[:closeBracket.start()]

        #
        # If this is the first inner element, set it up. If not, add to existing one
        #

        if elementInnerNum == 0:
          constraintElements[elementNum].append([element])

        else:
          constraintElements[elementNum][-1].append(element)

        assiLine = assiLine[closeBracket.end():]

        #
        # Is there an OR? TODO: Get OR/XOR whatever in between opening next bracket...
        #

        innerOr = self.patt[self.format + 'RestrInnerOr'].match(assiLine)

        if innerOr:
          assiLine = assiLine[innerOr.end():]
          elementInnerNum += 1

        #
        # Check if there's still brackets (shouldn't be if no innerOr really...)
        #

        openInnerBracket = self.patt['bracketOpen'].match(assiLine)

      #
      # No inner elements, just handle this element...
      #

      closeBracket = self.patt['bracketClose'].search(assiLine)

      element = assiLine[:closeBracket.start()]

      #
      # Set up element, reset elementInnerNum
      #

      if not self.patt['emptyline'].search(element):

        constraintElements[elementNum].append([element])

      elementNum += 1
      elementInnerNum = 0

      assiLine = assiLine[closeBracket.end():]

      openBracket = self.patt['bracketOpen'].search(assiLine)

      #
      # If there is an opening bracket, check if there's an OR
      #

      if openBracket:

        #
        # Is there now an OR? TODO: Should OR/XOR whatever in between opening next bracket really...
        #

        outerOr = self.patt[self.format + 'AssignOr'].search(assiLine[:openBracket.start()])

        if outerOr:
          origRestrNumString = self.patt[self.format + 'RestrNum'].search(assiLine[:openBracket.start()])

          if origRestrNumString:
            if constraint.origId != returnInt(origRestrNumString.group(1)):
              print "  Problem parsing file %s: can't handle OR statements" % self.name

          assiLine = assiLine[outerOr.end():]
          elementNum = 0
          openBracket = self.patt['bracketOpen'].search(assiLine)
    
    #
    # Loop over all elements of this constraint and set it up...
    #
 
    for i in range(0,len(constraintElements[0])):

      elementList = []

      for el in range(numMembers):
        elementList.append(constraintElements[el][i])


      self.setConstraintItem(constraint,elementList,[])
    """

  def setConstraintItem(self,constraint,singleElement):
    
    constraint.items.append(CnsConstraintItem())

    for member in singleElement:

      #
      # Set up new constraint item for each combination within element as well!
      #

      (chainCode,seqCode,atomName) = self.interpreteElement(member)
      constraint.items[-1].members.append(CnsConstraintMember(chainCode,seqCode,atomName))

  def interpreteElement(self,element):

    constraintItemSearch = self.patt[self.format + 'SeqCode'].search(element)

    if constraintItemSearch:

      seqCode = constraintItemSearch.group(1)
      
      # Look for atom name - if none found is interresidue constraint
      constraintItemAtomSearch = self.patt[self.format + 'AtomName'].search(element)      
      if constraintItemAtomSearch:     
        atomName = constraintItemAtomSearch.group(1)
      else:
        atomName = None

      chainCodeString = self.patt[self.format + 'ChainCode'].search(element)

      if chainCodeString:
        chainCode = chainCodeString.group(1)
        chainCode = chainCode.strip('"')
      else:
        chainCode = self.defaultMolCode

      return (chainCode,seqCode,atomName)

    else:

      print " Problems parsing constraint element ('%s')" % element

      return (None,None,None)


  def readGeneric(self,CnsConstraintClass):

    elements = self.constraintElements

    fin = open(self.name, 'rU')

    restrNum = -1  # Or start at 0?

    assiLine = ''

    lineErrors = []

    lines = fin.readlines()
    numLines = len(lines)
    
    isComment = False

    for line in lines:
    
      if self.patt['exclamation'].search(line) or self.patt['emptyline'].search(line) or self.patt[self.format + 'Class'].search(line):
        continue
      
      # NBNB Rasmus 27/6/13
      if not assiLine and self.patt['hash'].search(line):
        # Skip hash comments at the start of the file only. 
        # It is against the format rules, but harmless and not unlikely in practice.
        continue
      
      # Cut line at exclamation mark - are comments (TODO track these!)  
      if line.count("!"):
        line = line[:line.index("!")]
        
      # Get rid of {} comments.

      origLine = line 
      commentStart = self.patt[self.format + 'LongCommentStart'].search(line)
      commentStartIndex = None

      if commentStart:
        isComment = True
        commentStartIndex = commentStart.start()
        
        # Check if not peak assignment info
        if not self.patt['emptyline'].search(line[:commentStartIndex]):
          if self.patt[self.format + 'RestrNum'].search(line):
            isComment = False
          else:
            line = line[:commentStartIndex]
        else:
          line = ""
              
      if isComment:
      
        commentEnd = self.patt[self.format + 'LongCommentEnd'].search(origLine)
 
        if commentEnd:
          isComment = False
          commentEndIndex = commentEnd.end()

          # Check if not peak assignment info
          if not self.patt['emptyline'].search(origLine[commentEndIndex:]):
            if not self.patt[self.format + 'RestrNum'].search(origLine):
              if commentStartIndex != None:
                line = origLine[:commentStartIndex] + origLine[commentEndIndex:]
              else:
                line = origLine[commentEndIndex:]
          else:
            line = ""
       
          #if commentStartIndex != None:
          #  print "ONE LINE",  origLine[commentStartIndex:commentEndIndex]
          #else:
          #  print "END", origLine[:commentEndIndex]
            
        # If it is a comment and origLine not set means it's a standalone line
        #                    and line not set means it's been handled by start or end comment (nothing left)
        elif commentStartIndex == None:
          #print "MIDDLE",line
          continue       
        elif not line:
          #print "START", origLine[commentStartIndex:]
          continue
      
      # Doublecheck if line now empty
      if self.patt['emptyline'].search(line):
        continue

      checkLine = self.checkLinePattern(line)

      if not checkLine:
        return None

      elif checkLine == -1:
        pass

      elif self.patt[self.format + 'Assign'].search(line):

        if assiLine != '':
          if self.constraints[-1].getSpecificInfo(assiLine):
            self.interpreteAssignment(assiLine,self.constraints[-1],elements)
          else:
            lineErrors.append("  Error in cns line (not read):\n%s" % assiLine)
            del self.constraints[-1]

          assiLine = ''

        restrNum += 1

        origRestrNumString = self.patt[self.format + 'RestrNum'].search(line)

        if origRestrNumString:
          origRestrNum = origRestrNumString.group(1)
        else:
          origRestrNum = None

        #
        # Set up new constraint
        #

        self.constraints.append(CnsConstraintClass(restrNum,origRestrNum,self.patt,self.format))

        assiLine += line

      else:
        assiLine += line

    #
    # Handle final line
    #
    
    if self.constraints[-1].getSpecificInfo(assiLine):
      self.interpreteAssignment(assiLine,self.constraints[-1],elements)
    else:
      lineErrors.append("  Error in cns line (not read):\n%s" % assiLine)
      del self.constraints[-1]

    fin.close()
    
    #
    # Check if crap was read
    #
    
    fileReadOk = True
    
    if len(lineErrors) > numLines * 0.1:
      print "  Error: CNS file not valid, too many errors"
      fileReadOk = False
    else:
      for lineError in lineErrors:
        print lineError
    
    return fileReadOk
    
class CnsConstraintItem:
  def __init__(self):
    
    self.members = []
    
class CnsConstraintMember:
  def __init__(self,chainCode,seqCode,atomName):
    
    self.chainCode = chainCode
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.atomName = atomName
