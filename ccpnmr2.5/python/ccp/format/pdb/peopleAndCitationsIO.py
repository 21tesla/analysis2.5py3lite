"""
======================COPYRIGHT/LICENSE START==========================

peopleAndCitationsIO.py: I/O for PDB citation and people information

Copyright (C) 2006 Wim Vranken (European Bioinformatics Institute)

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

import os, string

# Import general functions
from memops.universal.Util import returnInt
from memops.universal.Io import getTopDirectory

from ccp.format.general.formatIO import Person, Citation

# Import header reader
from ccp.format.pdb.generalIO import PdbGenericFile, PdbFile

#####################
# Class definitions #
#####################  

class PdbPeopleAndCitationsFile(PdbGenericFile):
  """
  Information on file level - is dummy file for generic PdbFile object really
  """
  def initialize(self):

    self.citations = []
    self.persons = []
    self.authors = []
    
  def read(self, pdbFile = None, getAuthors = True, verbose = 0):
    
    if not pdbFile:
    
      #
      # Read in header if not passed in.
      #

      pdbFile= PdbFile(self.name)
      pdbFile.readHeader()
   
    self.pdbFile = pdbFile
    
    #
    # PDB entry authors
    #

    if getAuthors and self.pdbFile.headerVars.has_key('Authors'):
      authorText = self.pdbFile.headerVars['Authors'].strip()
      authorList = authorText.split(',')
      for author in authorList:
        authorBits = author.split('.')
        
        if len(authorBits) <= 1:
          continue
        elif not authorBits[-1]:
          # Problem with string, occurs quite often, family name first
          firstBits = authorBits[0].split()
          if len(firstBits) > 1:
            authorBits[0] = firstBits[-1]
            authorBits[-1] = ' '.join(firstBits[:-1])
            print "  Warning: fixed author name string %s!" % author
        
        for i in range(len(authorBits)-1):
          authorBits[i] = authorBits[i].strip()
      
        # Possible that family-type name has . in it, apparently.
        while (1):
          if len(authorBits) == 1:
            print "  Warning: setting family name as '%s'" % authorBits[-1]
            break
          if ' ' in authorBits[-2]:
            authorBits[-1] = authorBits.pop(-2) + '.' + authorBits[-1]
            print "  Warning: setting family name as '%s'" % authorBits[-1]
          else:
            break
          
        if authorBits[-1][-4:] in (' JR.',' SR.'):
          authorBits[-1] = authorBits[-1][:-4].strip()
          familyTitle = authorBits[-1][-4:].strip()
        else:
          familyTitle = None

        firstName = initials = None

        if len(authorBits) >= 2:
          firstName = authorBits[0]
          if len(authorBits) > 2:
            initials = authorBits[1:-1]

        self.authors.append(Person(self,firstName,initials,authorBits[-1],familyTitle = familyTitle))
    
    #
    # PDB  entry citation (if any)
    #
    
    journalInfoList = []
    isPrimary = False
    
    if self.pdbFile.headerVars.has_key('Journal'):
            
      journalInfoList.append(self.pdbFile.headerVars['Journal'])
      isPrimary = True

    # Additional reference info from REMARK 1
    if self.pdbFile.otherCitations:
      refIds = self.pdbFile.otherCitations.keys()
      refIds.sort()
      for refId in refIds:
        journalInfoList.append(self.pdbFile.otherCitations[refId])
    
    for journalInfo in journalInfoList:
    
      journalType = 'journal'
      
      if journalInfo.has_key('PUBL'):
        if journalInfo.has_key('EDIT'):
          journalType = 'book'
        else:
          journalType = 'thesis'
        
        # In principle, place of publication is listed first, followed by
        # a space, a colon, another space, and then the name of the publisher/issuer.
        # The alternative way is publisher, comma, then publication place.
        if journalInfo['PUBL'].count(' : '):
          (publisherPlace,publisher) = journalInfo['PUBL'].split(' : ')
        elif journalInfo['PUBL'].count(','):
          publisherInfo = journalInfo['PUBL'].split(',')
          publisher = publisherInfo[0]
          publisherPlace = ', '.join(publisherInfo[1:])
        else:
          publisher = journalInfo['PUBL']
          publisherPlace = None
          
        if publisherPlace:
          publisherPlace = publisherPlace.strip()
      
      else:
        publisherPlace = publisher = None
      
      # Again possible that REF not present...

      if not journalInfo.has_key('REF') or journalInfo['REF'] == 'TO BE PUBLISHED':
        status = 'to be published'
        pubName = None
        volume = None
        issue = None
        firstPage = None
        year = None
      else:
        refInfo = journalInfo['REF']

        pubName = refInfo[0:28].strip()

        year = firstPage = volume = issue = None

        # TODO - not sure how to deal with status in these cases.

        #print 'TYPE: [%s]' % journalType

        if journalType == 'journal' and refInfo[30:32] == 'V.':
          status = 'published'
          year = returnInt(refInfo[43:47].strip())
          firstPage = refInfo[37:42].strip()
          volume = refInfo[32:36].strip()
          #issue is not in PDB file
        elif refInfo[43:47].strip():
          status = 'published'
          year = returnInt(refInfo[43:47].strip())
        else:
          status = 'to be published'

        #year = returnInt(refInfo[-4:].strip())
        #firstPage = refInfo[-10:-4].strip()
        #volume = refInfo[-15:-10].strip()

        #pubName = refInfo[:-18].strip()

      # Necessary to check 'cos of crappy old PDB files...
      if journalInfo.has_key('TITL'):
        title = journalInfo['TITL'].strip()
      else:
        title = None

      citation = Citation(self,journalType,title,status,pubName,volume,issue,firstPage,year,publisherPlace, publisher, country = None, isPrimary = isPrimary)
      self.citations.append(citation)
      
      # Only first citation is primary...
      isPrimary = False
      
      for (attrName,pdbName) in (('authors','AUTH'),('editors','EDIT')):
        if journalInfo.has_key(pdbName):
          authorText = journalInfo[pdbName].strip()
          authorList = authorText.split(',')

          for author in authorList:
            authorBits = author.split('.')

            if len(authorBits) <= 1:
              continue
            elif not authorBits[-1]:
              # Problem with string, occurs quite often, family name first
              firstBits = authorBits[0].split()
              if len(firstBits) > 1:
                authorBits[0] = firstBits[-1]
                authorBits[-1] = ' '.join(firstBits[:-1])
                print "  Warning: fixed author name string %s!" % author
            
            for i in range(len(authorBits)-1):
              authorBits[i] = authorBits[i].strip()
            
            # Possible that family-type name has . in it, apparently.
            while (1):
              if len(authorBits) == 1:
                print "  Warning: setting family name as '%s'" % authorBits[-1]
                break

              if ' ' in authorBits[-2]:
                authorBits[-1] = authorBits.pop(-2) + '.' + authorBits[-1]
                print "  Warning: setting family name as '%s'" % authorBits[-1]
              else:
                break

            if authorBits[-1][-4:] in (' JR.',' SR.'):
              authorBits[-1] = authorBits[-1][:-4].strip()
              familyTitle = authorBits[-1][-4:].strip()
            else:
              familyTitle = None

            personsFound = []
            for tPerson in self.persons:
              if tPerson.familyName == authorBits[-1] and tPerson.firstName == authorBits[0]:
                personsFound.append(tPerson)
            
            person = None
            if personsFound:
              for tPerson in personsFound:
                if tPerson.initials == authorBits[1:-1] and (not familyTitle or tPerson.familyTitle == familyTitle):            
                  person = tPerson
              if not person:
                if len(personsFound) > 1:
                  print "  Warning: picking first person from found list as default (no full initial match)" 
                person = personsFound[0]

            # Create new person if no match...
            if not person:
              firstName = initials = None
              if len(authorBits) >= 2:
                firstName = authorBits[0]
                if len(authorBits) > 2:
                  initials = authorBits[1:-1]

              person = Person(self,firstName,initials,authorBits[-1],familyTitle = familyTitle)
              self.persons.append(person)
              
            getattr(citation,attrName).append(person)
