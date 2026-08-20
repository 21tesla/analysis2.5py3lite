"""
======================COPYRIGHT/LICENSE START==========================

peopleAndCitationsIO.py: I/O for NMR-STAR citation and people information

Copyright (C) 2008 Wim Vranken (European Bioinformatics Institute)

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
#from memops.universal.Util import returnInt

from ccp.format.general.formatIO import Person, Citation, Keyword

from ccp.format.nmrStar.generalIO import NmrStarFile, NmrStarGenericFile

#####################
# Class definitions #
#####################

# Can not be read in separately, only as part of projectIO
   
class NmrStarPeopleAndCitationsFile(NmrStarGenericFile):

  """
  Information on file level.
  """
  
  def initialize(self,parent,saveFrame,infoType):

    self.persons = []
    self.personInGroups = []
    self.citations = []
    self.keywords = []
    
    self.saveFrame = saveFrame
    self.infoType = infoType
    
    self.parent = parent
    self.version = parent.version

    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    if self.version in ('3.0','3.1'): 
      
      personTagNames = ['Given_name','Middle_initials','Family_name','Family_title','Name_salutation']
      
      #
      # Entry level information
      #
                                
      if self.saveFrame.name == 'entry_information':
        
        if self.infoType == 'contactPersons':

          personInGroupTags = ['Email_address', 'Department_and_institution', 'Mailing_address', 'Address_1',
                               'Address_2', 'Address_3', 'City', 'State_province', 'Country', 'Postal_code',
                               'Phone_number', 'FAX_number', 'Role', 'Organization_type']
        
          #
          # Contact person information
          #
  
          tableName = '_Contact_person'
          if self.saveFrame.tables.has_key(tableName):
            
            tableTags = self.saveFrame.tables[tableName].tags

            for i in range(0,len(tableTags['ID'])):

              personInfo = []
              for personTagName in personTagNames:
                personInfo.append(tableTags[personTagName][i])           
              
              self.persons.append(Person(self,serial = tableTags['ID'][i],*personInfo))

              personInGroupInfo = []
              for persInGrpTag in personInGroupTags:
                personInGroupInfo.append(tableTags[persInGrpTag][i])

              self.personInGroups.append(PersonInGroup(self,serial = tableTags['ID'][i],*personInGroupInfo))


        elif self.infoType == 'authors':
                
          #
          # Entry author information
          #
          
          tableName = '_Entry_author'
          if self.saveFrame.tables.has_key(tableName):
            tableTags = self.saveFrame.tables[tableName].tags
    
            for i in range(0,len(tableTags['Ordinal'])):
    
              personInfo = []
              for personTagName in personTagNames[:-1]:
                personInfo.append(tableTags[personTagName][i])           
              
              self.persons.append(Person(self,serial = tableTags['Ordinal'][i],*personInfo))
      
      #
      # Citations
      #
         
      elif self.saveFrame.name == 'citations':
                
        #
        # Create citation information...
        #

        citationType = self.saveFrame.tags['Type']
        
        if self.saveFrame.tags['Class'] ==  "entry citation":
          self.isPrimary = True
        else:
          self.isPrimary = False
        
        citationTagNames = ['Type','Title','Status']
        
        if citationType == 'journal':
          citationTagNames += ['Journal_abbrev','Journal_volume','Journal_issue','Page_first','Year',None,None,None,'isPrimary','Page_last','Journal_name_full','PubMed_ID','Details',None,None,'MEDLINE_UI_code','DOI']
        elif citationType == 'book':
          citationTagNames += [None,'Book_volume',None,'Page_first','Year','Book_publisher_city','Book_publisher',None,'isPrimary','Page_last',None,'PubMed_ID','Details','Book_title','Book_series','MEDLINE_UI_code','DOI']
        elif citationType == 'conference' or citationType == 'abstract':
          citationTagNames += [None,None,None,'Page_first','Year','Conference_site',None,'Conference_country','isPrimary','Page_last',None,'PubMed_ID','Details','Conference_title',None,'MEDLINE_UI_code','DOI']
        elif citationType == 'thesis':
          citationTagNames += [None,None,None,'Page_first','Year','Thesis_institution_city','Thesis_institution','Thesis_institution_country','isPrimary','Page_last',None,'PubMed_ID','Details',None,None,'MEDLINE_UI_code','DOI']
        elif citationType == 'BMRB only':
          citationTagNames += [None,None,None,'Page_first','Year',None,None,None,'isPrimary','Page_last',None,'PubMed_ID','Details',None,None,'MEDLINE_UI_code','DOI']

        # TODO
        #'Journal_ASTM': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Journal_ISSN': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Journal_CSD': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Book_chapter_title': [None,lambda x = value: returnStarString(x,length = 255),None,False],
        #'Book_ISBN': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Conference_state_province': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Conference_start_date': [None,returnStarDateTime,None,False],
        #'Conference_end_date': [None,returnStarDateTime,None,False],
        #'Conference_abstract_number': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'WWW_URL': [None,returnStarString,None,False],
        
        citationInfo = [self]
        for tagName in citationTagNames:
          value = None
          if tagName in ('isPrimary',):
            value = getattr(self,tagName)
          elif self.saveFrame.tags.has_key(tagName):
            value = self.saveFrame.tags[tagName]
          citationInfo.append(value)

        self.citations.append(Citation(*citationInfo))
        
        #
        # Citation author information
        #
        
        tableName = '_Citation_author'
        if self.saveFrame.tables.has_key(tableName):
          tableTags = self.saveFrame.tables[tableName].tags
  
          for i in range(0,len(tableTags['Ordinal'])):
  
            personInfo = []
            for personTagName in personTagNames[:-1]:
              personInfo.append(tableTags[personTagName][i])           
            
            self.persons.append(Person(self,serial = tableTags['Ordinal'][i],*personInfo))
            
            self.citations[-1].setAuthor(self.persons[-1])

        tableName = '_Citation_editor'
        if self.saveFrame.tables.has_key(tableName):
          tableTags = self.saveFrame.tables[tableName].tags
  
          for i in range(0,len(tableTags['Ordinal'])):
  
            personInfo = []
            for personTagName in personTagNames[:-1]:
              personInfo.append(tableTags[personTagName][i])           
            
            self.persons.append(Person(self,serial = tableTags['Ordinal'][i],*personInfo))
            
            self.citations[-1].setEditor(self.persons[-1])

        tableName = '_Citation_keyword'

        tagNames = ('Keyword',)

        if self.saveFrame.tables.has_key(tableName):
          tableTags = self.saveFrame.tables[tableName].tags
  
          for i in range(0,len(tableTags['Keyword'])):
  
            infoList = []
            for tagName in tagNames:
              infoList.append(tableTags[tagName][i])           
            
            self.keywords.append(Keyword(self,*infoList))

            self.citations[-1].setKeyword(self.keywords[-1])
            
    elif self.version == '2.1.1': 
      pass

class PersonInGroup:
  
  def __init__(self, parent, emailAddr, dept, mailAddr, addr1, addr2, addr3, city, stateProv, country, postCode, phoneNo, faxNo, role, orgType, serial = None):

    self.parent = parent
    
    self.emailAddr = emailAddr
    self.dept = dept
    self.mailAddr = mailAddr
    self.addr1 = addr1
    self.addr2 = addr2
    self.addr3 = addr3
    self.city = city
    self.stateProv = stateProv
    self.country = country
    self.postCode = postCode
    self.phoneNo = phoneNo
    self.faxNo = faxNo
    self.role = role
    self.orgType = orgType
    
    self.serial = serial
