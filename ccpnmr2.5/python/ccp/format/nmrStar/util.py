
"""
======================COPYRIGHT/LICENSE START==========================

util.py: Useful functions for scripts in this directory

Copyright (C) 2005-2010 Wim Vranken (European Bioinformatics Institute)

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
import re

from memops.universal.Util import returnInt, returnFloat

anySpacePatt = re.compile("\s")
multiSpacePatt = re.compile("\s{2,}")
bracketsPatt = re.compile("[\(\)]")
yesPatt = re.compile("^[yYtT]")
noPatt = re.compile("^[nNfF]")

#############
# Functions #
#############

def removeMultiSpaces(value):

  value = str(value.strip())

  if multiSpacePatt.search(value):
    (value,numSubs) = multiSpacePatt.subn(' ',value)

  return value

def getNmrStarValue(tagtable,tagname,valueIndex = 0):

  if tagtable.tagnames.count(tagname):
    tagIndex = tagtable.tagnames.index(tagname)
    return tagtable.tagvalues[tagIndex][valueIndex]
    
  else:
    return None
       
def returnStarInt(value):

  if value in ('.','?') or value is None:
    return None
  else:
    return returnInt(value)

def returnStarFloat(value):

  if value in ('.','?') or value is None:
    return None
  else:
    return returnFloat(value)

def returnStarString(value,strip = 0, length = None):

  value = str(value)
  
  if strip:
    value = value.strip()
  
  #if value.count("\n"):
  #  value = value.replace("\n"," ")

  value = removeMultiSpaces(value)

  if length:
    origValue = value
    value = value[:length]
    if len(origValue) > length:
      print "  Warning: shortened value '%s' to '%s' for NMR-STAR export!" % (origValue,value)

  if value in ('.','?') or not value:
    return None
  else:
    return value

def returnStarDateTime(value, length = None):
  # TODO need to do this correctly...
  return returnStarString(value,length = length)

def returnStarAtCode(value,length = None):
  # Basically accepts more or less anything, just do default.
  return returnStarString(value,length = length)

def returnStarCode(value,length = None):

  value = str(value)
  
  value = removeMultiSpaces(value)

  if anySpacePatt.search(value):
    (value,numSubs) = anySpacePatt.subn('_',value)

  return returnStarString(value,length = length)

def returnStarLine(value,length = None):

  value = str(value)

  #values = value.split("\n")
  #values = [ val.strip() for val in values ]
  #value = "\n".join(values)

  #if value.count("\n"):
  #  value = value.replace("\n"," ")

  value = removeMultiSpaces(value)

  value = value.strip()

  return returnStarString(value,length = length)

def returnStarName(value,length = None):

  if not value[0] == '_':
    value = '_' + value
  
  return returnStarCode(value,length = length)
    
def returnStarIdName(value,length = None):

  return returnStarCode(value,length = length)
    
def returnStarYesNo(value,length = None):

  if value:
    searchObj1 = yesPatt.search(str(value) )
    searchObj2 = noPatt.search(str(value) )
    if searchObj1:
      value = 'yes'
    elif searchObj2:
      value = 'no'
    else:
      value = 'yes'

  else:
    value = 'no'

  return value

def returnStarFaxPhoneEmail(value,length = None):

  value = str(value)

#  if bracketsPatt.search(value):
#    (value,numSubs) = bracketsPatt.subn('_',value)
    
  return returnStarString(value,length = length)

def returnStarLabel(value,length = None):

  if not value or value in ('.',):
    value =  ''

  elif not value[0] == '$':
    value = '$' + value

  return returnStarCode(value,length = length)
