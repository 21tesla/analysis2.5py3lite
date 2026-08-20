"""
======================COPYRIGHT/LICENSE START==========================

Util.py: Useful functions for scripts in this directory and its subdirectories

Copyright (C) 2005-2008 Wim Vranken (European Bioinformatics Institute)

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
from ccp.format.general.Constants import defaultSeqInsertCode
from memops.universal.ElementTree import ElementTree as ET

# characters that can be used as a one-digit index
indexChars = '123456789'

#######################
# Regular expressions #
#######################

def getRegularExpressions(format = None):

  patt = {}
  
  patt['emptyline'] = re.compile("^\s*$")
  patt['hash'] = re.compile("^\s*\#")
  patt['exclamation'] = re.compile("^\s*\!")
  patt['underscore'] = re.compile("^\s*\_")
  patt['colon'] = re.compile(":")
  patt['onlyDollar'] = re.compile("^\$$")
  patt['onlydigit'] = re.compile("^\d+$")
  patt['onlyFloat'] = re.compile("^\d+\.\d+$")
  patt['bracketEither'] = re.compile("\(|\)")
  patt['bracketOpen'] = re.compile("\s*\(")
  patt['bracketClose'] = re.compile("\)")
  patt['nucleusLetterDigit'] = re.compile("^([A-Za-z])(\d+)")
  patt['seqInsertCode'] = re.compile("^([A-Za-z]*)(\d+)([A-Za-z]*)")
  patt['anySpace'] = re.compile("\s")
  
  if format == 'amber':

    patt[format + 'RestraintStart'] = re.compile("(\&rst)")
    patt[format + 'RestraintEnd'] = re.compile("(\&end)")
  
  elif format == 'ansig':

    patt[format + 'CurlyBrace'] = re.compile("\{(.*)\}")

  elif format == 'auremol':

    patt[format + 'Section'] = re.compile("^\s*section_(.+)$")
    patt[format + 'SubSection'] = re.compile("^\s*([A-Z_]+)\:")

  elif format == 'autoAssign':
    patt[format + 'specPhase'] = re.compile("phase:\s*\{(.*)\}\s*\{\s*$")
    patt[format + 'onlyCurlyBraceEnd'] = re.compile("^\}$")   
    patt[format + 'CurlyBrace'] = re.compile("\{\s*(.*)\s*\}")
    patt[format + 'seqCode1Or3LetterAndCode'] = re.compile("([A-Z][a-z]{0,2})(\d+)")
    patt[format + 'spinSystemInfo'] = re.compile(".+(\(.+\))")

  elif format == 'aqua':
    patt[format + 'Count'] = re.compile("^count (\d+) type (.+)$")   
    patt[format + 'UpperDistance'] = re.compile("^NOEUPP")
    patt[format + 'LowerDistance'] = re.compile("^NOELOW")
    patt[format + 'EndComment']    = re.compile("\#(.+)$")

  elif format == 'bruker':
    patt[format + 'StartDoubleHash'] = re.compile("^\#\#")
    patt[format + 'StartDoubleHashDollar'] = re.compile("^\#\#\$")
    patt[format + 'StartDoubleDollar'] = re.compile("^\$\$")
    patt[format + 'SharpBracketOpen'] = re.compile("^\<")
    patt[format + 'SharpBracketClose'] = re.compile("\>$")
    patt[format + 'SharpBracketEither'] = re.compile("\<|\>")
    patt[format + 'BracketsListIndicator'] = re.compile("\(\d+\.\.\d+\)")
    patt[format + 'Date'] = re.compile("2\d\d\d-\d\d-\d\d")
    patt[format + 'DigitDotDigit'] = re.compile("^\-?\d*\.?\d+$")
    patt[format + 'Dot'] = re.compile("\.")
    patt[format + 'BracketCloseNumber'] = re.compile("\)\s*\d")
    patt[format + 'FinalReturn'] = re.compile("\n$")
    patt[format + 'BracketMultiValue'] = re.compile("\((\d+)\-(\d+)\)")
    patt[format + 'InOrDecrement'] = re.compile("\s*(i|d)d(\d+)")
    patt[format + 'HashComment'] = re.compile("^\s*\#(.+)")
    patt[format + 'SemiColonComment'] = re.compile("^\s*\;(.*)")
    patt[format + 'PulseName'] = re.compile("\/([^\/]+)\"")
    patt[format + 'FnMode'] = re.compile("FnMode?\s*(.+)\s*", re.I)
    patt[format + 'EchoAntiEcho'] = re.compile("igrad|\*EA", re.I)
    patt[format + 'PathwayElements'] = re.compile("(F\d+)\s*\(([^\)]+)\)")
    patt[format + 'DimIncrement'] = re.compile("in(\d+)\s*\=\s*inf(\d+)\s*\\s*/\s*(\d+)")
    patt[format + 'phaseSensitive'] = re.compile("\;\s*phase sensitive")
    patt[format + 'constantTime'] = re.compile("constant time")
    patt[format + 'semiconstantTime'] = re.compile("semi-constant time")
    patt[format + 'tSearch'] = re.compile("(t\d)")
    patt[format + 'Version'] = re.compile("Version\s*(\d+\.*\d*)")

  elif format == 'charmm':

    patt[format + 'NumberAtoms'] = re.compile("^\s*(\d+)\s*$")
    patt[format + 'AtomLine'] = re.compile("^\s*\d+\s+\d+\s+[A-Za-z]+\s+")

  # This is for CNS/xplor
  elif format == 'cns':
  
    patt[format + 'DistancePeakInfoLine'] = re.compile("peak",re.IGNORECASE)
    patt[format + 'MultiSign'] = re.compile("[\%\#\*\+]")
    patt[format + 'RestrNum'] = re.compile("\{\s*(\-?\d+)\s*\}")
    patt[format + 'Assign'] = re.compile("assi",re.IGNORECASE)
    patt[format + 'AssignOr'] = re.compile("or",re.IGNORECASE)
    patt[format + 'ChemShiftFormat'] = re.compile("attr\s+store",re.IGNORECASE)
    patt[format + 'ChemShiftStore'] = re.compile("do.+store\d+\s*\=\s*(\-?\d*\.?\d*)\s*\)",re.IGNORECASE)
    patt[format + 'ChainCode'] = re.compile("segid\s+(\"([A-Za-z0-9 ]+)\"|([A-Za-z0-9]+))\s*",re.IGNORECASE)
    patt[format + 'Class'] = re.compile("class\s+",re.IGNORECASE)
    patt[format + 'SeqCode'] = re.compile("resid?u?e?\s+(\-?\d+[a-zA-Z]?)\s*",re.IGNORECASE)
    patt[format + 'AtomName'] = re.compile("name\s+([A-Z0-9\%\#\*\+]+\'?)",re.IGNORECASE)
    patt[format + 'RestrDistances'] = re.compile("([\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.?]+)")
    patt[format + 'RestrAngles'] = re.compile("(\d+\.?\d*)\s+(-?\+?\d+\.?\d*)\s+(\+?\d+\.?\d*)\s+(\d)")
    patt[format + 'RestrCoupling'] = re.compile("(\d+\.?\d*)\s+(\d+\.?\d*)")
    patt[format + 'RestrCsa'] = re.compile("(-?\d+\.?\d*)\s+(\d+\.?\d*)\s*(\d+\.?\d*)?")
    patt[format + 'RestrRdc'] = re.compile("(-?\d+\.?\d*)\s+(\d+\.?\d*)\s*(\d+\.?\d*)?")
    patt[format + 'RestrInnerOr'] = re.compile("\s*or",re.IGNORECASE)
    patt[format + 'InnerElementPatt'] = re.compile("\(([^\)\(]+)\)",re.M)
    patt[format + 'LongCommentStart'] = re.compile("\{")
    patt[format + 'LongCommentEnd'] = re.compile("\}")

  elif format == 'csi':
  
    patt[format + 'Comment'] = re.compile("^\s*\>")

  elif format in ['dyana','cyana']:
  
    patt[format + 'CoordinateInfoLine'] = re.compile("\s*\d+\s+[A-Z0-9]+\s+\d+\s+")
    patt[format + 'CoordinateAtomLine'] = re.compile("^ATOM|^HETATM")
    patt[format + 'NewModel'] = re.compile("^MODEL\s+(\d+)")

  elif format == 'mars':
    
    patt[format + 'AtomNameHeader'] = re.compile("^\s*([\sNCHABO\-1])+\s*$")
    patt[format + 'AtomInfo'] = re.compile("([A-Z]+)(\-\d)?")

  elif format == 'mol':

    patt[format + 'Counts'] = re.compile("^([0-9 ][0-9 ][0-9]){6}")
    patt[format + 'Atoms'] = re.compile("^([0-9\- ][0-9\- ][0-9\- ][0-9\- ][0-9]\.[0-9][0-9 ][0-9 ][0-9 ]){3}")
    patt[format + 'Bonds'] = re.compile("^([0-9 ][0-9 ][0-9]){4}([0-9 ][0-9 ][0-9 ])([0-9 ][0-9 ][0-9])")

  elif format == 'mol2':
    
    patt[format + 'TriposTag'] = re.compile("\@\<TRIPOS\>(.+)")

  elif format == 'monte':

    patt[format + 'Comment'] = re.compile("^\s*(\#\#|\%\%)")
    patt[format + 'Assignment'] = re.compile("([A-Za-z]+)?([0-9]+|\?)?")
    patt[format + 'AtomInfo'] = re.compile("([A-Za-z0-9]+)(\(.+\)|\-1)?")

  elif format in ['nmrDraw','talos']:
  
    patt[format + 'Remark'] = re.compile("^REMARK")
    patt[format + 'Dataline'] = re.compile("^DATA") 
    patt[format + 'Vars'] = re.compile("^VARS")
    patt[format + 'Format'] = re.compile("^FORMAT")
    
  elif format == 'nmrView':

    patt[format + 'DigitSpace'] = re.compile("^\d+\s+")
    patt[format + 'CurlyBrace'] = re.compile("\{([^{}]*)\}")
    patt[format + 'CurlyBraceStart'] = re.compile("\{")    
    patt[format + 'CurlyBraceEnd'] = re.compile("\}\s+")    
    patt[format + 'NumbersNoBrace'] = re.compile("^(\d+\.\d+\s*){2}")

  elif format == 'nmrStar':

    patt[format + 'EndSaveTag'] = re.compile("save_$")

  elif format == 'pdb':

    patt[format + 'NewModel'] = re.compile("^MODEL\s+(\d+)")
    patt[format + 'AllAtom'] = re.compile("^ATOM|^HETATM")
    patt[format + 'HetAtom'] = re.compile("^HETATM")
    patt[format + 'Header'] = re.compile("^HEADER")
    patt[format + 'Title'] = re.compile("^TITLE")
    patt[format + 'Remark4'] = re.compile("^REMARK   4 ([a-zA-Z0-9]{4}) COMPLIES WITH FORMAT V.\s*(\d+\.\d+)\s*\,")
    patt[format + 'Compound'] = re.compile("^COMPND\s+\d*\s+")
    patt[format + 'DbReference'] = re.compile("^DBREF")
    patt[format + 'SequenceChange'] = re.compile("^SEQADV")
    patt[format + 'Sequence'] = re.compile("^SEQRES")
    patt[format + 'Source'] = re.compile("^SOURCE\s+\d*\s+(.+)")
    patt[format + 'Keywds'] = re.compile("^KEYWDS\s+\d*\s+(.+)")
    patt[format + 'ExpData'] = re.compile("^EXPDATA\s+\d*\s+(.+)")
    patt[format + 'Authors'] = re.compile("^AUTHOR\s+\d*\s+(.+)")
    patt[format + 'Journal'] = re.compile("^JRNL\s+(AUTH|TITL|REF|REFN|PUBL|EDIT)\s+\d*\s+(.+)")
    patt[format + 'Reference'] = re.compile("^REFERENCE\s+(\d+)")
    patt[format + 'ReferenceJournal'] = re.compile("\s*(AUTH|TITL|REF|REFN|PUBL|EDIT)\s+\d*\s+(.+)")
    patt[format + 'Remarks'] = re.compile("^REMARK\s+(\d+)\s+(.+)")
    patt[format + 'HetGroup'] = re.compile("^HET\s+")
    patt[format + 'HetName'] = re.compile("^HETNAM")
    patt[format + 'HetSynonym'] = re.compile("^HETSYN")
    patt[format + 'HetFormula'] = re.compile("^FORMUL")
    patt[format + 'Bonds'] = re.compile("^CONECT")
    patt[format + 'SsBond'] = re.compile("^SSBOND")
    patt[format + 'Link'] = re.compile("^LINK")
    patt[format + 'SecStruc'] = re.compile("^(HELIX|TURN|SHEET)")
    patt[format + 'Termination'] = re.compile("^TER ")

  elif format == 'pipp':
    patt[format + '?_AXIS'] = re.compile("^(.+)_AXIS")
    patt[format + 'Shift'] = re.compile("^\s*(.+)\s+([0-9]*\.?[0-9]+)\s+\((.+)\)\s*$")
    patt[format + 'ShiftNoAss'] = re.compile("^\s*(.+)\s+([0-9]*\.?[0-9]+)\s*$")

  elif format == 'pistachio':
    patt[format + 'Comment'] = re.compile("^\s*\%.+")
  
  elif format == 'sparky':
    patt[format + 'SharpBracketBetween'] = re.compile("^\<(.+)\>$")
    patt[format + 'LabelCodeName'] = re.compile("([A-Z]?)((\d+\,?)+)(.+)")
    patt[format + 'BracketBetween'] = re.compile("\((.+)\)")

  # Addition Maxim Mayzel, Gothenburg
  elif format == 'targetedAcquisition':
    patt[format + 'seqCode1Or1LetterAndCode'] = re.compile("([A-Z])(\d+)")
    #patt[format + 'spinSystemInfo'] = re.compile(".+(\(.+\))")
    patt[format + 'figOfMerit'] = re.compile("([B|G|M])\s*(\d+\.*\d*)")

  elif format == 'xeasy':
    patt[format + 'IName'] = re.compile("^\#\s*INAME\s+(\d)\s*(.+)$")
    patt[format + 'CyanaFormat'] = re.compile("^\#\s*CYANAFORMAT\s+(.+)\s*$")
    patt[format + 'PeakInfo'] = re.compile("^\s*([^#]+)\s*\#?.*$")


  """
  
  #
  # TODO: WHERE ARE THESE USED!!?!
  #
  patt['infoDirs'] = re.compile("^info\..+")
  patt['projectDirs'] = re.compile("^project\..+")
  patt['varFiles'] = re.compile("^vars\..+")
  patt['nmrStarFiles'] = re.compile("^bmr(\d+)\.str")
   # patt[format + 'DMXTypes'] = re.compile("^DMX|DRX|Avance$")
   # patt[format + 'AllTypes'] = re.compile("^DMX|DRX|Avance|AM|AMX$")

  # This is for nmrPipe
  patt['var2pipe'] = re.compile("var2pipe$")
  patt['bruk2pipe'] = re.compile("bruk2pipe$")
  patt['pipe'] = re.compile("\|")
  patt['backslash'] = re.compile('\\\\')
  patt['dashstart'] = re.compile("^\-")
  patt['x'] = re.compile("^x")
  patt['y'] = re.compile("^y")
  patt['z'] = re.compile("^z")
  patt['a'] = re.compile("^a[^(q2d)]")

  # This is for dataNavigator
  patt['pipe2Dproc'] = re.compile("^test\.ft$")
  patt['pipe3Dproc'] = re.compile("^test(\d+)\.ft$")

  # These are for nmrPipe and Azara exec... scripts.
  patt['pipeConv'] = re.compile("^autoconv\.([HCN]+)\.com$")
  patt['pipeProc'] = re.compile("^autoproc\.([HCN]+)\.com$")

  # These are for DataManager
  patt['Varian'] = re.compile("procpar$")
  patt['Bruker'] = re.compile("acqu\d?s?$")

  # For converters, exportNmrProcPars
  patt['sensitivityEnhanced'] = re.compile("(echoAntiEcho|ranceKay) dim (\d)")
  """
  return patt

patt = getRegularExpressions()

def getSeqAndInsertCode(seqCode):

  seqInsertCode = defaultSeqInsertCode
  
  if seqCode != None:
    try:
      seqCode = int(seqCode)
    except:
      searchObj = patt['seqInsertCode'].search(str(seqCode))
      if searchObj:
        try:
          seqCode = int(searchObj.group(2))       
          if searchObj.group(1) or searchObj.group(3):
            seqInsertCode = searchObj.group(1) + searchObj.group(3)

        except:
          seqCode = None
      else:
        seqCode = None
  
  return (seqCode,seqInsertCode)
  
def standardNucleusName(name):

  if (not name):
    name = '1H'
  elif ((name[0].upper() == 'H') or (name.upper() == '1H')):
    name = '1H'
  elif ((name[0].upper() == 'C') or (name.upper() == '13C')):
    name = '13C'
  elif ((name[0].upper() == 'N') or (name.upper() == '15N')):
    name = '15N'
  elif ((name[0].upper() == 'P') or (name.upper() == '31P')):
    name = '31P'
  else:
    name = '1H'

  return name


def parseProgramName(fname, parseInfo, specialCases=()):
  """
  Split fname into components using passed-in input
  
  Input: 
    fname: string to parse
    parseInfo: dictionary with a list for each of the following keys:
      skipNamesEnding, skipNamesStarting, 
      prefixes, startnames, midnames, endnames
    specialCases: list of (matchExpr, func) tuples
  
  Output:
    result = {
     'prefixes': [], 'names': [], 'endnames': [], 'suffix': '', 'unused': '',
     'status':None
    }
    #status is in None, 'skip', 'partial', or 'OK'
  
  Parsing rules:
  - Try to apply input lists in order given
  - midnames may appear multiple times - prefixes and endnames may not
  - Return first complete parsing of fname
  - Failing that return parsing that leaves shortest 'unused' string
  - 0) fname that match a matchExpr in specialCases is treated by relevant func
  For each expname in startnames:
  - 1) Look for match of expname anywhere in filename
  - 2) if found check if text before namestart fits concatenation of prefixes
  - 3) If OK look for secondary names at start of remaining text
  - 4) Parse remaining text for endnames:
    - 4.1) anything after first dot '.' considered suffix
    - 4.2) underscores split text into substrings parsed separately
    - 4.3) in each substring look for endnames starting at left.
    - 4.4) if a substring can be parsed except for a terminal single non-zero
           digit (123456789) the digit is treated as part of the suffix.
    - 4.5) Any left-over text is returned as unused, unused from separate
           substrings being joined with underscores ('_')
  - 5) If no startname is found, try to find suffix and endnames from end of 
       fname and return start as unused.
  """
  
  # default result
  result = {
   'names': [],
   'prefixes': [],
   'unused': fname,
   'endnames': [],
   'suffix': '',
   'status':None,
  }
  
  # skip files with undesirable endings
  for tag in parseInfo['skipNamesEnding']:
    if fname.endswith(tag):
      result['status'] = 'skip'
      return result
      
  # skip files with undesirable beginnings
  for tag in parseInfo['skipNamesStarting']:
    if fname.startswith(tag):
      result['status'] = 'skip'
      return result
  
  # teat special cases:
  for expr, func in specialCases:
    if re.match(expr, fname):
      specialResult = func(fname, parseInfo)
      if specialResult is not None:
        return specialResult
  
  # start matching on startnames
  numUnused = len(fname)
  foundprefs = []
  for testname in parseInfo['startnames']:
    
  
    nameend = 0
    foundprefs = []
    nameat = fname.find(testname)
    
    if nameat >= 0:
      # we have some match
      
      if nameat == 0:
        # name starts at beginning
        nameend = len(testname)
 
      elif nameat > 0:
        # name starts further in - check start of fname for prefixes
        text = fname[:nameat]
        nameend = len(testname) + nameat
        while True:
          for prefix in parseInfo['prefixes']:
            if text.startswith(prefix) and prefix not in foundprefs:
              foundprefs.append(prefix)
              text = text[len(prefix):]
              break
          else:
            break
      
        if text:
          # prefixes did not exhaust start of fname - try with other startName
          continue 
      
      # we have a startName
      foundnames = [testname]
      text = fname[nameend:]
      
      # now search for additional names:
      while True:
        for name2 in parseInfo['midnames']:
          if text.startswith(name2):
            foundnames.append(name2)
            text = text[len(name2):]
            break
        else:
          break
          
      
      unused, endnames, suffix = parseEndNames(parseInfo['endnames'], text)
      
      if len(unused) < numUnused:
        # improved match
        result['names'] = foundnames
        result['prefixes'] = foundprefs
        result['unused'] = unused
        result['endnames'] = endnames
        result['suffix'] = suffix
      
        if unused:
          # this is still a partial match
          numUnused = len(unused)
          result ['status'] = 'partial'
  

        else:
          # full match - stop looking
          result ['status'] = 'OK'
          break
  
  if not result['names']:
    # no match found - parse codes from string end
    unused, endnames, suffix = reverseParseEndNames(parseInfo['endnames'], 
                                                    fname)
    result['unused'] = unused
    result['endnames'] = endnames
    result['suffix'] = suffix
  
  
  #
  return result


def parseEndNames(endNameList, text):
  """ parse text for endnames
  """
  
  # set up
  endnames = []
  unuseds = []
  suffixes = []
  
  # split off dot-separated suffix
  ll = text.split('.',1)
  if len(ll) > 1:
    text, suffix = ll
  else:
    suffix = ''
  
  # split on underscore and treat components separately
  for text in text.split('_'):
  
    while text:
    
      for endName in endNameList:
        if endName not in endnames and text.startswith(endName):
          # found matching endName - chop it off
          endnames.append(endName)
          text = text[len(endName):]
          break
 
      else:
        # no matching endName found
        # still text left
        if len(text) == 1 and text in indexChars:
          # single digit - add to suffix
          suffixes.append(text)
        else:
          # parting incomplete - add to unuseds
          unuseds.append(text)
 
        break
  #
  unused = '_'.join(unuseds)
  if suffix:
    suffixes.append(suffix)
  #
  return (unused, endnames, '_'.join(suffixes))



def reverseParseEndNames(endNameList, text):
  """ parse endnames from text end 
  """
  
  endnames = []
  unused = ''
  suffix = ''
  
  # split dot-separated suffix
  ll = text.split('.',1)
  if len(ll) > 1:
    text, suffix = ll
  
  while text:
  
    if text.endswith('_'):
      # remove trailing underscore
      text = text[:-1]
    
    else:
      for endName in reversed(endNameList):
        if endName not in endnames and text.endswith(endName):
          # found matching endName - chop it off
          endnames.append(endName)
          text = text[:-len(endName)]
          break
      else:
        # no matching endName found
        unused = text
        break
        
  # reverse to normal order
  endnames.reverse()
  
  #
  return (unused, endnames, suffix)


def convertParseInfo(parseInfo):
  """ Convert parsing information to form needed for parsing
  """
  useTags = ['skipNamesEnding', 'skipNamesStarting', 
             'prefixes', 'midnames', 'endnames']
  
  result = {}
  for tag in useTags:
    result[tag] = [x[0] for x in parseInfo[tag]]
  
  result['nameMaps'] = dd = {}
  for name, ll in parseInfo['nameMaps']:
    for synonym in ll:
      dd[synonym] = name
  
  result ['expNameMaps'] = parseInfo['expNameMaps']
  
  result['startnames'] = [dd['namesStartWith'][0] 
                          for dd in parseInfo['expNameMaps']]
  #
  return result
 

def loadAllParseInfo(stdFilePath, localFilePath=None):
  """ load string parsing info from stdFilePath and from localFilePath
  If localFilePath is given its contents will be at the head of the lists
  overriding the standard values
  """
  
  stdData = loadParseInfo(stdFilePath)
  
  
  if localFilePath is None:
    result = stdData
  
  else:
    result = loadParseInfo(localFilePath)
    for tag in stdData:
      ll = result.get(tag)
      if ll:
        ll.extend(stdData[tag])
      else:
        result[tag] = stdData[tag]
  #
  return result

def loadParseInfo(filePath):
  """ load string parsing info from XML file
  
    Input: path to file
    
    Output: a dictionary that reflects the file information
            See code for details
  
  """
  
  # set up result
  result = {}
  result['skipNamesEnding'] = skipNamesEnding = []
  result['skipNamesStarting'] = skipNamesStarting = []
  result['nameMaps'] = nameMaps = []
  result['prefixes'] = prefixes = []
  result['midnames'] = midnames = []
  result['endnames'] = endnames = []
  result['expNameMaps'] = expNameMaps = []
  
  # parse file
  eTree = ET.parse(filePath)
  
  # read skipNamesEnding
  container = eTree.find('skipNamesEnding')
  for elem in container.findall('StringWithDoc'):
    skipNamesEnding.append((elem.get('name'), elem.get('text')))
  
  # read skipNamesStarting
  container = eTree.find('skipNamesStarting')
  for elem in container.findall('StringWithDoc'):
    skipNamesStarting.append((elem.get('name'), elem.get('text')))
  
  # read prefixes
  container = eTree.find('prefixes')
  for elem in container.findall('StringWithDoc'):
    prefixes.append((elem.get('name'),elem.get('text')))
  
  # read midnames
  container = eTree.find('midnames')
  for elem in container.findall('StringWithDoc'):
    midnames.append((elem.get('name'), elem.get('text')))
  
  # read endnames
  container = eTree.find('endnames')
  for elem in container.findall('StringWithDoc'):
    endnames.append((elem.get('name'), elem.get('text')))
  
  # read nameMappings
  for elem in eTree.findall('NameMap'):
    nameMaps.append((elem.get('name'), elem.get('synonyms').split()))
  
  # read name-prototype mapping
  for elem in eTree.findall('ExpNameMap'):
    dd = {}
    expNameMaps.append(dd)
    dd['namesStartWith'] = elem.get('namesStartWith', '').split()
    dd['hasnames'] = elem.get('hasnames', '').split()
    dd['hasprefixes'] = elem.get('hasprefixes', '').split()
    dd['doc'] = elem.get('text')
    dd['expPrototypes'] = expPrototypes = []
    for el2 in elem.findall('ExpPrototype'):
      dd2 = {}
      expPrototypes.append(dd2)
      dd2['name'] = el2.get('name')
      dd2['isReversed'] = el2.get('isReversed')
  #
  return result


#def testNameParsing(infoFile, altInfoFile=None, specialCases=None, *testnames):
#  """ Parse testnames using informatio in infoFile and altInfoFile
#  return list of dicts with results
#  """
#  
#  parseInfo = loadAllParseInfo(infoFile, altInfoFile)
#  
#  useInfo = convertParseInfo(parseInfo)
#  
#  result = {}
#  for testName in testnames:
#    result[testName] = parseProgramName(testName, useInfo, specialCases)
#  #
#  return result
  
  
  

