
"""
======================COPYRIGHT/LICENSE START==========================

Util.py: Utility code for CCPN code generation framework

Copyright (C) 2005  (CCPN Project)

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
# miscellaneous useful functions

import math, string, sys, types, os

######################################################################
# hack for Python 2.1 compatibility  NBNB                                                      #
######################################################################
try:
  junk = True
  junk = False
except:
  True = not 0
  False = not 1

def formatFloat(x, places = 3):
  """autoformat float to specified number of significant figures
  currently does not do scientific formatting at all
  """

  try:
    s = int(math.floor(math.log(abs(x)) / math.log(10)))
  except:
    s = 0
  d = max(0, places - s - 1)
  t = '%%.%sf' % d
  return t % x
 
# autoformat float to specified number of decimals
def formatDecimals(x, decimals = 0):
  """autoformat float to specified number of decimals
  """

  t = '%%.%df' % max(0, decimals)
  return t % x
 
def upperFirst(s):
  """uppercase first letter
  """
  return s[0].upper() + s[1:]

def lowerFirst(s):
  """lowercase first letter
  """
  return s[0].lower() + s[1:]

def capitalizeUnderscore(s):
  """capitalize entire string putting underscores in
  before existing capital letters (except first one)
  """

  t = ''
  for n in range(len(s)):

    if ((n > 0) and (s[n] in string.uppercase)):
      t = t + '_'

    t = t + s[n].upper()

  return t

def substituteName(string, name, n = 1):
  """substitute name n times in format string
  """

  return string % (n * (name,))

def toBoolean(x):
  """ Convert x to true/false. 
  Could be inlined, but function is more intelligible.
  """
  return x and True or False

def checkListIsSet(w):
  """ Checks if w (which could be list or tuple) is a set, i.e. no repeats.
  Creates a copy of the list or tuple to do the check.
  """
  z = list(w)
  for x in w:
    z.remove(x)
    assert x not in z, 'repeated element "%s"' % x

def isArray(x):
  """ Returns true if x is tuple or list, false otherwise.
  """

  if (type(x) in (types.TupleType, types.ListType)):
    return True
  else:
    return False

def isBigEndian():

  """ Returns true if platform is big endian, false if little endian.
  """

  if (sys.byteorder == 'big'):
    return True
  else:
    return False

def productArray(array):

  """ Returns product of entries of array
  """

  return reduce(lambda x, y: x*y, array, 1)

def cumulativeProductArray(array):

  """ Returns product of entries of array and also cumulative array of products
  """

  n = len(array)
  cumulative = n*[0]
  product = 1
  for i in range(n):
    cumulative[i] = product
    product = product * array[i]

  return (product, cumulative)

def indexOfArray(array, cumulative):

  index = 1
  for i in range(len(cumulative)):
    index = index + array[i] * cumulative[i]

  return index

def arrayOfIndex(index, cumulative):

  n = len(cumulative)
  array = n*[0]
  for i in range(n-1, -1, -1):
    array[i] = index / cumulative[i]
    index = index % cumulative[i]

  return array


def compactStringList(stringList, separator='', maxChars=80):
  """ compact stringList into shorter list of longer strings,
  each either made from a single start string, or no longer than maxChars
  
  From previous breakString function.
  Modified to speed up and add parameter defaults, Rasmus Fogh 28 Aug 2003
  Modified to split into two functions 
  and to add separator to end of each line, Rasmus Fogh 12 Sep 03
  Modified to separate string breaking from list modification
  Rasmus Fogh 29/6/06
  Modified to return single-element lists unchanged
  Rasmus Fogh 29/6/06
  """
  
  result = []
  
  if not stringList:
    return result
  elif len(stringList) ==1:
    return stringList[:]
  
  seplength = len(separator)
  
  nchars = len(stringList[0])
  start=0
  for n in range(1,len(stringList)):
    i = len(stringList[n])
    if nchars + i + (n-start)*seplength > maxChars:
      result.append(separator.join(stringList[start:n] + ['']))
      start = n
      nchars = i
    else:
      nchars = nchars + i
  result.append(separator.join(stringList[start:len(stringList)]))
  
  return result


def divideString(text, separator=' ', maxChars=72):
  """ divide string in series of substrings no longer than maxchars
  
      From previous breakString function.
      Modified to speed up and add parameter defaults, Rasmus Fogh 28 Aug 2003
      Modified to split into two functions
      and to add separator to end of each line, Rasmus Fogh 12 Sep 03
      Modified to separate string breaking from list modification
      Rasmus Fogh 29/6/06
      Added special case for text empty or None
      Rasmus Fogh 07/07/06
  """
  
  if not text:
    return ''
  
  return compactStringList(text.split(separator), separator=separator, 
                           maxChars=maxChars)

def breakString(text, separator=' ', joiner='\n', maxChars=72):
 
  """ Splits text on separator and then joins pieces back together using joiner
      so that each piece either single element or no longer than maxChars
      
      Modified to speed up and add parameter defaults, Rasmus Fogh 28 Aug 2003
      Modified to split into two functions 
      and to add separator to end of each line
      Modified to separate string breaking from list modification
      Rasmus Fogh 29/6/06
      Added special case for text empty or None
      Rasmus Fogh 07/07/06
  """
  
  if not text:
    return ''
  
  t = compactStringList(text.split(separator), separator=separator, 
                        maxChars=maxChars)
  
  return joiner.join(t)

def documentationFormat(text):
  """ Converts text to a multiline string, and returns a string literal
  expression with one line per line that evaluates to the multiline string
  Modified to separate string breaking from list modification
  Rasmus Fogh 29/6/06
  """

  if not text:
    return '""'
  
  ll = []
  
  for ss in text.splitlines(True):
    ll.extend(compactStringList(ss.split(' '), separator=' ',  maxChars=60))
  
  if len(ll) == 1:
    return repr(ll[0])
  
  if not ll[-1]:
    ll[-1] = '\n'
  elif ll[-1][-1] != '\n':
    ll[-1] = ll[-1] + '\n'
  
  return """(%s
)""" % '\n'.join(map(repr,ll))

def returnFloat(x,default = 0.0, verbose = True):

  """
  Returns a float, or default if fails
  Note that <number>e+n strings are converted OK!
  """

  try:
    x = float(x)
  except:
    if verbose:
      print "Error converting '" + str(x) + "' to float: set to %s" % str(default)
    x = default
  return x

def returnFloats(xlist, verbose = True):

  """
  Converts elements of a list to floats
  """

  for n in range(0,len(xlist)):
    xlist[n]=returnFloat(xlist[n], verbose = verbose)
  return xlist

def returnLong(x, default = 0.0, verbose = 1):

  """
  Returns a long, or default if fails
  """

  try:
    x = long(x)
  except:
    if verbose:
      print "Error converting '" + str(x) + "' to long: set to %s" % str(default)
    x = default
  return x

def returnLongs(xlist, verbose = True):

  """
  Converts elements of a list to longs
  """

  for n in range(0,len(xlist)):
    xlist[n]=returnLong(xlist[n],verbose = verbose)
  return xlist

def returnInt(x,default = 0, verbose = True):

  """
  Returns an int, or default if fails
  """

  try:
    x = int(x)
  except:
    if verbose:
      print "Error converting '" + str(x) + "' to integer: set to %s." % str(default)
    
    x = default
      
  return x

def returnInts(xlist,verbose = True):

  """
  Converts elements of a list to ints
  """

  for n in range(0,len(xlist)):
    xlist[n]=returnInt(xlist[n], verbose = verbose)
  return xlist

def returnStrings(xlist):

  """
  Converts elements of a list to strings
  """

  newList = xlist[:]
  for n in range(0,len(newList)):
    newList[n]=str(newList[n])
  return newList

def returnList(value):
  
  """
  Returns a list from a single value or tuple
  """
  
  if type(value) == type(()):
    value = list(value)
  elif value and type(value) != type([]):
    value = [value]
  
  return value

def unquote(cols,quote):
  
  """
  Remove quote for each element in list
  """
  
  for n in range(0,len(cols)):
    if cols[n][0] == quote:
      cols[n] = cols[n][1:]
    if cols[n][-1] == quote:
      cols[n] = cols[n][0:-1]
  return cols

def joinUnquote(cols,quote,joinString = " "):
  
  """
  Remove quote
  """
  
  if cols[0][0] == quote:
    cols[0] = cols[0][1:]
  if cols[-1][-1] == quote:
    cols[-1] = cols[-1][0:-1]
  # Rejoin with single spaces
  return joinString.join(cols)

def drawBox(text,indent = "",liner = "#"):

  """
  Draw a box around a string with an indent and a selected liner
  """

  box = indent
  box += liner * (len(text)+4)
  box += os.linesep + indent + liner + " " + text + " " + liner + os.linesep + indent
  box += liner * (len(text)+4)
  box += os.linesep
  return box 
  
def getUpperPowerTwo(value):
  
  """
  Get next or current power of two
  
  NOTE: if value = power of two, then this function returns this value
  """
  
  for i in range(0,25):

    t = 2 ** i
    
    if t >= value:
      return t

  return None

def getPowerTwoExp(value):
  
  """
  Gives the exponent (2 ** n) that leads to the value
  (or the first power of two above this value)
  """
  
  for i in range(0,25):

    t = 2 ** i
    
    if t >= value:
      return i
  
  return None


def getMeanValue(valueList, valueNumber = None, valueTotal = None):
  
  """
  Get the mean value for the values in valueList with total
  valueTotal and number of values valueNumber (if given, else is calculated)
  """

  if not valueTotal:
    
    valueTotal = 0
    
    for value in valueList:
      valueTotal += value
      
  if not valueNumber:
    valueNumber = len(valueList)

  if valueNumber < 1:
    valueMean = None
  else:
    valueMean = (valueTotal * 1.0) / valueNumber
  
  return valueMean

def getStandardDev(valueList,valueTotal = None):

  """
  Get the standard deviation for the values in valueList with total
  valueTotal (if given, else is calculated)
  """
  
  valueNumber = len(valueList)
  valueAverage = getMeanValue(valueList,valueNumber = valueNumber, valueTotal = valueTotal)  
  
  if valueAverage == None:
    return None

  squaredSum = 0.0
  
  for value in valueList:
  
    squaredSum += (value - valueAverage) ** 2
  
  if valueNumber > 1:
    standardDev = (squaredSum/(valueNumber - 1)) ** 0.5
  else:
    standardDev = 0.0
  
  return standardDev

def getRms(valueList, total = None):
  
  """
  Calculates the root mean square for a list of values.
  """

  sqSum = 0.0
  
  for value in valueList:
    sqSum += value ** 2
    
  if not total:
    total = len(valueList)
  
  rms = (sqSum / total) ** 0.5
  
  return rms

def stringToList(listString):
  
  """
  Converts a str([]) back to []
  NOTE: still have to convert back to integer/float!
  """
  
  listString = listString[1:-1]
  
  newList = listString.split(',')
  
  return newList
  
  
def nameToSqlName (name):
  """
  Function that returns sql name from a python name.
  Changes anExampleString51 to AN_EXAMPLE_STRING_51
  """

  sqlName = ""
  previousChar = ''
  for char in name:
    if char.isupper():
      sqlChar = "_" + char.lower()
      sqlName = sqlName + sqlChar
    elif char.isdigit() and not previousChar.isdigit():
      sqlChar = "_" + char
      sqlName = sqlName + sqlChar
    else:
      sqlName = sqlName + char
    previousChar = char

  if sqlName[0] == "_":
    sqlName = sqlName[1:]
    
  return sqlName.upper()

def frange(start, end=None, inc=None):

  "A simple range function that accepts float increments. Taken from ASPN: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66472"

  if end == None:
    end = start + 0.0
    start = 0.0
            
  if inc == None:
    inc = 1.0

  L = []

  while 1:
    next = start + len(L) * inc
    if inc > 0 and next >= end:
      break
    elif inc < 0 and next <= end:
      break
    L.append(next)
                    
  return L
      
def makePowerSet(valueList):

  """
  Returns the powerset for a list.
  
  Use as: powersetlist = makePowerSet(inputList)
  
  For [1,2,3] or [1,3,2] both will return 
  [[], [1], [2], [1, 2], [3], [1, 3], [2, 3], [1, 2, 3]]
  
  """
  
  # copy (to avoid side effects) and sort
  values = list(valueList)
  values.sort()
  
  # calculate result
  result = [[]]
  for value in values:
    length = len(result)
    for ll in result[:length]:
      result.append(list(ll))
      result[-1].append(value)
  #
  return result


def semideepcopy(dd, doneDict=None):
  """ does a semi-deep copy of a nested dictionary, for copying mappings.
  Dictionaries are copied recursively, .
  Lists are copied, but not recursively.
  In either case a single copy is made from a single object 
  no matter how many times it appears.
  Keys and other values are passed unchanged
  """
  
  if doneDict is None:
    doneDict = {}
  
  key = id(dd)
  result = doneDict.get(key)
  if result is None:
    result = {}
    doneDict[key] = result
 
    for kk,val in dd.items():
 
      if type(val) == types.DictType:
        result[kk] = semideepcopy(val, doneDict)
 
      elif type(val) == types.ListType:
        key2 = id(val)
        newval = doneDict.get(key2)
        if newval is None:
          newval = val[:]
          doneDict[key2] = newval
 
        result[kk] = newval
 
      else:
        result[kk] = val
  #
  return result

def isWindowsOS():

  return sys.platform[:3].lower() == 'win'

def isMacOS():

  return sys.platform.lower() == 'darwin'

def useWheelMouse():
  # NOTE:ED - small function so that all versions use the same mouse type
  #         left-mouse    -> selection box
  #         right-mouse   -> menu
  #         middle-mouse  -> drag
  #         scroll-wheel  -> zoom
  return isWindowsOS() or isMacOS()

def OSButton(button):
  MACOSBUTTONBINDING = (None, 1, 3, 2)

  if not (0 < button < 4):
    raise ValueError('Button must be 1-3')

  return int(MACOSBUTTONBINDING[button] if isMacOS() else button)

def buttonClick(button):
  return '<Button-%i>' % OSButton(button)

def buttonPress(button):
  return '<ButtonPress-%i>' % OSButton(button)

def buttonRelease(button):
  return '<ButtonRelease-%i>' % OSButton(button)

def buttonMotion(button):
  return '<B%i-Motion>' % OSButton(button)
