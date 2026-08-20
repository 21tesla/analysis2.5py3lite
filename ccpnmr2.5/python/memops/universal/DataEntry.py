"""
======================COPYRIGHT/LICENSE START==========================

DataEntry.py: Utility code for CCPN code generation framework

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
import string

# returns None (for 'cancel') if and only if initial_value = ''
# and user just responds to prompt with empty string
def askData(func_type, prompt, initial_value, min_value = None, max_value = None):

  assert func_type in ('string', 'integer', 'float')
  if (initial_value is not ''):
    s = ' [' + str(initial_value) + ']'
  else:
    s = ''
  prompt = prompt + s + ': '

  done = 0
  while (not done):

    done = 1 # true for most cases below, so set it here
    answer = raw_input(prompt)
    answer = answer.strip()
    if (initial_value is not '' and not answer):
      answer = initial_value

    if (answer):
      if (func_type == 'integer'):
        try:
          answer = int(answer)
        except:
          print 'Not an integer, please try again.'
          done = 0
      elif (func_type == 'float'):
        try:
          answer = float(answer)
        except:
          print 'Not a floating point value, please try again.'
          done = 0
    else:
      answer = None

  return answer

def askString(title, prompt, initial_value = '', *args, **kw):

  return askData('string', prompt, initial_value)

def askInteger(title, prompt, initial_value = '', min_value = None,
               max_value = None, *args, **kw):

  return askData('integer', prompt, initial_value, min_value, max_value)

def askFloat(title, prompt, initial_value = '', min_value = None,
               max_value = None, *args, **kw):

  return askData('float', prompt, initial_value, min_value, max_value)

class DataEntry:

  def askString(self, title, prompt, initial_value = '', *args, **kw):

    return askData('string', prompt, initial_value)

  def askInteger(self, title, prompt, initial_value = '', *args, **kw):

    return askData('integer', prompt, initial_value)

  def askFloat(self, title, prompt, initial_value = '', *args, **kw):

    return askData('float', prompt, initial_value)

  def askFile(self, title, prompt, initial_value = '', *args, **kw):

    return askData('string', prompt, initial_value)

  def askDir(self, title, prompt, initial_value = '', *args, **kw):

    return askData('string', prompt, initial_value)

dataEntry = DataEntry()

if (__name__ == '__main__'):

  print dataEntry.askString('ask string title', 'ask string prompt')
  print dataEntry.askInteger('ask integer title', 'ask integer prompt')
  print dataEntry.askFloat('ask float title', 'ask float prompt')
