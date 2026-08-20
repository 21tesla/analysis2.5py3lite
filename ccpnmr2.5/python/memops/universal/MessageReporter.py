"""
======================COPYRIGHT/LICENSE START==========================

MessageReporter.py: Utility code for CCPN code generation framework

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
def showMessage(message_type, message):

  print message_type + ': ' + message

def showInputMessage(input_types, message):

  assert len(input_types) == 2

  print message + ':'
  prompt = 'Choose one of: ' + ', '.join(input_types) + ' [' + input_types[0] + ']: '

  # danger: n = 0, 1 backwards below because want yes/ok to be 1, no/cancel to be 0
  n = -1
  while (n < 0):
    answer = raw_input(prompt)
    if (not answer):
      n = 1
    else:
      answer = answer.lower().strip()
      k = len(answer)
      if (answer == input_types[0][:k]):
        n = 1
      elif (answer == input_types[1][:k]):
        n = 0

  #print 'showInputMessage: returning', n
  return n

class MessageReporter:

  def showError(self, title, message, *args, **kw):
 
    showMessage('Error', message)

  def showInfo(self, title, message, *args, **kw):
 
    showMessage('Information', message)

  def showWarning(self, title, message, *args, **kw):
 
    showMessage('Warning', message)

  def showOkCancel(self, title, message, *args, **kw):
 
    return showInputMessage(('ok', 'cancel'), message)

  def showYesNo(self, title, message, *args, **kw):
 
    return showInputMessage(('yes', 'no'), message)

messageReporter = MessageReporter()

def showError(title, message, *args, **kw):
 
  messageReporter.showError(title, message, *args, **kw)

def showInfo(title, message, *args, **kw):
 
  messageReporter.showInfo(title, message, *args, **kw)

def showWarning(title, message, *args, **kw):
 
  messageReporter.showWarning(title, message, *args, **kw)

def showOkCancel(title, message, *args, **kw):
 
  return messageReporter.showOkCancel(title, message, *args, **kw)

def showYesNo(title, message, *args, **kw):
 
  return messageReporter.showYesNo(title, message, *args, **kw)

if (__name__ == '__main__'):

  messageReporter.showError('error title', 'error message')
  messageReporter.showInfo('info title', 'info message')
  messageReporter.showWarning('warning title', 'warning message')
  print messageReporter.showOkCancel('ok cancel title', 'ok cancel message')
  print messageReporter.showYesNo('yes no title', 'yes no message')
