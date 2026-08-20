
"""
======================COPYRIGHT/LICENSE START==========================

Ticker.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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
import os, time, threading

import Tkinter

from memops.gui.Label import Label
from memops.gui.Canvas import Canvas

class Ticker(Label):

  def __init__(self, parent, message=None, nRepeats=1, length=40, **kw):
  
    self.font     = '-schumacher-clean-medium-r-normal--14-140-75-75-c-70-iso646.1991-irv',
    self.length   = length
    self.message  = message or ' '
    self.buffer   = (' ',1)
    self.display  = self.length * ' '
    self.nRepeats = nRepeats
       
    Label.__init__(self, parent=parent, font=self.font, text=self.display,  fg='#000060', **kw)
       
    self.offset = 0
    self.thread = None
  
  def newThread(self):
  
    return threading.Thread(name='Ticker',target=self.daemon,args=(),kwargs={})
    
  def start(self):
    
    if not self.thread:
      self.thread = self.newThread()
      self.thread.start()

  def stop(self):
    
    #print "thread done"
    #del self.thread
    self.thread = None
    self.offset = 0
    self.message = ' '
  
  def daemon(self, *args, **kwargs):
    
    while self.thread:
      self.tick()
      time.sleep(0.08)
        
  def tick(self):
  
    self.display = self.display[1:] + self.message[self.offset]
    self.set(self.display)
    
    self.offset += 1
    if self.offset >= len(self.message):
      self.offset = 0
      self.nRepeats -= 1
      if self.nRepeats < 0:
        # just in case setText has been called in the mean time
        if self.buffer[1]:
          self.message,self.nRepeats  = self.buffer
          self.buffer = (self.length  * ' ', 0)
        else:
          self.stop()
     
      # lastly fill with blanks
      if self.nRepeats == 0:
        self.message,self.nRepeats  = self.buffer
        self.buffer = (self.length  * ' ', 0)
      
  
  def setMessage(self, text=None, repeats=1):
  
    self.buffer  = (text or self.length  * ' ', repeats)
    if not self.thread:
      self.start()
  
  def destroy(self):
    
   if self.thread:
     self.stop()
     
   Label.destroy(self)
  
