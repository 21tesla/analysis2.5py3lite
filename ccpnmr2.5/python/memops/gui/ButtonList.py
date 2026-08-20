
"""
======================COPYRIGHT/LICENSE START==========================

ButtonList.py: <write function here>

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
import types, os
import Tkinter

from memops.gui.BasePopup import BasePopup
from memops.gui.Button import Button
from memops.gui.Frame import Frame
from memops.gui.MessageReporter import showInfo
from memops.gui.Spacer import Spacer
from memops.gui.Base import getPopup
from memops.gui.WebBrowser import WebBrowser

from memops.gui import Documentation as Doc

from memops.universal.Io import getTopDirectory

GFX_DIR = os.path.join(getTopDirectory(),'python','memops','gui','graphics')

STRING_TYPE   = type('')
UNICODE_TYPE  = type(u'')
HORIZONTAL = Tkinter.HORIZONTAL
VERTICAL = Tkinter.VERTICAL

def notYetImplemented(text = ''):

  if (text):
    text = '"' + text + '" '

  showInfo('Not Yet Implemented', text + 'not yet implemented')

class ButtonList(Frame):

  def __init__(self, parent, texts, commands, images=None,
               buttonBorderwidth=1, spacers=1, expands=True, 
               docKeys=None, tipTexts=None, font=None, direction=HORIZONTAL,
               sticky='ew', *args, **kw):

    assert direction in (HORIZONTAL, VERTICAL), 'direction = %s' % direction

    if (type(commands) is types.DictType):
      have_dict = True
      for key in commands.keys():
        assert key in texts, 'key = %s' % key
    else:
      have_dict = False
      assert len(texts) == len(commands), 'len(texts) = %s, len(commands) = %s' % (len(texts), len(commands))

    if not tipTexts:
      tipTexts = [None] * max([len(commands) or [], len(texts) or []])

    if not docKeys:
      docKeys = [None] * max([len(commands) or [], len(texts) or []])

    if not kw.has_key('bg'):
      kw['bg'] = 'grey82'
      
    Frame.__init__(self, parent, sticky=sticky, *args, **kw)

    if (direction == Tkinter.HORIZONTAL):
      self.grid_rowconfigure(0, weight=1)
    else:
      self.grid_columnconfigure(0, weight=1)

    if images is None:
      images = [None] * len(texts)

    self.font = font
    self.spacers = spacers
    r = 0
    c = 0

    if expands:
      self.spacers = spacers = 0
      for i in range(len(texts)):
        if texts[i]:
          if (direction == Tkinter.HORIZONTAL):
            self.grid_columnconfigure(i, weight=1)
          else:
            self.grid_rowconfigure(i, weight=1)
          
    else:
      if spacers and (direction == Tkinter.HORIZONTAL):
        self.grid_columnconfigure(c, weight=1)
        self.grid_columnconfigure(len(texts)+1, weight=1)
        spacer1 = Spacer(self)
        spacer2 = Spacer(self)
        spacer1.grid(row=r, column=c, sticky=Tkinter.NSEW)
        spacer2.grid(row=r, column=len(texts)+1, sticky=Tkinter.NSEW)
        c += 1 
 
    self.buttons = []
    bg = kw['bg']

    for n, text in enumerate(texts):
      
      if n < len(images):
        image = images[n]
      else:
        image = None  
      
      if type(text) not in (STRING_TYPE, UNICODE_TYPE):
        text = str(text)
      
      if (have_dict):
        try:
          command = commands[text]
        except:
          command = lambda x=text: notYetImplemented(x)
      else:
        command = commands[n]
        if not command:
          command = lambda x=text: notYetImplemented(x)
      
      tipText = tipTexts[n]
      docKey = docKeys[n]
      button = Button(self, text=text, command=command, bg=bg, image=image,
                      borderwidth=buttonBorderwidth, docKey=docKey, tipText=tipText,
                      font=self.font)
      button.grid(row=r, column=c, sticky=Tkinter.NSEW)
      self.buttons.append(button)

      if direction == Tkinter.HORIZONTAL:
        c = c + 1
      else:
        r = r + 1

    # make buttons of equal width
    # for vertical buttons sticky=EW does this

    # below removed because it wastes space
    #if (direction == Tkinter.HORIZONTAL) and not spacers:
    #  # below so this happens after mapped
    #  self.after_idle(self.makeEqualWidths)

  def makeEqualWidths(self):

    self.update_idletasks()
    total =  int(self.winfo_width())

    width = total/float(len(self.buttons))
    for n in range(0,len(self.buttons)):
      self.grid_columnconfigure(n, minsize=max(30,width))

class UtilityButtonList(ButtonList):

  def __init__(self, parent,
               webBrowser=None, helpUrl=None, helpMsg=None,
               doClone=True, doHelp=True, doClose=True,
               cloneText=None, helpText=None, closeText=None,
               cloneCmd=None, helpCmd=None, closeCmd=None,
               cloneTip=None, helpTip=None, closeTip=None,
               *args, **kw):

    self.helpUrl = helpUrl
    self.helpMsg = helpMsg
    self.popup = getPopup(parent)
 
    cloneTip = cloneTip or Doc.cloneDoc.documentation
    helpTip = helpTip or Doc.helpDoc.documentation

    if not closeTip: # a bit of a hack, this
      if closeText == 'Cancel':
        closeTip = Doc.cancelDoc.documentation
      else:
        closeTip = Doc.closeDoc.documentation

    cloneKey = cloneText or Doc.cloneDoc.name
    helpKey = helpText or Doc.helpDoc.name
    closeKey = closeText or Doc.closeDoc.name

    if not isinstance(self.popup, BasePopup):
      self.popup = None
   
    if self.popup and not webBrowser:
    
      webBrowser = WebBrowser(self.popup.parent, popup=self.popup)
 
    self.webBrowser = webBrowser

    _commands = []
    _texts    = []
    _images   = []
    _docKeys = []
    _tipTexts = []
    
    _doActions = [(doClone, cloneCmd, self.clonePopup, cloneText, 'twoWindows.gif', cloneTip, cloneKey),
                  (doHelp, helpCmd, self.launchHelp, helpText, 'help.gif', helpTip, helpKey),
                  (doClose, closeCmd, self.closePopup, closeText, 'cancel.gif', closeTip, closeKey),]
    
    for doAction, userCmd, defaultCmd, text, imageFile, tipText, docKey in _doActions:
      if doAction:
        _commands.append(userCmd or defaultCmd)
        _docKeys.append(docKey)
        _tipTexts.append(tipText)
        
        if text:
          _texts.append(text)
          _images.append(None)
        else:
          icon = Tkinter.PhotoImage(file=os.path.join(GFX_DIR,imageFile))
          _texts.append('')
          _images.append(icon)
       
       
    # Needed for photoimage reference count/persistence
    self.icons = [img for img in _images if img]
  
    images = kw.get('images')
    if images:
      kw['images'] = images + _images
    else:
      kw['images'] = [None] * len(kw.get('texts', [])) + _images
  
    texts = kw.get('texts', [])
    commands = kw.get('commands', [])
    docKeys = kw.get('docKeys') or [None]*max([len(texts), len(commands)])
    tipTexts = kw.get('tipTexts') or [None]*max([len(texts), len(commands)])
    
    kw['docKeys'] = docKeys + _docKeys
    kw['tipTexts'] = tipTexts + _tipTexts
    kw['commands'] = commands + _commands
    kw['texts'] = texts + _texts
  
    ButtonList.__init__(self, parent, *args, **kw)

  def clonePopup(self):
     
    if self.popup:
      try:
        newPopup = self.popup.__class__(self.popup.parent)
        wh,x,y = self.popup.geometry().split('+')
        newPopup.geometry('%s+%d+%d' % (wh,int(x)+25,int(y)+25))
      except:
        pass
  
  def launchHelp(self):
  
    if self.helpUrl and self.webBrowser:
      self.webBrowser.open(self.helpUrl)
      
    elif self.popup:
      from memops.gui import HelpPopup
      HelpPopup.showHelpText(self.popup.parent, self.helpMsg, popup=self.popup)  
  
  def closePopup(self):
  
    if self.popup:
      self.popup.close()
      
    else:
      self.destroy()  
  
      
if __name__ == '__main__':

  import sys

  root = Tkinter.Tk()

  texts = ('quit1', 'quit2', 'tbd')
  commands = {}
  commands['quit1'] = commands['quit2'] = sys.exit

  button_list = ButtonList(root, texts, commands, direction=Tkinter.HORIZONTAL, tipTexts=texts)
  button_list.grid()

  texts = ('quit3', 'quit4', 'tbd')
  commands = {}
  commands['quit3'] = commands['quit4'] = sys.exit

  button_list = ButtonList(root, texts, commands, direction=Tkinter.VERTICAL, tipTexts=texts)
  button_list.grid()
  
  button_list = UtilityButtonList(root)
  button_list.grid()
  
  root.mainloop()
