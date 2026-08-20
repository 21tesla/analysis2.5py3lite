
"""
======================COPYRIGHT/LICENSE START==========================

HelpPopup.py: <write function here>

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
import Tkinter


import memops.gui.BasePopup as BasePopup
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.Entry import Entry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from MessageReporter import showError, showInfo
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledHtml import ScrolledHtml
from memops.gui.ToggleArrow import ToggleArrow
import memops.gui.Util

class HelpFrame(Frame):

  def __init__(self, parent, width = 70, height = 20,
               xscroll = False, yscroll = True, *args, **kw):

    apply(Frame.__init__, (self, parent) + args, kw)

    self.grid_columnconfigure(1, weight=1)

    row = 0
    texts = ('back', 'forward', 'load')
    commands = (self.prevUrl, self.nextUrl, self.loadUrl)
    self.buttons = ButtonList(self, texts=texts, commands=commands)
    self.buttons.grid(row=row, column=0)
    self.url_entry = Entry(self, returnCallback=self.loadUrl)
    self.url_entry.grid(row=row, column=1, columnspan=2, sticky=Tkinter.EW)
    row = row + 1

    frame = Frame(self)
    frame.grid(row=row, column=0, columnspan=3, sticky=Tkinter.W)
    self.createFindFrame(frame)
    row = row + 1

    self.grid_rowconfigure(row, weight=1)
    self.help_text = ScrolledHtml(self, width=width, height=height,
                                  xscroll=xscroll, yscroll=yscroll,
                                  startUrlCallback = self.startOpenUrl,
                                  endUrlCallback = self.endOpenUrl,
                                  enterLinkCallback=self.setLinkLabel,
                                  leaveLinkCallback=self.setLinkLabel)
    self.help_text.grid(row=row, column=0, columnspan=3, sticky=Tkinter.NSEW)
    row = row + 1

    self.link_label = Label(self)
    self.link_label.grid(row=row, column=0, columnspan=2, sticky=Tkinter.W)
    button = memops.gui.Util.createDismissButton(self, dismiss_cmd=self.close)
    button.grid(row=row, column=2)

    self.urls = []
    self.current = -1
    self.updateUrlList = True

    self.setButtonState()

  def close(self):

    if (self.popup):
      if (self.popup.modal):
        self.popup.do_grab()

    popup = self.parent
    while (not hasattr(popup, 'top')):
      popup = popup.parent
    popup.close()

  def createFindFrame(self, master):

    frame = Frame(master)
    frame.grid(row=0, column=0, sticky=Tkinter.W)

    arrow = ToggleArrow(frame, callback = self.toggleFindFrame)
    arrow.grid(row=0, column=0)
    button = Button(frame, text='find:', command=self.findPhrase)
    button.grid(row=0, column=1)
    self.find_entry = Entry(frame, width=20, returnCallback=self.findPhrase)
    self.find_entry.grid(row=0, column=2)

    self.find_frame = frame = Frame(master)

    entries = ('search forwards', 'search backwards')
    self.direction_buttons = PulldownMenu(frame, entries=entries)
    self.direction_buttons.grid(row=0, column=0, sticky=Tkinter.W, padx=5)

    self.forwards_entries = entries = ('wrap search', 'stop at end')
    self.backwards_entries = ('wrap search', 'stop at beginning')
    self.wrap_buttons = PulldownMenu(frame, entries=entries)
    self.wrap_buttons.grid(row=0, column=1, sticky=Tkinter.W, padx=5)

    self.direction_buttons.callback=self.setWrapText

    entries = ('case insensitive', 'case sensitive')
    self.case_buttons = PulldownMenu(frame, entries=entries)
    self.case_buttons.grid(row=0, column=2, sticky=Tkinter.W, padx=5)

    entries = ('literal search', 'regex search')
    self.pattern_buttons = PulldownMenu(frame, entries=entries)
    self.pattern_buttons.grid(row=0, column=3, sticky=Tkinter.W, padx=5)

    self.countVar = Tkinter.IntVar()

  def setWrapText(self, entry_ind, entry):

    if (entry_ind == 0):
      entries = self.forwards_entries
    else:
      entries = self.backwards_entries

    self.wrap_buttons.replace(entries, selected_index=self.wrap_buttons.getSelectedIndex())

  def findPhrase(self, *extra):

    phrase = self.find_entry.get()
    if (not phrase):
      showError('No find phrase', 'Enter phrase in entry box.')
      return

    if (self.direction_buttons.getSelectedIndex() == 0):
      forwards = 1
      backwards = 0
    else:
      forwards = 0
      backwards = 1

    if (self.case_buttons.getSelectedIndex() == 0):
      nocase = 1
    else:
      nocase = 0

    if (self.pattern_buttons.getSelectedIndex() == 0):
      exact = 1
      regexp = 0
    else:
      exact = 0
      regexp = 1

    start = self.help_text.tag_nextrange(Tkinter.SEL, '1.0')
    if (start):
      start = start[0]
      if (forwards):
        start = '%s+1c' % start
      else:
        start = '%s-1c' % start
    else:
      start = Tkinter.CURRENT

    if (self.wrap_buttons.getSelectedIndex() == 0):
      match = self.help_text.search(phrase, start, count=self.countVar,
                forwards=forwards, backwards=backwards,
                nocase=nocase, exact=exact, regexp=regexp)
    elif (forwards):
      match = self.help_text.search(phrase, start, count=self.countVar,
                forwards=forwards, backwards=backwards,
                nocase=nocase, exact=exact, regexp=regexp, stopindex=Tkinter.END)
    else:
      match = self.help_text.search(phrase, start, count=self.countVar,
                forwards=forwards, backwards=backwards,
                nocase=nocase, exact=exact, regexp=regexp, stopindex='1.0')

    if (match):
      self.help_text.see(match)
      self.help_text.setSelection(match, '%s+%dc' % (match, self.countVar.get()))
      self.help_text.tag_config(Tkinter.SEL, background='gray70')
    else:
      showInfo('No match', 'No match found for phrase.')

  def toggleFindFrame(self, isClosed):

    if (isClosed):
      self.find_frame.grid_forget()
    else:
      self.find_frame.grid(row=1, column=0, sticky=Tkinter.W)

  def setButtonState(self):

    n = len(self.urls)

    if ((n >= 2) and (self.current > 0)):
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED
    self.buttons.buttons[0].config(state=state)

    if ((n >= 2) and (self.current < (n-1))):
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED
    self.buttons.buttons[1].config(state=state)

  def setLinkLabel(self, url = ''):

    self.link_label.set(url)

  def loadUrl(self, *extra):

    url = self.url_entry.get()
    if (url):
      if (url.find(':') == -1):
        if (url[:2] != '//'):
          url = '//' + url
        url = 'http:' + url
      self.showUrl(url, forceLoad = True)

  def newUrl(self, url):

    self.current = self.current + 1
    self.urls[self.current:] = [[url, 0.0]]
    self.setButtonState()
    #print 'newUrl', self.current, self.urls

  def prevUrl(self):

    if (self.current > 0): # True if called via button
      (url, offset) = self.urls[self.current-1]
      self.updateUrlList = False
      self.showUrl(url, offset=offset, allowModifyPath=False)
      self.updateUrlList = True
      self.current = self.current - 1
      self.setButtonState()
    #print 'prevUrl', self.current, self.urls
      
  def nextUrl(self):

    if (self.current < (len(self.urls) - 1)): # True if called via button
      (url, offset) = self.urls[self.current+1]
      self.updateUrlList = False
      self.showUrl(url, offset=offset, allowModifyPath=False)
      self.updateUrlList = True
      self.current = self.current + 1
      self.setButtonState()
    #print 'nextUrl', self.current, self.urls

  def showText(self, message, popup = None):

    self.popup = popup
    if (popup):
      self.parent.modal = popup.modal
    self.help_text.setState(Tkinter.NORMAL) # so that can add text
    self.help_text.clear() # clear existing text
    self.help_text.append(message)
    self.help_text.setState(Tkinter.DISABLED) # make text read-only
    self.open()

  def startOpenUrl(self, yview):

    if (self.urls):
      self.urls[self.current][1] = yview[0]

  def endOpenUrl(self, url):

    #print 'endOpenUrl', url, self.updateUrlList
    self.url_entry.set(url)
    if (self.updateUrlList):
      self.newUrl(url)
    self.setLinkLabel()
    self.open()

  def showUrl(self, url, offset=None, forceLoad = False,
              allowModifyPath = True, popup = None):

    self.popup = popup
    if (popup):
      self.parent.modal = popup.modal
    self.help_text.openUrl(url, forceLoad=forceLoad, allowModifyPath=allowModifyPath)
    if (offset is not None):
      self.help_text.yview(Tkinter.MOVETO, str(offset))

def createHelpPopup(parent):

  if (not hasattr(parent, 'help_popup') or not parent.help_popup):

    kw = { 'transient': True }
    help_popup = BasePopup.createPopup(parent, frameClass=HelpFrame,
                                       title='Help Module', popup_kw=kw)
    help_popup.close()
    parent.help_popup = help_popup

  return parent.help_popup

def showHelpText(parent, message, popup = None):

  help_popup = createHelpPopup(parent)
  help_popup.showText(message, popup=popup)

def showHelpUrl(parent, url, popup = None):

  help_popup = createHelpPopup(parent)
  help_popup.showUrl(url, forceLoad=True, popup=popup)

if __name__ == '__main__':

  n = 0

  def showHelp():

    global n, root
    n = n + 1
    showHelpText(root, 'sample help text ' + str(n))

  root = Tkinter.Tk()

  button = Button(root, text='show help', command = showHelp)
  button.pack(side=Tkinter.LEFT)

  button = memops.gui.Util.createQuitButton(root)
  button.pack(side=Tkinter.LEFT)

  root.mainloop()
