
"""
======================COPYRIGHT/LICENSE START==========================

WebBrowser.py: <write function here>

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
import tempfile
import webbrowser as wb
from memops.gui.MessageReporter import showWarning
from memops.gui.PulldownList import PulldownList

browserNames = ['firefox','netscape','mozilla','konqueror','kfm','mosaic',
                'grail','w3m','windows-default','internet-config']

class WebBrowser:

  def __init__(self, parent, name=None, popup=None):
  
    names = getBrowserList()
    if names and (not name):
      name = names[0]
    
    self.name   = name
    self.popup  = popup
    self.parent = parent
  
  def open(self, url):
    
    try:
      browser = wb.get(self.name)
      browser.open(url)
      
    except:
      from memops.gui.HelpPopup import showHelpUrl
      showHelpUrl(self.parent, url, popup=self.popup)
  
  def openHtml(self, htmlString):

    tmpfile = tempfile.NamedTemporaryFile()
    tmpfile.write(htmlString)
    tmpfile.flush()
    self.open('file://%s' % tmpfile.name)

class WebBrowserPulldown(PulldownList):

  def __init__(self, parent, browser=None, **kw):

    self.browserList = getBrowserList()

    if not browser:
      browser = getDefaultBrowser()

    if self.browserList:
      if (not browser) or (browser not in self.browserList):
        browser = self.browserList[0]

    self.browser = browser

    PulldownList.__init__(self, parent, **kw)

    if self.browserList:
      self.setup(self.browserList, self.browserList, self.browserList.index(self.browser))
    
    self.callback = self.setWebBrowser

  def setWebBrowser(self, name):

    if name != self.browser:
      self.browser = name

  def destroy(self):

    pass

def getBrowserList():

  browsers = []
  default  = getDefaultBrowser()
  if default:
    browsers = [default,]
  
  for name in browserNames:
    if name == default:
      continue
  
    try:
      wb.get(name)
      browsers.append(name)
    except:
      
      try:
        if wb._iscommand(name):
          wb.register(name, None, wb.Netscape(name))
          wb.get(name)
          browsers.append(name)
      except:
        continue

  return browsers

def getDefaultBrowser():

  try:
    br = wb.get()
  except:
    return
  
  if not hasattr(br, 'name'):
    # Max OS X
    return
  
  try:
    wb.get(br.name)
  except:
    wb.register(br.name, None, br)
  
  return br.name
