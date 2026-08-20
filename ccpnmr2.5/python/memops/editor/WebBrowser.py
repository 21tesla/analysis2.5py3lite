
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
from memops.gui.WebBrowser import WebBrowser, WebBrowserPulldown, getDefaultBrowser

class ProjectWebBrowser(WebBrowser):

  def __init__(self, parent, name=None, popup=None, project=None):
  
    self.project = project
    if project and (not name):
      if hasattr(project,'application'):
        name = getProjectWebBrowser(project)
    
    WebBrowser.__init__(self, parent, name=name, popup=popup)
  
  def open(self, url):

    if self.project:
      name = getProjectWebBrowser(self.project)
      if name:
        self.name = name

    WebBrowser.open(self, url)
    
class ProjectWebBrowserPulldown(WebBrowserPulldown):

  def __init__(self, parent, project=None, **kw):

    self.project = project
    browser = None

    if project:
      browser = getProjectWebBrowser(project)
            
    WebBrowserPulldown.__init__(self, parent, browser=browser, **kw)

  def setWebBrowser(self, index, name):
    # overwrites superclass function

    if name != self.browser:
      self.browser = name
      if self.project:
        setProjectWebBrowser(self.project, name)

def getProjectWebBrowser(project):

  browser = getDefaultBrowser()
  if hasattr(project, 'application'):
    application = project.application
    ### TBD v2
    ###browser = application.getValue(project, keyword='projectWebBrowser', defaultValue=browser)

  return browser

def setProjectWebBrowser(project, name):

  if hasattr(project, 'application'):
    application = project.application
    ### TBD v2
    ###application.setValue(project, keyword='projectWebBrowser', value=name)
  
