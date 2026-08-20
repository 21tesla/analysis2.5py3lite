
"""
======================COPYRIGHT/LICENSE START==========================

BasePopup.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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
import os

from memops.gui.Base import getPopup

from memops.universal.Io import splitPath
from memops.universal.Io import getTopDirectory

import memops.editor.BasePopup

from ccpnmr.analysis.Analysis import LOCAL_HELP_DOC_DIR

class BasePopup(memops.editor.BasePopup.BasePopup):

  def __init__(self, parent, *args, **kw):

    helpUrl = kw.get('help_url')
    if not helpUrl:
      helpUrl = determineHelpUrl(self.__class__)
      
    self.help_url = helpUrl

    project = kw.get('project')
    parentPopup = getPopup(parent)
    if not project and hasattr(parentPopup, 'project'):
      project = parentPopup.project

    self.project = project

    if project:
      self.nmrProject = project.currentNmrProject
      self.analysisProject = project.currentAnalysisProject
      self.analysisProfile = project.currentAnalysisProfile
    else:
      self.nmrProject = None
      self.analysisProject = None
      self.analysisProfile = None

    application = kw.get('application')
    if not application:
      if hasattr(self, 'application'):
        application = self.application
      elif hasattr(parent, 'application'):
        application = parent.application

    self.project = project
    self.application = application

    memops.editor.BasePopup.BasePopup.__init__(self, parent, *args, **kw)

def determineHelpUrl(clazz = None, name = None):

  assert clazz or name

  topDir = getTopDirectory()

  #dir = splitPath(__file__)[0]
  #if (dir):
  #  dir = dir + '/'

  if (not name):
    name = clazz.__name__
  help_url = 'file:' + topDir + LOCAL_HELP_DOC_DIR + 'popups/' + name + '.html'

  return help_url
