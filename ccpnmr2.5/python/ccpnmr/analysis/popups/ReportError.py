
"""
======================COPYRIGHT/LICENSE START==========================

OpenMacro.py: Part of the CcpNmr Analysis program

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
import re
import os.path
import traceback
import os
import platform
import sys

import string, Tkinter, time
from memops.universal.Io import splitPath
from Tkinter import *

from memops.gui.Label               import Label
from memops.gui.ScrolledText        import ScrolledText
from memops.gui.ToggleLabel         import ToggleLabel
from memops.gui.ScrolledFrame       import ScrolledFrame
from memops.gui.LabelFrame          import LabelFrame
from memops.gui.ButtonList          import ButtonList
from memops.universal.Url           import *

from ccpnmr.analysis                    import Copyright
from ccpnmr.analysis.popups.BasePopup   import BasePopup                      

class ReportErrorPopup(BasePopup):

  def __init__(self, parent, traceBack, exc, val, tb, *args, **kw):

    self.parent = parent
    self.error = val[0]
    self.messageToReport = "".join(traceBack[0:])
    
    self.platform = sys.platform
    self.pythonVersion = sys.version
    
    version = Copyright.version
    versionData = [version.major, version.minor, version.release]
    self.analysisVesion = ' '.join([Copyright.suite, Copyright.program, '.'.join([str(x) for x in versionData])])
    
    if self.parent.analysisProfile.sendBugReports == 'maybe':
      BasePopup.__init__(self, parent=parent, title='Report Error', **kw)
    elif self.parent.analysisProfile.sendBugReports == 'yes':
      self.send()
    else:
      BasePopup.destroy(self)

  def body(self, guiParent):

    self.geometry('800x700')
    
    guiParent.expandGrid(2, 0)


    frame = LabelFrame(guiParent, text='Oops...', grid=(0,0))
    
    msg = 'An error was encountered when running CCPN software.'
    Label(frame, text=msg, grid=(0,0), sticky='w')
    msg = 'The CCPN developers would be grateful if this error was reported.'
    Label(frame, text=msg, grid=(1,0), sticky='w')
    msg = 'Please select one of the below options to automatically send a report or suppress further nagging.'
    Label(frame, text=msg, grid=(2,0), sticky='w')
        
    texts = ['Never Send Reports','Always Send Reports', 'Send This Time Only']
    commands = [self.neverSend, self.alwaysSend, self.send]
    
    self.chainButtons = ButtonList(frame, texts=texts, grid=(3,0),
                                   commands=commands, sticky='ew')
    
    frame1 = LabelFrame(guiParent, text='Data to send', grid=(1,0))
    
    frame1.grid_columnconfigure(0, weight=1)
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_rowconfigure(1, weight=1)
    frame1.grid_rowconfigure(2, weight=1)
    frame1.grid_rowconfigure(3, weight=1)
    frame1.grid_rowconfigure(4, weight=1)
    frame1.grid_rowconfigure(5, weight=1)

    uname = platform.uname()
    Label(frame1, text='Platform: '+self.platform+' '+' '.join(uname[:-1]), sticky = 'nw', grid=(0,0))
    Label(frame1, text='Python version: '+sys.version, sticky = 'nw', grid=(1,0))
    Label(frame1, text='Analysis version: '+self.analysisVesion, sticky = 'nw', grid=(2,0))
    
    #IntEntry(guiParent, grid=(2,1), width=5)
    Label(frame1, text='\nComment: ', sticky = 'sw', grid=(3,0))
    self.comment = ScrolledText(frame1, height=7, grid=(4, 0), xscroll=False, sticky='new')
    
    
    frame15 = LabelFrame(guiParent, text='Error Traceback', grid=(2,0))
    frame15.grid_columnconfigure(0, weight=1)
    frame15.grid_rowconfigure(0, weight=1)    
    sFrame = ScrolledFrame(frame15, sticky='nsew',  grid=(0,0), doExtraConfig=False)
    Label(sFrame.frame, text=self.messageToReport, justify='left', sticky = 'nsew', grid=(0,0))
    

  def neverSend(self):
    self.parent.analysisProfile.sendBugReports = 'no'
    self.parent.analysisProfile.save()
    self.close()
  
  def alwaysSend(self):
    self.send()
    self.parent.analysisProfile.sendBugReports = 'yes'
    self.parent.analysisProfile.save()
    
  def send(self):
    try:
      if self.parent.analysisProfile.sendBugReports == 'maybe':
        userComment = self.comment.text_area.getText()
      else:
        userComment = ""
      values = {'ErrorMessage' : self.error,
        'TraceBack' : self.messageToReport,
        'userName' : self.parent.analysisProfile.userName,
        'userOrganization' : self.parent.analysisProfile.userOrganisation,
        'userEmail' : self.parent.analysisProfile.userEmail,
        'platform' : self.platform,
        'analysisVersion' : self.analysisVesion,
        'pythonVersion' : self.pythonVersion,
        'comment' : userComment }

      print fetchUrl('http://www2.ccpn.ac.uk/cgi-bin/karolis/SubmitBug.py', values, timeout=3)
      try:
        self.close();
      except:
        print '"Always send" option'
        
    except:
      print "Report failed ", sys.exc_info()[0]
      try:
        self.close();
      except:
        print '"Always send" option'
  def close(self):
    BasePopup.close(self)
    self.destroy()
  
  def destroy(self):
    BasePopup.destroy(self)
