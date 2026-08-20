
"""
======================COPYRIGHT/LICENSE START==========================

Register.py: Part of the CcpNmr Analysis program

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
from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showInfo

from ccpnmr.analysis.popups.BasePopup import BasePopup

from ccpnmr.analysis.core.Register import updateRegister

class RegisterPopup(BasePopup):

  """
  **Register User with CCPN**

  The purpose of this dialog is to allow the user to register with CCPN.
  This is mainly to allow CCPN to keep track of the number of users,
  which is important for grant applications.  This information is saved
  both in a file on your computer but also in a private CCPN server
  database if your computer is connected to the internet.  (So unfortunately
  you will need to register once for each computer you use.)

  The required information: the user name, organisation and email address.
  The other information stored on the CCPN server: the version number, the
  first time and most recent time the server has been notifed.

  If you are registered, then the server is notified every time a project is
  opened or saved (if connected to the internet).

  Click "Register Now" to register, or "Register Later" not to register.
"""

  def __init__(self, parent, isModal=False, *args, **kw):

    title='Project : Register with CCPN'
    BasePopup.__init__(self, parent=parent, title=title, modal=isModal, **kw)

  def body(self, guiFrame):

    self.geometry('600x250+600+250')

    analysisProfile = self.analysisProfile
    userName = analysisProfile.userName
    userOrganisation = analysisProfile.userOrganisation
    userEmail = analysisProfile.userEmail

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(1, weight=1)

    explainText = 'To keep track of our users and which versions are being used\n' \
                      'we would like you to register your details with us.\n' \
                'Collating the number of users is important for grant applications.\n' \
                'Please do not use accents in any of the information\n' \
                'because Python does not handle it gracefully.'
    row = 0
    label = Label(guiFrame, text=explainText, grid=(row, 0), gridSpan=(1,2), sticky='ew')
    row += 1

    licenseAgreeText = 'I agree to abide by the rules of the CCPN licensing agreement.'

    self.agreeButton = CheckButton(guiFrame, licenseAgreeText, tipText=licenseAgreeText)
    self.agreeButton.grid(row=row, column=1, columnspan=1, sticky='nsew')

    row += 1

    self.entryWidgets = []
    for (text, value) in (('Name', userName), ('Organisation', userOrganisation), ('Email', userEmail)):
      label = Label(guiFrame, text=text+':', grid=(row,0))
      entry = Entry(guiFrame, text=value or '', grid=(row,1),
                    sticky='ew', tipText='Your ' + text)
      self.entryWidgets.append(entry)
      row += 1

    texts = [ 'Register Now', 'Read License', 'Register Later' ]
    tipTexts = [ 'Register now', 'Read License', 'Register later' ]
    commands = [ self.register, self.openLicense, self.close ]
    buttons = UtilityButtonList(guiFrame, helpUrl=self.help_url, grid=(row,0), gridSpan=(1,2),
                                commands=commands, texts=texts, tipTexts=tipTexts)
    self.buttons = buttons

  def openLicense(self):
    # licenseUrl = 'http://www.ccpn.ac.uk/license'
    licenseUrl = 'https://www.ccpn.ac.uk/industry/licencefees'
    try:
      import webbrowser
    except ImportError:
      showInfo('License Agreement', 'The CCPN License Agreement is available at %s' % licenseUrl,  parent=self)
      return

    wb = webbrowser.get()
    wb.open( licenseUrl )

  def register(self):

    if not self.agreeButton.get():
      showError('License Agreement', 'Please tick the box indicating that you agree to the CCPN licensing terms and conditions.', parent=self)
      return

    analysisProfile = self.analysisProfile
    attrs = ('userName', 'userOrganisation', 'userEmail')
    for n, entry in enumerate(self.entryWidgets):
      attr = attrs[n]
      value = entry.get().strip()
      if not value:
        showError('Blank value', '%s is blank, value required' % attr[4:], parent=self)
        return
      if attr == 'userEmail' and '@' not in value:
        showError('Illegal email', 'Email must have "@" in it', parent=self)
        return
      try:
        setattr(analysisProfile, attr, value)
      except Exception, e:
        showError('Attribute setting', 'Error setting %s: %s' % (attr[4:], e), parent=self)
        return
    analysisProfile.save()
    try:
      result = updateRegister(analysisProfile)
      showInfo('Registering', result, parent=self)
      self.close()
    except Exception, e:
      showError('Registering', str(e), parent=self)

