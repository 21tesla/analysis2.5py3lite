
"""
======================COPYRIGHT/LICENSE START==========================

EditProfiles.py: Part of the CcpNmr Analysis program

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
# Popup to edit analysis profiles

# Find all that are available

# Edit mutable ones

# Copy from immutable

from ccpnmr.analysis.core import Util
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.frames.KeysymList import KeysymList, no_keysym
from ccpnmr.analysis.popups.OpenMacro import OpenMacroPopup
from ccpnmr.analysis.core.Util import runMacro, reloadMacro
from ccpnmr.analysis.core.PeakBasic     import refreshPeakAnnotations

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.DataEntry import askString
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel, showWarning, showYesNo
from memops.gui.MultiWidget import MultiWidget
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.WebBrowser import getBrowserList, getDefaultBrowser

from ccp.general.Command import Command

import memops.gui.Color as Color
from memops.general.Util  import copySubTree

# TBD
#
# Better RefExperiment name selection

# Default, immutable : user location : 

# Moved/repoint macros

# Mark & Ruler Functionality

# PulldownList categories index bug


DEFAULT_COLOR = Color.black.hex

molTypes = ['protein','DNA','RNA','carbohydrate']

bugReportsDict = {
  'yes': 'Always Send',
  'no': 'Never Send',
  'maybe': 'Always Ask',
}

def testEditProfiles(argServer):


  p = EditProfilesPopup(argServer.parent)
  p.open()
  

class EditProfilesPopup(BasePopup):
  """
  **Edit User-Specific Profiles and Macros**
  
  This popup window is used to control various settings that are stored as a
  user-specific profile. Typically this data is stored in the .ccpn/
  sub-directory of the users home directory and not in any particular CCPN
  project, although the information is saved whenever a project is. The same
  profile data will be used when any CCPN project is opened by the user
  as long as they have the same home directory.
  
  Several kinds of information are stored in the user profile and currently these
  can be grouped into four categories: general options for control and graphical
  display in Analysis; "Macros", which represent Python scripts to add extra
  functionality; customisable residue codes and customisable colour schemes.

  **Main Options**
  
  This tab controls the selection of the current user-specific profile for
  Analysis; one will always be present, but the user can make as many as
  required and swap between them using the pulldown menu. It should be noted that
  if the selected profile is changed the effect will not be apparent unless the 
  CCPN project is saved and Analysis is restarted.

  Various profile parameters are listed in the main table and the user can
  adjust these to suit personal preference by double clicking in the "value"
  column. The details of the settings are described in the last column. Perhaps
  the most notable setting in this panel is the "Graphics Handler", which
  controls what technology is used for the drawing of contours in  spectrum
  windows. Although the "Tk" setting  is perhaps more reliable overall the
  "OpenGL" option can be much faster for drawing spectra, especially if the
  computer has a good graphics card and hardware linked OpenGL driver software
  (typically from NVidia or ATI).

  **Macros**
  
  The Macros tab lists all of the Python scripts that are available to run as
  in-program commands. These will include number of pre-defined, inbuilt
  commands and any user-specific scripts that extend the functionality of
  Analysis. All of the inbuilt macro scripts are associated with a keyboard
  shortcut (editable via the "Shortcut" column). Indeed, listing these commonly
  used commands as macros provides the mechanism by which the editable Analysis
  keyboard shortcuts are handled in general.

  Analysis macro scripts are functions written the the Python programming
  language. They are normal Python functions, but they must be able to operate
  with one mandatory input argument "argServer". This input is a Python object
  of the ccpnmr.analysis.macros.ArgumentServer class and serves to link the
  Python function with the Analysis program, to give access to the CCPN
  project data and the Analysis graphical interface.

  All macros are available for execution (being run) via the main "Macros" table
  but they may also be entered in to the "Macro" sections of the main Analysis
  menu and the right-mouse spectrum window menu. To place a macro script in the
  menus the user needs to toggle the "In ... menu?" columns to "Yes". Any script
  that requires access to a particular position in a spectrum, recorded when the
  right-mouse menu opens, will only be able to use that point if the macro is in
  the spectrum window menu; the main menu and macro table have no concept of
  current position."

  New macros can be added to the current profile (and hence project) via the
  [Add Macro] function. This will open a file browser that the user can use to
  find a Python file, and then the Python function within that file to load as a
  macro script. See the `Open Macro`_ documentation for more information. Note
  that the script itself is not stored inside the CCPN project or profile. By
  adding a macro the user is merely locating the required Python function on
  disk; the actual code remains at its original file location. If a macro's
  Python code is moved to a different location on disk the script will no longer
  be executable by CCPN. Under such circumstances however, the old macro (and
  hence recorded location) can be removed and a fresh one put in its place. If
  the contents of a macros code have changed on disk since it was last executed
  (in the current Analysis session) then the user may reload the macro so that
  it is using the newest version. This mechanism provides a convenient way of
  developing and debugging CCPN Python scripts without having to restart
  Analysis each time.

  Note that macros with names that are identical to existing, inbuilt Python
  functions should be avoided, e.g. "math", "test", "copy" etc. In any case using
  more verbose but descriptive function names is recommended. In in doubt about
  whether a function name is in use the user can issue an "import" command at the
  Python prompt to check.
  
  **Residue Codes**
  
  This section allows the user to customise how residue (chemical compound) "Ccp"
  codes are displayed in a Analysis. For example if a residue code like
  "dglc-hex-1-5" (glucose) is too long a shorter alternative can be provided
  e.g. "Gluc". By default all profiles map the DNA and RNA codes to three-letter
  forms like "Ade", "Cyt", even though the actual underlying Ccp codes are
  single letters ("A", "C").

  To change how a particular kind of residue is represented the use selects the
  appropriate biopolymer type in the upper pulldown menu, then double clicks to
  edit the "User GUI Code" column of the table. Note that the existing textual
  peak annotations in the spectrum displays will not change until the [Update
  Peak Annotation] button is pressed; this can be slow if there are many peaks
  in the project.

  **Colour Schemes**

  This section allows the user to define new named colours and colour schemes.
  In this system a single colour like "red" is still stored as a colour scheme
  although it only has a single colour specification. The existing colour
  schemes  are listed in the upper table and they may be copied and deleted. 
  Clicking on a scheme in the upper table lists its component colours in the
  lower table. It should be noted that deleting a colour scheme that is still in
  use will not effect anything that uses those colours in Analysis, but that
  colour scheme name will not appear when setting any colours in the future.

  New colour schemes may be created via duplication or [Create]: to start this
  makes a simple scheme with three colours. Colours may be added to a scheme by
  clicking [Add] and then double clicking to set the "Colour" column. Colours
  may be changed and removed at any time, but such changes will not effect
  anything existing in Analysis; a scheme is merely used to set colours, and no
  active link is maintained.

  To assist with arranging multi-colour schemes, there are several functions 
  below the lower table. These allow the user to move the selected colour
  within  its scheme, i.e to the top, bottom, up one place and down one place.
  To create schemes with a gradient of colours the [Smooth] function so that
  there is a smooth transition between the first and last colours of a selection.
  Selection multiple colours in a scheme can be done using left mouse-click with 
  the <Shift> key.

  .. _`Open Macro`: OpenMacroPopup.html

  """

  def __init__(self, parent, *args, **kw):
  
    self.profile = None
    self.macro = None
    self.expProfile = None
    self.molType  = 'protein'
    self.resProfile = None
    self.rowObj = None
    self.openMacroPopup = None
    self.colorScheme = None
    self.colorIndex = 0
    self.color = None
    
    BasePopup.__init__(self, parent=parent, title='Project : Preferences : User Options', **kw)

  def body(self, guiFrame):

    self.geometry('700x600')

    guiFrame.expandGrid(0,0)
    
    tipTexts = ['A table of the main program-wide user options affecting behaviour and appearance',
                'A table of all of the Python macro scripts known to the project',
                'A table of user-specific mappings of residue codes, e.g. DNA "A" -> "Ade"',
                'A system to edit and create the colour schemes used at various points in Analysis']
    #'A system to set the default contour and peak colours for spectra according to experiment type',
    options = ['Main Options','Macros', #'Experiment Colours',
               'Residue Codes','Colour Schemes']
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              grid=(0,0), tipTexts=tipTexts)
    self.tabbedFrame = tabbedFrame
    
    #frameA, frameB, frameC, frameD, frameE = tabbedFrame.frames
    frameA, frameB, frameD, frameE = tabbedFrame.frames

    # Profiles
    self.pulldownList = PulldownList(self, callback=self.setValue)
    
    frameA.expandGrid(1,1)
    
    label = Label(frameA, text='Current Analysis Profile: ', grid=(0,0))
    
    tipText = 'Selects which set of profile settings to use'
    self.profilePulldown = PulldownList(frameA, callback=self.selectProfile,
                                        grid=(0,1), tipText=tipText)
    
    tipTexts = ['Copy the currently displayed profile settings into a differently named set',
                'Make an entirely new group of profile settings, with default values',
                'Delete the currently viewed set of profiles, if there is more than one available']
    texts = ['Copy','New','Delete']
    commands = [self.copyProfile, self.newProfile, self.deleteProfile]
    self.profileButtons = ButtonList(frameA, texts=texts, commands=commands, 
                                     grid=(0,2), tipTexts=tipTexts)

    tipTexts = ['The short name of the parameter that the user may set a value for',
                'The current value associated with a parameter in the selected profile',
                'A description of what the setting is for']
    editWidgets      = [ None, True, None ]
    editGetCallbacks = [ None, self.getValue, None  ]
    editSetCallbacks = [ None, self.setValue, None  ] 
    headingList = ['Parameter','Value','Description']
    self.profileMatrix = ScrolledMatrix(frameA, headingList=headingList,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks, 
                                        editWidgets=editWidgets,
                                        multiSelect=False, highlightType=0,
                                        callback=self.selectRowObj,
                                        grid=(1,0), gridSpan=(1,3), tipTexts=tipTexts)
                                            
    # Macros
    
    # Sub-tabs? Available Macros, File Details, Add Macro
    frameB.grid_columnconfigure(0, weight=1)
    frameB.grid_rowconfigure(0, weight=1)
     
    tipTexts = ['The serial number of the macro, within the current profile',
                'A short textual name for the macro script which will appear in menus',
                'The keyboard shortcut used to invoke this macro function',
                'Whether the macro is visible and callable from the "Macro" section of the main menu',
                'Whether the macro is visible and callable from the right-mouse-button spectrum window menu; required for macros that use spectru, coordinates etc.',
                'A priority based raking of macros to enable use-configurable sorting and grouping',
                'The name of the Python module in which the macro (a Python function) resides',
                'The name of the Python function, as found in its module, which the macro specification loads and calls',
                'The location on disk of the Python module that contains the macro function']
    headingList = ('#','Name','Shortcut','In main\nmenu?',
                   'In mouse\nmenu?','Priority','Module','Function','Path')

    self.editMacroNameEntry     = Entry(self, text='', returnCallback=self.setMacroName, width=20)
    self.editMacroPriorityEntry = IntEntry(self, text='', returnCallback=self.setMacroName, width=5)
    self.editMacroShortcutMenu  = KeysymList(self, callback=self.setMacroShortcut)
    
    editWidgets      = [None,
                        self.editMacroNameEntry,
                        self.editMacroShortcutMenu,
                        None, None,
                        self.editMacroPriorityEntry,
                        None, None, None]
    editGetCallbacks = [None,
                        self.getMacroName,
                        self.getMacroShortcut,
                        self.toggleMacroInMenu,
                        self.toggleMacroInMouseMenu,
                        self.getMacroPriority,
                        None, None, None]
    editSetCallbacks = [None,
                        self.setMacroName,
                        self.setMacroShortcut,
                        None, None,
                        self.setPriority,
                        None, None, None]
    
    self.macroMatrix = ScrolledMatrix(frameB,
                                      headingList=headingList,
                                      editWidgets=editWidgets,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks,
                                      callback=self.selectMacroCell,
                                      deleteFunc=self.removeMacro,
                                      grid=(0,0), gridSpan=(1,2), tipTexts=tipTexts)
 
    guiFrame.expandGrid(0,0)

    tipTexts = ['Reload the selected macro function; attempts to re-import the Python module',
                'Execute the selected macro',
                'Removes the record of the macro function from the current profile; any Python code on disk is preserved',
                'Add a new macro function to the current profile by finding a Python module and file on disk']
    texts    = ['Reload','  Run  ','Remove Macro','Add Macro']
    commands = [self.reloadSelectedMacro,self.runSelectedMacro,
                self.removeMacro,self.addMacro]
    self.bottomMacroButtons = ButtonList(frameB, texts=texts, tipTexts=tipTexts,
                                         commands=commands, grid=(1,0))

    
    # Exp Colours
    
    #frameC.expandGrid(0,0)
    #self.refExpNameSelect = MultiWidget(self, Entry, callback=self.setRefExpNames, 
    #                                   minRows=1, useImages=False)
 
    #tipTexts = ['A descriptive name for the group of experiments which are to be similarly colored e.g. "HSQC" or "NOESY"',
    #            'A list of the full CCPN experiment type definitions within this group of experiments',
    #            'The number of separate colour profiles which will be applied in a cyclic manner to experiments of the class']
    #headingList = ['Class','Experiment\nTypes',
    #               'Num Colour\nProfiles']
    #editWidgets      = [None, self.refExpNameSelect, None]
    #editGetCallbacks = [None, self.getRefExpNames, None]
    #editSetCallbacks = [None, self.setRefExpNames, None]
    #
    #self.expProfileMatrix = ScrolledMatrix(frameC,
    #                                       headingList=headingList,
    #                                       editWidgets=editWidgets,
    #                                       editSetCallbacks=editSetCallbacks,
    #                                       editGetCallbacks=editGetCallbacks,
    #                                       callback=self.selectExpProfile,
    #                                       grid=(0,0), tipTexts=tipTexts)
    
    #tipTexts = ['Add a new, named experiment group for which a default colour profile will be applied',
    #            'delete the selected experiment profile grouping']
    #texts = ['Add Class','Delete Class']
    #commands = [self.addExpProfile, self.removeExpProfile]
    #self.expProfileButtons = ButtonList(frameC, texts=texts, commands=commands,
    #                                    grid=(1,0), tipTexts=tipTexts)

    #self.posSchemePulldown = PulldownList(self, callback=self.setPosScheme)
    #self.negSchemePulldown = PulldownList(self, callback=self.setNegScheme)
    #self.symColorPulldown = PulldownList(self, callback=self.setSymColor)
    #self.annColorPulldown = PulldownList(self, callback=self.setAnnColor)

    
    #tipTexts = ['The number of the colour specification within the selected experiment group; occasions will be used in order for spectra that match the experiment type',
    #            'The positive spectrum contour colour for this specification, to be used for experiments of the stated type in a cyclic manner',
    #            'The negative spectrum contour colour for this specification, to be used for experiments of the stated type in a cyclic manner',
    #            'The colour of the peak cross/symbol for this specification, to be used for experiments of the stated type in a cyclic manner',
    #            'The colour of th textual peak annotations for this specification, to be used for experiments of the stated type in a cyclic manner']
    #headingList = ['Occasion','Pos Contours','Neg Contours','Peak Symbol','Peak Text']
    #editWidgets      = [None, self.posSchemePulldown,
    #                    self.negSchemePulldown,
    #                    self.symColorPulldown,
    #                    self.annColorPulldown]
    #editGetCallbacks = [None, self.getPosScheme,
    #                    self.getNegScheme,
    #                    self.getSymColor,
    #                    self.getAnnColor]
    #editSetCallbacks = [None, self.setPosScheme,
    #                    self.setNegScheme,
    #                    self.setSymColor,
    #                    self.setAnnColor]
    
    #self.expColorMatrix = ScrolledMatrix(frameC, highlightType=0,
    #                                     headingList=headingList,
    #                                     editWidgets=editWidgets,
    #                                     editSetCallbacks=editSetCallbacks,
    #                                     editGetCallbacks=editGetCallbacks,
    #                                     callback=self.selectExpColorIndex,
    #                                     grid=(3,0), tipTexts=tipTexts)
    
    #tipTexts = ['Add a new colour specification for this experiment profile group; each specification will be used in turn for members of the group',
    #            'Delete the selected colour specification, within the selected experiment profile group']
    #texts = ['Add Specification','Delete Specification']
    #commands = [self.addExpColor, self.removeExpColor]
    #self.expColorButtons = ButtonList(frameC, texts=texts, commands=commands,
    #                                  grid=(4,0), tipTexts=tipTexts)
    
    
    # Res Codes

    
    row = 0

    label = Label(frameD, text='Molecule Type:', grid=(row,0))
    
    tipText = 'Selects which kind of bio-polymer residues to display codes for'
    self.molTypePulldown = PulldownList(frameD, texts=molTypes, callback=self.changeMolType,
                                        grid=(row,1), tipText=tipText)

    row +=1
    frameD.expandGrid(row, 1)
    
    self.userCodeEntry = Entry(self,text='', returnCallback=self.setResidueCode, width=6)
    tipTexts = ['The official CCPN code for the residue ',
                'A user-configurable alternative residue code, if required']
    headingList      = ['Internal Ccp Code','User GUI Code']
    editWidgets      = [None, self.userCodeEntry]
    editGetCallbacks = [None, self.getResidueCode]
    editSetCallbacks = [None, self.setResidueCode]
    self.resProfileMatrix = ScrolledMatrix(frameD, headingList=headingList,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         callback=self.selectResProfile,
                                         grid=(row,0), gridSpan=(1,2), tipTexts=tipTexts)
                                        

    row +=1
    tipTexts = ['Manually trigger a complete update of all on-spectrum peak annotations in light of any residue code changes',]
    commands = [self.refreshPeakAnnotations,]
    texts    = ['Update Peak Annotations',]
    buttons = ButtonList(frameD, texts=texts, commands=commands, grid=(row,0),
                         gridSpan=(1,2), tipTexts=tipTexts)

    # Colour Schemes

    frameE.grid_rowconfigure(2, weight=1)
    frameE.grid_columnconfigure(0, weight=1)
    
    tipTexts = ['Row number',
                'Name of the colour scheme, for selection within menus etc.',
                'Number of colours in the scheme; single colours are a scheme of one']
    headings = ('#', 'Scheme Name', 'Number of Colours')
    self.colorSchemeMatrix = ScrolledMatrix(frameE, headingList=headings,
                                      multiSelect=True, callback=self.selectScheme,
                                      deleteFunc=self.deleteScheme,
                                      grid=(0,0), tipTexts=tipTexts)
 
    tipTexts = ['Create a new color scheme; the user is prompted for a name',
                'Create a new color scheme by copying the selected scheme; the user is prompted for a name for the copy',
                'Delete the currently selected colour scheme; does not effect anything that used the scheme, but the scheme will not appear as an option']
    texts = [ 'Create', 'Copy', 'Delete' ]
    commands = [ self.createScheme, self.copyScheme, self.deleteScheme ]
    self.colorSchemeButtons = ButtonList(frameE, texts=texts, grid=(1,0),
                                         commands=commands, tipTexts=tipTexts)
                          
    labelFrame = LabelFrame(frameE, text='Scheme Colours:', grid=(2,0))
    labelFrame.expandGrid(0,0)

    tipTexts = ['The number of the colour within the scheme',
                'The hexadecimal red-green-blue color name of the form "#RRGGBB", selectable using pre-defined colour names',
                'The fractional red component of the colour, rounded to give the nearest whole number in hexadecimal representation',
                'The fractional green component of the colour, rounded to give the nearest whole number in hexadecimal representation',
                'The fractional blue component of the colour, rounded to give the nearest whole number in hexadecimal representation',
                'The hue value of the colour, rounded to give the nearest whole number in hexadecimal representation',
                'The saturation level of the colour, rounded to give the nearest whole number in hexadecimal representation',
                'The brightness value of the colour, rounded to give the nearest whole number in hexadecimal representation']
    headings = ('#','Colour','R','G','B','H','S','V')
    self.editColorWidget = PulldownList(self, callback=self.setColor)
    self.editRedEntry = FloatEntry(self, text='', returnCallback=self.setRed, width=6)
    self.editGreenEntry = FloatEntry(self, text='', returnCallback=self.setGreen, width=6)
    self.editBlueEntry = FloatEntry(self, text='', returnCallback=self.setBlue, width=6)
    self.editHueEntry = FloatEntry(self, text='', returnCallback=self.setHue, width=6)
    self.editSaturationEntry = FloatEntry(self, text='', returnCallback=self.setSaturation, width=6)
    self.editBrightnessEntry = FloatEntry(self, text='', returnCallback=self.setBrightness, width=6)

    editWidgets      = [ None, self.editColorWidget,
                         self.editRedEntry,
                         self.editGreenEntry,
                         self.editBlueEntry,
                         self.editHueEntry,
                         self.editSaturationEntry,
                         self.editBrightnessEntry,
                       ]
                        
    editGetCallbacks = [ None, self.getColor,
                         self.getRed,
                         self.getGreen,
                         self.getBlue,
                         self.getHue,
                         self.getSaturation,
                         self.getBrightness,
                        ]
    editSetCallbacks = [ None, self.setColor,
                         self.setRed,
                         self.setGreen,
                         self.setBlue,
                         self.setHue,
                         self.setSaturation,
                         self.setBrightness,
                        ]
                        
    self.colorTable = ScrolledMatrix(labelFrame,
                                     headingList=headings,
                                     callback=self.selectColor,
                                     multiSelect=True, highlightType=0,
                                     editWidgets=editWidgets,
                                     editGetCallbacks=editGetCallbacks,
                                     editSetCallbacks=editSetCallbacks,
                                     deleteFunc=self.removeColor,
                                     grid=(0,0), tipTexts=tipTexts)
                                     
    
    tipTexts = ['Add a new colour element to the current scheme, defaults to black but can be changed',
                'Delete the selected colour element from the current scheme',
                'Make a smooth colour transition between the first and last colour elements selected by the user in the table',
                'Reverse the order of the colour elements in the scheme',
                'Move the selects colour element to be first in the scheme',
                'Move the selects colour element to be last in the scheme',
                'Move the selects colour element one place earlier in the scheme',
                'Move the selects colour element one place later in the scheme']
    texts = [ 'Add', 'Delete', 'Smooth', 'Reverse',
              'Top', 'Bottom', 'Up', 'Down',]
    commands = [ self.addColor, self.removeColor,
                 self.smoothColors, self.reverseColors,
                 self.colorToTop,  self.colorToBottom,
                 self.colorUp,  self.colorDown]
    self.colorButtons = ButtonList(labelFrame, texts=texts, commands=commands,
                                   grid=(1,0), tipTexts=tipTexts)

    # Main 
    
    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame,
                                           helpUrl=self.help_url)
    self.bottomButtons.grid(sticky='e')    
    self.administerNotifiers(self.registerNotify)
    
    self.updateAfter()
    self.updateButtons()
    
  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete','setBgColor','setFgColor',
                 'setFont, setGraphicsHandler','setWebBrowser',
                 'setUseGlobalShortcuts','setUseCrosshair',
                 'setTwoCharShortcuts','setTransientWindows',
                 'setTransientDialogs','setPanView',
                 'setRulersColor','setMarksColor','setSendBugReports'):
      notifyFunc(self.updateAfter, 'ccpnmr.AnalysisProfile.AnalysisProfile', func)
    
    notifyFunc(self.updateAnalyisProj, 'ccpnmr.Analysis.AnalysisProject', 'setContourToUnaliased')

    for func in ('__init__', 'delete', 'setName', 'setOrdering', 'setShortcut'):
      notifyFunc(self.updateMacrosAfter, 'ccpnmr.AnalysisProfile.Macro', func)
    
    for func in ('delete','__init__','setGuiName'):
      notifyFunc(self.updateResCodes, 'ccpnmr.AnalysisProfile.ResidueProfile', func)

    for func in ('__init__', 'delete', 'addColor', 'removeColor', 'setColors'):
      notifyFunc(self.updateColorsAfter, 'ccpnmr.AnalysisProfile.ColorScheme', func)
    
    #for func in ('__init__', 'delete','setPosColorSchemes','setNegColorSchemes',
    #             'setPeakSymbolColors','setPeakTextColors',
    #             'setRefExpNames','addRefExpName','removeRefExpName'):
    #  notifyFunc(self.updateExpProfiles, 'ccpnmr.AnalysisProfile.RefExpProfile', func)
      
    #for func in ('__init__', 'delete','setPosColorSchemes','setNegColorSchemes',
    #             'setPeakSymbolColors','setPeakTextColors'):
    #  notifyFunc(self.updateExpColors, 'ccpnmr.AnalysisProfile.RefExpProfile', func)
    
  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)

  def destroy(self):
 
    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)
 
  def updateButtons(self):
    
    buttons = self.colorSchemeButtons.buttons
      
    if self.colorScheme:
      buttons[1].enable()
      buttons[2].enable()
    
    else:
      buttons[1].disable()
      buttons[2].disable()
  
  def copyProfile(self, profile=None):
  
    if not profile:
      profile = self.profile
  
    if profile:
      name = None
      
      while not name:
        name = askString('Query','Enter profile name',parent=self)
        if name is None:
          return        
        
        name = name.strip()
        
        if name:
          if self.project.findFirstAnalysisProfile(name=name):
            showWarning('Retry','Name already in use',parent=self)
            name = None
          elif ' ' in name:
            showWarning('Retry','Name cannot contain spaces',parent=self)
            name = None
  
      profile = copySubTree(profile, self.project, 
                            topObjectParameters={'name':name},
                            objectMap={})
      
      self.selectProfile(profile)
  
  def newProfile(self):
  
    name = 'default'
    
    i = 1
    while self.project.findFirstAnalysisProfile(name=name):
      name = 'default%d' % i  
      i += 1
    
    profile = self.project.newAnalysisProfile(name=name)
    self.selectProfile(profile)

  def deleteProfile(self):
  
    profiles = self.project.analysisProfiles
    
    if len(profiles) == 1:
      msg = 'Cannot delete the only Analysis settings profile'
      showWarning('Warning', msg, parent=self)
      return
      
    msg = 'Really remove profile "%s" and all its macro entries, ' % self.profile.name
    msg += 'residue names and colour schemes?'
    if not showYesNo('Confirm', msg, parent=self):
      return
  
    if self.profile:
      self.profile.delete()
    
    profile = self.project.findFirstAnalysisProfile()
    if not profile:
      profile = self.project.newAnalysisProfile(name='default')
      
    self.selectProfile(profile)
    
      
  def updateProfiles(self):
  
    index = 0
    names = []
    profiles = self.project.sortedAnalysisProfiles()
    profile = self.profile
        
    if profiles:
      if profile not in profiles:
        profile = profiles[0]
        
      names = [p.name for p in profiles]
      index = profiles.index(profile)  
    
    else:
      profile = None
    
    if profile is not self.profile:
      self.selectProfile(profile)
    
    self.profilePulldown.setup(names, profiles, index)

  def updateAnalyisProj(self, analysisProj=None):
  
      self.after_idle(self.updateProfile)

  def updateAfter(self, profile=None):
  
    self.updateProfiles() # Pulldown
  
    if not profile or (profile is self.profile):
      self.after_idle(self.updateProfile)
      self.after_idle(self.updateResCodes)
      self.updateMacrosAfter()
      self.updateColorsAfter()
      #self.after_idle(self.updateExpProfiles)

  
  def selectProfile(self, profile):
  
    if profile is not self.profile:
      self.profile = profile
      self.macro = None
      self.expProfile = None
      self.resProfile = None
      self.rowObj = None
      self.colorScheme = None
      self.colorIndex = 0
      self.color = None
      self.expProfile = None
      self.expColorIndex = 0 

      self.updateProfiles()
      self.updateProfile()
      self.updateMacros()
      self.updateResCodes()
      self.updateColorsAfter()
      #self.updateExpProfiles()

  def updateProfile(self):
  
    textMatrix  = []
    objectList  = []
    colorMatrix = []
    
    nullColors = [None] * 3
      
    profile = self.profile
    if profile:
      if self.project.currentAnalysisProfile is not profile:
        self.project.currentAnalysisProfile = profile
     
      textMatrix.append(['Name', profile.name,'Name of the currently active profile'])
      objectList.append([profile,'name', None, None, None])
      colorMatrix.append(nullColors)
     
      textMatrix.append(['Window Background', profile.bgColor,'Default colour for the background of spectrum windows'])
      objectList.append([profile,'bgColor', self.pulldownList, self.getWindowBg, self.setPulldownMenu])
      colorMatrix.append([None, profile.bgColor, None])

      textMatrix.append(['Window Foreground', profile.fgColor,'Default colour for the foreground of spectrum windows'])
      objectList.append([profile,'fgColor', self.pulldownList, self.getWindowFg, self.setPulldownMenu])
      colorMatrix.append([None, profile.fgColor, None])
     
      textMatrix.append(['Main Font', profile.font,'The typeface for the Analysis menus and popup windows'])
      objectList.append([profile,'font', self.pulldownList, self.getFont, self.setPulldownMenu])
      colorMatrix.append(nullColors)
     
      scheme = None
      if profile.marksColor:
        scheme = profile.marksColor.name
      textMatrix.append(['Mark Colour Scheme', scheme,'Colour scheme for multidimensional marker lines'])
      objectList.append([profile,'marksColor', self.pulldownList, self.getProfileColorScheme, self.setPulldownMenu])
      colorMatrix.append(nullColors)
     
      scheme = None
      if profile.rulersColor:
        scheme = profile.rulersColor.name
      textMatrix.append(['Ruler Colour Scheme', scheme,'Colour scheme for single-dimension ruler lines'])
      objectList.append([profile,'rulersColor', self.pulldownList, self.getProfileColorScheme, self.setPulldownMenu])
      colorMatrix.append(nullColors)

      msg  = 'Whether to automatically tile the contour display to'
      msg += '\ncover any peaks outside the normal spectrum width.'
      textMatrix.append(['Contour To Unaliased?', self.analysisProject.contourToUnaliased, msg])
      objectList.append([self.analysisProject, 'contourToUnaliased', None, self.toggleBoolean, None])
      colorMatrix.append(nullColors)

      textMatrix.append(['Graphics Handler', profile.graphicsHandler,'Whether to use OpenGL or Tk for contour graphics'])
      objectList.append([profile,'graphicsHandler', self.pulldownList, self.getGraphicsHandler, self.setPulldownMenu])
      colorMatrix.append(nullColors)

      msg  = 'If true, navigation keys move the spectrum window view,'
      msg += '\notherwise the keys move the spectrum background'
      textMatrix.append(['Keys Pan Window View?', profile.panView,msg])
      objectList.append([profile,'panView', None, self.toggleBoolean, None])
      colorMatrix.append(nullColors)
     
      textMatrix.append(['Transient Dialogs?', profile.transientDialogs,'Transient dialogues have no toolbar representation'])
      objectList.append([profile,'transientDialogs', None, self.toggleBoolean, None])
      colorMatrix.append(nullColors)

      textMatrix.append(['Transient Windows?', profile.transientWindows,'Transient windows have no toolbar representation'])
      objectList.append([profile,'transientWindows', None, self.toggleBoolean, None])
      colorMatrix.append(nullColors)
     
      textMatrix.append(['Two Character Shortcuts?', profile.twoCharShortcuts,'Whether to use two-character shortcuts for macros etc.'])
      objectList.append([profile,'twoCharShortcuts', None, self.toggleBoolean, None])
      colorMatrix.append(nullColors)

      textMatrix.append(['Use Crosshairs?', profile.useCrosshair,'Whether to have crosshairs in the spectrum displays'])
      objectList.append([profile,'useCrosshair', None, self.toggleBoolean, None])
      colorMatrix.append(nullColors)
     
      textMatrix.append(['Use Global Shortcuts?', profile.useGlobalShortcuts,'When not global, keyboard shortcuts\nonly act in the active window'])
      objectList.append([profile,'useGlobalShortcuts', None, self.toggleBoolean, None])
      colorMatrix.append(nullColors)

      application = profile.root.application
      value = application.getValue(profile, 'yAxisDecimalsExpand', defaultValue=True)
      textMatrix.append(['Y axis decimals expand?', value,'Whether decimal places in y-axis tick labels expand as you zoom in\n(you must save and reopen project for this to take effect)'])
      objectList.append([profile,'yAxisDecimalsExpand', None, self.toggleAppBoolean, None])
      colorMatrix.append(nullColors)

      textMatrix.append(['Web Browser', profile.webBrowser,'The default web browser to use for help information etc.'])
      objectList.append([profile,'webBrowser', self.pulldownList, self.getWebBrowser, self.setPulldownMenu])
      colorMatrix.append(nullColors)
    
      textMatrix.append(['Send Bug Reports?', bugReportsDict[profile.sendBugReports],'Whether bug reports should be sent.'])
      objectList.append([profile,'sendBugReports', self.pulldownList, self.getBugReports, self.setPulldownMenu])
      colorMatrix.append(nullColors)
    
    self.profileMatrix.update(objectList=objectList,
                              textMatrix=textMatrix,
                              colorMatrix=colorMatrix)

  def getProfileColorScheme(self, widget, obj, attrName):
  
   
    value = getattr(obj, attrName)
    schemes = self.profile.sortedColorSchemes()
    
    names  = [s.name for s in schemes]
    colors = [list(s.colors) for s in schemes]
    
    if value in schemes:
      index = schemes.index(value)
    else:
      index = 0
  
    widget.setup(names,schemes,index,colors=colors) 
 

  def getWindowBg(self, widget, obj, attrName):
  
    value = getattr(obj, attrName) 
    names = [c.name for c in Color.standardColors]
    colors = [c.hex for c in Color.standardColors]
    
    if value in colors:
      index = colors.index(value)
    else:
      index = 0
  
    widget.setup(names,colors,index,colors=colors) 

  
  def getWindowFg(self, widget, obj, attrName):
    
    value = getattr(obj, attrName) 
    names = [c.name for c in Color.standardColors]
    colors = [c.hex for c in Color.standardColors]
    
    if value in colors:
      index = colors.index(value)
    else:
      index = 5
  
    widget.setup(names,colors,index,colors=colors)
     
  
  def getFont(self, widget, obj, attrName):
  
    value = getattr(obj, attrName) 
   
    names = []
    fonts= []
    categories = []
    
    for face in ('Courier','Helvetica','Lucida','System','Times'):
      for size in (8,10,12,14):

        # NOTE:ED added size 14 and the italic modes for testing
        for mod in ('', ' bold'):                                   #, ' italic', ' bold italic'):
          name = '%s %d%s' % (face,size,mod)
          names.append(name)
          fontName = (face, size)
          if mod:
              for mm in mod.split():
                fontName += (mm,)
          fonts.append(fontName)
          categories.append(face)
          
    if value in names:
      index = names.index(value)
    else:
      index = 0
    
    widget.setup(names,names,index,categories=categories,fonts=fonts)
  
  def getGraphicsHandler(self, widget, obj, attrName):
  
    value = getattr(obj, attrName) 
    names = Util.graphics_handlers
    
    if value in names:
      index = names.index(value)
    else:
      index = 0
  
    widget.setup(names,names,index) 
     
  def getWebBrowser(self, widget, obj, attrName):
  
    value = getattr(obj, attrName) 
    names = getBrowserList()
    
    if value not in names:
      value = getDefaultBrowser()
    
    if value in names:
      index = names.index(value)
    else:
      index = -1
  
    widget.setup(names,names,index) 

  def getBugReports(self, widget, obj, attrName):
  
    value = getattr(obj, attrName) 
    names = ['yes', 'no', 'maybe']
    
    if value not in names:
      value = names[-1]
    
    index = names.index(value)
    texts = [bugReportsDict[nn] for nn in names]
  
    widget.setup(texts,names,index) 

  #
  # Generic widget calls
  #
    
  def selectRowObj(self, obj, row, col):
  
    self.rowObj = obj
    
  def getValue(self, rowObj):
 
    # get correct object, widget and get/set functions
    obj, attrName, widget, getter, setter = rowObj

    # Bit of a hack because widgets are normally set at construction per column
    # Now setting it per row, and hence on-the-fly
    self.profileMatrix.editWidget = widget
        
    # get current value & setup widget
    if getter:
      getter(widget, obj, attrName)
  
  def setValue(self, null): 
  
    # get correct widget and get/set functions
    obj, attrName, widget, getter, setter = self.rowObj
        
    # set and check the appropriate parameter value from current edit widget
    
    # no setter for boolean toggles - the getter will have done all the toggling already
    # no setter for file selects - the file popup gets and sets the value and cannot be interrupted
    if setter:
      setter(widget, obj, attrName)
  
  def setPulldownMenu(self, widget, obj, attrName):

    value = widget.getObject()
    setattr(obj, attrName, value)
    self.updateProfile()
  
  def toggleBoolean(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    setattr(obj, attrName, not value)

  def toggleAppBoolean(self, widget, obj, attrName):

    application = obj.root.application
    value = application.getValue(obj, attrName, defaultValue=True)
    value = not value
    application.setValue(obj, attrName, value)
    self.updateAfter()

  # Macro functions
  
  def selectMacroCell(self, macro, row, col):
    
    if macro:
      self.macro = macro
      for button in self.bottomMacroButtons.buttons[0:3]:
        button.enable()
          
  def runSelectedMacro(self):

    if self.macro:
      runMacro(self.macro,self.parent.argumentServer)
  
  def reloadSelectedMacro(self):

    if self.macro:
      reloadMacro(self.macro,self.parent.argumentServer)
  
  def addMacro(self):
  
    if self.openMacroPopup:
      self.openMacroPopup.open()
    else:
      self.openMacroPopup = OpenMacroPopup(self)
    
  def removeMacro(self, *event):
  
    if self.macro:
      if showOkCancel('Remove macro','Are you sure?', parent=self):
        self.macro.delete()
        self.macro = None
       
  def toggleMacroInMenu(self, macro):
  
    if macro:
      macro.isInMenu = not macro.isInMenu
      self.parent.setMacroMenu()
      self.updateMacros()

  def toggleMacroInMouseMenu(self, macro):
  
    if macro:
      macro.isInMouseMenu = not macro.isInMouseMenu
      self.updateMacros()
    
  def setMacroName(self, event):
  
    text = self.editMacroNameEntry.get()
    if text:
      self.macro.name = text
      self.updateMacros()
    
  def getMacroName(self, macro):
  
    if macro:
      self.editMacroNameEntry.set(macro.name)
    
  def setPriority(self, event):
  
    try:
      integer = self.editMacroPriorityEntry.get()
    except:
      showWarning('Error', 'Invalid integer entry', parent=self)
      return
    
    if integer:
      self.macro.ordering = integer
      self.updateMacros()
    
  def getMacroPriority(self, macro):
  
    if macro:
      self.editMacroPriorityEntry.set(macro.ordering)

  def setMacroShortcut(self, *obj):

    shortcut= self.editMacroShortcutMenu.getObject()

    if self.macro.shortcut == shortcut:
      return

    if shortcut:
      # unset shortcut for macro(s) with this shortcut
      macros = self.profile.findAllMacros(shortcut=shortcut)
      for macro in macros:
        macro.shortcut = None

    self.macro.shortcut = shortcut

  def getMacroShortcut(self, macro):

    if macro:
      shortcut = macro.shortcut
      if not shortcut:
        shortcut = no_keysym
      self.editMacroShortcutMenu.setSelected(shortcut)

  def updateMacrosAfter(self, obj=None):
  
    if obj and (obj.analysisProfile is not self.profile):
      return
  
    self.after_idle(self.updateMacros)
  
  def updateMacros(self):
  
    if not self.macro:
      for button in self.bottomMacroButtons.buttons[0:3]:
        button.disable()
  
    textMatrix = []
    colorMatrix = []
    macros  = self.profile.sortedMacros()
    blankColors = [None] * 9
    
    if self.profile:
      for macro in macros:
 
        text = [ macro.serial,
                 macro.name,
                 macro.shortcut,
                 macro.isInMenu and 'Yes' or 'No',
                 macro.isInMouseMenu and 'Yes' or 'No',
                 macro.ordering,
                 macro.module,
                 macro.function,
                 macro.path]
        
        colors = blankColors[:]
        if macro.isInMenu:
          colors[3] = '#B0FFB0'
        if macro.isInMouseMenu:
          colors[4] = '#B0B0FF'
 
        textMatrix.append(text)
        colorMatrix.append(colors)

    self.macroMatrix.update(objectList=macros,colorMatrix=colorMatrix,
                            textMatrix=textMatrix)

  # Residue Profile Functions
  
  def refreshPeakAnnotations(self):
  
    refreshPeakAnnotations(self.project)

  def getResidueCode(self, resProfile):
  
    if resProfile.className == 'ResidueProfile':
      userCode = resProfile.guiName
    else:
      userCode = resProfile.ccpCode
    
    self.userCodeEntry.set(userCode)

  
  def setResidueCode(self, event=None):
  
    newCode = self.userCodeEntry.get().strip()

    currentGuiName = None
    if self.resProfile.className == 'ResidueProfile':
      currentGuiName = self.resProfile.guiName

    if newCode and (newCode != currentGuiName):
      existing = self.profile.findFirstResidueProfile(molType=self.molType, guiName=newCode)
      if existing:
        msg = 'User code %s already in use for %s' % (newCode, existing.ccpCode)
        showWarning('Failure', msg, parent=self)
      
      else:
        if self.resProfile.className != 'ResidueProfile':
          func = self.profile.newResidueProfile
          molType = self.molType
          ccpCode = self.resProfile.ccpCode # Actually a chemComp stand-in
          self.resProfile = func(molType=molType, ccpCode=ccpCode, guiName=newCode)
      
        else:
          self.resProfile.guiName = newCode

  def changeMolType(self, molType):
    
    if molType != self.molType:
      self.molType = molType
      self.updateResCodes()
     
  def selectResProfile(self, obj, row, col):

    self.resProfile = obj
    
  def updateResCodes(self, obj=None):
  
    textMatrix = []
    objectList = []
    
    if self.profile:
      resProfiles = self.getResProfiles()
    
      for ccpCode, obj in resProfiles:
        if obj.className != 'ResidueProfile':
          datum = [ccpCode, None]
        else:
          datum = [ccpCode, obj.guiName]
       
        textMatrix.append(datum)
        objectList.append(obj)
        
    self.resProfileMatrix.update(textMatrix=textMatrix,
                                 objectList=objectList)
      
  def getResProfiles(self):
  
    doneDict = {}
    resProfiles = []
    for resProfile in self.profile.residueProfiles:
      if resProfile.molType == self.molType:
        ccpCode = resProfile.ccpCode
        doneDict[ccpCode] = True
        resProfiles.append((ccpCode, resProfile))

    for chemComp in self.project.chemComps:
      if chemComp.molType == self.molType:
        ccpCode = chemComp.ccpCode
        if not doneDict.get(ccpCode):
          resProfiles.append((ccpCode, chemComp))
      
    resProfiles.sort()
    
    return resProfiles
      
  # Color scheme functions
  
      
  def updateColorsAfter(self, scheme=None):
  
    self.updateSchemeList()
    self.updateColorList()

  def updateSchemeList(self, *extra):

    profile = self.profile
    schemes = [s for s in profile.sortedColorSchemes()]
    textMatrix = []
    n = 0
    for scheme in schemes:
      n = n + 1
      text = [n,scheme.name,len(scheme.colors)]
      textMatrix.append(text)
 
    self.colorSchemeMatrix.update(objectList=schemes,
                             textMatrix=textMatrix)

    if schemes and not self.colorScheme:
      self.colorSchemeMatrix.selectObject(schemes[0])
    
  def selectScheme(self, scheme, row, col):

    self.colorScheme = scheme
    self.updateColorList()
    self.updateButtons()

  def copyScheme(self):

    if self.colorScheme:
      self.createScheme(colors=self.colorScheme.colors)

  def createScheme(self, colors=None):

    used = [s.name for s in Color.standardColors]
    name = None

    while not name:
      msg  = 'Enter name for new color scheme:'
      name = askString('Name Query', msg, parent=self)
      if name is None:
        return
      
      name.strip()
      
      if ' ' in name:
        msg = 'Name cannot contain spaces'
        showWarning('Retry', msg, parent=self)
        name = None
         
      elif name in used:
        name = None
        msg = 'Cannot use a standard colour name'
        showWarning('Retry', msg, parent=self)
      
      elif self.profile.findFirstColorScheme(name=name):
        name = None
        msg = 'Colour/scheme name already in use'
        showWarning('Retry', msg, parent=self)
      
    if name:
      if not colors:
        colors = [Color.red.hex, Color.green.hex, Color.blue.hex]

      scheme = self.profile.newColorScheme(name=name, colors=colors)
      self.colorSchemeMatrix.selectObject(scheme)
        

  def deleteScheme(self, *event):

    schemes = self.colorSchemeMatrix.currentObjects
    if len(schemes) == 1:
      msg = 'Delete color scheme "%s"?' % schemes[0].name
      if (showYesNo('Delete color scheme', msg, parent=self)):
        schemes[0].delete()
    
    elif len(schemes) > 1:
      names = ' '.join([s.name for s in schemes])
      msg = 'Delete %d color schemes [%s]?' % (len(schemes),names)
      if (showYesNo('Delete color schemes', msg, parent=self)):
        for scheme in schemes:
          scheme.delete()

  def updateColorList(self, *extra):

    rgb = Color.hexToRgb
    hsb = Color.hexToHsb

    scheme = self.colorScheme

    textMatrix  = []
    colorMatrix = []
    objectList  = []

    if scheme:
      n = 0
      for color in scheme.colors:
        r, g, b = rgb(color)
        h, s, v = hsb(color)
      
        datum = [n+1,
                 color,
                 r, g, b,
                 h, s, v,
                 ]
 
        textMatrix.append(datum)
        colorMatrix.append([None, color, None])
        objectList.append((n, color))

        n += 1
 
    self.colorTable.update(objectList=objectList,
                           textMatrix=textMatrix,
                           colorMatrix=colorMatrix) 
    
  def selectColor(self, schemeColor, row, col):

    index, color = schemeColor
    self.colorIndex = index
    self.color = color

  def getRed(self, schemeColor):
    
    index, color = schemeColor
    r, g, b = Color.hexToRgb(color) 
    self.editRedEntry.set(r)
  
  def getGreen(self, schemeColor):

    index, color = schemeColor
    r, g, b = Color.hexToRgb(color) 
    self.editGreenEntry.set(g)

  def getBlue(self, schemeColor):

    index, color = schemeColor
    r, g, b = Color.hexToRgb(color) 
    self.editBlueEntry.set(b)

  def getHue(self, schemeColor):

    index, color = schemeColor
    h, s, v = Color.hexToHsb(color) 
    self.editHueEntry.set(h)

  def getSaturation(self, schemeColor):

    index, color = schemeColor
    h, s, v = Color.hexToHsb(color) 
    self.editSaturationEntry.set(s)

  def getBrightness(self, schemeColor):

    index, color = schemeColor
    h, s, v = Color.hexToHsb(color) 
    self.editBrightnessEntry.set(v)

  
  def setRed(self, event):
  
    r, g, b = Color.hexToRgb(self.color)
    r = self.editRedEntry.get()
    r = max(0.0, min(1.0, r))

    colors = list(self.colorScheme.colors)
    colors[self.colorIndex] = Color.hexRepr(r,g,b)
    self.colorScheme.colors = colors
  
  def setGreen(self, event):
  
    r, g, b = Color.hexToRgb(self.color)
    g = self.editGreenEntry.get()
    g = max(0.0, min(1.0, g))
    
    colors = list(self.colorScheme.colors)
    colors[self.colorIndex] = Color.hexRepr(r,g,b)
    self.colorScheme.colors = colors
  
  def setBlue(self, event):
  
    r, g, b = Color.hexToRgb(self.color)
    b = self.editBlueEntry.get()
    b = max(0.0, min(1.0, b))
    
    colors = list(self.colorScheme.colors)
    colors[self.colorIndex] = Color.hexRepr(r,g,b)
    self.colorScheme.colors = colors
  
  def setHue(self, event):
  
    h, s, v = Color.hexToHsb(self.color)
    h = self.editHueEntry.get()
    h = max(0.0, min(1.0, h))
    
    r, g, b = Color.hsbToRgb(h, s, v)
    colors = list(self.colorScheme.colors)
    colors[self.colorIndex] = Color.hexRepr(r,g,b)
    self.colorScheme.colors = colors
    
  def setSaturation(self, event):
  
    h, s, v = Color.hexToHsb(self.color)
    s = self.editSaturationEntry.get()
    s = max(0.0, min(1.0, s))
    
    r, g, b = Color.hsbToRgb(h, s, v)
    colors = list(self.colorScheme.colors)
    colors[self.colorIndex] = Color.hexRepr(r,g,b)
    self.colorScheme.colors = colors
  
  def setBrightness(self, event):
  
    h, s, v = Color.hexToHsb(self.color)
    v = self.editBrightnessEntry.get()
    v = max(0.0, min(1.0, v))
    
    r, g, b = Color.hsbToRgb(h, s, v)
    colors = list(self.colorScheme.colors)
    colors[self.colorIndex] = Color.hexRepr(r,g,b)
    self.colorScheme.colors = colors
                         
  def getColor(self, schemeColor):

    index = 0
    names = []
    colors = []
    
    standards = Color.standardColors
    
    if self.colorScheme:
      color   = self.colorScheme.colors[self.colorIndex]
      colors  = [s.hex for s in standards]
      names   = [s.name for s in standards]
      
      if color in colors:
        index = colors.index(color) 
    
    self.editColorWidget.setup(names, colors, index, colors=colors)


  def setColor(self, event):

    if self.colorScheme:
      colors = list(self.colorScheme.colors)
      i = self.colorIndex
      color = self.editColorWidget.getObject()
      
      if i < len(colors):
        colors[i] = color
        self.colorScheme.colors = colors


  def addColor(self):

    # TBD: To be really cunning set color based on existing ones
    # in scheme - intermediate or follow trend

    scheme = self.colorScheme
    if not scheme:
      return

    scheme.addColor(DEFAULT_COLOR)
    
    lastObj = self.colorTable.objectList[-1]
    self.colorTable.selectObject(lastObj)


  def removeColor(self, *event):

    scheme = self.colorScheme
    if not scheme:
      return

    colors = list(scheme.colors)
    
    if len(colors) > 1:
      i = self.colorIndex
      
      if i < len(colors):
        colors.pop(i)
        scheme.colors = colors
 
      
  def smoothColors(self):
   
    objs = self.colorTable.currentObjects
    n = len(objs)
   
    if n > 2:
      i, start = objs[0]
      j, end = objs[-1]
      newColors = []
      r1, g1, b1 = Color.hexToRgb(start)
      r2, g2, b2 = Color.hexToRgb(end)
      
      m = n-1.0
      dr = (r2-r1)/m
      dg = (g2-g1)/m
      db = (b2-b1)/m
   
      for f in range(n):
        r = r1+(dr*f)
        g = g1+(dg*f)
        b = b1+(db*f)
        
        #r,g,b = Color.hsbToRgb(h,s,v)
        newColors.append(Color.hexRepr(r,g,b))
   
      colors = list(self.colorScheme.colors)
      colors = colors[:i] + newColors + colors[j+1:]
 
      self.colorScheme.colors = colors
   
  def reverseColors(self):
  
    if self.colorScheme:
      colors = list(self.colorScheme.colors)
      colors.reverse()
      self.colorScheme.colors = colors
  
  
  def colorUp(self):
  
    if self.colorScheme:
      i = self.colorIndex
      if i > 0:
        colors = list(self.colorScheme.colors)
        a, b = colors[i-1:i+1]
        colors[i-1:i+1] = b, a
        self.colorScheme.colors = colors
        self.colorTable.selectNthObject(i-1)
  
  def colorDown(self):

    if self.colorScheme:
      i = self.colorIndex
      colors = list(self.colorScheme.colors)
      if i < len(colors)-1:
        a, b = colors[i:i+2]
        colors[i:i+2] = b, a
        self.colorScheme.colors = colors
        self.colorTable.selectNthObject(i+1)
  
  def colorToTop(self):

    if self.colorScheme:
      i = self.colorIndex
      colors = list(self.colorScheme.colors)
      color = colors.pop(i)
      colors = [color,] + colors
      self.colorScheme.colors = colors
      self.colorTable.selectNthObject(0)
  
  def colorToBottom(self):

    if self.colorScheme:
      i = self.colorIndex
      colors = list(self.colorScheme.colors)
      color = colors.pop(i)
      colors = colors + [color,]
      self.colorScheme.colors = colors
      self.colorTable.selectNthObject(len(colors)-1)
            
  # Experiment Colors functions
  
  def selectExpProfile(self, obj, row, col):
  
    self.expProfile = obj
    self.updateExpColors()

  def removeExpProfile(self):
  
    if self.expProfile:
      self.expProfile.delete()
      self.expProfile = None
      
  def addExpProfile(self):
  
     if self.profile:
       name = None
       
       while name is None:
         msg = 'Name for class of experiments'
         name = askString('Name Query', msg, parent=self)
         if name is None:
           return

         if ' ' in name:
            name = None
            msg = 'Name cannot have spaces'
            showWarning('Retry', msg, parent=self)           

         elif self.profile.findFirstRefExpProfile(name=name):
            name = None
            msg = 'Name already in use within current Analysis profile'
            showWarning('Retry', msg, parent=self)           

       expProfile = self.profile.newRefExpProfile(name=name)
       self.expProfile = expProfile
       self.addExpColor()

  def updateExpProfiles(self, obj=None):

    if obj and (obj.analysisProfile is not self.profile):
      return

    textMatrix  = []
    expProfiles = []

    if self.profile:
      expProfiles = self.profile.sortedRefExpProfiles()

      for refExpProfile in expProfiles:
        refExpText = ''
        for name in refExpProfile.refExpNames:
          if len(refExpText) > 50:
            refExpText += '\n'
          elif refExpText:
            refExpText += ' '
          
          refExpText += name    
        
        n1 = len(refExpProfile.peakSymbolColors)
        n2 = len(refExpProfile.peakTextColors)
        n3 = len(refExpProfile.posColorSchemes)
        n4 = len(refExpProfile.negColorSchemes)
      
        nColorProfiles = max(n1, n2, n3, n4)
      
        datum = [refExpProfile.name,
                 refExpText,
                 nColorProfiles]
 
        textMatrix.append(datum)
 
    self.expProfileMatrix.update(objectList=expProfiles,
                                 textMatrix=textMatrix)
  
  def selectExpColorIndex(self, obj, row, col):

    self.expColorIndex = obj

  

  def updateExpColors(self, obj=None): 
  
    if obj and (obj is not self.expProfile):
      return
  
    textMatrix  = []
    colorMatrix = []
    objectList  = []

    if self.expProfile:
      peakSymbolColors = self.expProfile.peakSymbolColors
      peakTextColors   = self.expProfile.peakTextColors
      posColorSchemes  = self.expProfile.posColorSchemes
      negColorSchemes  = self.expProfile.negColorSchemes
      
      n1 = len(peakSymbolColors)
      n2 = len(peakTextColors)
      n3 = len(posColorSchemes)
      n4 = len(negColorSchemes)
      n  = max(n1, n2, n3, n4)
    
      for i in range(n):
        
        posText = negText = syMtext = annText = None
        colors = [None, None, None, None, None]
        
        if i < len(peakSymbolColors):
          symText = peakSymbolColors[i]
          colors[3] = symText
        
        if i < len(peakTextColors):
          annText = peakTextColors[i]
          colors[4] = annText
        
        if i < len(posColorSchemes):
          posText = posColorSchemes[i].name  
        
        if i < len(negColorSchemes):
          negText = negColorSchemes[i].name
        
        datum = [i+1,posText,negText,symText,annText]
        
        textMatrix.append(datum)
        objectList.append(i)
        colorMatrix.append(colors)
         


    self.expColorMatrix.update(objectList=objectList,
                               textMatrix=textMatrix,
                               colorMatrix=colorMatrix)
    
    
    
  def getPosScheme(self, i):
  
    texts = []
    objects = []
    colors = []
    index = 0
  
    if self.expProfile:
      schemes = self.profile.sortedColorSchemes()
    
      if i <len(self.expProfile.posColorSchemes):
        scheme = self.expProfile.posColorSchemes[i]
        index = schemes.index(scheme)
  
      for scheme in schemes:
        texts.append(scheme.name)
        objects.append(scheme)
        colors.append(list(scheme.colors))
  
    self.posSchemePulldown.setup(texts,objects,index, colors=colors)
   
  def getNegScheme(self, i):
  
    texts = []
    objects = []
    colors = []
    index = 0
    
    if self.expProfile:
      schemes = self.profile.sortedColorSchemes()
    
      if i <len(self.expProfile.negColorSchemes):
        scheme = self.expProfile.negColorSchemes[i]
        index = schemes.index(scheme)
  
      for scheme in schemes:
        texts.append(scheme.name)
        objects.append(scheme)
        colors.append(list(scheme.colors))
  
    self.negSchemePulldown.setup(texts,objects,index, colors=colors)
   
  def getSymColor(self, i):
  
    texts = []
    objects = []
    colors = []
    index = 0
   
    if self.expProfile:
      objects = Color.standardColors
      colors = [c.hex for c in objects]
      texts = [c.name for c in objects]
      
      if i < len(self.expProfile.peakSymbolColors):
        color = self.expProfile.peakSymbolColors[i]
        
        if color in colors:
          index = colors.index(color)
   
    self.symColorPulldown.setup(texts,objects,index, colors=colors)
   
  def getAnnColor(self, i):
  
    texts = []
    objects = []
    colors = []
    index = 0

    if self.expProfile:
      objects = Color.standardColors
      colors = [c.hex for c in objects]
      texts = [c.name for c in objects]
      
      if i < len(self.expProfile.peakTextColors):
        color = self.expProfile.peakTextColors[i]
        
        if color in colors:
          index = colors.index(color)
      
    self.annColorPulldown.setup(texts,objects,index, colors=colors)
   
  def setPosScheme(self, null):
 
    scheme = self.posSchemePulldown.getObject()
    
    if scheme:
      i = self.expColorIndex
      schemes = list(self.expProfile.posColorSchemes)
 
      if i < len(schemes):
        schemes[i] = scheme
      else:
        schemes.append(scheme)
 
      self.expProfile.posColorSchemes = schemes


  def setNegScheme(self, null):
 
    scheme = self.negSchemePulldown.getObject()
 
    if scheme:
      i = self.expColorIndex
      schemes = list(self.expProfile.negColorSchemes)
 
      if i < len(schemes):
        schemes[i] = scheme
      else:
        schemes.append(scheme)
 
      self.expProfile.negColorSchemes = schemes
 
  def setSymColor(self, color):
 
    i = self.expColorIndex
    colors = list(self.expProfile.peakSymbolColors)
  
    if i < len(colors):
      colors[i] = color.hex
    else:
      colors.append(color.hex)  

    self.expProfile.peakSymbolColors = colors
 
  def setAnnColor(self, color):
 
    i = self.expColorIndex
    colors = list(self.expProfile.peakTextColors)
  
    if i < len(colors):
      colors[i] = color.hex
    else:
      colors.append(color.hex)  

    self.expProfile.peakTextColors = colors

  def addExpColor(self):
  
    if self.expProfile:
  
      sym = list(self.expProfile.peakSymbolColors)
      ann = list(self.expProfile.peakTextColors)
      pos = list(self.expProfile.posColorSchemes)
      neg = list(self.expProfile.negColorSchemes)
   
      color = Color.inverseGrey(self.profile.bgColor)
      scheme = self.profile.findFirstColorScheme()
   
      sym.append(color)
      ann.append(color)
      pos.append(scheme)
      neg.append(scheme)
  
      self.expProfile.peakSymbolColors = sym
      self.expProfile.peakTextColors = ann
      self.expProfile.posColorSchemes = pos
      self.expProfile.negColorSchemes = neg
  
  def removeExpColor(self):
  
    i = self.expColorIndex
    
    if i is not None:
      sym = list(self.expProfile.peakSymbolColors)
      ann = list(self.expProfile.peakTextColors)
      pos = list(self.expProfile.posColorSchemes)
      neg = list(self.expProfile.negColorSchemes)
      
      if i < len(sym):
        del sym[i]

      if i < len(ann):
        del ann[i]
        
      if i < len(pos):
        del pos[i]
        
      if i < len(neg):
        del neg[i]    

      self.expProfile.peakSymbolColors = sym
      self.expProfile.peakTextColors = ann
      self.expProfile.posColorSchemes = pos
      self.expProfile.negColorSchemes = neg

  def setRefExpNames(self, opt):
  
    values = self.refExpNameSelect.get()
    self.expProfileMatrix.keyPressEscape()
  
    if values is not None:
      names = [v for v in values if v]
      self.expProfile.refExpNames = names
    
  def getRefExpNames(self, expProfile):
  
    values  = list(expProfile.refExpNames)
    
    self.refExpNameSelect.set(values=values)
