
"""
======================COPYRIGHT/LICENSE START==========================

WindowPopup.py: Part of the CcpNmr Analysis program

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

import Tkinter

from memops.universal.Io import getTopDirectory

from memops.gui.ButtonList          import ButtonList, UtilityButtonList
from memops.gui.CheckButton         import CheckButton
from memops.gui.Label               import Label
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.ScrolledFrame       import ScrolledFrame

from ccpnmr.analysis.popups.BasePopup import BasePopup

from ccpnmr.analysis.core.ExperimentBasic import changeSpectrumNContours
from ccpnmr.analysis.popups.GetCenterCoords import GetCenterCoordsPopup
 
from ccpnmr.analysis.frames.WindowFrame import WindowFrame
from ccpnmr.analysis.core import Util
from ccpnmr.analysis.core.WindowBasic import setWindowPeakListState, \
                  toggleSpectrum, swapStripAxis, windowPaneHasValueAxis

TOP_BAR_TIP_TEXTS = {'Spectra':'Opens a panel toggle the display of individual spectra on or off',
                     'Contours':'Opens a panel to change contour levels for the currently displayed spectra',
                     'Peaks':'Opens a panel toggle the display of individual peak lists (cross marks and annotations) on or off',
                     'Strips':'Opens a panel to control the number, order and orientation of strips (parallel window sub-divisions)',
                     'x,y':'The current cursor location for the window axes, click to center the view or change cursor position units'}

class WindowPopup(BasePopup):
  """
  **Spectrum Window Contour Display**
  
  The spectrum display in Analysis displays NMR spectra of arbitrary
  dimensionality and allows the user to identify signal peaks and assign them
  to atomic resonances. A CCPN project may contain as many spectrum windows as
  memory allows with whichever axis types the user requires.
  
  Simple 1-dimensional spectra are displayed as intensity graphs, but higher
  dimensionality spectra are displayed as contour plots. A contour plot
  represents the intensities in the plane of an NMR spectrum and each contour
  line represents points of equal data value (height). For a 2-dimensional
  spectrum there is only one displayable spectrum plane, but for 3D and higher
  there are many planes, which can be imagined as being stacked on top of one
  another along depth axes. The user can navigate to view any  position in a
  spectrum by changing the region of the plane displayed in the screen and
  independently adjust what the current depth position is. Also, in windows with
  depth axes Analysis has the ability to display several planes at once in a
  specified depth region, i.e. a given thickness. A spectrum window always has
  the same kind of spectrum dimensions (e.g. in terms of isotope) represented
  by the on-screen (X,Y) and depth (Z) axes. If an alternative, orthogonal view
  is required, for example to show different NMR axes in the plane of the
  screen then a new window must be made. However, if the user  simply wishes to
  swap spectrum data axes that have the *same* kind of isotope then this can be
  done without having to make a new window; the mapping of spectra into windows
  can be adjusted via the "Spectrum & Peak List Mappings" section of the main
  `Windows`_ option.

  There are three basic components to a spectrum window popup: there is the
  main graphical spectrum display with axes and scrollbars, which may be
  sub-divided into several strips/panels; a toolbar of buttons at the top,
  which control what is displayed in the window; and also a menu accessible
  with a right mouse click over the main display, which has a large amount of
  important functionality.

  **Top Toolbar**

  The toolbar of buttons above the main spectrum display area allows the user
  to control various aspects of what is displayed. After clicking on one of the
  top buttons an extra panel will appear to give access to further options, and
  clicking on the same top button again will make the options panel disappear.

  *Spectra*

  The [Spectra] option allows the user to toggle individual spectra within the
  main display on or off. The spectra that can be displayed are naturally
  those that have the right kinds of isotopes on their axes for the spectrum
  window. Although, a spectrum may be displayed at least the X & Y axes match,
  i.e. a 2D spectrum can be displayed in the X-Y plane of a 3D window. Spectra
  may be prevented from appearing within this toolbar, to avoid overcrowding,
  by changing the "Spectrum In Toolbar" column within the "Spectrum & Peak List
  Mappings" section of the main `Windows`_ option.

  Spectra may also be toggled on an off using keyboard shortcuts, which are set
  in the main `Spectra`_ display options. Also, the order in which spectra
  are drawn is reflected in the order that they appear in the top panel;
  this is set via the "Rank" of the `Spectra`_ display options.

  *Contours*

  The [Contours] option allows the user to adjust the contour base level (green
  arrows), the number of contour lines and flip between positive and negative
  levels. These settings will be applied to *all of the visible spectra* in
  the window. Hence, it is common for the user to toggle certain spectra off
  first if they should not be affected. These contour settings are fairly
  simple, but more detailed control of contours is given by the [More..]
  button which opens the main `Spectrum Contour Levels`_ popup.

  *Peaks*
  
  The [Peaks] option allows the user to control which peak lists, in terms of
  the peak markers and annotations, will be visible in the main graphical
  display. The user simply clicks on the button of a peak list to turn it on or
  off. By default the "Multi-list" option is set so that any number of peak 
  list will be displayed, otherwise the user can uncheck this so that only one
  peak list is shown at a time. Having only one peak list shown is useful to
  avoid picking peaks in the wrong spectra; peaks will only be picked in
  displayed peak lists.

  It should be noted that when a spectrum has more than one peak list only one
  of these will be active; the active peak list is the one that the user can
  select in the spectrum and the one newly picked signal maxima are placed in.
  To change the active peak list for a spectrum and how peak lists are shown
  (symbol type and colour) the user should edit the main `Peak Lists`_ table.

  *Strips*
  
  The [Strips] option controls how the graphical spectrum display is
  sub-divided. These sub-divisions are called "strips" and allow the user to
  have different panels in the window showing different regions on the spectra.
  All strips will be tied together along one axis; with vertical strips the
  Y-axis is the same for all sub-divisions, and for horizontal strips the
  X-axis is common. Strips can be added and removed by the user, although
  many systems in Analysis will automatically curate strips (e.g, for protein
  sequence assignment). The strips' panel also indicates which 
  strip in the window is deemed to be active (its button will be green);
  it is the active strip that can be moved with the green arrow and is
  the target of any navigation options into the window.

  Spectrum windows in Analysis have another way of defining sub-divisions which
  are called "separators". These are distinct from strips and are not generally
  administered from the toolbar (except the ability to remove them all). A
  separator is a sub-division that goes at right angles to the main strip
  divisions and serves to split the axis which is not affected by the strips.
  For example a window with vertical strips can have horizontal separators
  which split the Y-axis into separate regions to display separate PPM ranges
  (e.g. CA, CB, C'). Separators may be added via the right mouse `Window
  Menu`_. 

  *Coordinate Location*
  
  The last button in the top toolbar displays the current location of the mouse
  cursor for the window. Clicking on the button  opens the `Center
  coordinates`_ panel which  allows the user to specify a particular location
  for the center of the window view and also allows the cursor location to be
  displayed in Hz units, rather than PPM.

  **Spectrum Display Area**

  The main graphical display area is where spectrum lines, contours and peak
  markers are displayed. At the edge of the area are tick axes with numbers
  indicating the display region (usually in PPM). The display region may be
  changed by using the scrollbars or the various mouse controls listed below.
  Note that the scrollbars may be removed to save space using the right mouse
  `Window Menu`_.  The box at the bottom right corner where the tick axes meet
  shows the position of the current view, according to zoom level and
  position, relative to the maximum possible extent of the window axes. 

  For 3D and higher dimensionality windows there will be extra scrollbars at
  the bottom representing the depth axes. The current position of these is
  displayed to the left of the scrollbar, and this may be set by typing in the
  entry box. These depth (or "orthogonal") scrollbars are different to the X-Y
  scrollbars because they also control the thickness of the depth view (i.e.
  the number of planes deep to show contours for). To change the thickness of a
  view the left or right edge of the depth slider can be expanded or shrunk
  using the middle mouse button.
  
  *Window Navigation*
  
  The cursor may change the region of the spectrum (or spectra) displayed in
  the graphical area by using both keyboard and mouse controls that allow the
  view to pan, zoom and change orthogonal depth position. These controls are
  summarised as follows, noting that the middle mouse button is specialised
  for these functions:
  
  :Pan Window: Arrow keys
               Middle-click + drag

  :Zoom In/Out: <PageDown>/<PageDown>
                Scroll wheel forward/back
                <Shift> + Middle-click + drag up/down

  :Zoom To Region: <Control> + Middle-click + drag region

  :Depth Scroll: <Control> + scroll wheel forward/back
                 <Shift> + scroll wheel forward/back

  *Peak Picking*

  The intensity maxima and minima in NMR spectra which represent the locations
  of resonance signals may be marked as being a peak, i.e. "picked". Picked
  peaks are marked with crosses at their centre and carry an assignment
  annotation to indicate the resonances from which the peak derives. For an
  unassigned peak this annotation will just be dashes (one for each spectrum
  dimension).

  Spectrum peaks may be picked in three basic ways in Analysis: by specifying
  an exact peak location with the mouse, by dragging a rectangular peak-find
  region with the mouse to find extrema, and by using the "Region Peak Find"
  functionality of the `Peak Finding`_ section to find extrema in a large
  spectrum region. In all cases newly picked peaks are placed in the active
  peak list for a given spectrum (set in the `Peak Lists` popup). The mouse
  and keyboard controls used to pick peaks are listed below, noting that the
  <Control> key and left mouse button are specialised for this purpose:

  :Pick Peak at Mouse Position: <Control> + left-click
  
  :Peak Peaks in Region: <Control> + <shift> + left-click & drag
  
  It should be noted that when peaks are picked all visible spectra are
  considered for picking, but only if their active peak list is toggled on (via
  the top tab). By default, if an exact peak position is specified then new
  peaks will be picked in all visible spectra. Similarly, when dragging a box
  to find intensity minima/maxima in a region, peaks are picked for the extrema
  found in all visible spectra. If the user wants to avoid picking peaks in
  specific spectra then a spectrum can be toggled off in the top panel, or the
  active peak list for the spectrum can be toggled off (leaving any contours
  visible). For the region search peaks will be located according to the "Find
  Parameters" set in the `Peak Finding`_ popup; by default extrema are only
  found above the same threshold that contours are displayed for.

  Note that if a spectrum is tiled so that its contours are visible outside the
  normal bounds any peak picked in a tiled region will be picked at the
  apparent PPM value; the peak is recorded internally at the equivalent
  position inside the spectrum, but its position is automatically "unaliased"
  to move it outside the spectrum width. Here the underlying, aliased peak
  inside the spectrum width will also be displayed, albeit with a dashed cross.

  *Peak Selection*
  
  Specific visible peaks may be selected to that specific operations may be
  performed on them (assign, delete...). Peaks that have been selected are
  displayed with a solid box around their marker crosses. Peaks can be selected
  in several different spectra and in several different windows, even if the
  peak locations are not currently visible; so caution should be used when
  deleting peaks. Beak selections can be defined and added to according to the
  mouse controls listed below, noting common use of the left mouse button.
  Also, a number of operations may be performed on the selected peak(s) by
  using keyboard shortcuts.

  :Select Single Peak: Left-click at center

  :Select Several Peaks: Left-click & drag box over peaks
  
  :Add Single Peak To Selection: <Shift> + Left-click
  
  :Add Several Peaks To Selection: <Shift> + Left-click & drag box over peaks

  :Clear Peak Selection: Left-click in a region without peaks

  :Move Selected Peak: <p> + drag mouse

  :Re-center Selected Peaks: <P>

  :Assign Peak Under Mouse: <a> (opens `Assignment Panel`)
  
  :Show Selected Peaks Table: <s> (opens `Selected Peaks`)
  
  :Delete Selected Peaks: <Delete>

  :Move Selected Peaks Label: <q> + drag mouse
  
  Note that the keyboard shortcuts listed are the default settings
  and could actually be different. All keyboard shortcuts are listed
  (and can be changed) in the "Macros" section of the `User Options`_
  popup.
  
  *Marker Lines & Rulers*
  
  The user may place dashed, multi-dimensional marker lines (i.e. a big cross)
  at the current mouse location using <m>. Also a 1-dimensional dashed ruler
  line may be placed at the mouse location using <h> for a horizontal line and
  <v> for a vertical line. Any multi-dimensional mark that is suitably close to
  a peak will automatically lock on to the exact peak location. The <n> key is
  used to remove all marks and rulers. 
  
  For a detailed description of marker lines and rulers and their various
  options see the `Marks and Rulers`_ documentation.

  *1D Traces*
  
  An spectrum window can additionally carry a number of 1-dimensional (i.e.
  intensity graph) slice traces, which show a cross section through the
  displayed spectra at the current mouse location. Such slice traces may be
  added via the right mouse menu "Window" option, selection either crosshair
  traces that are superimposed in the main graphical display or side traces
  that are shown as extra panels to the side of the window.

  **Right Mouse Menu**

  The right mouse button, when pressed over the graphical display area, will
  open up a menu that allow the user to perform various operations. These
  functions are sub-divided into  the main categories listed below. For more
  detailed descriptions see the `Window Menu`_ section.

  *Assign*

  These options control the assignment of the selected peak or peaks to
  resonances and spin systems (intra residue groups). For example the user can
  assign a peak, propagate (spread) assignments over several peaks, remove
  assignments etc.

  *Peak*
  
  These options relate to general operations performed on peaks, including
  deletion, setting aliasing, setting details and arranging labels. Many of
  these functions have keyboard shortcuts (see "Macros" section of the `User
  Options`_).

  *Locate Peaks*
  
  These functions are used to find the selected peaks in the spectrum window or
  use the locations of peaks to navigate to matching peaks (or equivalent
  positions). Some of these options include navigation to other windows.

  *View*
  
  The view functions control the current spectrum region being
  displayed and the contour levels.
  
  *Navigate*
 
  The navigation options allow the user to locate equivalent positions in other
  spectrum windows based upon the mouse location. Such locations are found in
  windows that have at least some of the same kinds of axes as the current
  window. For example the user may navigate from a 2D HN window displaying an
  HSQC spectrum to the equivalent HN (amide) position in a 3D window.

  *Strip*
  
  These options are used to control strips and separators; the two kinds
  of spectrum window sub-division.
 
  *Markers*

  The marker options control the presence of dashed position marker lines; 
  both  multi-dimensional markers and 1-dimensional rulers.

  *Window*

  This section provides various general functions for the spectrum window which
  include: making a duplicate "clone" of a window, printing a window to file
  (e.g. PDF) and controlling the 1D slice traces.

  *Macros*
 
  This menu contains lists any Python macro scripts that have been added to
  Analysis (see `User Options`_) to operate in the context of a spectrum
  window.

  .. _`Center coordinates`: GetCenterCoordsPopup.html
  .. _`Windows`: EditWindowPopup.html
  .. _`Spectra`: EditSpectrumPopup.html
  .. _`Spectrum Contour Levels`: EditContourLevelsPopup.html
  .. _`Peak Lists`: EditPeakListsPopup.html
  .. _`Window Menu`: ../menu/WindowMenus.html
  .. _`Peak Finding`: EditPeakFindParamsPopup.html
  .. _`User Options`: EditProfilesPopup.html
  .. _`Marks and Rulers`: EditMarksPopup.html
  .. _`Assignment Panel`: EditAssignmentPopup.html
  .. _`Selected Peaks`: EditAssignmentPopup.html
  """
  
  UNIT_PPM = 'ppm'
  UNIT_HZ = 'Hz'
  UNITS = [UNIT_PPM, UNIT_HZ]

  def __init__(self, parent, window, *args, **kw):

    self.window = window

    self.fontChanged = False
    self.activeButtonFrame = None
    self.waitDraw = False
    self.waitStrips   = False
    self.waitPeakList = False
    self.activeButtonFrame = None
    self.isBeingDestroyed = False
    self.monoPeakList = None

    self.positionUnit = self.UNIT_PPM
    self.positionView = None
    self.positionKeepView = False

    # TBD: must remember and be able to modify activeWindowPane and activeWindowFrame
    self.activeWindowPane = self.window.findFirstSpectrumWindowPane()
    self.activeWindowFrame = None

    BasePopup.__init__(self, parent=parent, **kw)

  def body(self, master):
    
    path = os.path
    gfxDir = path.join(getTopDirectory(),'python','memops','gui','graphics')
    MakeImage = Tkinter.PhotoImage
    self.upArrow    = MakeImage(file=path.join(gfxDir,'arrowUp.gif'))
    self.downArrow  = MakeImage(file=path.join(gfxDir,'arrowDown.gif'))
    self.leftArrow  = MakeImage(file=path.join(gfxDir,'arrowLeft.gif'))
    self.rightArrow = MakeImage(file=path.join(gfxDir,'arrowRight.gif'))
    self.vertLines  = MakeImage(file=path.join(gfxDir,'verticalBars.gif'))
    self.horizLines = MakeImage(file=path.join(gfxDir,'horizontalBars.gif'))
    self.switchHV   = MakeImage(file=path.join(gfxDir,'switchHV.gif'))
    self.plusIcon   = MakeImage(file=path.join(gfxDir,'plus.gif'))
    self.minusIcon  = MakeImage(file=path.join(gfxDir,'minus.gif'))
    self.toolIcon   = MakeImage(file=path.join(gfxDir,'tools.gif'))

    self.font = self.analysisProfile.font or 'Helvetica 10'

    activeWindowPane = self.activeWindowPane

    master.grid_rowconfigure(2, weight=1)
    master.grid_columnconfigure(0, weight=1)

    # TBD: need to hide Contour button if activeWindowPane has value axis
    # for now just look at findFirstSpectrumWindowPane()
    if windowPaneHasValueAxis(self.window.findFirstSpectrumWindowPane()):
      # TBD: a bit of a hack
      if len(self.window.spectrumWindowPanes) > 1:
        options = ['Spectra','Peaks','Shapes','x,y']
        commands = [self.selSpectra,
                    self.selPeaks, self.selShapes,
                    self.getCenterCoords]
        tipTexts = [TOP_BAR_TIP_TEXTS.get(x) for x in options]
      
      else:
        options = ['Spectra','Peaks','Strips','x,y']
        commands = [self.selSpectra,
                    self.selPeaks, self.selStrips,
                    self.getCenterCoords]
        tipTexts = [TOP_BAR_TIP_TEXTS.get(x) for x in options]

    else:
      options = ['Spectra','Contours','Peaks','Strips','x,y']
      commands = [self.selSpectra, self.selContours,
                  self.selPeaks, self.selStrips,
                  self.getCenterCoords]
      tipTexts = [TOP_BAR_TIP_TEXTS.get(x) for x in options]

    scrolledFrame = ScrolledFrame(master, bd=1, relief='raised', yscroll=False)
    scrolledFrame.grid(row=0, column=0, sticky='nsew')
    frame = scrolledFrame.frame
    
    self.buttonList = UtilityButtonList(frame, texts=options, commands=commands, 
                                        helpUrl=self.help_url, doClone=False,
                                        doClose=False, tipTexts=tipTexts, grid=(0,0))
    self.locationButton = self.buttonList.buttons[-2]
      
    if activeWindowPane:
      xname = activeWindowPane.findFirstAxisPanel(label='x').axisType.name
      yname = activeWindowPane.findFirstAxisPanel(label='y').axisType.name
    else:
      xname = yname = '1H'
    
    xyName = '(%s,%s)' % (xname,yname)
    self.locationButton.config(text=xyName)
    
    self.updatePositionView()

    # Contour frame
    
    self.contourFrame = ScrolledFrame(master, bd=1, relief='raised', yscroll=False)
    frame = self.contourFrame.frame
    #frame.grid_columnconfigure(0, weight=1)

    tipTexts = ['Move the contour base level, for displayed spectra, up (away from zero)',
                'Move the contour base level, for displayed spectra, down (towards zero)',
                'Add another of contour level to the displayed spectra (to the top)',
                'Remove a contour level from the displayed spectra (from the top)',
                'Cycle the contour display through positive & negative contour options; negative, positive both',
                'Open a popup window to give greater control of contouring for individual spectra']
    texts = ['','','+1','-1','Pos/Neg','More..']
    commands = [self.contourLevelsUp,
                self.contourLevelsDown,
                self.ncontoursIncr,
                self.ncontoursDecr,
                self.posNegContoursToggle,
                self.generalContours]
      
    buttonList = ButtonList(frame, texts=texts, tipTexts=tipTexts,
                              commands=commands, grid=(0,0))
    buttons = buttonList.buttons
    buttons[0].config(image=self.upArrow)
    buttons[1].config(image=self.downArrow)

    # Spectrum Frame

    self.specFrame = ScrolledFrame(master, bd=1, relief='raised', yscroll=False)
    frame = self.specFrame.frame
    #frame.grid_columnconfigure(1, weight=1)
    
    tipText = 'Selects which spectra are currently displayed in the window (assuming they have contours)'
    self.spectrumToggleSelector = PartitionedSelector(frame, self.toggleSpectrum,
                                                      font=self.font, tipText=tipText,
                                                      grid=(0,0))
    
    # Peak Frame
    
    self.peakFrame = ScrolledFrame(master, bd=1, relief='raised', yscroll=False)
    frame = self.peakFrame.frame
    frame.grid_columnconfigure(2, weight=1)
    
    tipText = 'Selects which peak lists (cross marks and text annotations) are displayed in the window'
    self.peaksToggleSelector = PartitionedSelector(frame, self.togglePeaks,
                                                   font=self.font, tipText=tipText,
                                                   grid=(0,0), sticky='ew')
                                                    
    tipText = 'Sets whether multiple peak lists are displayed; the peak list buttons work independently or in "radio" mode, where only one is selected at a time'
    self.peakListCheck = CheckButton(frame, text='Multi-list', font=self.font,
                                     callback=self.setMultiPeakLists, tipText=tipText,
                                     selected=self.window.useMultiplePeakLists, grid=(0,1))
    self.setMultiPeakLists(self.window.useMultiplePeakLists)

    # Strip Frame
   
    self.stripFrame = ScrolledFrame(master, bd=1, relief='raised', yscroll=False)
    frame = self.stripFrame.frame
    frame.grid_columnconfigure(3, weight=1)
 
    tipTexts = ['Add a new strip to the window (parallel sub-divisions)',
                'Remove the active strip from the window (double click to select active)',
                'Shuffle the active strip one position to the left (double click to select active)',
                'Shuffle the active strip one position to the right (double click to select active)',
                'Switch between vertical strips (parallel to Y-axis) and horizontal strips (parallel to X-axis)']
    texts    = ['','','','','']
    commands = [self.addStrip,
                self.deleteStrip,
                self.moveStripsLeft,
                self.moveStripsRight,
                self.toggleStripDir]
                
    self.stripsButtons = ButtonList(frame, texts=texts,
                                    commands=commands,
                                    tipTexts=tipTexts)
                                    
    buttons = self.stripsButtons.buttons
    buttons[0].config(image=self.plusIcon)
    buttons[1].config(image=self.minusIcon)
    if self.window.stripAxis == 'x':
       buttons[2].config(image=self.leftArrow)
       buttons[3].config(image=self.rightArrow)
    else:
       buttons[2].config(image=self.upArrow)
       buttons[3].config(image=self.downArrow)
    buttons[4].config(image=self.switchHV)
    
    self.stripsButtons.grid(row=0, column=0, sticky='w')

    label = Label(frame, text='Clear:', grid=(0,1))
    
    tipTexts = ['Delete all vertical strips, or sub-divisions, from the window',
                'Delete all horizontal strips, or sub-divisions, from the window']
    texts    = ['','']
    commands = [self.deleteVert,
                self.deleteHoriz]
    
    clearButtons = ButtonList(frame, texts=texts, grid=(0,2),
                              commands=commands, tipTexts=tipTexts)
    
    buttons = clearButtons.buttons
    buttons[0].config(image=self.vertLines)
    buttons[1].config(image=self.horizLines)
   
    # self.stripsLabel = Label(frame, text='Active:', font=self.font, grid=(0,3))
    self.stripsLabel = Label(frame, text='Active:', grid=(0,3))

    tipText = 'Selects which strips are active; their locations are moved with the Z-axis scrollbar and the first active strip is the one shuffled or deleted'
    self.stripsToggleSelector = PartitionedSelector(frame, self.toggleStrips,
                                                    radio=True, font=self.font,
                                                    grid=(0,4), tipText=tipText)
    
    if activeWindowPane:
      ncols = len(activeWindowPane.findFirstAxisPanel(label='x').axisRegions)
      nrows = len(activeWindowPane.findFirstAxisPanel(label='y').axisRegions)
    else:
      ncols = nrows = 1
    if (self.window.stripAxis == 'x' and ncols < 2) or (self.window.stripAxis == 'y' and nrows < 2):
      self.stripsButtons.buttons[1].disable()

    # Shapes Frame
   
    self.shapeFrame = ScrolledFrame(master, bd=1, relief='raised', yscroll=False)
    frame = self.shapeFrame.frame
    frame.grid_columnconfigure(3, weight=1)
 
    # Window Panes

    self.windowFrames = []
    for windowPane in self.window.sortedSpectrumWindowPanes():
      windowFrame = WindowFrame(master, windowPane)
      self.windowFrames.append(windowFrame)
      if windowPane == self.activeWindowPane:
        self.activeWindowFrame = windowFrame

    self.gridAll()

    self.setWindowTitle()

    #self.firstPass = True
    self.bind('<Configure>', self.getWindowLocation)
    #self.bind('<Expose>', self.checkIconified)

    self.setSpectrumSelector()
    self.setStripsSelectorAfter()

    #if (self.window.isIconified):
    #  self.after_idle(self.close)

    self.curateNotifiers(self.registerNotify)

  # following five functions just used to keep isIconified up-to-date

  def updateAll(self):

    self.setSpectrumSelector()
    self.setPeaksSelector()
    self.setStripsSelector()

    self.drawAllAfter()

  def open(self):

    #print 'WindowPopup open'
    self.window.isIconified = False
    self.update_idletasks()
    self.updateAll()
    BasePopup.open(self)

  def deiconify(self):

    #print 'WindowPopup deiconify'
    self.window.isIconified = False
    BasePopup.deiconify(self)

  def iconify(self):

    #print 'WindowPopup iconify'
    for windowFrame in self.windowFrames:
      windowFrame.unpostMenu()
    
    self.window.isIconified = True
    BasePopup.iconify(self)

  def close(self):

    #print 'WindowPopup close'
    for windowFrame in self.windowFrames:
      windowFrame.unpostMenu()
    
    self.window.isIconified = True
    BasePopup.close(self)

  def withdraw(self):

    #print 'WindowPopup withdraw'
    for windowFrame in self.windowFrames:
      windowFrame.unpostMenu()
    
    self.window.isIconified = True
    BasePopup.withdraw(self)

  """
  def checkIconified(self, *event):

    # if you do this too early then window location is trampled when uniconified
    # a bit of a pain because it leads to a flash of the window
    if (self.firstPass):
      self.firstPass = False
      if (self.window.isIconified):
        self.close()
  """

  def changedIsIconified(self, window):

    if self.window == window:
      if window.isIconified:
        self.close()
      else:
        self.open()

  def getPositionUnits(self):

    return self.UNITS

  def getPositionUnit(self):

    ##analysisProject = self.window.analysisProject
    # TBD: put into data model
    ##unit = self.parent.application.getValue(analysisProject, keyword='positionUnit')
    ##if not unit:
    ##  unit = self.UNIT_PPM
    unit = self.positionUnit

    return unit

  def setPositionUnit(self, unit):

    ##analysisProject = self.window.analysisProject
    ##self.parent.application.setValue(analysisProject, keyword='positionUnit', value=unit)
    self.positionUnit = unit

  def getPositionViews(self):

    def isOkView(view, includePeakListCheck=True):

      if not view:
        return False

      analysisSpectrum = view.analysisSpectrum
      spectrum = analysisSpectrum.dataSource

      # need spectrum check because changing app date on delete calls this update before the view is removed
      if not (spectrum and self.isViewVisible(view) and view.isInToolbar):
        return False

      if includePeakListCheck:
        return [x for x in view.windowPeakLists if x.isSymbolDrawn]
      else:
        return True

    views = self.getSpectrumViews()
    vv = [x for x in views if isOkView(x, includePeakListCheck=True)]
    if not vv:
      vv = [x for x in views if isOkView(x, includePeakListCheck=False)]
    if not vv:
      vv = views

    return vv

  def updatePositionView(self):

    views = self.getPositionViews()
    if views and self.positionView not in views:
      self.positionView = views[0]

  def getPositionView(self):

    return self.positionView

  def setPositionView(self, view):

    self.positionView = view

  # TBD: should we get rid of this? (it needs a windowPane)
  # 14 Oct 09: now also being used to determine position unit
  def getCenterCoords(self, *event):

    self.updatePositionView()
    popup = GetCenterCoordsPopup(self)
    center = popup.center
    region_number = popup.region_number
    popup.destroy()

    if center:
      (cx, cy) = center
      (rx, ry) = region_number
      windowPane = self.activeWindowPane
      findFirstAxisPanel = windowPane.findFirstAxisPanel
      xaxisPanel = findFirstAxisPanel(label='x')
      yaxisPanel = findFirstAxisPanel(label='y')
      (x0, x1) = xaxisPanel.sortedAxisRegions()[rx].region
      (y0, y1) = yaxisPanel.sortedAxisRegions()[ry].region
      tx = (cx - 0.5*(x0+x1)) / (x1 - x0)
      ty = (cy - 0.5*(y0+y1)) / (y1 - y0)
      if xaxisPanel.axisUnit.isBackwards:
        tx = - tx
      if yaxisPanel.axisUnit and yaxisPanel.axisUnit.isBackwards:
        ty = - ty
      
      windowFrame = self.activeWindowFrame
      canvas = windowFrame.scrolled_window.canvases[ry][rx]
      windowFrame.scrolled_window.translate(canvas, tx, ty)

  def toggleStripButtons(self):

    buttons = self.stripsButtons.buttons
    if self.window.stripAxis == 'y':
       buttons[2].config(image=self.leftArrow)
       buttons[3].config(image=self.rightArrow)
    else:
       buttons[2].config(image=self.upArrow)
       buttons[3].config(image=self.downArrow)

  def toggleStripDir(self):

    self.toggleStripButtons()
    window  = self.window
    swapStripAxis(self.parent, window)

  def setSpectrumSelector(self, *notifyObj):

    spectra = []
    colors  = []
    labels  = []
    views   = []
    refViews = self.getSpectrumViews()

    for view in refViews:
      analysisSpectrum = view.analysisSpectrum
      spectrum = analysisSpectrum.dataSource
      # need below check because changing app date on delete calls this update before the view is removed
      if spectrum and view.isInToolbar:
        views.append(view)
        
        posColors = list(analysisSpectrum.posColors)
        if posColors:
          color = posColors[int(len(posColors)*0.5)]
        else:
          color = '#808080'
          analysisSpectrum.posColors = [color,]
          print 'Warning %s missing positive color scheme' % analysisSpectrum

        spectra.append( spectrum )
        colors.append( color )
        labels.append('%s:%s' %  (spectrum.experiment.name,spectrum.name))

    self.spectrumToggleSelector.update(objects=spectra,labels=labels,colors=colors)
    n = 0
    for view in views:
      self.spectrumToggleSelector.setButtonState(n, self.isViewVisible(view))
      n = n + 1
    self.setPeaksSelectorAfter()

  def getWindowXySf(self, *extra):

    view = self.positionView

    if view:
      sfs = []
      for label in ('x', 'y'):
        axisMapping = view.findFirstAxisMapping(label=label)
        sf = 1.0
        if axisMapping:
          expDimRef = axisMapping.analysisDataDim.dataDim.expDim.findFirstExpDimRef()
          if expDimRef:
            sf = expDimRef.sf
        sfs.append(sf)
    else:
      sfs = [1.0, 1.0]

    return sfs

  def setMultiPeakLists(self, allowMultiple):

    if not allowMultiple:
      winPeakLists = self.peaksToggleSelector.objects
      if winPeakLists:
        active = None
        for winPeakList in winPeakLists:
          if winPeakList.isSymbolDrawn:
            self.monoPeakList = active = winPeakList
            setWindowPeakListState(active, True)
            break

        for winPeakList in winPeakLists:
          if winPeakList is not active:
            setWindowPeakListState(winPeakList, False)

    if allowMultiple:
      boolean = True
      self.monoPeakList = None
    else:
      boolean = False

    self.window.useMultiplePeakLists = boolean

  def setPeaksSelectorAfter(self, *notifyObj):

    if self.waitPeakList:
      return

    else:
     self.waitPeakList = 1
     self.after_idle(self.setPeaksSelector)

  def setPeaksSelector(self, *notifyObj):

    if self.window.isDeleted:
      return

    peakLists = []
    colors  = []
    labels  = []
    winPeakLists = []
    allViews = self.getSpectrumViews()

    for view in allViews:
      analysisSpectrum = view.analysisSpectrum
      spectrum = analysisSpectrum.dataSource
      # need below check because changing app date on delete calls this update before the view is removed
      if spectrum and self.isViewVisible(view) and view.isInToolbar:
        
        if view.windowPeakLists:
          posColors = list(analysisSpectrum.posColors)
          color = posColors[int(len(posColors)*0.8)]

          sortList = []
          for winPeakList in view.windowPeakLists:
            try:
              peakList = winPeakList.analysisPeakList.peakList
            except:
              peakList = None
            if not peakList or peakList.isDeleted:
              continue
            
            sortList.append((peakList.serial, peakList, winPeakList))
          
          sortList.sort()
          for serial, peakList, winPeakList in sortList:
            winPeakLists.append(winPeakList)
            peakLists.append(peakList)
            colors.append( color )
            data = (spectrum.experiment.name,spectrum.name,serial)
            labels.append('%s:%s:%d' % data )
            
    if not self.window.useMultiplePeakLists:
      active = [wpl for wpl in winPeakLists if wpl.isSymbolDrawn]

      if self.monoPeakList not in active:
        if active:
          self.monoPeakList = active[0]
        
        else:
          self.monoPeakList = None # A conscious decision
        
      for winPeakList in winPeakLists:
        if winPeakList is not self.monoPeakList:
          setWindowPeakListState(winPeakList, False)

    self.peaksToggleSelector.update(objects=winPeakLists,labels=labels,colors=colors)

    for n, winPeakList in enumerate(winPeakLists):
      self.peaksToggleSelector.setButtonState(n, winPeakList.isSymbolDrawn)

    self.waitPeakList = False
    self.updatePositionView()

  def setStripsSelectorAfter(self, *notifyObj):

    if notifyObj:
      axisRegion = notifyObj[0]
      if axisRegion.axisPanel.spectrumWindowPane != self.activeWindowPane:
        return

    if self.waitStrips:
      return
    else:
     self.waitStrips = 1
     self.after_idle(self.setStripsSelector)

  def setStripsSelector(self):

    if self.window.isDeleted:
      return

    activeWindowPane = self.activeWindowPane
    if not activeWindowPane:
      return

    axisRegions = activeWindowPane.findFirstAxisPanel(label=self.window.stripAxis).sortedAxisRegions()
    N = len(axisRegions)

    if N == 0:
      self.stripsToggleSelector.grid_forget()

    elif N ==1:
      self.stripsToggleSelector.update(objects=[axisRegions[0],],labels=['1',],colors=['#80A080',])
      self.stripsToggleSelector.selectButton(0, doCallback=False)

    else:
      objects = ['all',]
      labels  = ['All',]
      colors  = ['#80A080',]
      for i in range(N):
        objects.append(axisRegions[i])
        labels.append('%d' % (i+1))
        colors.append('#80A080')

      self.stripsToggleSelector.update(objects=objects,labels=labels,colors=colors)

      for i in range(N):
        axisRegion = axisRegions[i]
        if axisRegion.isActive:
          axisPanels = [panel for panel in activeWindowPane.sortedAxisPanels()[2:] if not panel.isDeleted]
          for axisPanel in axisPanels:
            axisRegions2 = axisPanel.sortedAxisRegions()

            if i < len(axisRegions2): # For safety during notified updates
              if axisPanel.axisType.isSampled:
                if self.activeWindowFrame:
                  self.activeWindowFrame.setPseudoState(axisPanel)
              else:
                (r0, r1) = Util.checkSwapRegion(axisRegions2[i].region, axisPanel.axisUnit)
                axisPanel.region_selector.setViewRegion(r0, r1)

          break

      M = 0
      I = 0
      for i in range(N):
        if axisRegions[i].isActive:
          I  = i
          M += 1

      if M == N:
        for i in range(N):
          self.stripsToggleSelector.selectButton(0, doCallback=False)

      else:
        self.stripsToggleSelector.selectButton(I+1, doCallback=False)

    self.waitStrips = 0

  def toggleSpectrum(self, spectrum):

    toggleSpectrum(self.window, spectrum=spectrum)

  def togglePeaks(self, winPeakList):

    for windowFrame in self.windowFrames:
      windowFrame.togglePeaks(winPeakList)

  def toggleStrips(self, object):

    axisRegions = self.activeWindowPane.findFirstAxisPanel(label=self.window.stripAxis).sortedAxisRegions()

    if object == 'all':
      for axisRegion in axisRegions:
        axisRegion.isActive = True
    else:
      for axisRegion in axisRegions:
        axisRegion.isActive = False
      
      if object in axisRegions:
        axisRegion = object
        
      axisRegion.isActive = True

      num = axisRegions.index(axisRegion)
      self.activeWindowFrame.updateRegionSelectors(num)

  def configTopBar(self, frame, buttonIndex):
  
    #x = 0
    #y = int(self.buttonList.winfo_x()) + int(self.buttonList.winfo_height())
    #w = int(self.scrolled_window.corner_canvas.winfo_x())
    
    buttons = self.buttonList.buttons
    for button in buttons:
      button.config(bg='#E6E6E6')
    
    for frameOther in (self.contourFrame, self.peakFrame,
                       self.stripFrame, self.specFrame):
                       
      if frameOther is not frame:
        frameOther.grid_forget()

    if self.activeButtonFrame is frame:
      self.activeButtonFrame = None
      frame.grid_forget()
    
    else:
      frame.grid(row=1, column=0, sticky='ew')
    
      #wf = int(frame.winfo_reqwidth())
      #if wf > w:
      #  frame.place(x=x, y=y, width=w)
      #else:
      #  frame.place(x=x, y=y, width=wf)
        
      self.activeButtonFrame = frame
      buttons[buttonIndex].config(bg='#D0B0A0')
          
  def selSpectra(self):

    self.configTopBar(self.specFrame, 0)
   
  def selContours(self):
  
    self.configTopBar(self.contourFrame, 1)

  def selPeaks(self):
  
    self.configTopBar(self.peakFrame, -4)
     
  def selStrips(self):
  
    self.configTopBar(self.stripFrame, -3)
  
  def selShapes(self):
  
    self.configTopBar(self.shapeFrame, -3)
  
  def gridAll(self):

    row = 2  # rows 0, 1 are for scrolledFrame
    for windowFrame in self.windowFrames:
      windowFrame.parent.grid_rowconfigure(row, weight=1)
      windowFrame.grid(row=row, column=0, columnspan=2, sticky='nsew')
      row += 1

  def destroy(self):

    self.isBeingDestroyed = True
    self.curateNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  # TBD: clean up
  def curateNotifiers(self, notify):

    notify(self.setWindowTitle, 'ccpnmr.Analysis.SpectrumWindow', 'setName')
    notify(self.setWindowLocation, 'ccpnmr.Analysis.SpectrumWindow', 'setLocation')
    notify(self.changedIsIconified, 'ccpnmr.Analysis.SpectrumWindow', 'setIsIconified')

    notify(self.setSpectrumSelector, 'ccpnmr.Analysis.SpectrumWindowView', 'setIsInToolbar')
    notify(self.setSpectrumSelector, 'ccp.nmr.Nmr.DataSource', 'setName')
    notify(self.setSpectrumSelector, 'ccp.nmr.Nmr.Experiment', 'setName')

    notify(self.setPeaksSelectorAfter, 'ccpnmr.Analysis.WindowPeakList', '__init__')
    notify(self.setPeaksSelectorAfter, 'ccpnmr.Analysis.WindowPeakList', 'delete')

    notify(self.setStripsSelectorAfter, 'ccpnmr.Analysis.AxisRegion', 'setIsActive')

  def setWindowTitle(self, window = None):

    if ((window == None) or (window == self.window)):
      self.setTitle(self.window.name)

  def setWindowLocation(self, window):

    if window is self.window:
      #print 'setWindowLocation1', self.geometry(), window.location
      location = '+%s+%s' % (window.location[0], window.location[1])
      self.geometry(location)
      #print 'setWindowLocation2', self.geometry()

  def getWindowLocation(self, *event):

    location = self.geometry()
    # below causes lots of redraws on OSX10.5 because motion events (so drawing
    # crosshair) also seems to come with configure events, which calls this function
    ###self.drawAllAfter()

    #print 'getWindowLocation1', location
    # below only so that geometry manager does not resize without being told to
    #w = int(location.split('x')[0])
    #if ((w > 1) and self.firstPass):
    #  self.geometry(location)
    #  self.firstPass = False

    # TBD: below supposedly does not give location on all platforms
    # (when called because of binding to 'Configure')
    n = location.find('+')
    xy = location[n+1:]
    (x, y) = map(int, xy.split('+'))
    #print 'getWindowLocation2', (x, y)
    if not self.window.isDeleted and not self.isBeingDestroyed:
      self.unregisterNotify(self.setWindowLocation, 'ccpnmr.Analysis.SpectrumWindow', 'setLocation')
      self.window.location = (x, y)
      self.registerNotify(self.setWindowLocation, 'ccpnmr.Analysis.SpectrumWindow', 'setLocation')
      # the above would have in turn called setWindowLocation() if (x, y) had changed
      # but there is no point in doing that and can cause windows to judder back to previous location

    # Deal with the inner option frame geometry
    
    """frame = self.activeButtonFrame
    if frame:      
      x = int(frame.winfo_x())
      y = int(frame.winfo_y())
      w = int(self.scrolled_window.corner_canvas.winfo_x())
      wf = int(frame.winfo_reqwidth())
      if wf > w:
        frame.place(x=x, y=y, width=w)
      else:
        frame.place(x=x, y=y, width=wf)"""

  def setLocationLabel(self, text):

    self.locationButton.config(text=text)

  def ncontoursIncr(self):

    views = self.getActiveSpectrumViews()
    for view in views:
      spectrum = view.analysisSpectrum.dataSource
      changeSpectrumNContours(spectrum, 1)

  def ncontoursDecr(self):

    views = self.getActiveSpectrumViews()
    for view in views:
      spectrum = view.analysisSpectrum.dataSource
      changeSpectrumNContours(spectrum, -1)

  def posNegContoursToggle(self):

    views = self.getActiveSpectrumViews()
    havePos = False
    haveNeg = False
    
    for view in views:
      if view.analysisSpectrum.posLevels:
        havePos = True
      if view.analysisSpectrum.negLevels:
        haveNeg = True
      
    if havePos and haveNeg:
      for view in views:
        view.analysisSpectrum.negLevels = []
       
    elif havePos:
      for view in views:
        posLevels = view.analysisSpectrum.posLevels
        view.analysisSpectrum.negLevels = [-x for x in posLevels]
        view.analysisSpectrum.posLevels = []
        
    elif haveNeg:
      for view in views:
        negLevels = view.analysisSpectrum.negLevels
        view.analysisSpectrum.posLevels = [-x for x in negLevels]

  def contourLevelsUp(self):

    views = self.getActiveSpectrumViews()
    for view in views:
      analysisSpectrum = view.analysisSpectrum
      levelChanger = analysisSpectrum.autoLevelChanger
      changeMode = analysisSpectrum.autoLevelMode

      if changeMode == 'add':
        if levelChanger <= 0:
          continue
      else:
        if levelChanger <= 1:
          continue

      Util.changeSpectrumContourLevels(analysisSpectrum, levelChanger, changeMode)

  def contourLevelsDown(self):

    views = self.getActiveSpectrumViews()
    for view in views:
      analysisSpectrum = view.analysisSpectrum
      levelChanger = analysisSpectrum.autoLevelChanger
      changeMode = analysisSpectrum.autoLevelMode

      if levelChanger > 0:
        if changeMode == 'add':
          levelChanger = - levelChanger
        else:
          levelChanger = 1.0 / levelChanger

        Util.changeSpectrumContourLevels(analysisSpectrum, levelChanger, changeMode)

  def generalContours(self):

    views = self.getActiveSpectrumViews()
    if (views):
      spectrum = views[0].analysisSpectrum.dataSource
    else:
      spectrum = None

    self.parent.editContourLevels(spectrum)

  def addStrip(self):

    if self.activeWindowFrame:
      self.activeWindowFrame.addStrip()

  def deleteStrip(self):

    if self.activeWindowFrame:
      self.activeWindowFrame.deleteStrip()

  def moveStripsLeft(self):

    if self.activeWindowFrame:
      self.activeWindowFrame.moveStripsLeft()

  def moveStripsRight(self):

    if self.activeWindowFrame:
      self.activeWindowFrame.moveStripsRight()

  def deleteVert(self):
  
    if self.window.stripAxis == 'x':
      self.deleteStrips()
    else:
      self.deleteSeparators()
  
  def deleteHoriz(self):
  
    if self.window.stripAxis == 'x':
      self.deleteSeparators()
      
    else:
      self.deleteStrips()            

  def deleteSeparators(self):
  
    if self.activeWindowFrame:
      self.activeWindowFrame.deleteSeparators()
      
  def deleteStrips(self):

    if self.activeWindowFrame:
      self.activeWindowFrame.deleteStrips()

  def getSpectrumViews(self):

    if self.activeWindowFrame:
      views = self.activeWindowFrame.getSpectrumViews()
      views.sort(self.activeWindowFrame.compareViewOrder)
    else:
      views = []

    return views

  def isViewVisible(self, view):

    return self.activeWindowFrame.isViewVisible(view)

  def getActiveSpectrumViews(self):

    if self.activeWindowFrame:
      views = self.activeWindowFrame.getActiveSpectrumViews()
    else:
      views = []

    return views

  def drawCrosshairs(self, typeLocation, originatingWindowPopup):

    for windowFrame in self.windowFrames:
      windowFrame.drawCrosshairs(typeLocation, originatingWindowPopup)

  def endCrosshair(self):

    for windowFrame in self.windowFrames:
      windowFrame.endCrosshair()

  def drawAllAfter(self, *extra):

    for windowFrame in self.windowFrames:
      windowFrame.drawAllAfter()

  def drawAll(self, *extra):

    for windowFrame in self.windowFrames:
      windowFrame.drawAll()

  #
  # turnDrawRequestsOff / turnDrawRequestsOn is for cases where you change a lot of
  # data model attributes but don't want a re-drawing to occur until at the end
  # turnDrawRequestsOff returns False if someone else has already turned drawing
  # requests off, otherwise it returns True
  # You should only call turnDrawRequestsOn if turnDrawRequestsOff returned True
  # Code should look like:
  #
  # turnedOff = self.turnDrawRequestsOff()
  # try:
  #   ... (lots of code with API model changes)
  # finally:
  #   if turnedOff: self.turnDrawRequestsOn()
  #

  def turnDrawRequestsOff(self):

    if self.waitDraw:
      return False

    self.waitDraw = True

    for windowFrame in self.windowFrames:
      windowFrame.turnDrawRequestsOff()

    return True

  def turnDrawRequestsOn(self, doDraw=True, doLift=False):

    self.waitDraw = False

    for windowFrame in self.windowFrames:
      windowFrame.turnDrawRequestsOn(doDraw)

    if doLift:
      self.lift()

