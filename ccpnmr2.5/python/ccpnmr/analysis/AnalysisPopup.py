
"""
======================COPYRIGHT/LICENSE START==========================

AnalysisPopup.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

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
import shutil
import sys
import traceback

import Tkinter

from ccp.general.Command import Command

from memops.api import Implementation as Impl

from memops.universal.Io import joinPath, getTopDirectory
from memops.universal.Util import isWindowsOS

from memops.general      import Implementation
from memops.general.Io   import saveProject
from memops.general.Util import isProjectModified

from memops.gui.Button          import Button
from memops.gui.DataEntry       import askString
from memops.gui.Label           import Label
from memops.gui.FileSelect      import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FontMenu        import FontMenu
from memops.gui.Menu            import Menu
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.WebBrowser      import WebBrowser

from memops.editor.ArchiveProjectPopup   import ArchiveProjectPopup
from memops.editor.OpenProjectPopup      import OpenProjectPopup
from memops.editor.SaveProjectPopup      import SaveProjectPopup
from memops.editor.WebBrowser            import ProjectWebBrowser
from memops.editor.WidgetCountPopup      import WidgetCountPopup

from ccp.api.nmr import Nmr

from ccp.gui.Io import loadProject
from ccp.gui.NmrExpPrototypeEditor import NmrExpPrototypePopup

from ccpnmr.analysis import Copyright

from ccpnmr.analysis.Analysis                        import Analysis, LOCAL_HELP_DOC_DIR

from ccpnmr.analysis.core.ExperimentBasic            import getSpectra, isSpectrum
from ccpnmr.analysis.core.Register                   import updateRegister, isRegistered
from ccpnmr.analysis.core.WindowBasic                import isActiveWindow, initWindowMacros, toggleSpectrum
from ccpnmr.analysis.core.WindowBasic                import togglePeakList, cloneSpectrumWindow, copySpectrumWindowViewProperties, copySpectrumWindowRegions
from ccpnmr.analysis.core import Util
from ccpnmr.analysis.core.ErrorHandler               import ErrorHandler

from ccpnmr.analysis.macros.ArgumentServer           import ArgumentServer


from ccpnmr.analysis.popups.AddContourFile           import AddContourFilePopup
from ccpnmr.analysis.popups.BackupProject            import BackupProjectPopup
from ccpnmr.analysis.popups.BasePopup                import BasePopup, determineHelpUrl
from ccpnmr.analysis.popups.BrowseAtoms              import BrowseAtomsPopup
from ccpnmr.analysis.popups.BrowseConstraints        import BrowseConstraintsPopup
from ccpnmr.analysis.popups.BrowsePeakGroups         import BrowsePeakGroupsPopup
from ccpnmr.analysis.popups.BrowseReferenceShifts    import BrowseReferenceShiftsPopup
from ccpnmr.analysis.popups.BrowseResonances         import BrowseResonancesPopup, ResonanceInfoPopup
from ccpnmr.analysis.popups.CalcHeteroNoe            import CalcHeteroNoePopup
from ccpnmr.analysis.popups.CalcHnHaCoupling         import CalcHnHaCouplingPopup
from ccpnmr.analysis.popups.CalcDistConstraints      import CalcDistConstraintsPopup
from ccpnmr.analysis.popups.CalcRates                import CalcRatesPopup
from ccpnmr.analysis.popups.CalcShiftDifference      import CalcShiftDifferencePopup
from ccpnmr.analysis.popups.ConfirmSeqSpinSystems    import ConfirmSeqSpinSystemsPopup
from ccpnmr.analysis.popups.CopyAssignments          import CopyAssignmentsPopup
from ccpnmr.analysis.popups.CreateContourFile        import CreateContourFilePopup
from ccpnmr.analysis.popups.EditAssignment           import EditAssignmentPopup
from ccpnmr.analysis.popups.EditAxisPanel            import EditAxisPanelPopup
from ccpnmr.analysis.popups.EditCalculation          import EditCalculationPopup
from ccpnmr.analysis.popups.EditContourFiles         import EditContourFilesPopup
from ccpnmr.analysis.popups.EditContourLevels        import EditContourLevelsPopup
from ccpnmr.analysis.popups.EditExperiment           import EditExperimentPopup
from ccpnmr.analysis.popups.EditExperimentSeries     import EditExperimentSeriesPopup
from ccpnmr.analysis.popups.EditMarks                import EditMarksPopup
from ccpnmr.analysis.popups.EditMeasurementLists     import EditMeasurementListsPopup
from ccpnmr.analysis.popups.EditMolLabelling         import EditMolLabellingPopup
from ccpnmr.analysis.popups.EditPeakAliasing         import EditPeakAliasingPopup
from ccpnmr.analysis.popups.EditPeakDrawParams       import EditPeakDrawParamsPopup
from ccpnmr.analysis.popups.EditPeakFindParams       import EditPeakFindParamsPopup
from ccpnmr.analysis.popups.EditPeak                 import EditPeakPopup
from ccpnmr.analysis.popups.EditPeakLists            import EditPeakListsPopup, SelectedPeaksPopup
from ccpnmr.analysis.popups.EditProfiles             import EditProfilesPopup
from ccpnmr.analysis.popups.EditMolecules            import EditMoleculesPopup
from ccpnmr.analysis.popups.EditSpectrum             import EditSpectrumPopup
from ccpnmr.analysis.popups.EditSpinSystem           import EditSpinSystemPopup
from ccpnmr.analysis.popups.EditStructures           import EditStructuresPopup
from ccpnmr.analysis.popups.EditWindow               import EditWindowPopup
from ccpnmr.analysis.popups.FollowShiftChanges       import FollowShiftChangesPopup
from ccpnmr.analysis.popups.InitRootAssignments      import InitRootAssignmentsPopup
from ccpnmr.analysis.popups.IsotopeSchemeEditor      import IsotopeSchemeEditor
from ccpnmr.analysis.popups.LinkPeakLists            import LinkPeakListsPopup
from ccpnmr.analysis.popups.LinkSeqSpinSystems       import LinkSeqSpinSystemsPopup
from ccpnmr.analysis.popups.LinkNoeResonances        import LinkNoeResonancesPopup
from ccpnmr.analysis.popups.MakeHbondRestraints      import MakeHbondRestraintsPopup
from ccpnmr.analysis.popups.NewWindow                import NewWindowPopup
from ccpnmr.analysis.popups.OpenSpectrum             import OpenSpectrumPopup
from ccpnmr.analysis.popups.PrintWindow              import PrintWindowPopup
from ccpnmr.analysis.popups.ProjectSummary           import ProjectSummaryPopup
from ccpnmr.analysis.popups.Register                 import RegisterPopup
from ccpnmr.analysis.popups.ReportError              import ReportErrorPopup  
from ccpnmr.analysis.popups.ResidueInfo              import ResidueInfoPopup
from ccpnmr.analysis.popups.SecStructureGraph        import SecStructureGraphPopup
from ccpnmr.analysis.popups.SecStructurePredict      import SecStructurePredictPopup
from ccpnmr.analysis.popups.SequenceShiftPredict     import SequenceShiftPredictPopup
from ccpnmr.analysis.popups.SpinSystemTyping         import SpinSystemTypeScoresPopup, SpinSystemTypingPopup
from ccpnmr.analysis.popups.ViewAssignment           import ViewAssignmentPopup
from ccpnmr.analysis.popups.ViewChemicalShifts       import ViewChemicalShiftsPopup
from ccpnmr.analysis.popups.ViewNoeMatrix            import ViewNoeMatrix
from ccpnmr.analysis.popups.ViewQualityReports       import ViewQualityReportsPopup
from ccpnmr.analysis.popups.ViewRamachandran         import ViewRamachandranPopup
from ccpnmr.analysis.popups.ViewStructure            import ViewStructurePopup 
from ccpnmr.analysis.popups.WindowPopup              import WindowPopup



from ccpnmr.format.converters.NmrStarFormat import NmrStarFormat
from ccpnmr.eci.ReadPdb import ReadPdb

from cambridge.dangle.DangleGui import DangleGui
from nijmegen.cing.CingPopup import CingPopup
from utrecht.haddock.HaddockPopup import HaddockPopup
from rutgers.rpf.PyRPF import PyRpfPopup

from ccpnmr.nexus.AutoBackbonePopup import AutoBackbonePopup
from paris.aria.AriaExtendNmrFrame import AriaPopup

# NB new
from gothenburg.prodecomp.ProdecompFrame import ProdecompPopup


try:
  import numpy
  HAVE_NUMPY = True
  
except:
  HAVE_NUMPY = False
  print ''
  print 'WARNING: Python NumPy module not installed or accessible.'
  print 'NumPy is required for CcpNmr automatic assignment and'
  print 'peak separator routines.'
  print ''


from ccpnmr.update.UpdatePopup import UpdatePopup

ProjectMenu    = 'Project'
ExperimentMenu = 'Experiment'
WindowMenu     = 'Window'
DataMenu       = 'Data Analysis'
OtherMenu      = 'Other'
PeaksMenu      = 'Peak'
MoleculeMenu   = 'Molecule'
MacroMenu      = 'Macro'
AssignMenu     = 'Assignment'
ResonanceMenu  = 'Resonance'
ChartMenu      = 'Chart'
StructureMenu  = 'Structure'

window_popup_prefix = 'window_'

after_delay = 10 # msec

DEFAULT_FONT = 'Helvetica 10'
GFX_DIR = os.path.join(getTopDirectory(),'python','memops','gui',
                       'graphics','16x16')
ICON_FILES = ['preferences-system','emblem-system',
              'help-browser','document-properties',
              'chart','view-refresh',
              'weather-overcast','window-new',
              'printer','font-x-generic',
              'system-log-out','document-save',
              'document-save-as','process-stop',
              'document-open','go-jump',
              'applications-system','folder-open']

class AnalysisPopup(BasePopup, Analysis):

  #def __init__(self, root, project = None, mem_cache = None):

  allowedCallbackFuncnames = ('init', 'save', 'close')

  def __init__(self, root, cache_size = 64, glDirect = None):

    cache_size = 1024 * 1024 * cache_size

    if glDirect is not None:
      self.glDirect = not not glDirect # makes this True or False rather than integer
    # leave out otherwise graphics windows will appear before this popup on toolbar
    #self.project = project

    self.popupFuncDict = {
      'open_project': self.openProject,
      'open_spectrum': self.openSpectrum,
      'edit_spectrum': self.editSpectrum,
      'add_contour_file': self.addSpectrumContourFile,
      'create_contour_file': self.createSpectrumContourFile,
      'edit_experiment': self.editExperiment,
      'edit_experiment_series': self.editExpSeries,
      'edit_experiment_prototypes': self.editExpPrototype,
      'shift_lists': self.editSampleConditionSets,
      'edit_measurement_lists': self.editMeasurementLists,
      'edit_peak': self.editPeak,
      'edit_peak_aliasing': self.editPeakAliasing,
      'selected_peaks': self.viewSelectedPeaks,
      'view_peaks': self.viewPeaks,
      'view_peak_groups': self.viewPeakGroups,
      'edit_peak_lists': self.editPeakLists,
      'copy_assignments': self.copyAssignments,
      'browse_atoms': self.browseAtoms,
      'browse_reference_shifts': self.browseReferenceShifts,
      'isotope_scheme_editor': self.isotopomerEditor,
      'isotope_labelling': self.editIsotopeLabelling,
      'browse_residue_info': self.browseResidueInfo,
      'browse_resonance_info': self.browseResonanceInfo,
      'browse_resonances': self.browseResonances,
      'selected_resonances': self.viewSelectedResonances,
      'edit_spin_system': self.editSpinSystems,
      'initialise_root_spectra': self.initialiseRootSpectra,
      'link_peaklists': self.linkPeakLists,
      'link_seq_spin_systems': self.linkSeqSpinSystems,
      'link_noe_resonances': self.linkNoeResonances,
      'confirm_seq_spin_systems': self.confirmSeqSpinSystems,
      'type_spin_systems': self.typeSpinSystems,
      'type_spin_system': self.typeSpinSystem,
      'edit_assignment': self.assignmentPanel,
      'quality_reports': self.qualityReports,
      'auto_backbone_assign': self.autoBackboneAssign,
      #'auto_backbone_assign': self.activateMars,
      'view_chem_shifts': self.chemShiftsTable,
      #'sequence_shift_predict': self.sequenceShiftPredict,
      'view_assignment': self.assignmentGraph,
      'sec_struc_graph': self.secStructureGraph,
      'view_noe_matrix': self.viewNoeMatrix,
      'edit_structures': self.editStructures,
      'plot_ramachandran': self.plotRamachandran,
      'view_structure': self.viewStructure,
      'calc_dist_constraints': self.calDistConstraints,
      'browse_restraints': self.browseConstraints,
      'make_hydrogen_bonds': self.makeHbonds,
      'dangle': self.startDangle,
      #'d2d': self.startD2D,
      'edit_peak_draw_params': self.peakDrawParams,
      'edit_peak_find_params': self.peakFindParams,
      'edit_molecules': self.editMolecules,
      'calc_rates': self.calcRates,
      'calc_shift_difference': self.calcShiftDifferences,
      'follow_shift_changes': self.followShiftChanges,
      'calc_hnha_coupling': self.calcHnHaCoupling,
      'calc_hetero_noe': self.calcHeteroNoe,
      'edit_calculation': self.editCalculation,
      'pales': self.pales,
      'blackledge_module': self.blackledge_module,
      'view_widget_count': self.viewWidgetCount,
      'backup_project': self.backupProject,
      'save_project': self.askSaveFile,
      'new_window': self.newWindow,
      'edit_window': self.editWindow,
      'print_window': self.printWindow,
      'edit_profiles': self.popupEditProfiles,
      'edit_axis_panel': self.editAxisPanel,
      'edit_marks': self.editMarks,
      'setup_clouds': self.setupClouds,
      'setup_bacus': self.setupBacus,
      'setup_midge': self.setupMidge,
      'setup_hcloudsmd': self.setupHcloudsMd,
      'setup_filter_clouds': self.setupFilterClouds,
      'setup_cloud_threader': self.setupCloudThreading,
      'setup_cloud_homologue': self.setupCloudHomologue,
      'run_format_converter': self.runFormatConverter,
      'entry_completion_interface': self.startECI,
      'edit_contour_levels': self.editContourLevels,
      'edit_contour_files': self.editContourFiles,
      'register_analysis': self.registerAnalysis,
      'project_summary': self.projectSummary,
    }

    # Needed for photoimage reference count/persistence
    self.icons = {}
    for iconFile in ICON_FILES:
      img = Tkinter.PhotoImage(file=os.path.join(GFX_DIR, iconFile + '.gif'))
      self.icons[iconFile] = img

    self.iconTool = self.icons['emblem-system']
    self.iconSpecialTool = self.icons['applications-system']
    self.iconPrefs = self.icons['preferences-system']
    self.iconHelp  = self.icons['help-browser']
    self.iconTable = self.icons['document-properties']
    self.iconChart = self.icons['chart']
    self.iconRefresh = self.icons['view-refresh']
    self.iconClouds = self.icons['weather-overcast']
    self.iconPrint = self.icons['printer']
    self.iconNewWindow = self.icons['window-new']
    self.iconFont = self.icons['font-x-generic']
    self.iconQuit = self.icons['system-log-out']
    self.iconSave = self.icons['document-save']
    self.iconSaveAs = self.icons['document-save-as']
    self.iconClose = self.icons['process-stop']
    self.iconOpen = self.icons['folder-open']
    self.iconOpenFile = self.icons['document-open']
    self.iconImport = self.icons['go-jump']

    self.printedRegistration = False

    version = Copyright.version
    versionData = [version.major, version.minor, version.release]
    self.program_name = ' '.join([Copyright.suite,
                                  Copyright.program,
                                  '.'.join([str(x) for x in versionData])])

    Analysis.__init__(self, cache_size)
    
    BasePopup.__init__(self, parent=root, title=self.program_name,
                       location='+100+100', class_=self.application.name,
                       quitFunc=self.quit, bg='grey90')
    self.showYesNo = showYesNo
    self.webBrowser = self.defaultWebBrowser = WebBrowser(self, popup=self)
    
    #Error handling section
    self.errorHandler = ErrorHandler()
    
    def show_error(*args):
      try:
        self._root().report_callback_exception_tk_native(*args)
        if self.analysisProfile.sendBugReports != 'no':
          formatedTb= traceback.format_exception(*args)
          if self.errorHandler.reportNeeded(formatedTb, *args):  
            popup = ReportErrorPopup(self, formatedTb, *args) 
      except:
         print 'Automated report failed' 
    
    self._root().report_callback_exception_tk_native = self._root().report_callback_exception
    self._root().report_callback_exception = show_error
    #END:Error handling section

  def after_idle(self, func):

    # ordinary after_idle does not seem to always wait until idle so hack it
    #BasePopup.after_idle(func)
    self.after(after_delay, func)
    
  def body(self, master):

    self.menus = {}
    self.menu_items = {}
    
    self.fixedActiveMenus = {}

    self.popups = {}

    self.callbacksDict = {}

    self.selected_objects = []

    self.menubar = Menu(master)

    self.font = DEFAULT_FONT
    
    self.setProjectMenu()
    self.initProject()
    self.setPeaksMenu()
    self.setMoleculeMenu()
    self.setAssignMenu()
    self.setResonanceMenu()
    self.setDataMenu()
    self.setStructureMenu()
    self.setChartMenu()
    self.setMacroMenu()
    self.setOtherMenu()
    self.setMenuState() # need to do it again because of OtherMenu state

    self.config(menu=self.menubar)
    
    master.grid_rowconfigure(0, weight=1)
    master.grid_columnconfigure(0, weight=1)
    
    c = Tkinter.Canvas(master)
    item = c.create_text(0,0, text='   '.join([k for k in self.menus]) + 'A'*6 , font=self.font)
    bbox = c.bbox(item)    
    del c   
    
    self.geometry('%dx1' % min(int(self.winfo_screenheight()), bbox[2]-bbox[0]))


  def changeFont(self, analysisProfile):
  
    if analysisProfile is self.analysisProfile:
      if self.font != analysisProfile.font:
        self.font = analysisProfile.font or DEFAULT_FONT 
        self.setFont()
      
  def selectFont(self, font): # Only called from main menu option
    
    if self.project:
      self.analysisProfile.font = font    

  def setFont(self, font=None):

    if font is None:
      font = self.font or DEFAULT_FONT 
    
    BasePopup.setFont(self, font)


  def drawWindows(self):

    windows = self.getActiveWindows()
    for window in windows:
      if not window.isIconified:
        popup = self.getWindowPopup(window.name)
        
        if popup:
          popup.drawAll()

  def isCrosshairVisible(self):

    if self.project:
      isVisible = self.analysisProfile.useCrosshair
    else:
      isVisible = True

    return isVisible

  def drawCrosshairs(self, typeLocation, originatingWindowFrame=None):

    if not self.isCrosshairVisible():
      return

    popup = self.popups.get('browse_reference_shifts')
    if popup and popup.state() == 'normal':
      popup.drawCrosshairs(typeLocation)

    self.drawWindowCrosshairs(typeLocation, originatingWindowFrame)

  def drawWindowCrosshairs(self, typeLocation, originatingWindowFrame=None):

    windows = self.getActiveWindows()
    for window in windows:
      popup = self.getWindowPopup(window.name)
      if popup and popup.state() == 'normal': # be safe in case called while window being constructed
        popup.drawCrosshairs(typeLocation, originatingWindowFrame)

  def endCrosshair(self):

    if not self.isCrosshairVisible():
      return

    windows = self.getActiveWindows()
    for window in windows:
      popup = self.getWindowPopup(window.name)
      
      if popup:
        popup.endCrosshair()

  def startSelection(self):

    self.selection_changed = False

  def endSelection(self, redraw = False):

    if redraw or self.selection_changed:
      self.selection_changed = False
      self.drawWindows()
      
    #self.updateAfter()

  def updateSelectionPopups(self, peaks):

    popup = self.popups.get('selected_peaks')
    if popup:
      popup.update(peaks=self.currentPeaks)
     
    peak = self.currentPeak
    popup = self.popups.get('edit_assignment')
    if peak and popup:
      self.after_idle(lambda: popup.update(peak))
    
  def clearSelected(self):
 
    if self.selected_objects:
      for object in self.selected_objects:
        if isinstance(object, Nmr.Peak):
          try:
            object.cPeak.setIsSelected(False)
          except:
            pass
 
      self.currentPeaks = []
      self.currentPeak  = None
      self.updateSelectionPopups(self.currentPeaks)
      self.selected_objects = []
      self.selection_changed = True

  def isSelected(self, object):

    return object in self.selected_objects

  def setCurrentPeaks(self, peaks):

    self.clearSelected()
    for peak in peaks:
      self.selected_objects.append(peak)
      self.currentPeaks.append(peak)
      if self.currentPeak not in self.currentPeaks:
        self.currentPeak = self.currentPeaks[0]

      peak.cPeak.setIsSelected(True)
        
    self.updateSelectionPopups(self.currentPeaks)
    self.selection_changed = True
    self.drawWindowPopups()

  def addSelected(self, object):

    if not self.isSelected(object):
      self.selected_objects.append(object)

      if isinstance(object, Nmr.Peak):
        self.currentPeaks.append(object)
        if self.currentPeak not in self.currentPeaks:
          self.currentPeak = self.currentPeaks[0]

        object.cPeak.setIsSelected(True)
        
      self.updateSelectionPopups(self.currentPeaks)
      self.selection_changed = True

  def removeSelected(self, object):

    if self.isSelected(object):
      self.selected_objects.remove(object)

      if isinstance(object, Nmr.Peak):
        self.currentPeaks.remove(object)
        if object is self.currentPeak:
          self.currentPeak = None

        object.cPeak.setIsSelected(False)
      
      self.updateSelectionPopups(self.currentPeaks)
      self.selection_changed = True

  def deleteSelected(self):

    if self.selected_objects:
      if self.currentPeak in self.currentPeaks:
        self.currentPeak = None
        
      for object in self.selected_objects:
        #print 'deleteSelected', object.serial
        object.delete()
        # automatically redrawn by notify

      # Hopefully blank already from C peak notifier
      self.currentPeaks = []
      self.updateSelectionPopups(self.currentPeaks)
      self.selected_objects = []

  def queryDeleteSelected(self, originator = None):

    if not originator:
      originator = self

    if self.selected_objects:
      n = len(self.selected_objects)
      
      if (n == 1):
        s1 = s2 = ''
      else:
        s1 = 's'
        s2 = str(n) + ' '

      className = self.selected_objects[0].className
        
      if showYesNo('Delete object%s' % s1,
                   'Delete %sselected %s%s?' % (s2, className, s1), parent=originator):
        self.deleteSelected()

  def deletePeaks(self, peaks):

    peaks = set(peaks)
    self.selected_objects = list(set(self.selected_objects) - peaks)
    # the below should not be done because the C deletion notifier does this
    #self.currentPeaks = list(set(self.currentPeaks) - peaks)
    #if self.currentPeak in peaks:
    #  self.currentPeak = None
      
    for peak in peaks:
      peak.delete()
      
    self.updateSelectionPopups(self.currentPeaks)
    
  def queryDeletePeaks(self, peaks, originator=None):

    if not originator:
      originator = self

    if peaks:
      n = len(peaks)
      
      if n == 1:
        s1 = s2 = ''
      else:
        s1 = 's'
        s2 = str(n) + ' '
        
      if showYesNo('Delete Peak%s' % s1,
                   'Delete %sPeak%s?' % (s2, s1), parent=originator):
        self.deletePeaks(peaks)

  def toggleSpectrum(self, window, shortcut=None, spectrum=None):

    assert shortcut or spectrum
    
    isGlobal = 0
    if self.analysisProfile.useGlobalShortcuts:
      isGlobal = 1

    return toggleSpectrum(window, shortcut=shortcut, spectrum=spectrum, isGlobal=isGlobal)

  def togglePeakList(self, window, peakList):
    
    togglePeakList(window, peakList)

  def curatePopupNotifiers(self, notify):
  
    notify(self.changeFont, 'ccpnmr.AnalysisProfile.AnalysisProfile', 'setFont')
    notify(self.setMacroMenu, 'ccpnmr.AnalysisProfile.Macro', 'delete')
    notify(self.setMacroMenu, 'ccpnmr.AnalysisProfile.Macro', 'setName')
    notify(self.changedSpectrum, 'ccp.nmr.Nmr.DataSource', '__init__')
    notify(self.deletedSpectrum, 'ccp.nmr.Nmr.DataSource', 'delete')
    notify(self.initSpectrumWindow, 'ccpnmr.Analysis.SpectrumWindow', '__init__')
    notify(self.deleteSpectrumWindow, 'ccpnmr.Analysis.SpectrumWindow', 'delete')
    notify(self.changedSpectrumWindowName, 'ccpnmr.Analysis.SpectrumWindow', 'setName')
    notify(self.deleteSpectrumWindowGroup, 'ccpnmr.Analysis.SpectrumWindowGroup', 'delete')
    notify(self.changedWindowGroupWindows, 'ccpnmr.Analysis.SpectrumWindowGroup', 'addSpectrumWindow')
    notify(self.changedWindowGroupWindows, 'ccpnmr.Analysis.SpectrumWindowGroup', 'removeSpectrumWindow')
    notify(self.changedWindowGroupWindows, 'ccpnmr.Analysis.SpectrumWindow', 'addSpectrumWindowGroup')
    notify(self.changedWindowGroupWindows, 'ccpnmr.Analysis.SpectrumWindow', 'removeSpectrumWindowGroup')
    notify(self.changedNumPoints, 'ccp.nmr.Nmr.FreqDataDim', 'setNumPoints')
  
  def setProjectMenu(self):

    """
    from ccpnmr.update.UpdateAgent import getNumUninstalledUpdates
    numUpdates = getNumUninstalledUpdates()
    """
    numUpdates = None  # TBD

    if numUpdates:
      updateText = 'Updates * %d available *' % numUpdates
    else:
      updateText = 'Updates'

    # Imports submenu
    
    
    importsMenu = Menu(self.menubar, tearoff=False)
    importsMenu.add_command(label='Via Format Converter', shortcut='F', command=self.runFormatConverter)
    importsMenu.add_command(label='NMR-STAR 2.1.1', command=self.importNmrStar211 )
    importsMenu.add_command(label='NMR-STAR 3.1', shortcut='N', command=self.importNmrStar31 )
    importsMenu.add_command(label='PDB 3.20', shortcut='P',command=self.importPdb )
    importsMenu.add_command(label='Coordinates (PDB-style)', shortcut='C',command=self.importCoordinates )

    # # NOTE:ED not needed
    # importsMenu.add_command(label='Nef', command=self.importNefFile)
    # exportsMenu = Menu(self.menubar, tearoff=False)
    # exportsMenu.add_command(label='Nef', command=self.exportNefFile)

    # Preferences submenu

    fontsMenu = FontMenu(self.menubar, self.selectFont, sizes=(8,10,12,14),
                         doItalic=False, doBoldItalic=False, tearoff=0)
    
    prefsMenu = Menu(self.menubar, tearoff=False)
    prefsMenu.add_cascade(label='Fonts', shortcut='F', 
                          image=self.iconFont, compound='left',
                          menu=fontsMenu,
                          tipText='Select font to use in the graphical interface')
    prefsMenu.add_command(label='Colour Schemes', 
                          image=self.iconTable, compound='left',
                          shortcut='C',
                          tipText='Edit and create colour schemes',
                          command=self.editColorSchemes)
    prefsMenu.add_command(label='Residue Codes',  
                          image=self.iconTable, compound='left',
                          shortcut='R', tipText='User-specified codes that override the standard residue names',
                          command=self.editResidueCodes)
    prefsMenu.add_command(label='User Options',  
                          image=self.iconTable, compound='left',
                          shortcut='U', tipText='General options for Analysis behaviour',
                          command=self.editProfiles)
    
    # Help Submenu
    
    helpMenu = Menu(self.menubar, tearoff=0)
    helpMenu.add_command(label='Version', shortcut='V',
                         command=self.showVersion)
    helpMenu.add_command(label='About',   shortcut='A',
                         command=self.showAbout)
    helpMenu.add_command(label='Help',    shortcut='H',
                         command=self.showHelp)

    # 

    menu = Menu(self.menubar, tearoff=0)
    menu.add_command(label='New', shortcut='N',   
                     image=self.iconNewWindow, compound='left',
                     command=self.newProject,
                     tipText='Create a new, blank CCPN project (closes any open project)')
    menu.add_command(label='Open Project', shortcut='O',   
                     image=self.iconOpen, compound='left',
                     command=self.openProject,
                     tipText='Open a new CCPN project by selecting a project directory on disk')
    menu.add_command(label='Open Spectra', shortcut='t',   
                     image=self.iconOpenFile, compound='left',
                     command=self.openSpectrum,
                     tipText='Open spectrum data from disk, creating a default CCPN project if needed')

    menu.add_command(label='Load Nef',
                     image=self.iconOpenFile, compound='left',
                     command=self.loadNefFile,
                     tipText='Open Nef file from disk, creating a default CCPN project if needed')

    menu.add_command(label='Save', shortcut='S',
                     image=self.iconSave, compound='left',
                     command=self.saveProject,
                     tipText='Save the current CCPN project on disk')
    menu.add_command(label='Save As', shortcut='A',   
                     image=self.iconSaveAs, compound='left',
                     command=self.saveAsProject,
                     tipText='Save the current CCPN project under a different name (project directory)')

    menu.add_cascade(label='Import', shortcut='I',
                     image=self.iconImport, compound='left',
                     menu=importsMenu)

    # # NOTE:ED not needed
    # menu.add_cascade(label='Export',
    #                  image=self.iconImport, compound='left',
    #                  menu=exportsMenu)

    menu.add_command(label='Close', shortcut='C',
                     image=self.iconClose, compound='left',
                     command=self.closeProject,
                     tipText='Close the current CCPN project')
    menu.add_command(label='Quit', shortcut='Q',   
                     image=self.iconQuit, compound='left',
                     command=self.quit,
                     tipText='Quit Analysis, closing any open CCPN project')
    menu.add_separator()
    menu.add_command(label='Summary', shortcut='m',
                     image=self.iconTable, compound='left',
                     command=self.projectSummary,
                     tipText='Provide executive summary of project')
    menu.add_cascade(label='Preferences', shortcut='P',
                     image=self.iconPrefs, compound='left',
                     menu=prefsMenu)
    menu.add_command(label='Register', shortcut='R',
                     image=self.iconTool, compound='left',
                     command=self.registerAnalysis,
                     tipText='Register yourself with CCPN')
    menu.add_command(label='Validate',shortcut='V',
                     image=self.iconTool, compound='left',
                     command=self.validateProject,
                     tipText='Check the current CCPN project for data model consistency errors')
    menu.add_command(label='Backup', shortcut='B',
                     image=self.iconTool, compound='left',
                     command=self.backupProject,
                     tipText='Setup options for automated backup of CCPN project data')
    menu.add_command(label='Archive', shortcut='v',
                     image=self.iconTool, compound='left',
                     command=self.archiveProject,
                     tipText='Save the current CCPN project in an archived form, e.g. tar gzipped')
    menu.add_command(label=updateText, shortcut='U',
                     image=self.iconRefresh, compound='left',
                     command=self.updateAnalysis,
                     tipText='Get any new patches and updates to the Analysis program')
    menu.add_separator()
    menu.add_cascade(label='Help', shortcut='H',
                     image=self.iconHelp, compound='left',
                     menu=helpMenu)
    
    self.menubar.add_cascade(label=ProjectMenu, shortcut='j', menu=menu)
    self.menus[ProjectMenu] = menu
    self.menu_items[ProjectMenu] = ['New', 'Open Project', 'Open Spectra', 'Load Nef',
                                    'Save', 'Save As', 'Import', 'Close',
				                    'Quit', 'Summary', 'Preferences', 'Register', 'Validate',
                                    'Backup', 'Archive', updateText, 'Help']
    
    # Menus that area active in absence of a project
    # for ii in (0,1,2,7,15,17):
    for ii in (0,1,2,3,8,16,18):
      self.fixedActiveMenus[(ProjectMenu,ii)] = True

  def openWindowGroup(self, spectrumWindowGroup=None):

    if spectrumWindowGroup:
      windows = spectrumWindowGroup.spectrumWindows
    else:
      windows = []
      
    self.destroyWindows(windows)
    self.analysisProject.activeWindowGroup = spectrumWindowGroup
    self.openActiveWindows()
    self.setWindowMenu()

  def destroyWindows(self, windows):

    for key in self.popups.keys():
      if key.startswith(window_popup_prefix): # bit of a hack
        popup = self.popups[key]
        
        if popup.window not in windows:
          del self.popups[key]
          popup.destroy()

  def deleteSpectrumWindowGroup(self, spectrumWindowGroup):

    analysisProject = self.project.currentAnalysisProject
    group = analysisProject.activeWindowGroup
    if not group and spectrumWindowGroup not in analysisProject.spectrumWindowGroups:
      self.openWindowGroup()
      
  def changedWindowGroupWindows(self, *extra):

    windows = self.getActiveWindows()
    for key in self.popups.keys():
      if key.startswith(window_popup_prefix): # bit of a hack
        popup = self.popups[key]
        
        if popup.window not in windows:
          popup = self.popups[key]
          del self.popups[key]
          popup.destroy()
    
    #print 'changedWindowGroupWindows', [window.name for window in windows]
    for window in windows:
      if not window.isIconified:
        popup = self.getWindowPopup(window.name)
        if not popup:
          # can be in the middle of creating window, so use after_idle
          # after_idle below does not seem to work (for some unknown reason):
          # when the EditWindowPopup is open then openWindow executes immediately
          self.after_idle(lambda window=window: self.openWindow(window))
          #self.after(after_delay, lambda window=window: self.openWindow(window))
          # put back to after_idle now that all after_idles use after_delay

    self.setWindowMenu()

  def changedNumPoints(self, dataDim):

    dataSource = dataDim.dataSource
    
    for peakList in dataSource.peakLists:
      if peakList.peaks:
        break
    else:
      return

    msg = 'Changed number of points for dimension %d in spectrum %s %s,'
    msg += ' peaks will no longer be in correct location'
    data = (dataDim.dim, dataSource.experiment.name, dataSource.name)
    showWarning('Changed number of points', msg % data, parent=self)

  def openActiveWindows(self):

    windows = self.getActiveWindows()
    for window in windows:
      if not window.isIconified:
        self.openWindow(window)

  def openActivePopups(self):

    popup_names = self.getPopupsOpen()
    for popup_name in popup_names:
      func = self.popupFuncDict.get(popup_name)
      if func:
        func()

  def drawWindowPopups(self):

    popups = self.getWindowPopups()
    for popup in popups:
      popup.turnDrawRequestsOn()
      popup.drawAll()

  def changedSpectrum(self, dataSource):

    if not isSpectrum(dataSource):
      return

    self.setExperimentMenu()

  def deletedSpectrum(self, dataSource):

    #Analysis.deletedSpectrum(self, dataSource)
    self.changedSpectrum(dataSource)

    #self.after_idle(lambda: self.checkExperiment(dataSource.experiment))

  def setExperimentMenu(self):

    # below was done when idea was to list all spectra, not done currently
    menu = self.menus.get(ExperimentMenu)
    if menu:
      menu.delete(0, Tkinter.END)
      
    else:
      menu = Menu(self.menubar, tearoff=0)
      self.menubar.add_cascade(label=ExperimentMenu, shortcut='E', menu=menu)
      self.menus[ExperimentMenu] = menu

    # TBD
    #
    # EditType should be modal callable.
    #

    menu.add_command(label='Open Spectra' , shortcut='O', 
                     image=self.iconOpenFile, compound='left',
                     command=self.openSpectrum,
                     tipText='Open spectrum data files from disk')
    menu.add_separator()
    menu.add_command(label='Spectra', shortcut='S',  
                     image=self.iconTable, compound='left',
                     command=self.editSpectrum,
                     tipText='Spectrum information including display options, referencing & file details')
    menu.add_command(label='Experiments', shortcut='E',   
                     image=self.iconTable, compound='left',
                     command=self.editExperiment,
                     tipText='Experiment information including type and working shift list')
    menu.add_command(label='NMR Series', shortcut='N',   
                     image=self.iconTable, compound='left',
                     command=self.editExpSeries,
                     tipText='Setup of experiment series for relaxation and titration analysis etc.')
    menu.add_command(label='Experiment Prototypes',shortcut='P',  
                     image=self.iconTable, compound='left',
                     command=self.editExpPrototype,
                     tipText='Curate and manage the reference experiment types')

    menu_items = self.menu_items[ExperimentMenu] = ['Open Spectra', 'Spectra','Experiments',
                                                    'NMR Series','Experiment Prototypes']

    self.fixedActiveMenus[(ExperimentMenu,0)] = True

  def setPeaksMenu(self):

    # TBD
    #
    # Extra Peak Table
    # 
    # Duplicate popup option - view several peak lists
    #
    # Move Copy Peak Assignments to Assignments menu (a merged Copy Assignments)

    menu = Menu(self.menubar, tearoff=0)
    menu.add_command(label='Peak Lists', shortcut='P',   
                     image=self.iconTable, compound='left',
                     command=self.editPeakLists,
                     tipText='Tables of peaks and peak lists, including display options and synthetic peaks')
    menu.add_separator()
    menu.add_command(label='Selected Peaks',  shortcut='S',   
                     image=self.iconTool, compound='left',
                     command=self.viewSelectedPeaks,
                     tipText='A table of the peaks currently selected in spectrum windows')
    ## dans menu option
    #menu.add_command(label='Peak Separator', 
    #                 image=self.iconTool, compound='left',
    #                 command=self.peakSeparatorParams,
    #                 tipText='Separate merged peaks using peak models')
    menu.add_separator()
    menu.add_command(label='Peak Finding',    shortcut='F',  
                     image=self.iconPrefs, compound='left',
                     command=self.peakFindParams,
                     tipText='Options to determine how peak extrema are found during picking')
    menu.add_command(label='Draw Parameters', shortcut='D',   
                     image=self.iconPrefs, compound='left',
                     command=self.peakDrawParams,
                     tipText='Options on how to render peak crosses and peak/resonance annotations')

    self.menubar.add_cascade(label=PeaksMenu, shortcut='P', menu=menu)

    self.menus[PeaksMenu] = menu
    self.menu_items[PeaksMenu] = ['Peak Lists', 'Selected Peaks',
                                  'Peak Separator', 'Peak Finding', 'Draw Parameters']

  def setMoleculeMenu(self):

    # TBD
    #
    # Unified Molecule option - Systems and Templates - reinstate covalent links?
    # Isotopomer scheme tidy
    

    menu = Menu(self.menubar, tearoff=0)
    menu.add_command(label='Molecules', shortcut='M',   
                     image=self.iconTable, compound='left',
                     command=self.editMolSystems,
                     tipText='Setup the polymer chains, small molecules and sequences within the project')
    menu.add_command(label='Isotope Labelling',  shortcut='L',   
                     image=self.iconTable, compound='left',
                     command=self.editIsotopeLabelling,
                     tipText='Setup the isotope labelling patterns of molecules')
    menu.add_command(label='Reference Isotope Schemes',  shortcut='I',   
                     image=self.iconTable, compound='left',
                     command=self.isotopomerEditor,
                     tipText='Curate and manage per-residue isotope labelling schemes')
    menu.add_separator()
    menu.add_command(label='Atom Browser',   shortcut='A', 
                     image=self.iconTool, compound='left',
                     command=self.browseAtoms,
                     tipText='A table of all the atoms that are available for assignment')
    menu.add_command(label='Add Sequence',   shortcut='S', 
                     image=self.iconTool, compound='left',
                     command=self.addSequence,
                     tipText='Quickly add a new polymer sequence and assignment chain')
    menu.add_separator()
    menu.add_command(label='Residue Information', shortcut='R',    
                     image=self.iconChart, compound='left',
                     command=self.browseResidueInfo,
                     tipText='Graphical display of assignment status of atoms within residues')
    self.menubar.add_cascade(label=MoleculeMenu, shortcut='l', menu=menu)

    self.menus[MoleculeMenu] = menu
    self.menu_items[MoleculeMenu] = ['Molecules', 'Isotope Labelling',
                                     'Reference Isotope Schemes', 
                                     'Atom Browser', 'Add Sequence', 
                                     'Residue Information']

    # NOTE:ED - works now :)
    # if isWindowsOS():
    #   # Always inactivate 'Residue Information' under Windows (does not work)
    #   ind = self.menu_items[MoleculeMenu].index('Residue Information')
    #   # ind+2 because of two separators
    #   self.fixedActiveMenus[(MoleculeMenu, ind+2)] = False

  def setAssignMenu(self):
 
    # TBD
    #
    # Merge Initialise & Link?
    #auto_menu = Menu(self.menubar, tearoff=False)
    #menu.add_cascade(label='Automation',  shortcut='u', menu=auto_menu)
    #auto_menu.add_command(label='Setup Spin Systems',  shortcut='S', command=self.setupSeqSpinSystems)
    
    # Refactored by Rasmus 5/4/10 to ensure that menu_items remains in sync.
    # Too much work to extend to all the other cases.
    
    menu = Menu(self.menubar, tearoff=0)
    
    self.menu_items[AssignMenu] = menuNames = []
    
    label='Assignment Panel'
    menu.add_command(label=label, image=self.iconTool, compound='left',
                     shortcut='A', command=self.assignmentPanel,
                     tipText='The popup dialog used to edit peak dimension assignments')
    menuNames.append(label)
    label='Copy Assignments'
    menu.add_command(label=label, shortcut='C',
                     image=self.iconTool, compound='left',
                     command=self.copyAssignments,
                     tipText='Tools to copy assignments between molecules and peak lists')
    menuNames.append(label)
    label='Spin System Typing'
    menu.add_command(label=label, shortcut='T',
                     image=self.iconTool, compound='left',
                     command=self.typeSpinSystems,
                     tipText='Tools to predict the residue types of resonance spin systems')
    menuNames.append(label)
    menu.add_separator()
    label='Initialise Root Resonances'
    menu.add_command(label=label, shortcut='I',
                     image=self.iconSpecialTool, compound='left',
                     command=self.initialiseRootSpectra,
                     tipText='Add initial resonance and spin system numbers to HQSC & HNCO peaks')
    menuNames.append(label)
    label='Pick & Assign From Roots'
    menu.add_command(label=label, shortcut='P',
                     image=self.iconSpecialTool, compound='left',
                     command=self.linkPeakLists,
                     tipText='Pick and assign peaks based on HSQC or HNCO peak positions')
    menuNames.append(label)
    label='Protein Sequence Assignment'
    menu.add_command(label=label, shortcut='S',
                     image=self.iconSpecialTool, compound='left',
                     command=self.linkSeqSpinSystems,
                     tipText='Find, connect and assign sequentially related spin systems using peak matching')
    menuNames.append(label)
    label='Automated Seq. Assignment'
    menu.add_command(label=label, shortcut='u',
                     image=self.iconSpecialTool, compound='left',
                     command=self.autoBackboneAssign,
                     tipText='Automatic protein sequence assignment')
    menuNames.append(label)
    label='NOE Contributions'
    menu.add_command(label=label, shortcut='N', image=self.iconSpecialTool, compound='left',
                     command=self.linkNoeResonances,
                     tipText='Assign the contributions to NOE peaks using chemical shits and structure distances')
    menuNames.append(label)
    menu.add_separator()
    label='Assignment Graph'
    menu.add_command(label=label,  shortcut='G',     
                     image=self.iconChart, compound='left',
                     command=self.assignmentGraph,
                     tipText='A graphical display of the assigned atom groups and their connections')
    menuNames.append(label)
    label='Quality Reports'
    menu.add_command(label=label, shortcut='Q',     
                     image=self.iconChart, compound='left',
                     command=self.qualityReports,
                     tipText='Show reports on assignment completeness, NOE links, assignment errors etc.')
    menuNames.append(label)
    
    
    self.menubar.add_cascade(label=AssignMenu, shortcut='A', menu=menu)

    self.menus[AssignMenu] = menu


  def setResonanceMenu(self):

    # TBD Check for tabbification

    menu = Menu(self.menubar, tearoff=0)
    #menu.add_command(label='Selected Resonances', shortcut='e',   
    #                 image=self.iconTable, compound='left',
    #                 command=self.viewSelectedResonances)
    menu.add_command(label='Resonances',   shortcut='R',   
                     image=self.iconTable, compound='left',
                     command=self.browseResonances,
                     tipText='A table of resonances and their chemical shifts, assignment peak links etc.')
    menu.add_command(label='Spin Systems', shortcut='S',   
                     image=self.iconTable, compound='left',
                     command=self.editSpinSystems,
                     tipText='Tables of spin systems showing the resonances and shifts they contain')
    menu.add_command(label='Measurement Lists', shortcut='M',   
                     image=self.iconTable, compound='left',
                     command=self.editMeasurementLists,
                     tipText='The NMR measurements including chemical shifts, T1, T2 etc.')
    menu.add_separator()
    menu.add_command(label='Reference Chemical Shifts', shortcut='C',     
                     image=self.iconChart, compound='left',
                     command=self.browseReferenceShifts,
                     tipText='Tables and graphs of reference chemical shift information from BMRB & RefDB')
    menu.add_command(label='Chemical Shifts Table', shortcut='T',     
                     image=self.iconChart, compound='left',
                     command=self.chemShiftsTable,
                     tipText='A neat graphical representation of chemical shifts in tabular form')
    #menu.add_separator()
    #menu.add_command(label='CamCoil: Shift Predictions', shortcut='P', 
    #                 image=self.iconSpecialTool, compound='left',
    #                 command=self.sequenceShiftPredict,
    #                 tipText='Use CamCoil to get shift predictions from sequence')
    self.menubar.add_cascade(label=ResonanceMenu, shortcut='R', menu=menu)

    self.menus[ResonanceMenu] = menu
    self.menu_items[ResonanceMenu] = ['Resonances',
                                      'Spin Systems',
                                      'Measurement Lists', 
                                      'Reference Chemical Shifts',
                                      'Chemical Shifts Table']

  def setChartMenu(self):

    # TBD Check for tabbification

    menu = Menu(self.menubar, tearoff=0)
    menu.add_command(label='Assignment Graph',      shortcut='A',     
                     image=self.iconChart, compound='left',
                     command=self.assignmentGraph,
                     tipText='A graphical display of the assigned atom groups and their connections')
    menu.add_command(label='Chemical Shifts Table', shortcut='C',     
                     image=self.iconChart, compound='left',
                     command=self.chemShiftsTable,
                     tipText='A neat graphical representation of chemical shifts in tabular form')
    menu.add_command(label='Reference Chemical Shifts', shortcut='R',     
                     image=self.iconChart, compound='left',
                     command=self.browseReferenceShifts,
                     tipText='Tables and graphs of reference chemical shift information from BMRB & RefDB')
    menu.add_command(label='Secondary Structure Chart', shortcut='S',     
                     image=self.iconChart, compound='left',
                     command=self.secStructureGraph,
                     tipText='Tables and charts giving evidence for protein secondary structure')
    menu.add_command(label='Residue Interaction Matrix',shortcut='I',     
                     image=self.iconChart, compound='left',
                     command=self.viewNoeMatrix,
                     tipText='A density matrix plot of residue-residue interactions from peaks & restraints')
    menu.add_command(label='Ramachandran Plot',     shortcut='P',     
                     image=self.iconChart, compound='left',
                     command=self.plotRamachandran,
                     tipText='A plot of phi & psi protein backbone angles in structures')
    self.menubar.add_cascade(label=ChartMenu, shortcut='C', menu=menu)

    self.menus[ChartMenu] = menu
    self.menu_items[ChartMenu] = ['Assignment Graph',
                                  'Chemical Shifts Table',
                                  'Reference Chemical Shifts', 
                                  'Secondary Structure Chart',
                                  'Residue Interaction Matrix',
                                  'Ramachandran Plot']
  
  def setDataMenu(self):

    # TBD
    #
    # Merge Measurements and lists
    # Tabbify other options as required
    #

    menu = Menu(self.menubar, tearoff=0)
    menu.add_command(label='Measurement Lists',   shortcut='M',    
                     image=self.iconTable, compound='left',
                     command=self.editMeasurementLists,
                     tipText='The NMR measurements including chemical shifts, T1, T2 etc.')
    menu.add_command(label='NMR Series',     shortcut='N',    
                     image=self.iconTable, compound='left',
                     command=self.editExpSeries,                     
                     tipText='Setup of experiment series for relaxation and titration analysis etc.')
    menu.add_separator()
    menu.add_command(label='Shift Differences',   shortcut='D', 
                     image=self.iconTool, compound='left',
                     command=self.calcShiftDifferences,
                     tipText='Tools to calculate chemical shift differences between peak lists & shift lists')
    menu.add_command(label='Heteronuclear NOE',   shortcut='H', 
                     image=self.iconTool, compound='left',
                     command=self.calcHeteroNoe,
                     tipText='A tool to quickly calculate heteronuclear NOE values by peak intensity comparison')
    menu.add_command(label=u'3J H-H\u03B1 Coupling', shortcut='C', 
                     image=self.iconTool, compound='left',
                     command=self.calcHnHaCoupling,
                     tipText='A tool to extract amide H to alpha H 3J coupling and predict phi angles using HNHA experiments')
    menu.add_separator()
    menu.add_command(label='Follow Intensity Changes', shortcut='I', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.calcRates,
                     tipText='Fit peak intensities to graphs to estimate relaxation rates etc.')
    menu.add_command(label='Follow Shift Changes',shortcut='S', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.followShiftChanges,
                     tipText='Fit peak position changes to graphs to estimate binding constants etc.')
    menu.add_command(label='PALES: Alignment and RDCs', shortcut='P', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.pales,
                     tipText='Use PALES to determine alignment and analyse RDCs')
    menu.add_command(label='MODULE: Alignment and RDCs', shortcut='M', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.blackledge_module,
                     tipText='Use MODULE to determine alignment and analyse RDCs')
  
    self.menubar.add_cascade(label=DataMenu, shortcut='D', menu=menu)

    self.menus[DataMenu] = menu
    self.menu_items[DataMenu] = ['Measurement Lists',
                                 'NMR Series',
                                 'Shift Differences',
                                 'Heteronuclear NOE',
                                 u'3J H-H\u03B1 Coupling',
                                 'Follow Intensity Changes',
                                 'Follow Shift Changes',
                                 'PALES: Alignment and RDCs',
                                 'MODULE: Alignment and RDCs',]
 

  def setStructureMenu(self):

    menu = Menu(self.menubar, tearoff=0)
    cyanaSubmenu = Menu(self.menubar, tearoff=0)
    cyanaSubmenu.add_command(label='Setup CYANA calculation',
                     image=self.iconSpecialTool, compound='left',
                     command=self.setupCyanaCalculation,
                     tipText='Setup CYANA calculation to iteratively assign NOEs and calculate structures')
    cyanaSubmenu.add_command(label='Import CYANA calculation results',
                     image=self.iconSpecialTool, compound='left',
                     command=self.importCyanaData,
                     tipText='Import output data from a CYANA calculation into project')
    cyanaSubmenu.add_command(label='Run CYANA calculation',
                     image=self.iconSpecialTool, compound='left',
                     command=self.runCyana2Ccpn,
                     tipText='Setup and Run CYANA calculation and import calculation results')

    menu.add_command(label='Restraints and Violations', shortcut='R',
                     image=self.iconTable, compound='left',
                     command=self.browseConstraints,
                     tipText='Tables listing restraints, restraint lists, restraint violations etc.')
    menu.add_command(label='Structures', shortcut='S',   
                     image=self.iconTable, compound='left',
                     command=self.editStructures,
                     tipText='Tables of coordinate structures, RMSDs. Options to load & save PDB info.')
    menu.add_separator()
    menu.add_command(label='Structure Viewer',shortcut='V',   
                     image=self.iconTool, compound='left',
                     command=self.viewStructure,
                     tipText='A simple graphical view of structures to display NMR restraints etc.')
    menu.add_command(label='Make Distance Restraints', shortcut='M', 
                     image=self.iconTool, compound='left',
                     command=self.calDistConstraints,
                     tipText='Tools to make distance restraints from NOEs etc. Includes peak-shift matching.')
    menu.add_command(label='Make H Bond Restraints', shortcut='H', 
                     image=self.iconTool, compound='left',
                     command=self.makeHbonds,
                     tipText='Tools to make H-bond distance restraints')
    menu.add_separator()
    menu.add_command(label='DANGLE: Predict Dihedrals', shortcut='D', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.startDangle,
                     tipText='Predict protein backbone dihedral angles using chemical shifts')
    #menu.add_command(label='D2D: Predict Secondary Structure', shortcut='2', 
    #                 image=self.iconSpecialTool, compound='left',
    #                 command=self.startD2D,
    #                 tipText='Predict protein secondary structure using chemical shifts')
    menu.add_command(label='ARIA: Structure calculation', shortcut='A', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.startAria,
                     tipText='Use ARIA to iteratively assign NOEs and calculate structures')
    menu.add_cascade(label='Cyana', shortcut='y',
                     image=self.iconSpecialTool, compound='left',
                     menu=cyanaSubmenu)

    menu.add_command(label='HADDOCK: Structure Docking', shortcut='K', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.startHaddock,
                     tipText='Generate structures of complexes using high-ambiguity driven biomolecular docking')
    menu.add_command(label='MECCANO: Structures from RDCs', shortcut='M', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.meccano,
                     tipText='Use MECCANO to determine structures from backbone RDCs')
    menu.add_command(label='PyRPF: Validate Peaks vs Structure', shortcut='F', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.startPyRPF,
                     tipText='Compare through-space peak lists with structure distances to find missing and unexplained peak locations')
    menu.add_command(label='CING: Validate Structures', shortcut='C', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.submitCing,
                     tipText='Check structures and NMR data for quality and errors')
    menu.add_command(label='ECI: Database Deposition', shortcut='E', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.startECI,
                     tipText='Collate information form CCPN project for deposition to PDB & BMRB databases')
    menu.add_separator()
    menu.add_command(label='Secondary Structure Chart', shortcut='e',     
                     image=self.iconChart, compound='left',
                     command=self.secStructureGraph,
                     tipText='Tables and charts giving evidence for protein secondary structure')
    menu.add_command(label='Ramachandran Plot',     shortcut='P',     
                     image=self.iconChart, compound='left',
                     command=self.plotRamachandran,
                     tipText='A plot of phi & psi protein backbone angles in structures')

    self.menubar.add_cascade(label=StructureMenu, shortcut='S', menu=menu)

    self.menus[StructureMenu] = menu
    self.menu_items[StructureMenu] = ['Restraints and Violations',
                                      'Structures',
                                      'Structure Viewer',
                                      'Make Distance Restraints',
                                      'Make H Bond Restraints',
                                      'DANGLE: Predict Dihedrals',
                                      'HADDOCK: Structure Docking',
                                      'MECCANO: Structures from RDCs',
                                      'CING: Validate Structures',
                                      'ECI: Database Deposition',
                                      'Secondary Structure Chart',
                                      'Ramachandran Plot']

  def changedSpectrumWindowName(self, spectrumWindow):

    new_name = spectrumWindow.name
    if hasattr(spectrumWindow, 'old_name'):
      old_name = spectrumWindow.old_name
      old_popup_name = self.getWindowPopupName(old_name)
      popup = self.popups.get(old_popup_name)
      if popup:
        del self.popups[old_popup_name]
      else:
        print 'Warning: popup %s was expected to exist but did not, in popups list' % old_popup_name
      new_popup_name = self.getWindowPopupName(new_name)
      self.popups[new_popup_name] = popup
      
    spectrumWindow.old_name = new_name
    self.setWindowMenu()

  def initSpectrumWindow(self, spectrumWindow):

    Analysis.initSpectrumWindow(self, spectrumWindow)

    #self.after_idle(lambda: self.editWindow(spectrumWindow))
    # after_idle below does not seem to work (for some unknown reason):
    # when the EditWindowPopup is open then openWindow executes immediately
    #self.after_idle(lambda: self.openWindow(spectrumWindow))
    self.after(after_delay, lambda window=spectrumWindow: self.openWindow(window))
    self.after(after_delay, self.setWindowMenu)


  def initSpectrumWindowPane(self, windowPane, windowFrame=None):

    # Overwrites below (blank for now)
    #Analysis.initSpectrumWindowPane(self, windowPane)
    
    def getWindowPaneFrame(windowPane, windowFrame=None):
 
      if hasattr(windowPane, 'windowFrame'):
        return windowPane.windowFrame
 
      elif windowFrame:
        windowPane.windowFrame = windowFrame
        return windowFrame
 
      else:
        window = windowPane.spectrumWindow
 
        popupName = self.getWindowPopupName(window.name)
        windowPopup = self.popups.get(popupName)
 
        if not windowPopup:
          windowPopup = self.openWindow(window)
 
        windowPanes = window.sortedSpectrumWindowPanes()
        index = windowPanes.index(windowPane)
        windowFrame = windowPopup.windowFrames[index]
        windowPane.windowFrame = windowFrame
 
        return windowFrame

    windowPane.getWindowFrame = lambda wp=windowPane, wf=windowFrame: getWindowPaneFrame(wp, wf)


  def deleteSpectrumWindow(self, spectrumWindow):

    self.destroyWindow(spectrumWindow.name)

    self.setWindowMenu()

  def setWindowMenu(self):

    menu = self.menus.get(WindowMenu)
    if menu:
      menu.delete(0, Tkinter.END)
      
    else:
      menu = Menu(self.menubar, tearoff=0)
      self.menubar.add_cascade(label=WindowMenu, shortcut='W', menu=menu)
      self.menus[WindowMenu] = menu

    menu.add_command(label='Windows', shortcut='W',    
                     image=self.iconTable, compound='left',
                     command=self.editWindow,
                     tipText='Tables to setup spectrum windows, their display, axis mappings and groups')
    menu.add_command(label='Axes', shortcut='A',    
                     image=self.iconTable, compound='left',
                     command=self.editAxisPanel,
                     tipText='Edit and create new window axis')
    menu.add_separator()
    menu.add_command(label='New Window', shortcut='N',     
                     image=self.iconNewWindow, compound='left',
                     command=self.newWindow,
                     tipText='Tool to make a new spectrum window with required axes')
    menu.add_command(label='Print Window', shortcut='P',      
                     image=self.iconPrint, compound='left',
                     command=self.printWindow,
                     tipText='Render a spectrum window display as PS, EPS or PDF')
    menu.add_separator()
    menu.add_command(label='Marks and Rulers', shortcut='M',   
                     image=self.iconPrefs, compound='left',
                     command=self.editMarks,
                     tipText='Setup options for 1D ruler lines and multidimensional cross-marks')

    # for now do not put windows themselves in menu_items
    menu_items = self.menu_items[WindowMenu] = ['Windows',
                                                'Axes',
                                                'New Window',
                                                'Print Window',
                                                'Marks and Rulers']

    # TBD window shortcut
    windows = self.getActiveWindows()

    if windows:
      menu.add_separator()

      windowNames = []
      for window in windows:
        windowPanes = window.spectrumWindowPanes
        if len(windowPanes) == 1:
          windowPane = tuple(windowPanes)[0]
          axes = []
          axisPanels = windowPane.sortedAxisPanels()
 
          for axisPanel in axisPanels:
            axisType = axisPanel.axisType
 
            if axisType:
              aName = axisType.name
              if aName in axisType.isotopeCodes:
                while aName[0].upper() == aName[0].lower():
                  aName = aName[1:]
 
              else:
                aName = '(%s)' % aName
 
              axes.append(aName)
            else:
              axes.append('?')
          data = (len(axisPanels), ''.join(axes), window.name)
          key = '%dD %s %s' % data
        else:
          axes = ['*']
          data = (None, ''.join(axes), window.name)
          key = '%s %s %s' % data
 
        name = '%s: %s' % data[1:]
        windowNames.append((key, name, window))
 
      windowNames.sort()
      for key, name, window in windowNames:
        command = lambda selected=window: self.openWindow(selected)
        menu.add_command(label=name, command=command)

  def setMacroMenu(self, *opt):

    menu = self.menus.get(MacroMenu)
    if menu:
      menu.delete(0, Tkinter.END)
      
    else:
      menu = Menu(self.menubar, tearoff=1)
      self.menubar.add_cascade(label=MacroMenu, shortcut='M', menu=menu)
      self.menus[MacroMenu] = menu

    menu.add_command(label='Organise Macros',  shortcut='O',     
                     image=self.iconTable, compound='left',
                     command=self.editMacros,
                     tipText='Curate and manage Python macro scripts')
    menu.add_command(label='Reload Menu Macros', shortcut='R',     
                     image=self.iconRefresh, compound='left',
                     command=self.reloadMenuMacros,
                     tipText='Refresh the macro scripts listed below from disk versions')
    self.menu_items[MacroMenu] = ['Organise  Macros','Reload Menu Macros']

    macros = []
    if self.project:
      for macro in self.analysisProfile.macros:
        if macro.isInMenu:
          macros.append( (macro.name, macro) )
    
      if len(macros) > 0:
        menu.add_separator()
 
      macros.sort()
      for name, macro in macros:
        menu.add_command(label=name, command=lambda m=macro: self.runMacro(m))
        self.menu_items[MacroMenu].append(name)

  def reloadMenuMacros(self):
  
    for macro in self.project.currentAnalysisProfile.macros:
      if macro.isInMenu:
        Util.reloadMacro(macro,self.argServer)

  def runMacro(self,macro):

    sys.path.append(macro.path)
    try:
      command = Command(self.argServer,macro.name,macro.module,macro.function)
    finally:
      del sys.path[-1]
    command.run()

  def setOtherMenu(self):

    cloudsMenu = Menu(self.menubar, tearoff=1)
    cloudsMenu.add_command(label='Resonance Disambiguation',      command=self.setupBacus,
                           tipText='CCPN implementation of the program BACUS')
    cloudsMenu.add_command(label='Relaxation Matrix Optimisation',command=self.setupMidge,
                           tipText='CCPN implementation of the program MIDGE')
    cloudsMenu.add_command(label='Proton Molecular Dynamics',     command=self.setupHcloudsMd,
                           tipText='Create anonymous proton cloud structures using distance restraints')
    cloudsMenu.add_command(label='Filter Clouds',                 command=self.setupFilterClouds,
                           tipText='Quality control for proton cloud structures')
    cloudsMenu.add_command(label='Cloud Threading',               command=self.setupCloudThreading,
                           tipText='Thread a sequence through a proton cloud to assign')
    cloudsMenu.add_command(label='Cloud Homologue Assign',        command=self.setupCloudHomologue,
                           tipText='Assign a sequence by aligning a proton cloud with an homologous structure ')

    #for i in (1,):
    #  cloudsMenu.entryconfig(i, state=Tkinter.DISABLED)

    menu = Menu(self.menubar, tearoff=0)
    menu.add_command(label='NMR Calculations', shortcut='N',      
                     image=self.iconTable, compound='left',
                     command=self.editCalculation,
                     tipText='Curate and manage calculation jobs sent o external programs like, CING or ARIA')
    menu.add_command(label='Widget Counter',  shortcut='W',     
                     image=self.iconTool, compound='left',
                     command=self.viewWidgetCount,
                     tipText='A temporary developer module to keep track of graphical object creation/destruction')
    menu.add_command(label='Format Converter', shortcut='F',      
                     image=self.iconSpecialTool, compound='left',
                     command=self.runFormatConverter,
                     tipText='Export and import CCPN project data to and from a multitude of textual NMR formats')

    menu.add_command(label='Prodecomp', shortcut='P', 
                     image=self.iconSpecialTool, compound='left',
                     command=self.startProdecomp,
                     tipText='PRODECOMP: Process sets of projection spectra by decomposition')



    menu.add_cascade(label='CLOUDS',           shortcut='C',      
                     image=self.iconClouds, compound='left',
                     menu=cloudsMenu)

    self.menubar.add_cascade(label=OtherMenu, shortcut='O', menu=menu)

    self.menus[OtherMenu] = menu
    self.menu_items[OtherMenu] = ['NMR Calculations',
                                  'Widget Counter',
                                  'Format Converter',
                                  'Prodecomp',
                                  'CLOUDS',
                                  ]
    
    # FormatConverter always active
    self.fixedActiveMenus[(OtherMenu, 2)] = True
    
  def initProject(self, project=None, resizeTop=True):

    self.project = project
    self.projectName = None
    self.setTitle(self.program_name)
    self.argumentServer = self.argServer = None

    if project:

      self.projectName = project.name
      self.curatePopupNotifiers(self.registerNotify)
      self.argumentServer = ArgumentServer(self, inGui=True)
      self.argServer = self.argumentServer
      errText = 'Project invalid, please quit, fix and re-start: '
      
      try:
        Analysis.initProject(self, project)
      except Implementation.ApiError, e:
        showError('Project invalid', errText + e.error_msg, parent=self)
        self.curatePopupNotifiers(self.unregisterNotify)
        self.project = None
        #return
        raise
      
      except:
        exc = sys.exc_info()
        ee = exc[0]
        if type(ee) == type(type):
          ee = ee.__name__
        error_msg = '%s: %s' % (ee, exc[1])
        showError('Project invalid', errText + error_msg, parent=self)
        self.curatePopupNotifiers(self.unregisterNotify)
        self.project = None
        #return
        raise

      geom = self.getPopupGeometry('top')
      if geom and resizeTop:
        try:
          self.geometry(geom)
        except:
          # 27 Aug 2010: had case when geom was saved as "683x-1+236+0"
          pass

      self.setFont(self.analysisProfile.font)

      self.initPopupLocation('top')
      self.openActiveWindows()
      self.openActivePopups()
      #self.after_idle(self.drawWindowPopups)
      
      self.setMacroMenu()
      
      # First time reference experiments
      experiments = []
      found = 0
      for experiment in self.project.currentNmrProject.sortedExperiments():
        if not experiment.refExperiment:
          simulated = 0
          for spectrum in experiment.dataSources:
            if spectrum.isSimulated:
              simulated +=1
          
          if not simulated or (len(experiment.dataSources) != simulated):
            experiments.append(experiment)

        else:
          found = 1

      if experiments and not found:
        self.initRefExperiments(experiments)

      self.webBrowser = ProjectWebBrowser(self.top, popup=self, project=project)

      self.checkRegistration()

    self.setExperimentMenu()
    self.setWindowMenu()

    try:
      self.setMenuState()
    except: # fails first time around because OtherMenu not set up yet
      pass

    if project:
      callbacks = self.callbacksDict.get('init', ())
      for callback in callbacks:
        callback(project)

  def checkRegistration(self):

    analysisProfile = self.analysisProfile
    # TBD: remove check when attributes in released API

    # # test to unregister - need to load a project though
    # analysisProfile.userName = None
    # analysisProfile.userOrganisation = None
    # analysisProfile.userEmail= None

    if hasattr(analysisProfile, 'userName'):
      if isRegistered(analysisProfile):
        if not self.printedRegistration:
          print 'Registered (%s, %s, %s)' % (analysisProfile.userName, analysisProfile.userOrganisation, analysisProfile.userEmail)
          self.printedRegistration = True
        try:
          updateRegister(analysisProfile)
        except:
          pass  # not that important if above failed
      else:
        self.registerAnalysis(isModal=True)

  def registerCallback(self, funcname, callback):

    if funcname not in self.allowedCallbackFuncnames:
      raise Exception('funcname = %s, must be in %s' % (funcname, self.allowedCallbackFuncnames))

    if not self.callbacksDict.has_key(funcname):
      self.callbacksDict[funcname] = []
    self.callbacksDict[funcname].append(callback)

  def unregisterCallback(self, funcname, callback):

    self.callbacksDict[funcname].remove(callback)

  def initRefExperiments(self, experiments=None):
  
    popup = EditExperimentPopup(self, isModal=True)
  
    if experiments:
      popup.updateExpTypes(experiments)

  def setMenuState(self):
    
    # set default state, according to existence of project
    if self.project:
      state = Tkinter.NORMAL
    else:
      state = Tkinter.DISABLED
    
    for name,menu in self.menus.items():
      menu_items = self.menu_items[name]
      for n in range(len(menu_items)+4): # +4 in case of separators
        try:
          menu.entryconfig(n, state=state)
        except:
          pass
    
    # reset state for special cases
    for key, val in self.fixedActiveMenus.items():
      if val:
        state = Tkinter.NORMAL
      else:
        state = Tkinter.DISABLED
      name, index = key
      self.menus[name].entryconfig(index, state=state)
      
      
  def newProject(self, name=''):

    def validName(name):
      if len(name) > 32:
        return False
      if len(name) < 1:
        return False
      for cc in name:
        if cc != '_' and not cc.isalnum():
          return False
      return True

    if self.project:
      if not self.closeProject():
        return

    prompt = 'Enter project name:'
    while not name:
      name = askString(title='Project : New', prompt=prompt, parent=self)

      if name is None:
        break

      elif not validName(name):
        name = ''
        prompt = 'Name invalid.\nEnter project name\n(between 1 and 32 chars; alphanumeric and "_" only):'

    if name:
      project = Impl.MemopsRoot(name=name)
      self.initProject(project)

  def openPopup(self, popup_name, clazz, oldStyle=False, *args, **kw):

    popup = self.popups.get(popup_name)
  
    if popup:
      popup.open()
      
    else:
      if self.project:
        analysisProfile = self.analysisProfile
        geom = self.getPopupGeometry(popup_name)
      else:
        analysisProfile = None
        geom = None
        
      if popup_name.startswith(window_popup_prefix): # bit of a hack
        transient = analysisProfile.transientWindows
        name = None # location currently handled directly by data model
      else:
        if analysisProfile:
          transient = analysisProfile.transientDialogs
        else:
          transient = True
        name = popup_name
      
      if oldStyle:
        popup = self.popups[popup_name] = clazz(self, transient=transient, *args, **kw)
      else:
        popup = self.popups[popup_name] = clazz(self, project=self.project, popup_name=name,
                                                transient=transient, *args, **kw)
      # above automatically opens popup

      if geom:
        popup.geometry(geom)

    return popup

  def openProject(self):

    if self.project:
      if not self.closeProject():
        return

    load_project = lambda path: loadProject(self, path)
    self.openPopup('open_project', OpenProjectPopup,
                   callback=self.initProject,
                   load_project=load_project,
                   help_url=determineHelpUrl(OpenProjectPopup))

  def openSpectrum(self):

    self.openPopup('open_spectrum', OpenSpectrumPopup)

  def loadNefFile(self):

    from ccpnmr.v2io.NefIo import loadNefFile as _loadNefFile

    if self.project:
      if not self.closeProject():
        return

    fileTypes = [FileType('STAR', ['*.nef']),
                 FileType('All', ['*'])]

    fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                                      title='Load Nef file', dismiss_text='Cancel',
                                      selected_file_must_exist=True, multiSelect=False, )

    fileName = fileSelectPopup.getFile()

    if fileName:

      try:
        project = _loadNefFile(fileName, overwriteExisting=True)
        self.initProject(project)

      except Exception as es:
        print str(es)

  # # NOTE:ED not needed
  # def importNefFile(self):
  #   print '>>>Import Nef file'
  #
  # def exportNefFile(self):
  #   print '>>>Export Nef file'

  def editSpectrum(self, spectrum=None):

    popup = self.openPopup('edit_spectrum', EditSpectrumPopup)

    if spectrum:
      popup.spectrum = spectrum
      popup.updateAfter()
      popup.tabbedFrame.select(0)

    return popup

  def addSpectrumContourFile(self, spectrum = None):

    popup = self.openPopup('add_contour_file', AddContourFilePopup)

    if (spectrum):
      popup.setSpectrum(spectrum)

  def createSpectrumContourFile(self, spectrum = None):

    popup = self.openPopup('create_contour_file', CreateContourFilePopup)

    if (spectrum):
      popup.setSpectrum(spectrum)

  def editExperiment(self):

    popup = self.openPopup('edit_experiment', EditExperimentPopup)
    return popup

  def editExpType(self):

    popup = self.editExperiment()
    popup.tabbedFrame.select(1)
    return popup

  def editExpSeries(self, experiment=None):

    popup = self.openPopup('edit_experiment_series', EditExperimentSeriesPopup)
    
    if experiment:
      popup.update(experiment=experiment)

  def editExpPrototype(self):

    self.openPopup('edit_experiment_prototypes', NmrExpPrototypePopup,
                   title='Experiment: Experiment Prototypes')

  def editSampleConditionSets(self, sampleConditionSet=None):

    showWarning('Warning','Not implemented', parent=self)
    #self.openPopup('shift_lists', Popup)

  def editMeasurementLists(self, experiment=None):

    popup = self.openPopup('edit_measurement_lists', EditMeasurementListsPopup)
    
    if experiment:
      popup.setExperiment(experiment=experiment)

    return popup

  def editMeasurements(self, measurementList=None):

    popup = self.editMeasurementLists()
    popup.tabbedFrame.select(1)
    
    if measurementList:
      popup.updateMeasurements(measurementList=measurementList)


  def editPeak(self, peak = None, peakList = None):

    popup = self.openPopup('edit_peak', EditPeakPopup)

    if (peak or peakList):
      popup.update(peak=peak, peakList=peakList)

  def editPeakAliasing(self, peak=None):
  
    popup = self.openPopup('edit_peak_aliasing', EditPeakAliasingPopup)
   
    if peak:
      popup.updateAfter(object=peak)

  def viewSelectedPeaks(self):
  
    peaks = self.currentPeaks or []
    popup = self.openPopup('selected_peaks', SelectedPeaksPopup)
    popup.update(peaks)
   
  def viewPeaks(self, peaks=None):

    popup = self.openPopup('view_peaks', SelectedPeaksPopup)
    popup.update(peaks)

  def viewPeakGroups(self, groups=None, variableDim=None, peakLists=None,searchPeaks=None):
  
    popup = self.openPopup('view_peak_groups', BrowsePeakGroupsPopup)

    if groups or (variableDim is not None):
      popup.update(groups,variableDim,peakLists,searchPeaks)
 
  def editPeakList(self, peakList=None):

    popup = self.editPeakLists()

    if peakList:
       popup.peakTableFrame.update(peakList)
       
    popup.tabbedFrame.select(1)

  def editPeakLists(self):

    popup = self.openPopup('edit_peak_lists', EditPeakListsPopup)
    return popup

  def copyAssignments(self):
  
    self.openPopup('copy_assignments', CopyAssignmentsPopup)


  def browseAtoms(self, chain=None, requestor=None):

    popup = self.openPopup('browse_atoms', BrowseAtomsPopup)
    if chain:
      popup.setChain(chain)

    if requestor:
      popup.setRequestor(requestor)
      
    return popup  
  
  def browseReferenceShifts(self):

    self.openPopup('browse_reference_shifts', BrowseReferenceShiftsPopup)
  
  def isotopomerEditor(self):
  
    self.openPopup('isotope_scheme_editor', IsotopeSchemeEditor)

  def editIsotopeLabelling(self):
  
    self.openPopup('isotope_labelling', EditMolLabellingPopup)
    
  def browseResidueInfo(self):

    self.openPopup('browse_residue_info', ResidueInfoPopup)

  def browseResonanceInfo(self, resonance=None):

    popup = self.openPopup('browse_resonance_info', ResonanceInfoPopup,
                           resonance=resonance)
    popup.update(resonance)
    

  def browseResonances(self):

    popup = self.openPopup('browse_resonances', BrowseResonancesPopup)
    popup.update() # TBD: not sure if this is really needed
    
  
  def viewSelectedResonances(self, resonances=None, shiftList=None):
  
    # Could be none if restore prev windows
    # not storing prev selected resonances though
  
    if resonances is not None:
      def SelectedResonances(parent, *args, **kw):
        return BrowseResonancesPopup(parent, resonances, shiftList, *args, **kw)
      
      popup = self.openPopup('selected_resonances', SelectedResonances)
      popup.update(resonances, shiftList)
     
  def editSpinSystems(self, doClear=False):

    popup = self.openPopup('edit_spin_system', EditSpinSystemPopup)
    if doClear:
      popup.update(doClear)

  def initialiseRootSpectra(self):
  
    popup = self.openPopup('initialise_root_spectra', InitRootAssignmentsPopup)      

  def linkPeakLists(self, doClear=False):

    popup = self.openPopup('link_peaklists', LinkPeakListsPopup)  

  def linkSeqSpinSystems(self, doClear=False):

    popup = self.openPopup('link_seq_spin_systems', LinkSeqSpinSystemsPopup)  

  def linkNoeResonances(self, doClear=False):

    popup = self.openPopup('link_noe_resonances', LinkNoeResonancesPopup)  

  def confirmSeqSpinSystems(self):

    popup = self.openPopup('confirm_seq_spin_systems',ConfirmSeqSpinSystemsPopup)

  def typeSpinSystems(self):

    popup = self.openPopup('type_spin_systems', SpinSystemTypingPopup)

  def typeSpinSystem(self, spinSystem=None, shiftList=None):

    popup = self.openPopup('type_spin_system', SpinSystemTypeScoresPopup)

    if spinSystem:
      popup.update(spinSystem, shiftList=shiftList)

  def assignmentPanel(self):

    popup = self.openPopup('edit_assignment', EditAssignmentPopup)
    #popup.update() # TBD: not sure if this is really needed

  def qualityReports(self):

    self.openPopup('quality_reports', ViewQualityReportsPopup)

  def autoBackboneAssign(self):
  
    if not HAVE_NUMPY:
      msg  = 'You must have the NumPy Python module installed to run'
      msg += ' the CcpNmr automatic assignment routines.' 
      showWarning('Cannot launch', msg, parent=self)
    
    self.openPopup('auto_backbone_assign', AutoBackbonePopup,
                   title='Assignment : Automated Seq. Assignment')
    
  
  #def activateMars(self):
  
  #  self.openPopup('auto_backbone_assign', AutoBackbonePopup,
  #                 title='Assignment : Automated Seq. Assignment')
                   
  #  popup = self.popups['auto_backbone_assign']
  #  popup.activateMars()
  
  #def activatePales(self):
  #  #if self.project:
  #  menu = self.menus[DataMenu]
  #  menu.add_command(label='PALES: Reduced Dipolar Couplings', shortcut='P', 
  #                   image=self.iconSpecialTool, compound='left',
  #                   command=self.pales)
  #  #self.pales()
    
  #  #else:
  #  #  print 'PALES cannot be started without an open project'

  #def activateModule(self):

  #  #if self.project:
  #  menu = self.menus[DataMenu]
  #  menu.add_command(label='MODULE: Rigid Body Modelling using RDCs', shortcut='M', 
  #                   image=self.iconSpecialTool, compound='left',
  #                   command=self.blackledge_module)
  #  #self.blackledge_module()

  #  #else:
  #  #  print 'MODULE cannot be started without an open project'

  def chemShiftsTable(self):

    self.openPopup('view_chem_shifts', ViewChemicalShiftsPopup)

  #def sequenceShiftPredict(self):

  #  self.openPopup('sequence_shift_predict', SequenceShiftPredictPopup)

  def assignmentGraph(self):

    self.openPopup('view_assignment', ViewAssignmentPopup)
    
  def secStructureGraph(self):
  
    self.openPopup('sec_struc_graph', SecStructureGraphPopup)
          
    
  def viewNoeMatrix(self):
  
    self.openPopup('view_noe_matrix', ViewNoeMatrix)
        
  def editStructures(self):

    popup = self.openPopup('edit_structures', EditStructuresPopup)
    return popup

  def plotRamachandran(self):

    self.openPopup('plot_ramachandran', ViewRamachandranPopup)
 
  def startDangle(self):
  
    self.popups['dangle'] = DangleGui(self, project=self.project)

  #def startD2D(self):
  
  #  self.popups['d2d'] = SecStructurePredictPopup(self, project=self.project)

  def startHaddock(self):
  
    self.popups['haddock'] = HaddockPopup(self, self.project)

  def startPyRPF(self):
  
    self.popups['rpf'] = PyRpfPopup(self, self.project)

  def startProdecomp(self):
  
    self.popups['prodecomp'] = ProdecompPopup(self, self.project)

  def submitCing(self):
  
    self.popups['cing'] = CingPopup(self)
 
  def activateAriaSetup(self):
    # Leave for a bit so old tutorial script works.
  
    self.openPopup('aria_setup', AriaPopup)
 
  def startAria(self):
  
    self.openPopup('aria_setup', AriaPopup)

  def setupCyanaCalculation(self):

    from ccpnmr.analysis.macros.MultiStructure  import setupCyanaCalculationDialogue
    setupCyanaCalculationDialogue(self.argumentServer)


  def importCyanaData(self, calculationData=None):

    from ccpnmr.analysis.macros.MultiStructure  import importDataFromCyana
    if calculationData == None:
      dataSources = importDataFromCyana(self.argumentServer)
    else:
      from cyana2ccpn.cyana2ccpn import importFromCyana
      dataSources = importFromCyana(calculationData[0],calculationData[1])
    for dataSource in dataSources:
      # self.finishInitSpectrum(dataSource)
      Analysis.finishInitSpectrum(self, dataSource)
    print "DONE"

  def runCyana2Ccpn(self):

    from ccpnmr.analysis.macros.MultiStructure  import runCyana2CcpnDialogue
    calculationData = runCyana2CcpnDialogue(self.argumentServer)
    yy = self.argumentServer.askYesNo("Import Calculation Results")
    print 'calcData',calculationData
    if yy:
      print calculationData
      self.importCyanaData(calculationData=calculationData)


  def viewStructure(self, structure=None):

    popup = self.openPopup('view_structure', ViewStructurePopup)
    if structure:
      popup.update(structure)                   

  def calDistConstraints(self):
  
    self.openPopup('calc_dist_constraints',CalcDistConstraintsPopup)
            
  def browseConstraints(self, constraintList=None):

    popup = self.openPopup('browse_restraints', BrowseConstraintsPopup)
    
    if constraintList:
      popup.update(constraintList)
  
  def makeHbonds(self):
    
    popup = self.openPopup('make_hydrogen_bonds', MakeHbondRestraintsPopup)
     
  def peakDrawParams(self):

    self.openPopup('edit_peak_draw_params', EditPeakDrawParamsPopup)

  def peakFindParams(self):

    self.openPopup('edit_peak_find_params', EditPeakFindParamsPopup)

  def editMolecules(self, molecule=None):

    popup = self.openPopup('edit_molecules', EditMoleculesPopup)
    return popup
  
  def editMolSystems(self):

    popup = self.editMolecules()
    popup.tabbedFrame.select(0)

  
  def addSequence(self):
  
    popup = self.editMolecules()
    popup.tabbedFrame.select(5)
    popup.addSequence()
  
  def calcRates(self):
  
    self.openPopup('calc_rates', CalcRatesPopup)

  def calcShiftDifferences(self):
  
    self.openPopup('calc_shift_difference', CalcShiftDifferencePopup)

  def followShiftChanges(self):
  
    self.openPopup('follow_shift_changes', FollowShiftChangesPopup)

  def calcHnHaCoupling(self):
  
    self.openPopup('calc_hnha_coupling', CalcHnHaCouplingPopup)

  def calcHeteroNoe(self):
  
    self.openPopup('calc_hetero_noe', CalcHeteroNoePopup)

  def editCalculation(self):
  
    self.openPopup('edit_calculation', EditCalculationPopup)

  def pales(self):
    from gottingen.PalesPopup                             import PalesPopup
    self.openPopup('pales', PalesPopup)

  def blackledge_module(self):
    from grenoble.BlackledgeModule.BlackledgeModulePopup import BlackledgeModulePopup
    self.openPopup('blackledge_module', BlackledgeModulePopup)

  def meccano(self):
    try:
      from grenoble.meccano.MeccanoPopup import MeccanoPopup
    except Exception, e:
      showWarning('Meccano exception', str(e), parent=self)
      print e
      return
      
    self.openPopup('meccano', MeccanoPopup)

  def viewWidgetCount(self):
  
    self.openPopup('view_widget_count', WidgetCountPopup)

  def backupProject(self):

    self.openPopup('backup_project', BackupProjectPopup,
                   help_url=determineHelpUrl(name='BackupProjectPopup'))

  def archiveProject(self):

    self.openPopup('archive_project', ArchiveProjectPopup,
                   help_url=determineHelpUrl(name='ArchiveProjectPopup'))

  def saveCallback(self, *extra):

    if self.project:
      if self.project.name != self.projectName: # SaveAs can change project name
        self.setTitle(self.program_name)
        for key in self.popups.keys():
          popup = self.popups[key]
          try:
            title = popup.getTitle()
            popup.setTitle(title)
          except Exception, e:
            print 'saveCallback exception', str(e)
            continue
          except:
            continue
        self.projectName = self.project.name

    callbacks = self.callbacksDict.get('save', ())
    for callback in callbacks:
      callback(self.project)

  def setSaveState(self):

    self.timeStampTopObjects()
    self.setPopupsOpen()
    self.setPopupGeometries()

  def askSaveFile(self):

    self.setSaveState()
    self.openPopup('save_project', SaveProjectPopup, callback=self.saveCallback,
                   help_url=determineHelpUrl(name='SaveProjectPopup'))

  def modalAskSaveFile(self):

    self.setSaveState()
    popup = SaveProjectPopup(self, project=self.project, dismiss_text='Cancel Quit',
                             help_url=determineHelpUrl(name='SaveProjectPopup'), modal=True)
    did_save = popup.did_save
    popup.destroy()

    return did_save

  def saveProject(self):

    self.checkRegistration()

    if self.project.activeRepositories:
      self.saveFile()
    else:
      self.askSaveFile()

  def quitSaveProject(self):

    if self.project.activeRepositories:
      return self.saveFile()
    else:
      return self.modalAskSaveFile()

  def saveAsProject(self):

    self.askSaveFile()

  def importNmrStar211(self):
  
    if self.project:
    
      fileTypes = [  FileType('STAR', ['*.str']),
                     FileType('All', ['*'])]
      
      fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                          title='Import NMR-STAR 2.1.1 file', dismiss_text='Cancel',
                          selected_file_must_exist=True, multiSelect=False,)

      fileName = fileSelectPopup.getFile()

      if fileName:
        nmrStarObj = NmrStarFormat(self.project, self, verbose=True)
        nmrStarObj.readProject(fileName, minimalPrompts=True, version='2.1.1')
        
        msg = 'NMR-STAR 2.1 data loaded. If you imported resonance assignments '
        msg += 'consider running Other::FormatConverter::Process::Run linkResonances '
        msg += 'to see the original assignments in CCPN.'

        showInfo(self, msg, parent=self)
        
        return nmrStarObj

  def importNmrStar31(self):
  
    if self.project:
    
      fileTypes = [  FileType('STAR', ['*.str']),
                     FileType('All', ['*'])]
      
      fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                          title='Import NMR-STAR 3.1 file', dismiss_text='Cancel',
                          selected_file_must_exist=True, multiSelect=False,)

      fileName = fileSelectPopup.getFile()

      if fileName:
        nmrStarObj = NmrStarFormat(self.project, self, verbose=True)
        nmrStarObj.readProject(fileName, minimalPrompts=True, version='3.1')
        
        msg = 'NMR-STAR 3.1 data loaded. If you imported resonance assignments '
        msg += 'consider running Other::FormatConverter::Process::Run linkResonances '
        msg += 'to see the original assignments in CCPN.'

        showInfo(self, msg, parent=self)

        return nmrStarObj
 
  def importPdb(self):

    if self.project:
    
      fileTypes = [  FileType('PDB', ['*.pdb']),
                     FileType('PDB Entry', ['*.ent']),
                     FileType('All', ['*']) ]
                     
      fileSelectPopup = FileSelectPopup(self, file_types=fileTypes,
                          title='Import PDB 3.20 file', dismiss_text='Cancel',
                          selected_file_must_exist=True, multiSelect=False,)

      fileName = fileSelectPopup.getFile()

      if fileName:
        name = os.path.split(fileName)[1]
        if '.' in name:
          name = name.split('.')[0]
        
        return ReadPdb(fileName, self.project, 'PDB:'+ name)

  def importCoordinates(self):

    if self.project:
      from os.path import dirname
 
      dataRepository = self.project.findFirstRepository(name='userData')
      projPath = dirname(dataRepository.url.path)
 
      fileTypes = [  FileType('PDB', ['*.pdb']), FileType('All', ['*'])]
      fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
                 title = 'Import PDB-style Coordinates', dismiss_text = 'Cancel',
                 selected_file_must_exist = True, multiSelect=True,
                 directory=projPath)

      fileNames = fileSelectPopup.file_select.getFiles()

      if fileNames:
        popup = self.editStructures()
        popup.tabbedFrame.select(0)
        popup.importStructure(fileNames)
        popup.open()

  def validateProject(self):

    if not showYesNo('Project : Validate', 'Do a project validation?', parent=self):
      return

    if showYesNo('Project : Validate',
                    'Do a complete validation (can take some time)?', parent=self):
      complete = True
    else:
      complete = False

    try:
      self.project.checkAllValid(complete=complete)
      showInfo('Project Valid', 'Project is valid (with complete=%s)' % complete, parent=self)
    except:
      showWarning('Project Not Valid', 'Project is not valid, see console for reason', parent=self)
      raise

  def copyStorage(self, storage):

    path = joinPath(storage.url.path, storage.path)
    newPath = path + '.bak'

    try:
      shutil.copyfile(path, newPath)
    except:
      pass # not much else can do

  def copyModifiedStorages(self):

    modifiedStorages = [storage for storage in self.project.storages if storage.isModified]
    if (self.project.isModified):
      modifiedStorages.append(self.project)

    for storage in modifiedStorages:
      self.copyStorage(storage)

  def timeStampTopObjects(self):
  
    proj = self.project
    topObjects = list(proj.molSystems) + list(proj.nmrProjects) + list(proj.analysisProjects)
    
    for topObject in topObjects:
      if topObject.isLoaded and topObject.isModified:
        Util.setTopObjectAnalysisSaveTime(topObject)

  def getPopupGeometry(self, popup_name):

    key = 'popup_geometry:%s' % popup_name

    return self.application.getValue(self.analysisProject, key)

  def setPopupGeometry(self, popup, popup_name):

    key = 'popup_geometry:%s' % popup_name

    try:
      geometry = popup.geometry()
      self.application.setValue(self.analysisProject, key, str(geometry))
    except Exception, e:
      print 'setPopupGeometry exception', str(e)
    except:
      pass

  def setPopupGeometries(self):

    # TBD: what about other popups (i.e. created elsewhere)

    self.setPopupGeometry(self, 'top')

    for key in self.popups.keys():
      popup = self.popups[key]
      if not key.startswith(window_popup_prefix): # bit of a hack
        self.setPopupGeometry(popup, key)

  def getPopupsOpen(self):

    key = 'popups_open'
    value = self.application.getValue(self.analysisProject, key)
    if value:
      popup_names = value.split(':')
    else:
      popup_names = []

    return popup_names

  def setPopupsOpen(self):

    # TBD: what about other popups (i.e. created elsewhere)
    popup_names = []
    for key in self.popups.keys():
      popup = self.popups[key]
      try:
        state = popup.state()
      except Exception, e:
        print 'setPopupsOpen exception', str(e)
        continue
      except:
        continue
      if state == 'normal' and not key.startswith(window_popup_prefix): # bit of a hack
        popup_names.append(key)

    key = 'popups_open'
    if popup_names:
      value = ':'.join(popup_names)
    else:
      value = None

    try:
      self.application.setValue(self.analysisProject, key, value)
    except Exception, e:
      print 'setPopupsOpen setValue exception', str(e)

  def saveFile(self):

    if not self.project:
      return False

    try:
      self.setSaveState()
      if saveProject(self.project, createFallback=True, showWarning=showWarning):
        print 'successfully saved project'
        self.saveCallback()
        return True
      else:
        return False
    except IOError, e:
      showError('Saving file', str(e), parent=self)
      return False

  def checkSaving(self):

    if not self.project:
      return True

    if isProjectModified(self.project):
      if showYesNo('Save project', 'Save changes to project?', parent=self):
        if not self.quitSaveProject():
          return False

    return True

  def closeProject(self, queryClose = True, querySave = True):

    if queryClose:
      if (not showYesNo('Project : Close', 'Close current project?', parent=self)):
        return False

    if querySave:
      if (not self.checkSaving()):
        return False

    callbacks = self.callbacksDict.get('close', ())
    for callback in callbacks:
      callback(self.project)

    self.webBrowser = self.defaultWebBrowser
    self.curatePopupNotifiers(self.unregisterNotify)
    self.destroyPopups()
    # self.argServer cleaned up in below
    Analysis.closeProject(self)
    # self.closeProject()

    return True

  def destroyPopups(self):

    # TBD: what about other popups (i.e. created elsewhere)
    for key in self.popups.keys():
      popup = self.popups[key]
      popup.destroy()

    self.popups = {}

  def newWindow(self):

    self.openPopup('new_window', NewWindowPopup)

  def editWindow(self, window=None):

    popup = self.openPopup('edit_window', EditWindowPopup)

    if (window):
      popup.setWindow(window)
      

  def printWindow(self, window=None):

    popup = self.openPopup('print_window', PrintWindowPopup)

    if window:
      # TBD: setWindow currently takes WindowPane
      popup.setWindow(window.findFirstSpectrumWindowPane())

  def editSpectrumFile(self, spectrum=None):

    popup = self.editSpectrum()

    if spectrum:
      popup.fileSpectrum = spectrum
      popup.updateAfter()
      
    popup.tabbedFrame.select(4)  
  
  def editSpectrumReferencing(self, spectrum=None):

    popup = self.editSpectrum()

    if spectrum:
      popup.refSpectrum = spectrum
      popup.updateAfter()
      
    popup.tabbedFrame.select(2)

  def editSpectrumTolerances(self, spectrum=None):

    popup = self.editSpectrum()

    if spectrum:
      popup.tolSpectrum = spectrum
      popup.updateAfter()
      
    popup.tabbedFrame.select(3)

  def popupEditProfiles(self, tab=0):

    # change tab first where possible: Avoids flicker
    popup = self.popups.get('edit_profiles')
    if popup:
      popup.tabbedFrame.select(tab)
      popup = self.openPopup('edit_profiles', EditProfilesPopup)
    
    else:
      popup = self.openPopup('edit_profiles', EditProfilesPopup)
      popup.tabbedFrame.select(tab)

    return popup

  def editProfiles(self):
  
    popup = self.popupEditProfiles(tab=0)

    return popup

  def editMacros(self):
  
    self.popupEditProfiles(tab=1)

  def editColorSchemes(self):

    self.popupEditProfiles(tab=3)

  def editResidueCodes(self):

    self.popupEditProfiles(tab=2)

  def editAxisPanel(self):

    self.openPopup('edit_axis_panel', EditAxisPanelPopup)

  def editMarks(self):

    self.openPopup('edit_marks', EditMarksPopup)

  def peakSeparatorParams(self):

    if not HAVE_NUMPY:
      msg  = 'You must have the NumPy Python module installed to run'
      msg += ' the peak separator routines.' 
      showWarning('Cannot launch', msg, parent=self)
      return

    from cambridge.bayes.PeakSeparatorGui import PeakSeparatorGui
    self.openPopup('edit_peak_separator_params', PeakSeparatorGui)

  def setupClouds(self):
 
    from ccpnmr.clouds.CloudsPopup import CloudsPopup
    self.openPopup('setup_clouds', CloudsPopup)

  def setupBacus(self):
 
    from ccpnmr.clouds.BacusPopup import BacusPopup
    self.openPopup('setup_bacus', BacusPopup)

  def setupMidge(self):

    from ccpnmr.clouds.MidgePopup import MidgePopup
    self.openPopup('setup_midge', MidgePopup)

  def setupHcloudsMd(self):
  
    from ccpnmr.clouds.HcloudsMdPopup import HcloudsMdPopup
    self.openPopup('setup_hcloudsmd', HcloudsMdPopup)

  def setupFilterClouds(self):

    from ccpnmr.clouds.FilterCloudsPopup import FilterCloudsPopup
    self.openPopup('setup_filter_clouds', FilterCloudsPopup)

  def setupCloudThreading(self):

    from ccpnmr.clouds.CloudThreaderPopup import CloudThreaderPopup
    self.openPopup('setup_cloud_threader', CloudThreaderPopup)

  def setupCloudHomologue(self):

    from ccpnmr.clouds.CloudHomologueAssignPopup import CloudHomologueAssignPopup
    self.openPopup('setup_cloud_homologue', CloudHomologueAssignPopup)

  def runFormatConverter(self):

    from ccpnmr.format.gui.FormatConverter import FormatConverter

    threading = Util.getFormatConverterThreading(self.analysisProject)
    popup = self.openPopup('run_format_converter', FormatConverter, project=self.project, oldStyle=True, threading=threading)
    popup.initProject(self.project)
    popup.protocol('WM_DELETE_WINDOW', popup.close)

  def startECI(self):

    from ccpnmr.eci.EntryCompletionPopup import EntryCompletionPopup
    
    popup = self.openPopup('entry_completion_interface', EntryCompletionPopup)

  def editContourLevels(self, spectrum=None):

    popup = self.openPopup('edit_contour_levels', EditContourLevelsPopup)

    if (spectrum):
      popup.setSpectrum(spectrum)

  def editContourFiles(self):

    self.openPopup('edit_contour_files', EditContourFilesPopup)

  def getWindowPopups(self):

    popups = []
    for key in self.popups.keys():
      if key.startswith(window_popup_prefix): # bit of a hack
        popups.append( self.popups[key] )

    return popups

  def getWindowPopupName(self, window_name):

    popup_name = window_popup_prefix + window_name

    return popup_name

  def getWindowPopup(self, window_name, doOpen=False):

    popup_name = self.getWindowPopupName(window_name)
    popup      = self.popups.get(popup_name)
    
    if (not popup) and doOpen:
      windows = self.getActiveWindows()
      for window in windows:
        if window.name == window_name:
          popup = self.openWindow(window)
          break
     
    return popup

  def openWindow(self, window):

    if not self.project:
      return None

    if not isActiveWindow(window):
      return None

    #print 'openWindow', window.name
    window.old_name = window.name
    popup_name = self.getWindowPopupName(window.name)
    location = '+%d+%d' % window.location
    #print 'openWindow', window.name, window.location
    popup = self.openPopup(popup_name, WindowPopup, window=window, location=location)
    #self.update_idletasks() # TBD: not sure if this is needed

    #if (window.isIconified):
    #  # if not after_idle location is messed up
    #  self.after_idle(lambda: popup.close())

    window.isIconified = False

    return popup

  def destroyWindow(self, window_name):

    if (not self.project):
      return

    popup_name = self.getWindowPopupName(window_name)
    popup = self.popups.get(popup_name)
    if (popup):
      del self.popups[popup_name]
      popup.destroy()

  def destroy(self):
  
    # Analysis.destroy(self)
    BasePopup.destroy(self)

  def quit(self):

    if not showYesNo('Quit ' + self.program_name,
                     'Quit ' + self.program_name + '?',
                     parent=self):
      return

    if self.project:
      if (not self.checkSaving()):
        return

    # NOTE:ED - add a flag to show that the program has already exited
    self._hasQuitNormally = True

    self.destroy()

    if not isinstance(self.parent, Tkinter.Tk):  # Launched from somewhere else
      return

    # Need the below code to fix the problem with Bash
    # where the command line was getting screwed up on exit.
    if os.name == 'posix':
      os.system('stty sane')

    sys.exit(0)

  def projectSummary(self):

    self.openPopup('project_summary', ProjectSummaryPopup)

  def registerAnalysis(self, isModal=False):

    # TBD: remove check when when attributes in released API
    if hasattr(self.analysisProfile, 'userName'):
      if isModal:
        popup = RegisterPopup(self, isModal=True)
        popup.destroy()
      else:
        self.openPopup('register_analysis', RegisterPopup)
    else:
      showError('Registration not yet active', 'Registration is not active until next release of software', parent=self)

  def updateAnalysis(self):
  
    popup = self.popups.get('update_analysis')
  
    if popup:
      popup.open()
    
    else:  
      topDir = getTopDirectory()
      url = 'file:' + topDir + LOCAL_HELP_DOC_DIR + '/popups/UpdatePopup.html'
      version = Copyright.version
      serverDirectory = 'ccpNmrUpdate%d.%d' % (version.major, version.minor)
      popup = UpdatePopup(self, serverLocation='www2.ccpn.ac.uk',
                          serverDirectory=serverDirectory,
                          dataFile='__UpdateAgentData.db',
                          helpUrl=url)
      self.popups['update_analysis'] = popup

  def showVersion(self):

    showInfo('Version', self.versionInfo, parent=self)

  def showAbout(self):
    
    topDir = getTopDirectory()
    url = 'file:' + topDir + LOCAL_HELP_DOC_DIR + '/about.html'
    self.webBrowser.open(url)

  def showHelp(self):

    topDir = getTopDirectory()
    url = 'file:' + topDir + LOCAL_HELP_DOC_DIR + '/index.html'
    self.webBrowser.open(url)

  def gotoPeak(self, window_name, peak, row=0, col=0):
    # Function appears to be deprecated
    # All peak navigation is now per frame

    popup = self.getWindowPopup(window_name)
    
    windowFrame = popup.activeWindowPane
    
    if windowFrame:
      windowFrame.gotoPeak(peak, row=row, col=col)

  # Note: position is assumed to be in same units as axisPanels in window
  #       it is a dictionary keyed on label
  #       axes with labels as keys have region changed so that position is at center
  def gotoPosition(self, windowPane, position, row=None, col=None, doOpen=True):
    
    windowFrame = windowPane.getWindowFrame()
    windowFrame.gotoPosition(position, row, col, doLift=doOpen)

  
  """
  def gotoPopupPosition(self, popup, windowPane, position, row=None, col=None, doOpen=True):

    if popup:
      popup.gotoPosition(windowPane, position, row=row, col=col, doLift=doOpen)
      
      if doOpen:
        popup.open()
  """

  # override Analysis version
  def initMacros(self):

    Analysis.initMacros(self)
    
    module = __name__.split('.')
    path   = '/'.join(module[:-1])
    
    analysisProfile = self.project.currentAnalysisProfile
    
    function = 'browseAtomsMacro'
    if analysisProfile.findFirstMacro(function=function) is None:
      m = analysisProfile.newMacro(name='browseAtoms',
                                   path=path,
                                   function=function, module=module[-1],
                                   shortcut='b', ordering=1)

    
    path += '/core'
    function = 'removeAllMarksAndRulers'
    
    macro = analysisProfile.findFirstMacro(function=function)
    if macro and macro.path[-4:] != 'core':
      macro.delete()
      macro = None
      
    if macro is None:
      m = analysisProfile.newMacro(name='clearMarksAndRulers',
                                   path=path,
                                   function=function, module='MarkBasic',
                                   shortcut='n', ordering=1)
      
    initWindowMacros(self.project)

  def turnDrawRequestsOff(self, window_name):

    popup = self.getWindowPopup(window_name)
    if popup:
      popup.turnDrawRequestsOff()

  def turnDrawRequestsOn(self, window_name, doDraw = True, doLift = False):

    popup = self.getWindowPopup(window_name)
    if popup:
      popup.turnDrawRequestsOn(doDraw=doDraw, doLift=doLift)

  def cloneWindow(self, window):

    newWindow = cloneSpectrumWindow(window)
    self.after_idle(lambda: self.copyViewProperties(window, newWindow))

  def copyViewProperties(self, window, newWindow):

    popup = self.getWindowPopup(newWindow.name)
    turnedOff = popup.turnDrawRequestsOff()
    try:
      copySpectrumWindowViewProperties(window, newWindow)
      # TBD: below only needed if axisMapping not default one
      for windowPane in newWindow.spectrumWindowPanes:
        for view in windowPane.spectrumWindowViews:
          self.changedViewAxisMapping(view)
      copySpectrumWindowRegions(window, newWindow)
    finally:
      if turnedOff: popup.turnDrawRequestsOn()

def browseAtomsMacro(argServer):

  argServer.parent.browseAtoms()

