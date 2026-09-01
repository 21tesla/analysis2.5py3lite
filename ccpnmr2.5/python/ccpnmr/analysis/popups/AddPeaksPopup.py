"""
======================COPYRIGHT/LICENSE START==========================

AddPeaksPopup.py: Part of the CcpNmr Analysis program

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

from ccpnmr.analysis.core.AssignmentBasic import getShiftLists
from ccpnmr.analysis.core.PeakListImport import importTabPeaks
from ccpnmr.analysis.popups.BasePopup import BasePopup
from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showInfo
from memops.gui.PulldownList import PulldownList


class AddPeaksPopup(BasePopup):
    """
    **Add peaks from an NMRdraw .tab peak list file**

    The purpose of this dialog is to add the peaks of an NMRdraw .tab
    peak list file to a new peak list of the selected spectrum, creating
    and assigning resonances for the chemical shifts of the file (see
    str(Add Peaks)_.  The file's ASS column assigns each peak to an atom
    (e.g. W81-HE1); peaks that cannot be mapped to an atom in the
    molecule are ignored with a notification on stdout.  If a resonance
    assignment already exists for an atom, it is kept (and the new peak
    dimension left unassigned) unless the Overwrite resonance box is
    checked, in which case the existing resonance is re-used for the new
    peak dimension as well.

    The task cannot execute without a molecule (a MolSystem) from which
    to resolve the atoms - the Add Molecule button opens the Molecule
    dialog.

    .. _str(Add Peaks): (menu item Peaks > Add Peaks)
    """

    def __init__(self, parent, *args, **kw):

        self.spectrum = None
        self.shiftList = None
        self.fileChanged = False

        BasePopup.__init__(self, parent=parent, title="Add Peaks from NMRdraw .tab or NEF file", **kw)

    def body(self, master):

        self.geometry("640x220")
        master.grid_columnconfigure(1, weight=1)
        for n in range(5):
            master.grid_rowconfigure(n, weight=1)

        row = 0
        label = Label(master, text="Spectrum:")
        label.grid(row=row, column=0, sticky="e")
        tipText = "The spectrum whose new peak list will be populated with the file's peaks"
        self.spectrumPulldown = PulldownList(master, callback=self.setSpectrum, tipText=tipText)
        self.spectrumPulldown.grid(row=row, column=1, sticky="w")

        row = row + 1
        label = Label(master, text="File:")
        label.grid(row=row, column=0, sticky="e")
        tipText = "Choose an NMRdraw .tab or NEF peak list file"
        self.fileEntry = Entry(master, tipText=tipText)
        self.fileEntry.grid(row=row, column=1, sticky="ew")
        tipText = "Browse for an NMRdraw .tab or NEF peak list file"
        browseButton = Button(master, text="Browse...", command=self.selectFile, tipText=tipText)
        browseButton.grid(row=row, column=2, sticky="w")

        row = row + 1
        label = Label(master, text="Peak list name:")
        label.grid(row=row, column=0, sticky="e")
        tipText = "The name of the new peak list (defaults to the root of the file name)"
        self.nameEntry = Entry(master, tipText=tipText)
        self.nameEntry.grid(row=row, column=1, columnspan=2, sticky="ew")

        row = row + 1
        label = Label(master, text="Resonance list:")
        label.grid(row=row, column=0, sticky="e")
        tipText = "The resonance (shift) list to populate with the new resonances"
        self.shiftListPulldown = PulldownList(master, callback=self.setShiftList, tipText=tipText)
        self.shiftListPulldown.grid(row=row, column=1, sticky="w")
        tipText = "If an assignment already exists, overwrite it (re-using the existing resonance for the new peak) instead of keeping it"
        self.overwriteCheck = CheckButton(
            master,
            text="Overwrite resonance",
            selected=False,
            tipText=tipText,
        )
        self.overwriteCheck.grid(row=row, column=2, sticky="w")

        row = row + 1
        self.moleculeLabel = Label(master, text="", tipText="Molecules available for the atom assignments")
        self.moleculeLabel.grid(row=row, column=1, columnspan=2, sticky="w")

        texts = ["Add Peaks", "Add Molecule"]
        commands = [self.addPeaks, self.addMolecule]
        tipTexts = [
            "Add the file's peaks to a new peak list and assign the created resonances to the atoms",
            "Open the Molecule dialog - a molecule is required to resolve the file's atom assignments",
        ]
        self.buttons = UtilityButtonList(
            master, texts=texts, doClone=False, tipTexts=tipTexts, commands=commands, helpUrl=self.help_url
        )
        self.buttons.grid(row=row, column=0, columnspan=3, sticky="ew")

        self.curateNotifiers(self.registerNotify)
        self.updateSpectrum()
        self.updateMoleculeInfo()

    def destroy(self):

        self.curateNotifiers(self.unregisterNotify)

        BasePopup.destroy(self)

    def curateNotifiers(self, notifyFunc):

        for clazz in ("Experiment", "DataSource", "ShiftList"):
            for func in ("__init__", "delete", "setName"):
                notifyFunc(self.updateNotifier, "ccp.nmr.Nmr.%s" % clazz, func)
        for func in ("__init__", "delete", "setName"):
            notifyFunc(self.updateMoleculeNotifier, "ccp.molecule.MolSystem.MolSystem", func)

    def updateNotifier(self, *extra):

        self.updateSpectrum()
        self.updateShiftLists()

    def updateMoleculeNotifier(self, *extra):

        self.updateMoleculeInfo()

    def updateSpectrum(self, spectrum=None):

        if not spectrum and self.spectrum is not None and self.spectrum.isDeleted:
            self.spectrum = None

        spectra = self.parent.getSpectra()
        if spectra:
            if self.spectrum is None or self.spectrum not in spectra:
                self.spectrum = spectra[0]
            index = spectra.index(self.spectrum)
            names = ["%s:%s" % (s.experiment.name, s.name) for s in spectra]
        else:
            self.spectrum = None
            index = 0
            names = []

        self.spectrumPulldown.setup(names, spectra, index)

        self.updateShiftLists()

    def updateShiftLists(self):

        nmrProject = self.nmrProject
        if nmrProject is None:
            self.shiftListPulldown.setup([], [], 0)
            return

        shiftLists = getShiftLists(nmrProject)
        if self.shiftList is not None and (self.shiftList.isDeleted or self.shiftList not in shiftLists):
            self.shiftList = None
        if shiftLists:
            if self.shiftList is None:
                # prefer the experiment's own (working) shift list
                if self.spectrum is not None and self.spectrum.experiment.shiftList in shiftLists:
                    self.shiftList = self.spectrum.experiment.shiftList
                else:
                    self.shiftList = shiftLists[0]
            index = shiftLists.index(self.shiftList)
        else:
            self.shiftList = None
            index = 0

        names = ["%s" % (sl.name or "ShiftList %d" % sl.serial) for sl in shiftLists]
        self.shiftListPulldown.setup(names, shiftLists, index)

    def updateMoleculeInfo(self):

        molSystems = [ms for ms in self.project.sortedMolSystems() if not ms.isDeleted]
        if molSystems:
            names = [
                "%s (%d residues)" % (ms.name, sum(len(c.sortedResidues()) for c in ms.sortedChains())) for ms in molSystems
            ]
            self.moleculeLabel.set("Molecule(s): %s" % ", ".join(names))
        else:
            self.moleculeLabel.set("No molecule defined - use Add Molecule first")

    def setSpectrum(self, spectrum):

        self.spectrum = spectrum
        self.shiftList = None
        self.updateShiftLists()

    def setShiftList(self, shiftList):

        self.shiftList = shiftList

    def selectFile(self):

        directory = os.getcwd()
        if self.spectrum is not None and self.spectrum.dataStore is not None:
            location = self.spectrum.dataStore.dataLocation
            if location:
                directory = os.path.dirname(location)
        popup = FileSelectPopup(
            self,
            directory=directory,
            file_types=[
                FileType("NMRdraw .tab files", ["*.tab"]),
                FileType("NEF peak files", ["*.nef"]),
                FileType("All files", ["*"]),
            ],
        )
        fileName = popup.getFile()
        popup.destroy()
        if fileName:
            self.fileEntry.set(fileName)
            self.setDefaultName(os.path.basename(fileName))

    def setDefaultName(self, fileName):

        if not self.nameEntry.get():
            self.nameEntry.set(os.path.splitext(fileName)[0])

    def addPeaks(self):

        if self.spectrum is None:
            showError("No spectrum", "No spectrum available to take the new peak list", parent=self)
            return

        if not self.project.sortedMolSystems():
            showError("No molecule", "A molecule is required to resolve the atom assignments - use Add Molecule first", parent=self)
            return

        fileName = self.fileEntry.get()
        if not fileName or not os.path.exists(fileName):
            showError("No file", "Choose an existing NMRdraw .tab or NEF file first", parent=self)
            return

        listName = self.nameEntry.get().strip()
        if not listName:
            listName = os.path.splitext(os.path.basename(fileName))[0]
            self.nameEntry.set(listName)

        report = importTabPeaks(
            self.project,
            fileName,
            self.spectrum,
            listName=listName,
            shiftList=self.shiftList,
            overwrite=self.overwriteCheck.getSelected(),
        )

        if report["error"]:
            showError("Add Peaks failed", report["error"], parent=self)
            return

        showInfo(
            "Added peaks",
            "Peak list %r: %d peak(s) added (%d unassigned, %d row(s) skipped), %d resonance(s) created, "
            "%d assignment(s) applied (%d overwritten, %d kept)"
            % (
                listName,
                report["peaksAdded"],
                report["peaksUnassigned"],
                report["peaksSkipped"],
                report["resonancesCreated"],
                report["assignmentsApplied"],
                report["assignmentsOverwritten"],
                report["assignmentsKept"],
            ),
            parent=self,
        )

    def addMolecule(self):

        self.parent.editMolSystems()
