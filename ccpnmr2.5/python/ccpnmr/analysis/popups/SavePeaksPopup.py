import os

from ccpnmr.analysis.popups.BasePopup import BasePopup
from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showInfo
from memops.gui.PulldownList import PulldownList


class SavePeaksPopup(BasePopup):
    """
    **Save peaks to an external file (.tab, .nef, .peaks)**

    This popup allows exporting a selected peak list to standard external
    formats, including NMRdraw (.tab), NEF (.nef), or XEASY (.peaks).
    """

    def __init__(self, parent, *args, **kw):

        self.peakList = None

        BasePopup.__init__(self, parent=parent, title="Save Peaks to .tab, .nef, or .peaks", **kw)

    def body(self, master):

        self.geometry("640x180")
        master.grid_columnconfigure(1, weight=1)
        for n in range(3):
            master.grid_rowconfigure(n, weight=1)

        row = 0
        label = Label(master, text="Peak List:")
        label.grid(row=row, column=0, sticky="e")
        tipText = "The peak list to save/export"
        self.peakListPulldown = PulldownList(master, callback=self.setPeakList, tipText=tipText)
        self.peakListPulldown.grid(row=row, column=1, sticky="w")

        row = row + 1
        label = Label(master, text="Output File:")
        label.grid(row=row, column=0, sticky="e")
        tipText = "Choose the path and file format (.tab, .nef, .peaks) to save to"
        self.fileEntry = Entry(master, tipText=tipText)
        self.fileEntry.grid(row=row, column=1, sticky="ew")
        tipText = "Browse and choose output file path"
        browseButton = Button(master, text="Browse...", command=self.selectFile, tipText=tipText)
        browseButton.grid(row=row, column=2, sticky="w")

        row = row + 1
        texts = ["Save Peaks"]
        commands = [self.savePeaks]
        tipTexts = [
            "Write the selected peak list to the specified output file format",
        ]
        self.buttons = UtilityButtonList(
            master, texts=texts, doClone=False, tipTexts=tipTexts, commands=commands, helpUrl=self.help_url
        )
        self.buttons.grid(row=row, column=0, columnspan=3, sticky="ew")

        self.curateNotifiers(self.registerNotify)
        self.updatePeakLists()

    def destroy(self):

        self.curateNotifiers(self.unregisterNotify)
        BasePopup.destroy(self)

    def curateNotifiers(self, notifyFunc):

        for clazz in ("Experiment", "DataSource", "PeakList"):
            for func in ("__init__", "delete", "setName"):
                notifyFunc(self.updateNotifier, "ccp.nmr.Nmr.%s" % clazz, func)

    def updateNotifier(self, *extra):

        self.updatePeakLists()

    def updatePeakLists(self):

        nmrProject = self.nmrProject
        if nmrProject is None:
            self.peakListPulldown.setup([], [], 0)
            return

        peakLists = []
        for experiment in nmrProject.experiments:
            for dataSource in experiment.dataSources:
                for pl in dataSource.peakLists:
                    if not pl.isDeleted:
                        peakLists.append(pl)

        if self.peakList is not None and (self.peakList.isDeleted or self.peakList not in peakLists):
            self.peakList = None

        if peakLists:
            if self.peakList is None:
                self.peakList = peakLists[0]
            index = peakLists.index(self.peakList)
            names = [
                "%s:%s:%s (%d peaks)"
                % (
                    pl.dataSource.experiment.name,
                    pl.dataSource.name,
                    pl.name or "List %d" % pl.serial,
                    len(pl.peaks),
                )
                for pl in peakLists
            ]
        else:
            self.peakList = None
            index = 0
            names = []

        self.peakListPulldown.setup(names, peakLists, index)

    def setPeakList(self, peakList):

        self.peakList = peakList

    def selectFile(self):

        directory = os.getcwd()
        if self.peakList is not None and self.peakList.dataSource.dataStore is not None:
            location = self.peakList.dataSource.dataStore.fullPath
            if location:
                directory = os.path.dirname(location)
        popup = FileSelectPopup(
            self,
            directory=directory,
            file_types=[
                FileType("NMRdraw .tab files", ["*.tab"]),
                FileType("NEF peak files", ["*.nef"]),
                FileType("XEASY .peaks files", ["*.peaks"]),
                FileType("All files", ["*"]),
            ],
        )
        fileName = popup.getFile()
        popup.destroy()
        if fileName:
            self.fileEntry.set(fileName)

    def savePeaks(self):

        peakList = self.peakListPulldown.getObject()
        if peakList is None:
            showError("No peak list", "Select a peak list to save first", parent=self)
            return

        fileName = self.fileEntry.get().strip()
        if not fileName:
            showError("No file", "Choose an output file path first", parent=self)
            return

        directory = os.path.dirname(os.path.abspath(fileName))
        if not os.path.isdir(directory):
            try:
                os.makedirs(directory)
            except Exception as e:
                showError("Error", f"Failed to create directory {directory}: {e}", parent=self)
                return

        ext = os.path.splitext(fileName)[1].lower()
        try:
            if ext == ".tab":
                from ccpnmr.analysis.core.PeakListExport import exportTabPeaks
                exportTabPeaks(peakList, fileName)
            elif ext == ".nef":
                from ccpnmr.analysis.core.PeakListExport import exportNefPeaks
                exportNefPeaks(peakList, fileName)
            elif ext == ".peaks":
                from ccpnmr.analysis.core.PeakListExport import exportXeasyPeaks
                exportXeasyPeaks(peakList, fileName)
            else:
                showError("Unknown Format", "File extension must be .tab, .nef, or .peaks", parent=self)
                return
        except Exception as e:
            showError("Save Failed", f"Failed to save peak list: {e}", parent=self)
            return

        showInfo("Success", f"Successfully saved peak list to {os.path.basename(fileName)}", parent=self)
