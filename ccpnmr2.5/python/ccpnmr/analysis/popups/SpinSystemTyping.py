"""
======================COPYRIGHT/LICENSE START==========================

SpinSystemTyping.py: Part of the CcpNmr Analysis program

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

from ccpnmr.analysis.core.AssignmentBasic import (
    assignResonanceType,
    assignSpinSystemResidue,
    assignSpinSystemType,
    deassignResonance,
    getResonanceName,
    getShiftLists,
    removeSpinSystemResonance,
)
from ccpnmr.analysis.core.ChemicalShiftBasic import getShiftsChainProbabilities, lookupAtomProbability
from ccpnmr.analysis.core.MoleculeBasic import DEFAULT_ISOTOPES
from ccpnmr.analysis.popups.BasePopup import BasePopup
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showOkCancel
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix


class SpinSystemTypeScoresPopup(BasePopup):
    """
    **Predict Residue Type for a Spin System of Resonances**

    This tool aims to predict the residue type of a spin system based upon the
    chemical shifts of the resonances that it contains. The general principle is
    that different kinds of atoms in different kinds of residues have different
    observed distributions of chemical shifts. This system uses chemical shift
    distributions from the RefDB database, or otherwise from the BMRB where data
    is not available in RefDB. The observed chemical shifts of a spin system are
    compared to the per-atom distributions for each residue type and the residue
    types with the best matches are deemed to be more likely.

    This system can work with various levels of information, although the more
    information the better. Naturally, the more chemical shifts you have in a spin
    system then the better the prediction of type, and 13C resonances are more
    distinctive than 1H on the whole. Also, setting the atom type of a resonance
    can have a big influence on the type of residue predicted, for example knowing
    that a 13C resonance at 63 ppm is of type CB points  very strongly toward the
    residue being a serine. Atom type information can come from two sources: from
    a specific type assignment made by the user (via this popup or elsewhere) or
    by virtue of assignment in an experimental dimension that detects a
    restricted class of atom - e.g. 13C resonances in an HNCA experiment, assuming
    their shift matches, are of CA type as far as this prediction is concerned.
    Resonances that do not have a known atom type are compared with all of the
    unallocated types to find the combination that is most likely.

    The residue type prediction is based on the list of resonances displayed in the
    upper table. Here the user can see the chemical shifts (from the selected
    shift list) and any specific atom type setting. The user may set the atom type
    for any of the resonances, which would normally be done to reduce prediction
    ambiguity, by double-clicking in the "Atom Type" column.

    The lower table shows a ranked list of the probable residue types. All
    probability scores are normalised and represented as a percentage of the total
    of all scores, considering residue types in the selected chain. The type of a
    spin system may be set by clicking on a row of the lower table (hopefully a
    unique and high-scoring option) and  then selecting [Assign Spin System Type].
    If the user attempts to change the type of a spin system that is currently
    assigned to a specific residue then there is an opportunity to back out of the
    assignment, but otherwise any sequence specific information will be removed.

    **Caveats & Tips**

    It is assumed that the spectra from which the chemical shifts are derived are
    fairly well referenced.

    A type prediction will always be given, no matter how few resonances are
    present in a spin system. This system says which of the available types are
    most likely, *not how reliable* the prediction is; the latter depends largely
    on the amount of information present. The user should not for example make a
    judgement based only on amide resonances. Reliability scores will be added in
    the future.

    Rouge resonances in a spin system often adversely affect the prediction, if
    something is not genuinely in the spin system it should be removed.

    The system will never predict the residue type to be something that does not
    appear in the selected molecular chain. Thus, make sure the chain selection is
    appropriate for your prediction.

    **Reference**

    The residue type prediction method is not published independently but is very
    similar to the Bayesian method presented in: *Marin A, Malliavin TE, Nicolas P,
    Delsuc MA. From NMR chemical shifts to amino acid types: investigation of the
    predictive power carried by nuclei. J Biomol NMR. 2004 Sep;30(1):47-60.*

    One major difference however is that probabilities for resonances not being
    observed are not used. The CCPN prediction method is not only for complete
    spin systems and may be used at any time during the assignment process; here
    missing resonances are mostly due to the current assignment state and not such
    a useful indicator of residue type."""

    def __init__(self, parent, spinSystem=None, chain=None, *args, **kw):

        self.spinSystem = spinSystem
        self.shiftList = None
        self.resonance = None
        self.isotopes = ("1H", "13C", "15N")
        self.chain = chain
        self.ccpCode = None
        self.waiting = False
        self.atomTypes = {}

        self.project = parent.project
        self.guiParent = parent

        BasePopup.__init__(self, parent, title="Spin System Type Scores", **kw)

    def body(self, guiFrame):

        guiFrame.grid_columnconfigure(3, weight=1)

        row = 0
        label = Label(guiFrame, text="Spin System: ", grid=(row, 0))
        tipText = "Indicates which spin system the residue type prediction is done for"
        self.spinSystemLabel = Label(
            guiFrame, text="Serial:   Assignment:", grid=(row, 1), gridSpan=(1, 3), tipText=tipText
        )

        row += 1
        label = Label(guiFrame, text="Shift List: ", grid=(row, 0))
        tipText = (
            "Selects which shift list is the source of chemical shift information to make the residue type prediction"
        )
        self.shiftListPulldown = PulldownList(guiFrame, tipText=tipText, callback=self.setShiftList, grid=(row, 1))

        label = Label(guiFrame, text="Chain: ", grid=(row, 2))
        tipText = "Selects which molecular chain the prediction is for; sets prior probabilities for the various residue types"
        self.chainPulldown = PulldownList(guiFrame, self.changeChain, grid=(row, 3), tipText=tipText)

        row += 1
        labelFrame = LabelFrame(guiFrame, text="Resonances", grid=(row, 0), gridSpan=(1, 4))
        labelFrame.expandGrid(0, 0)

        self.atomTypePulldown = PulldownList(self, callback=self.setAtomType)

        editWidgets = [None, None, None, None, self.atomTypePulldown]
        editGetCallbacks = [None, None, None, None, self.getAtomType]
        editSetCallbacks = [None, None, None, None, self.setAtomType]

        tipTexts = [
            "The nuclear isotope type of the resonance within the current spin system",
            "The assignment annotation for the spin system resonance within the current spin system",
            "The chemical shift of the resonance in the stated shift list",
            "The weighted standard deviation of the resonance chemical shift",
            "The current atom type of the resonance; when set this helps refine residue type prediction",
        ]
        headingList = ["Isotope", "Name", "Shift\nValue", "Shift\nError", "Atom\nType"]
        self.resonanceMatrix = ScrolledMatrix(
            labelFrame,
            editWidgets=editWidgets,
            multiSelect=False,
            editGetCallbacks=editGetCallbacks,
            editSetCallbacks=editSetCallbacks,
            headingList=headingList,
            callback=self.selectResonance,
            grid=(0, 0),
            tipTexts=tipTexts,
        )

        tipTexts = [
            "Remove the selected resonance from the current spin system",
            "Remove residue type information from the current spin system",
            "Show a table of information for the  selected resonance, including a list of all peak dimension positions",
            "Show a table of the peaks to which the selected resonance is assigned",
        ]
        texts = ["Remove From\nSpin System", "Deassign\nResidue Type", "Resonance\nInfo", "Show\nPeaks"]
        commands = [self.removeResonance, self.deassignType, self.showResonanceInfo, self.showPeaks]
        buttonList = ButtonList(labelFrame, texts=texts, commands=commands, grid=(1, 0), tipTexts=tipTexts)
        self.resButtons = buttonList.buttons

        row += 1
        guiFrame.grid_rowconfigure(row, weight=1)
        labelFrame = LabelFrame(guiFrame, text="Type Scores", grid=(row, 0), gridSpan=(1, 4))
        labelFrame.expandGrid(0, 0)

        tipTexts = [
            "The ranking of the residue type possibility for the current spin system",
            "The CCPN residue code for the type",
            "The estimated percentage probability of the spin system being the residue type",
        ]
        headingList = ["Rank", "Ccp Code", "% Probability"]
        self.scoresMatrix = ScrolledMatrix(
            labelFrame, headingList=headingList, callback=self.selectCcpCode, grid=(0, 0), tipTexts=tipTexts
        )

        row += 1
        tipTexts = [
            "Assign the residue type of the current spin system to the kind selected in the lower table",
        ]
        texts = ["Assign Spin System Type"]
        commands = [self.assign]
        bottomButtons = UtilityButtonList(
            guiFrame,
            texts=texts,
            commands=commands,
            helpUrl=self.help_url,
            grid=(row, 0),
            gridSpan=(1, 4),
            tipTexts=tipTexts,
        )
        self.assignButton = bottomButtons.buttons[0]

        self.updateShiftLists()
        self.updateChains()
        self.getChainAtomTypes()
        self.update()

        self.curateNotifiers(self.registerNotify)

    def curateNotifiers(self, notifyFunc):

        for func in ("addResonance", "removeResonance", "setResonances", "delete", "setName"):
            notifyFunc(self.updateAfter, "ccp.nmr.Nmr.ResonanceGroup", func)

        for func in ("setResonanceSet", "addAssignName", "removeAssignName", "setAssignNames"):
            notifyFunc(self.updateAfter, "ccp.nmr.Nmr.Resonance", func)

        for func in ("__init__", "delete"):
            notifyFunc(self.updateAfter, "ccp.nmr.Nmr.ResonanceSet", func)
            notifyFunc(self.updateChains, "ccp.molecule.MolSystem.Chain", func)
            notifyFunc(self.updateShiftLists, "ccp.nmr.Nmr.ShiftList", func)

        for func in ("__init__", "setValue"):
            notifyFunc(self.updateShiftAfter, "ccp.nmr.Nmr.Shift", func)

    def getChainAtomTypes(self):

        doneResType = {}

        atomTypes = atomTypes = {}
        for isotope in DEFAULT_ISOTOPES.values():
            atomTypes[isotope] = set()

        if self.chain:
            for residue in self.chain.residues:
                molResidue = residue.molResidue
                ccpCode = molResidue.ccpCode
                molType = molResidue.molType
                key = "%s:%s:%s" % (ccpCode, molResidue.linking, molResidue.descriptor)

                if doneResType.get(key):
                    continue

                doneResType[key] = True

                for atom in residue.atoms:
                    chemAtom = atom.chemAtom
                    element = chemAtom.elementSymbol
                    isotope = DEFAULT_ISOTOPES.get(element)

                    if not isotope:
                        continue

                    atomTypes[isotope].add((ccpCode, atom.name, molType))

        self.atomTypes = atomTypes

    def getAtomType(self, resonance):

        index = 0
        atomNames = set(
            [
                "<None>",
            ]
        )

        assignNames = resonance.assignNames
        if assignNames:
            for atomName in assignNames:
                atomNames.add(atomName)

            if len(assignNames) > 1:
                orig = ",".join(assignNames)
                atomNames.add(orig)

        shift = resonance.findFirstShift(parentList=self.shiftList)

        if shift and self.chain:
            project = self.project
            atomTypes = self.atomTypes.get(resonance.isotopeCode, [])

            for ccpCode, atomName, molType in atomTypes:
                prob = lookupAtomProbability(project, ccpCode, atomName, shift.value, molType=molType)

                if prob >= 0.001:
                    atomNames.add(atomName)

        atomNames = list(atomNames)
        atomNames.sort()

        if resonance.assignNames:
            orig = ",".join(assignNames)
            index = atomNames.index(orig)

        atomNameObjs = atomNames[:]
        atomNameObjs[0] = None

        self.atomTypePulldown.setup(atomNames, atomNameObjs, index)

    def setAtomType(self, obj):

        atomNameStr = self.atomTypePulldown.getObject()

        if self.resonance:
            if atomNameStr:
                atomNames = atomNameStr.split(",")
                assignResonanceType(self.resonance, assignNames=atomNames)
            else:
                assignResonanceType(self.resonance, assignNames=None)

    def removeResonance(self):

        if self.resonance and self.spinSystem and (self.resonance in self.spinSystem.resonances):
            if showOkCancel("Confirm", "Really remove resonance from spin system?", parent=self):
                self.spinSystem.codeScoreDict = {}
                deassignResonance(self.resonance, clearAssignNames=False)
                removeSpinSystemResonance(self.spinSystem, self.resonance)

    def showResonanceInfo(self):

        if self.resonance:
            self.guiParent.browseResonanceInfo(self.resonance)

    def showPeaks(self):

        if self.resonance:
            peaksDict = {}
            for contrib in self.resonance.peakDimContribs:
                peaksDict[contrib.peakDim.peak] = 1

            peaks = peaksDict.keys()
            if len(peaks) > 0:
                self.guiParent.viewPeaks(peaks)

    def deassignType(self):

        if self.spinSystem:
            residue = self.spinSystem.residue

            if residue:
                resText = "%d%s" % (residue.seqCode, residue.ccpCode)
                msg = "Spin system assigned to %s. Continue and deassign residue?"
                if showOkCancel("Warning", msg % resText, parent=self):
                    assignSpinSystemResidue(self.spinSystem, None)
                    assignSpinSystemType(self.spinSystem, None)

            else:
                assignSpinSystemType(self.spinSystem, None)

    def assign(self):

        if self.spinSystem and self.ccpCode:
            if self.spinSystem.residue and (self.spinSystem.residue.ccpCode != self.ccpCode):
                resText = "%d%s" % (self.spinSystem.residue.seqCode, self.spinSystem.residue.ccpCode)
                msg = "Spin system is already assigned to %s. Continue?"
                if showOkCancel("Warning", msg % resText, parent=self):
                    assignSpinSystemResidue(self.spinSystem, residue=None)

                else:
                    return

            if self.spinSystem.ccpCode != self.ccpCode:
                assignSpinSystemType(self.spinSystem, self.ccpCode, "protein")
                self.update()

    def getChains(self):

        chains = []
        if self.project:
            for molSystem in self.project.sortedMolSystems():
                for chain in molSystem.sortedChains():
                    if chain.molecule.molType in ("protein", None):
                        text = "%s:%s" % (molSystem.code, chain.code)
                        chains.append([text, chain])

        return chains

    def changeChain(self, chain):

        if self.chain is not chain:
            self.chain = chain
            self.getChainAtomTypes()
            self.updateAfter()

    def updateChains(self, *chain):

        data = self.getChains()
        names = [x[0] for x in data]
        chains = [x[1] for x in data]
        chain = self.chain
        index = 0

        if chains:
            if chain not in chains:
                chain = chains[0]

            index = chains.index(chain)

        if chain is not self.chain:
            self.chain = chain
            self.getChainAtomTypes()
            self.updateAfter()

        self.chainPulldown.setup(names, chains, index)

    def updateShiftLists(self, *opt):

        shiftLists = getShiftLists(self.nmrProject)
        shiftList = self.shiftList
        names = ["%s [%d]" % (x.name or "<No name>", x.serial) for x in shiftLists]
        index = 0

        if names:
            if shiftList not in shiftLists:
                shiftList = shiftLists[0]

            index = shiftLists.index(shiftList)

        if shiftList is not self.shiftList:
            self.shiftList = shiftList
            self.updateAfter()

        self.shiftListPulldown.setup(names, shiftLists, index)

    def setShiftList(self, shiftList):

        if self.shiftList is not shiftList:
            self.shiftList = shiftList
            self.updateAfter()

    def updateButtons(self):

        if self.resonance:
            self.resButtons[0].enable()
            self.resButtons[2].enable()
            self.resButtons[3].enable()

        else:
            self.resButtons[0].disable()
            self.resButtons[2].disable()
            self.resButtons[3].disable()

        if self.spinSystem and (self.spinSystem.residue or self.spinSystem.ccpCode):
            self.resButtons[1].enable()
        else:
            self.resButtons[1].disable()

        if self.ccpCode and self.spinSystem:
            self.assignButton.enable()
        else:
            self.assignButton.disable()

    def selectResonance(self, resonance, row, col):

        self.resonance = resonance
        self.updateButtons()

    def selectCcpCode(self, ccpCode, row, col):

        self.ccpCode = ccpCode
        self.assignButton.enable()

    def clearSpinSystemCache(self):

        if self.spinSystem:
            self.spinSystem.sstTypes = []
            self.spinSystem.ssScore = None
            self.spinSystem.codeScoreDict = {}

    def updateShiftAfter(self, shift):

        if shift.parentList is not self.shiftList:
            return

        resonance = shift.resonance

        if resonance.resonanceGroup is not self.spinSystem:
            return

        self.clearSpinSystemCache()

        self.updateAfter()

    def updateAfter(self, obj=None):

        if obj and obj is self.spinSystem:
            self.clearSpinSystemCache()
            if obj.isDeleted:
                self.spinSystem = None

        if self.spinSystem:
            if obj is not None:
                if obj.className == "ResonanceSet":
                    for resonance in obj.resonances:
                        if resonance.resonanceGroup is self.spinSystem:
                            break

                    else:
                        return

                elif obj.className == "Resonance":
                    if obj.resonanceGroup is not self.spinSystem:
                        return

            self.clearSpinSystemCache()

        if self.waiting:
            return

        else:
            self.waiting = True
            self.after_idle(self.update)

    def update(self, spinSystem=None, chain=None, shiftList=None):

        if spinSystem is not None:
            if spinSystem is not self.spinSystem:
                self.ccpCode = None
                self.spinSystem = spinSystem

            if chain:
                self.chain = chain
            else:
                self.chain = None

            self.updateChains()
            self.getChainAtomTypes()

        if shiftList:
            self.shiftList = shiftList
            self.updateShiftLists()

        if self.resonance:
            if not self.spinSystem:
                self.resonance = None

            elif self.resonance.resonanceGroup is not self.spinSystem:
                self.resonance = None

        self.updateButtons()

        textMatrix = []
        objectList = []
        if self.spinSystem:
            if self.spinSystem.residue:
                resText = "Assignment: %d%s" % (self.spinSystem.residue.seqCode, self.spinSystem.residue.ccpCode)
            elif self.spinSystem.ccpCode:
                resText = "Type: %s" % self.spinSystem.ccpCode
            elif self.spinSystem.name:
                resText = "Name: %s" % self.spinSystem.name
            else:
                resText = "Unassigned"

            self.spinSystemLabel.set("Serial: %d  %s" % (self.spinSystem.serial, resText))

            for resonance in self.spinSystem.resonances:
                shift = resonance.findFirstShift(parentList=self.shiftList)
                if shift:
                    datum = [
                        resonance.isotopeCode,
                        getResonanceName(resonance),
                        shift.value,
                        shift.error,
                        "/".join(resonance.assignNames),
                    ]

                    objectList.append(resonance)
                    textMatrix.append(datum)

        self.resonanceMatrix.update(textMatrix=textMatrix, objectList=objectList)

        textMatrix = []
        objectList = []
        colorMatrix = []
        if self.spinSystem and self.chain and self.shiftList:
            shifts = []
            for resonance in self.spinSystem.resonances:
                if resonance.isotopeCode in self.isotopes:
                    shift = resonance.findFirstShift(parentList=self.shiftList)
                    if shift:
                        shifts.append(shift)

            scores = getShiftsChainProbabilities(shifts, self.chain)
            total = sum(scores.values())

            scoreList = []

            if total:
                ccpCodes = self.getCcpCodes(self.chain)
                baseLevel = 100.0 / len(ccpCodes)
                for ccpCode in ccpCodes:
                    scoreList.append((100.0 * scores[ccpCode] / total, ccpCode))

            scoreList.sort()
            scoreList.reverse()

            i = 0
            for score, ccpCode in scoreList:
                if not score:
                    continue

                i += 1
                datum = [i, ccpCode, score]

                if score >= min(100.0, 5 * baseLevel):
                    color = "#80ff80"
                elif score > 2 * baseLevel:
                    color = "#ffff80"
                elif score > baseLevel:
                    color = "#ffc080"
                else:
                    color = "#ff8080"

                colors = [color, color, color]
                objectList.append(ccpCode)
                textMatrix.append(datum)
                colorMatrix.append(colors)

        self.scoresMatrix.update(textMatrix=textMatrix, colorMatrix=colorMatrix, objectList=objectList)
        self.waiting = False

    def getCcpCodes(self, chain):

        ccpDict = {}
        for residue in chain.residues:
            ccpCode = residue.ccpCode

            # if (ccpCode == 'Cys') and (residue.descriptor == 'link:SG'):
            #  ccpCode = 'Cyss'

            ccpDict[ccpCode] = True

        ccpCodes = list(ccpDict.keys())
        ccpCodes.sort()

        return ccpCodes

    def destroy(self):

        self.curateNotifiers(self.unregisterNotify)

        BasePopup.destroy(self)

