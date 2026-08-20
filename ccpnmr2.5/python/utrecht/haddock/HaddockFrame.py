#!/usr/bin/env python

"""
=========================================================================
Package:    - Code for the Graphical User Interface fontend to the 
              Haddock model package in the CCPN data model. 
            - Code for the export of a Haddock compatible project. A 
              Haddock compatible project can either be a parameter file
              ready for submission to the Haddock webserver or a
              directory structure with necessary files for use with a 
              localy installed version of Haddock.

Dependencies: The CCPN Haddock package requires CCPN data model version
              2.0 or higher. The export of a webserver compatible 
              parameter file requires Haddock webserver version 2.1 or 
              higher and a valid user account. The export of a 'classic' 
              Haddock project requires Haddock version 2.0 or higher.

Copyright and License information:
              The Haddock data model as implemented in the CCPN data
              model as well as the use of CCPN GUI code elements is 
              licenced to the CCPN Projects (Copyright (C) 2008) and
              distributed under the terms of the GNU Lesser General
              Public License.
            
              The Haddock project export code as well as the use of 
              Haddock software is covert in the Haddock License
              agreement (Copyright (C) 2008 Haddock Project, Bijvoet
              Center for Biomolecular Research, Utrecht University,
              The Netherlands).

GNU LGPL:        This library is free software; you can redistribute it 
              and/or modify it under the terms of the GNU Lesser General 
              Public License as published by the Free Software 
              Foundation; either version 2.1 of the License, or (at 
              your option) any later version.
 
              A copy of this license can be found in LGPL.license
 
              This library is distributed in the hope that it will be 
              useful, but WITHOUT ANY WARRANTY; without even the implied 
              warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
              PURPOSE. See the GNU Lesser General Public License for 
              more details.
 
              You should have received a copy of the GNU Lesser General 
              Public License along with this library; if not, write to 
              the Free Software Foundation, Inc., 59 Temple Place, Suite 
              330, Boston, MA 02111-1307 USA.

Information:  For further information regarding CCPN, please contact:
              - CCPN website (http://www.ccpn.ac.uk/)
              - email: ccpn@bioc.cam.ac.uk
              
              For further information regarding Haddock, please contact
              Alexandre M.J.J. Bonvin:
              - http://haddock.chem.uu.nl
              - email: a.m.j.j.bonvin@uu.nl    

Citing:          If you are using this software for academic purposes, we 
                suggest quoting the following references:

              For CCPN:    
              Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
              Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. 
              Ionides and Ernest D. Laue (2005). A framework for 
              scientific data modeling and automated software development. 
              Bioinformatics 21, 1678-1684.
            
              For Haddock:
              Cyril Dominguez, Rolf Boelens and Alexandre M.J.J. Bonvin 
              (2003). HADDOCK: a protein-protein docking approach based 
              on biochemical and/or biophysical information. 
              J. Am. Chem. Soc. 125, 1731-1737.
            
              S.J. de Vries, A.D.J. van Dijk, M. Krzeminski, M. van Dijk, 
              A. Thureau, V. Hsu, T. Wassenaar and A.M.J.J. Bonvin (2007) 
              HADDOCK versus HADDOCK: New features and performance of 
              HADDOCK2.0 on the CAPRI targets. 
              Proteins: Struc. Funct. & Bioinformatic 69, 726-733.    
=========================================================================
"""

import os, sys, Tkinter

from memops.gui.Button          import Button
from memops.gui.ButtonList      import ButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.DataEntry       import askString
from memops.gui.Entry           import Entry
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FileSelect      import FileType
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.IntEntry        import IntEntry
from memops.gui.Label           import Label
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.MultiWidget     import MultiWidget
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.TabbedFrame     import TabbedFrame
from memops.gui.Util            import createDismissHelpButtonList

from memops.editor.BasePopup    import BasePopup

from HaddockBasic                 import setPartnerChains, setRunConstraintSet, copyRun
from HaddockBasic                import getStructureFromFile, addRdcParam, addDaniParam
from HaddockExportClassic         import exportClassic
from HaddockExportParam         import exportParam
from HaddockImportRunCns        import runCnsImporter
from HaddockLocal                 import *
from HaddockRdcPopup            import HaddockRdcPopup
from HaddockDaniPopup            import HaddockDaniPopup
from EditSymmetryPopup             import EditSymmetryPopup
from HaddockServerUpload        import HaddockServerUpload

activeColour  = '#B080F0'
passiveColour = '#A0E0A0'
flexibleColor = '#F0A0A0'
semiflexColor = '#F0F0A0'
headingColor  = '#80C080'

energyTermTypes = {'UNAMBIG':'DistanceConstraintList',
                   'AMBIG':'DistanceConstraintList',
                   'DIHEDRAL':'DihedralConstraintList',
                   'HBOND':'HBondConstraintList',
                   'RDC':'RdcConstraintList',
                   'DANI':None}

class HaddockFrame(Frame):
    
    def __init__(self, parent, ccpnProject=None):
   
        self.parent           = parent
        self.ccpnProject      = ccpnProject
        self.molPartner       = None
        self.residue          = None
        self.run              = None
        self.energyTerm       = None
        self.constraintSet      = None
        self.symmetrySet        = None
        self.symmetryPopup      = None
        self.haddockUpload     = None
        self.waiting          = False
        self.hProject = None
        
        Frame.__init__(self, parent=parent)

        # Generic widgets
        self.stringEntry  = Entry(self.parent, returnCallback=self.setValue)
        self.intEntry     = IntEntry(self.parent, returnCallback=self.setValue)
        self.floatEntry   = FloatEntry(self.parent, returnCallback=self.setValue)
        self.pulldownMenu = PulldownMenu(self.parent, callback=self.setValue, do_initial_callback=False)
    
        # Ensure that the first row and column in frame expand
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1, minsize=400)

        options = ['     General     ',' Import/Export  ','Molecular Partners','      Residues     ','       Runs       ','    Restraints   ']
        self.tabbedFrame = TabbedFrame(self, options=options, toggleOff=False, selected=0)
        self.tabbedFrame.grid(row=0, column=0, sticky='nsew')

        genFrame, ctrlFrame, molFrame, resFrame, runFrame, conFrame = self.tabbedFrame.frames    

        # GENERAL FRAME        
        genFrame.grid_columnconfigure(1, weight=1)
        genFrame.grid_rowconfigure(1, weight=1)

        label = Label(genFrame, text='HADDOCK Project: ')
        label.grid(row=0, column=0, sticky='nw')

        # Pulldown menu to choose projects
        self.hProjectPulldown = PulldownMenu(genFrame, callback=self.chooseHaddockProject)
        self.hProjectPulldown.grid(row=0, column=1, sticky='nw') 

        # setup generic table headings, justification and widget getters/setters
        headingList      = ['Parameter','Value','Description']
        justifyList      = ['center','center', 'left']
        editWidgets      = [None, True, None]
        editGetCallbacks = [None, self.getValue, None]
        editSetCallbacks = [None, self.setValue, None]
   
        self.generalMatrix = ScrolledMatrix(genFrame, headingList=headingList,
                                            justifyList=justifyList,
                                            editSetCallbacks=editSetCallbacks,
                                            editGetCallbacks=editGetCallbacks, 
                                            editWidgets=editWidgets,
                                            multiSelect=False, initialRows=10,
                                            passSelfToCallback=True,
                                            callback=self.selectRowObj)

        self.generalMatrix.grid(row=1, column=0, columnspan=2, sticky='nsew')
        self.generalMatrix.refreshFunc = self.updateGeneral

        # EXPORT FRAME
        ctrlFrame.grid_columnconfigure(0, weight=1)
        ctrlFrame.grid_rowconfigure(1, weight=1)
        
        Frame1 = LabelFrame(ctrlFrame, text='Project export options:')
        Frame1.grid(row=0,column=0,sticky='ew')
        Frame1.grid_columnconfigure(2, weight=1)
        
        classicExport = Button(Frame1,text='Export Classic Project',command=self.expClassic)
        classicExport.grid(row=0,column=0,sticky='w',pady=10)

        label2 = Label(Frame1, text='  Export a classical HADDOCK directory structure containing a\nnew.html,ensemble files, restraint files and pdb files')
        label2.grid(row=0, column=1, sticky='w',pady=10)

        paramExport = Button(Frame1,text='Export Parameter File',command=self.expParam)
        paramExport.grid(row=1, column=0, sticky='w')

        label3 = Label(Frame1, text='  Export a new-style HADDOCK parameter file for use with\nthe HADDOCK server. Supports 2 molecules only.')
        label3.grid(row=1, column=1, sticky='w')
        
        Frame2 = LabelFrame(ctrlFrame, text='Project import options:')
        Frame2.grid(row=1,column=0,sticky='nsew')
        Frame2.grid_columnconfigure(0, weight=1)
        Frame2.grid_rowconfigure(0, weight=1)
    
        cnsImport = Button(Frame2,text='  Import run.cns File   ',command=self.impRunCns)
        cnsImport.grid(row=0,column=0, sticky='nw', pady=10)
        
        label5 = Label(Frame2, text='  Import run specific settings from a Haddock run.cns parameter file')
        label5.grid(row=0,column=1, sticky='nw', pady=10)
        
        Frame3 = LabelFrame(ctrlFrame, text='Run on HADDOCK server:')
        Frame3.grid(row=2,column=0,sticky='ew')
        Frame3.grid_columnconfigure(2, weight=1)
    
        upload = Button(Frame3,text=' Run Haddock server ',command=self.uploadToServer)
        upload.grid(row=0, column=0, sticky='w', pady=10)
        
        label7 = Label(Frame3, text='  Requires a valid HADDOCK server username and password')
        label7.grid(row=0, column=1, sticky='w', pady=10)

        # MOLECULAR PARTNER FRAME
        molFrame.grid_columnconfigure(0, weight=1)
        molFrame.grid_rowconfigure(0, weight=1)

        self.molSysPulldown   = PulldownMenu(self, callback=self.setMolSystem, do_initial_callback=False)
        self.modelPulldown = PulldownMenu(self, callback=self.setModel, do_initial_callback=False)
        self.chainSelect = MultiWidget(self, CheckButton, callback=self.setChains, minRows=0, useImages=False)

        headingList = ['Partner','Mol System','Chains','Models','Force Field',
                       'Auto Set His\nProtonation State','Is nucleic acid?']
        editWidgets      = [None, self.molSysPulldown,  self.chainSelect,
                            self.modelPulldown, None,
                            None, None]    
        editGetCallbacks = [None, self.getMolSystem, self.getChains,
                            self.getModel, None,
                            self.toggleAutoHis, self.toggleIsDna]
        editSetCallbacks = [None, self.setMolSystem, self.setChains,
                            self.setModel, None,
                            None, None]
        self.molPartnerMatrix = ScrolledMatrix(molFrame, headingList=headingList,
                                            editSetCallbacks=editSetCallbacks,
                                            editGetCallbacks=editGetCallbacks, 
                                            editWidgets=editWidgets,
                                            multiSelect=False, initialRows=10,
                                            callback=self.selectMolPartner)
        self.molPartnerMatrix.grid(row=0, column=0, columnspan=1, sticky='nsew')

        texts = ['Add Mol Partner','Remove Mol Partner','Load From File']
        commands = [self.addMolPartner,self.removeMolPartner,self.loadPdbFromFile]

        self.molPartnerButtons = ButtonList(molFrame, texts=texts, expands=True,  commands=commands)
        self.molPartnerButtons.grid(row=1, column=0, sticky='ew')

        # RESIDUE FRAME
        resFrame.grid_columnconfigure(1, weight=1)
        resFrame.grid_rowconfigure(1, weight=1)

        label1 = Label(resFrame, text='Molecular Partner: ')
        label1.grid(row=0, column=0, sticky='nw')
        self.molPartnerPulldown = PulldownMenu(resFrame, callback=self.chooseMolPartner)
        self.molPartnerPulldown.grid(row=0, column=1, sticky='w') 

        label2 = Label(resFrame, text='Semi-flexible Mode: ')
        label2.grid(row=0, column=2, sticky='ne')
        self.setSemiflexMode = PulldownMenu(resFrame, callback=self.chooseSemiFlexMode)
        self.setSemiflexMode.grid(row=0, column=3, sticky='e')
        
        label3 = Label(resFrame, text='         AIR upper distance limit: ')
        label3.grid(row=0, column=4, sticky='ne')
        self.airUpperDistanceLimit = FloatEntry(resFrame, returnCallback=self.updateAirUpperDistanceLimit, width=6)
        self.airUpperDistanceLimit.bind('<Leave>', self.updateAirUpperDistanceLimit)
        self.airUpperDistanceLimit.grid(row=0, column=5, sticky='e')

        headingList = ['CCPN Chain','Residue','Active?','Passive?','Flexibility']
        editWidgets      = [None,None,None,None,None]
        editGetCallbacks = [None,None,self.toggleActiveResidue,self.togglePassiveResidue,None]
        editSetCallbacks = [None,None,None,None,None]
        self.residueMatrix = ScrolledMatrix(resFrame,
                                            headingList=headingList,
                                            multiSelect=True,
                                            editWidgets=editWidgets,
                                            editGetCallbacks=editGetCallbacks,
                                            editSetCallbacks=editSetCallbacks,
                                            callback=self.selectResidue)
        self.residueMatrix.grid(row=1, column=0, columnspan=6, sticky = 'nsew')

        texts = ['Set Active', 'Set Passive',
                 'Set Inactive', 'Set Fully Flexible',
                 'Set Semi-flexible', 'Set Inflexible']
        commands = [self.setResiduesActive, self.setResiduesPassive,
                    self.setResiduesInactive, self.setResiduesFullFlexible,
                    self.setResiduesSemiFlexible, self.setResiduesInflexible]
           
        self.residueButtons = ButtonList(resFrame, texts=texts, expands=True,  commands=commands)
        self.residueButtons.grid(row=2, column=0, columnspan=6, sticky='ew')

        # RUN FRAME
        runFrame.grid_columnconfigure(2, weight=1)
        runFrame.grid_rowconfigure(1, weight=1)

        label = Label(runFrame, text='Current Run: ')
        label.grid(row=0, column=0, sticky='w')
        self.runPulldown = PulldownMenu(runFrame, callback=self.chooseRun)
        self.runPulldown.grid(row=0, column=1, sticky='w') 

        symmButton = Button(runFrame,text='Edit DANI parameters',command=self.editDaniSets)
        symmButton.grid(row=0, column=2, sticky='e')

        symmButton = Button(runFrame,text='Edit RDC parameters',command=self.editRdcSets)
        symmButton.grid(row=0, column=3, sticky='e')

        symmButton = Button(runFrame,text='Edit Symmetries',command=self.editSymmetrySets)
        symmButton.grid(row=0, column=4, sticky='e') 

        headingList      = ['Parameter','Value','Description']
        justifyList      = ['center','center', 'left']
        editWidgets      = [None, True, None]    
        editGetCallbacks = [None, self.getValue, None]
        editSetCallbacks = [None, self.setValue, None]
        self.runMatrix = ScrolledMatrix(runFrame, headingList=headingList,
                                        justifyList=justifyList,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks, 
                                        editWidgets=editWidgets,
                                        multiSelect=False, initialRows=10,
                                        passSelfToCallback=True,
                                        callback=self.selectRowObj)

        self.runMatrix.grid(row=1, column=0, columnspan=5, sticky='nsew')
        self.runMatrix.refreshFunc = self.updateRuns
        self.runMatrix.doEditMarkExtraRules = self.runMatrixEditMark

        texts = ['Add New Default Run','Remove Current Run','Copy Current Run']
        commands = [self.addNewRun,self.removeRun,self.copyRun]

        self.runButtons = ButtonList(runFrame, texts=texts, expands=True,  commands=commands)
        self.runButtons.grid(row=2, column=0, columnspan=5, sticky='ew')

        # RESTRAINTS FRAME
        conFrame.grid_columnconfigure(5, weight=1)
        conFrame.grid_rowconfigure(1, weight=1)

        label = Label(conFrame, text='CCPN Constraint Set: ')
        label.grid(row=0, column=0, sticky='nw')
        self.constraintSetPulldown = PulldownMenu(conFrame, callback=self.chooseConstraintSet)
        self.constraintSetPulldown.grid(row=0, column=1, sticky='w') 

        label = Label(conFrame, text='Current Run: ')
        label.grid(row=0, column=2, sticky='w')

        self.runPulldown2 = PulldownMenu(conFrame, callback=self.chooseRun)
        self.runPulldown2.grid(row=0, column=3, sticky='w') 

        self.restrListPulldown   = PulldownMenu(self, callback=self.setRestrList, do_initial_callback=False)

        headingList = ['#','Restraint Type','Format','CCPN List','CNS File']
        editWidgets      = [None, None, None, self.restrListPulldown, None]
        editGetCallbacks = [None, None, None, self.getRestrList, self.getCnsFile]
        editSetCallbacks = [None, None, None, self.setRestrList, None]
                   
        self.restraintMatrix = ScrolledMatrix(conFrame,
                                            headingList=headingList,
                                            editWidgets=editWidgets,
                                            editGetCallbacks=editGetCallbacks,
                                            editSetCallbacks=editSetCallbacks,
                                            callback=self.selectEnergyTerm)
        self.restraintMatrix.grid(row=1, column=0, columnspan=6, sticky = 'nsew')

        texts = ['Remove Restraints','Add Restraints']
        commands = [self.removeEnergyTerm,self.addEnergyTerm]

        self.termButtons = ButtonList(conFrame, texts=texts, expands=True,  commands=commands)
        self.termButtons.grid(row=2, column=0, columnspan=4, sticky='ew')

        label = Label(conFrame, text='Restraint Type:')
        label.grid(row=2, column=4, sticky='w')

        entries = energyTermTypes.keys()
        entries.sort()
        self.termPullldown = PulldownMenu(conFrame, callback=None,
                                          entries=entries,
                                          selected_index=5,
                                          do_initial_callback=False)
        self.termPullldown.grid(row=2, column=5, sticky='w')

        self.initHaddockProject()
        self.updateConstraintSets()
        self.updateHaddockProjects()
        self.updateMolPartners()
        self.updateResRestraints()
        self.updateRestraintTerms()
        self.updateRuns()

        self.notify(self.parent.parent.registerNotify)

    def initHaddockProject(self):
            
        # Load current Haddock projects if none initialize default project with one run

        ccpnProject = self.ccpnProject
        if ccpnProject:
          if len(ccpnProject.sortedHaddockProjects()):
              self.hProject = ccpnProject.sortedHaddockProjects()[-1]
          else:
              self.hProject = ccpnProject.findFirstHaddockProject(name='Default')
          if not self.hProject:
              self.hProject = ccpnProject.newHaddockProject(name='Default',workingDir='.')
              self.hProject.newRun()
              
          self.run = self.hProject.findFirstRun()      

    def destroy(self):
  
        self.notify(self.parent.parent.unregisterNotify)
        Frame.destroy(self)

    def notify(self, notifyFunc):
        # respond to external data model events

        for func in ('__init__', 'delete', 'setSymmetryCode','setSegmentLength'):
            notifyFunc(self.updateAllAfter, 'molsim.Symmetry.Symmetry', func)

    def runMatrixEditMark(self, obj, row, col):
  
        if self.runMatrix.objectList[row]: return True
        else: return False  

    def expClassic(self):
    
        """Description:    Exports the project as a classic HADDOCK project directory. 
                           Calls the exportClassic class
           Input:        Haddock project instance, CCPN project instance
           Output:        Classic project on disk
        """
        classic = exportClassic(hProject=self.hProject,
                                latestRun=self.run,
                                ccpnProject=self.ccpnProject)    
    
    def expParam(self):

        """Exports the project as a new selfcontaining HADDOCK parameter file. Calls the exportParam class"""

        param = exportParam(hProject=self.hProject,
                            latestRun=self.run,
                            ccpnProject=self.ccpnProject)
        param.writeToFile()

    def uploadToServer(self):
        
        """Upload the HADDOCK project to the server. Requires a valid password and username"""
        
        if self.haddockUpload: self.haddockUpload.open()
        else: self.haddockUpload = HaddockServerUpload(self.parent,
                                                       hProject=self.hProject,
                                                       latestRun=self.run,
                                                       ccpnProject=self.ccpnProject)

    def impRunCns(self):    

        """Import the run specific settings from a Haddock run.cns file"""
        
        filetype = [ FileType("CNS file", ["*.cns"]) ]
        popup = FileSelectPopup(self, filetype)
        runcns = popup.getFile()
        popup.destroy()
        
        runCnsImporter(self.hProject,runcns)
        self.updateAllAfter()
        
    def selectEnergyTerm(self, obj, row, col):
        
        self.energyTerm = obj
  
    def editSymmetrySets(self):
  
        if self.symmetryPopup: self.symmetryPopup.open()
        else: self.symmetryPopup = EditSymmetryPopup(self.parent, self.ccpnProject)

    def editDaniSets(self):
        
        """Edit relaxation data energy terms.Popup can only be initiated if DANI restraint
           file is defined
        """    
        constraints = [c.code for c in self.run.sortedHaddockEnergyTerms()]
        
        if 'DANI' in constraints:
            self.haddockDaniPopup = HaddockDaniPopup(self.parent, self.run)
        else:
            showWarning('Warning', 'No DANI restraints present', parent=self)
            return

    def editRdcSets(self):

        """Edit Residual Dipolar Coupling related energy terms. Popup can only be initiated if RDC restraint
           file is defined
        """
        constraints = [c.code for c in self.run.sortedHaddockEnergyTerms()]
        
        if 'RDC' in constraints:
            self.haddockRdcPopup = HaddockRdcPopup(self.parent, self.run)
        else:
            showWarning('Warning', 'No RDC restraints present', parent=self)
            return
    
    def chooseConstraintSet(self, index, name):
  
        constraintSets = self.ccpnProject.sortedNmrConstraintStores()
        constraintSet = constraintSets[index]

        if constraintSet is not self.constraintSet:
            self.constraintSet = constraintSet
            self.updateRestraintTerms()  

    def updateConstraintSets(self):
  
        names = []
        index = -1

        if self.ccpnProject:
            constraintSets = self.ccpnProject.sortedNmrConstraintStores()

            if constraintSets:
                if self.constraintSet not in constraintSets: self.constraintSet = constraintSets[0]

                names = ['%s' % cs.serial for cs in constraintSets]
                index = constraintSets.index(self.constraintSet)

            else: self.constraintSet = None

        self.constraintSetPulldown.setup(names, index)
    
    def updateAirUpperDistanceLimit(self, event):
        
        if self.molPartner:
            self.molPartner.airUpperDistanceLimit = self.airUpperDistanceLimit.get()
        
    def chooseSemiFlexMode(self, index, name):
        
        """If their is a haddock partner then set it to the choosen semi-flexibility mode"""
        
        if self.molPartner:
            self.molPartner.semiFlexMode = name
    
    def updateSemiFlexModePulldown(self):    
        
        """Set the semi-flexibility mode pulldown menu for the Haddock partners. If no partner then their is
           nothing to choose. By default the mode is set to 'automatic'
        """
        names = []
        index = -1
        
        if self.hProject:
            modes = self.hProject.metaclass.container.getElement('SemiFlexMode').enumeration
    
            if modes and self.molPartner:
                names = modes
                index = modes.index(self.molPartner.semiFlexMode)
            else: self.molPartner = None
        
        self.setSemiflexMode.setup(names, index)
        
    def chooseMolPartner(self, index, name):

        """Return a list of all stored HADDOCK projects stored in the CCPN project.
           Called by self.molPartnerPulldown"""

        partners = self.hProject.sortedHaddockPartners()
        self.molPartner = partners[index]
        self.updateResRestraints()
    
    def updateMolPartnerPulldown(self):
  
        names = []
        index = -1

        if self.hProject:
            partners = self.hProject.sortedHaddockPartners()

            if partners:
                if self.molPartner not in partners: self.molPartner = partners[0]

                names = ['Partner %s' % p.code for p in partners]
                index = partners.index(self.molPartner)

            else: self.molPartner = None

        self.molPartnerPulldown.setup(names, index)

    def setResiduesInactive(self):
    
        """Set ambiquaty state of selected residue to inactive"""

        if self.residue:
            for residue in self.residueMatrix.currentObjects: residue.interaction = 'none'

        self.updateAllAfter()
        
    def setResiduesActive(self):

        """Set ambiquaty state of selected residue to active"""
        
        if self.residue:
            for residue in self.residueMatrix.currentObjects: residue.interaction = 'active'

        self.updateAllAfter()
        
    def setResiduesPassive(self):

        """Set ambiquaty state of selected residue to passive"""

        if self.residue:
            for residue in self.residueMatrix.currentObjects: residue.interaction = 'passive'       

        self.updateAllAfter()
        
    def setResiduesInflexible(self):
    
        """Set flexibility state of selected residue to inflexible"""

        if self.residue:
            for residue in self.residueMatrix.currentObjects: residue.flexibility = 'none'

        self.updateAllAfter()
    
    def setResiduesSemiFlexible(self):

        """Set flexibility state of selected residue to semi-flexible"""

        if self.residue:
            for residue in self.residueMatrix.currentObjects: residue.flexibility = 'semi'
    
        self.updateAllAfter()
        
    def setResiduesFullFlexible(self):

        """Set flexibility state of selected residue to fully-flexible"""

        if self.residue:
            for residue in self.residueMatrix.currentObjects: residue.flexibility = 'full'
     
        self.updateAllAfter()
        
    def toggleActiveResidue(self, *opt):

        """Toggle the residues of the selected molecular system to Active"""

        if self.residue: interaction = self.residue.interaction
    
        if interaction == 'active': self.residue.interaction = 'none'
        else: self.residue.interaction = 'active'

        self.updateAllAfter()
        
    def togglePassiveResidue(self, *opt):

        """Toggle the residues of the selected molecular system to Active"""

        if self.residue: interaction = self.residue.interaction
    
        if interaction == 'passive': self.residue.interaction = 'none'
        else: self.residue.interaction = 'passive'

        self.updateAllAfter()
        
    def selectResidue(self, object, row, col):
        
        self.residue = object
    
    def selectMolPartner(self, object, row, col):
  
        self.molPartner = object
        self.updateMolPartnerPulldown()
        self.updateSemiFlexModePulldown()
    
    def loadPdbFromFile(self):
        
        """Add structure from PDB file to new MolSystem with name of the pdb"""
        
        filetype = [ FileType("PDB structure file", ["*.pdb"]) ]
        popup = FileSelectPopup(self, filetype)
        pdb = popup.getFile()
        popup.destroy()
        
        if not pdb: return
                
        print("Load PDB structure from file: %s" % pdb)
        
        newMolsysName = (os.path.basename(pdb)).split('.')[0]
        
        molSystem = self.ccpnProject.findFirstMolSystem(code=newMolsysName)
        
        if molSystem:
            msg = 'Use existing molecular system "%s" for structure import?' % (newMolsysName)
          
            if not showYesNo('Query', msg):
                i = 1
                name = newMolsysName
                while self.ccpnProject.findFirstMolSystem(code=name):
                    name = '%s_%d' % (newMolsysName, i)
                    i += 1
            
                molSystem = self.ccpnProject.newMolSystem(code=name,name=name) 
        
        else:
            molSystem = self.ccpnProject.newMolSystem(code=newMolsysName,name=newMolsysName)
        
        structure = getStructureFromFile(molSystem, pdb, fileType='rough', doWarnings=True)
        self.addMolPartner()
        
    def chooseHaddockProject(self, index, name):
    
        """Is called from the 'General' tab, HADDOCK Project pulldown menu. If <new> is choosen then
           create a new HADDOCK project and update pulldown menu. 'updateHaddockProjects' refreshes
           complete list of all projects inclusing <new> and Default."""

        if name == '<New>':
            name   = ''; prompt = 'Enter name for HADDOCK project:'
            while not name:
                name = askString(title='Name Query', prompt=prompt, parent=self)
        
                if name is None: break
                elif len(name.split()) != 1:
                    name = ''
                    prompt = 'Name invalid.\nEnter project name\n(without spaces):'
                elif name:
                    self.hProject = self.ccpnProject.newHaddockProject(name=name,workingDir='.')
                    self.hProject.newRun()
                    self.updateHaddockProjects()
   
        else: 
            self.hProject = self.ccpnProject.findFirstHaddockProject(name=name)
     
        self.updateAllAfter()
        
    def updateHaddockProjects(self):
  
        """Updates the full list of available HADDOCK projects including the option <new> and Default."""

        index = -1
        names = []
        
        if self.ccpnProject:
            hProjects = self.ccpnProject.sortedHaddockProjects()
 
            if hProjects:
                names = [hProject.name for hProject in hProjects]
                if self.hProject not in hProjects: self.hProject = hProjects[0]
                index = hProjects.index(self.hProject)
            else:
                self.hProject = self.ccpnProject.newHaddockProject(name='Default',workingDir='.')
                names.append(self.hProject.name)
                index = 0
 
            names.append('<New>')
        
        self.hProjectPulldown.setup(names,index)

  #
  # Table update functions
  #

    def updateAllAfter(self, obj=None):
        
        """Function to funnel updates to one point bounce out if already updating to avoid unnecessary calls"""

        if self.waiting: return
        else:
          self.waiting = True
          self.after_idle(self.updateAll)
        
    def updateAll(self, project=None):
    
        if project and project is not self.ccpnProject:
            self.ccpnProject = project
            self.molPartner    = None
            self.residue  = None
            self.run = None
            self.energyTerm    = None
            self.constraintSet = None
            self.symmetrySet   = None
            self.initHaddockProject()
        
        self.updateHaddockProjects()
        self.updateConstraintSets()
        self.updateGeneral()
        self.updateMolPartners()
        self.updateResRestraints()
        self.updateRestraintTerms()
        self.updateRuns()
        self.waiting = False
        
    def updateGeneral(self):
    
        """Update the 'General' scrolled matrix. Set working directory calles 'getDirectory' to set directory. Default
           is current. Molecular partners displays currently set HADDOCK partners. Editing these partners calles the
           'Molecular Partners' tab. By default always show two molecular partners. If partners are allready in project
           or if they are set then display the partner, the chain, the number of active and passive restraints and the
           number of (semi)-flexible segments. Ambiques restraints and flexibility can be set by calling 'editResRestraints'
           function (Restraints tab). Runs displays the runs in the current HADDOCK project. A project will always be set
           with one default run (number 1). Edit runs calls 'editRuns' function (Runs tab).
        """
        table = self.generalMatrix
        
        if not self.hProject:
          table.update(objectList=[], textMatrix=[], colorMatrix=[])
          return
          
        blank = [None] * 3
 
        textMatrix = []; objectList = []; colorMatrix = []

        # Set working directory
        colorMatrix.append([headingColor]*3)
        textMatrix.append(['Project General',None,None])
        objectList.append(None)

        colorMatrix.append(blank)
        textMatrix.append(['Working Directory', self.hProject.workingDir, 'Location to store PDB files, restraint files and "new.html" file.'])
        objectList.append([self.hProject,'workingDir', None, self.getDirectory, None])

        # Set HADDOCK partners. 
        colorMatrix.append([headingColor]*3)
        textMatrix.append(['Molecules',None,None])
        objectList.append(None)
     
        if self.hProject.haddockPartners:
            for partner in self.hProject.sortedHaddockPartners():
                pCode = partner.code
                label = 'Mol Partner %s' % pCode
                molText = '%s: %s' % (partner.molSystem.code,','.join([c.chain.code for c in partner.chains]))

                colorMatrix.append(blank)
                textMatrix.append([label, molText, 'CCPN Mol System chains that make up haddock partner %s' % pCode])
                objectList.append([None,'', None, self.editMolPartners, None])

                active = 0; passive = 0; flexible = 0; semiflexible = 0

                for chain in partner.chains:
                    for residue in chain.residues:
                        interaction = residue.interaction
                        flexibility = residue.flexibility

                        if interaction == 'active': active  += 1
                        elif interaction == 'passive': passive += 1
                
                        if flexibility == 'full': flexible += 1 
                        elif flexibility == 'semi': semiflexible += 1 

                colorMatrix.append(blank)
                textMatrix.append(['Active Residues %s' % pCode, active,
                                   'Number of residues in partner %s making active ambiguous interactions' % pCode])
                objectList.append([partner,'', None, self.editResRestraints, None])

                colorMatrix.append(blank)
                textMatrix.append(['Passive Residues %s' % pCode, passive,
                                   'Number of residues in partner %s making passive ambiguous interactions ' % pCode])
                objectList.append([partner,'', None, self.editResRestraints, None])

                colorMatrix.append(blank)
                textMatrix.append(['Fully Flexible Residues %s' % pCode, flexible,
                                   'Number of fully flexible residues in partner %s' % pCode])
                objectList.append([partner,'', None, self.editResRestraints, None])

                colorMatrix.append(blank)
                textMatrix.append(['Semi-flexible Residues %s' % pCode, semiflexible,
                                   'Number of semi-flexible residues in partner %s' % pCode])
                objectList.append([partner,'', None, self.editResRestraints, None])

        else: # Show some empty fields
            for pCode in ('A','B'):
                label = 'Mol Partner %s' % pCode
        
                colorMatrix.append(blank)
                textMatrix.append([label, 'edit Partner', 'CCPN Mol System chains that make up haddock partner %s' % pCode]) 
                objectList.append([None,'', None, self.editMolPartners, None])

        # Set runs
        colorMatrix.append([headingColor]*3)
        textMatrix.append(['Runs',None,None])
        objectList.append(None)
 
        if self.hProject.runs:
            for run in self.hProject.sortedRuns():
                colorMatrix.append(blank)
                textMatrix.append(['Run %d' % run.serial, 'edit Run','HADDOCK docking parameter settings'])
                objectList.append([run,'', None, self.editRuns, None])

                i = 0
                for term in run.sortedHaddockEnergyTerms():  
                    if term.code in energyTermTypes:
                        colorMatrix.append(blank)
                        textMatrix.append(['Restraints %d:%d' % (run.serial,i), '%s:%d' % (term.code,term.termId),
                                           'Restraint term of type %s in run %s' % (term.code,run.serial)])
                        objectList.append([run,'', None, self.editRestraintTerms, None])
                        i += 1
        else:
            colorMatrix.append(blank)
            self.addNewRun()
            
            for run in self.hProject.sortedRuns():
                colorMatrix.append(blank)
                textMatrix.append(['Run %d' % run.serial, 'edit Run','HADDOCK docking parameter settings'])
                objectList.append([run,'', None, self.editRuns, None])

        # Finally update the ScrolledMatrix table
        table.update(objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)
        
    def updateMolPartners(self):
  
        """By default the 'useDnaRestraints' flag is set to True with the first occurance of a DNA partner"""

        textMatrix = []
        objectList = []
        
        if self.hProject:
            run = self.hProject.sortedRuns()[-1]
            
            for partner in self.hProject.sortedHaddockPartners():
                if partner.structureEnsemble.getDetails(): 
                    modelId = int(partner.structureEnsemble.getDetails())
                else: 
                    ensembles = self.ccpnProject.sortedStructureEnsembles()
                    if partner.molSystem: 
                        ensembles = [e for e in ensembles if e.molSystem is partner.molSystem]
                        models = []
                        for ensemble in ensembles:    
                            models += [model for model in ensemble.sortedModels()]    
                        if len(models): 
                            partner.structureEnsemble.setDetails(str(len(models)))
                            modelId = len(models)    
        
                datum = [partner.code,
                         partner.molSystem.code,
                         ','.join([c.chain.code for c in partner.chains]),
                         modelId,
                         partner.forceFieldCode,
                         partner.autoHistidinePstate and 'Yes' or 'No',
                         partner.isDna and 'Yes' or 'No']
        
                objectList.append(partner)
                textMatrix.append(datum)

                isDNA = False
                for partner in self.hProject.sortedHaddockPartners():
                    if partner.get('isDna') == True: isDNA = True
                run.useDnaRestraints = isDNA
    
        self.molPartnerMatrix.update(objectList=objectList, textMatrix=textMatrix)
        self.updateMolPartnerPulldown()
        self.updateSemiFlexModePulldown()
        
    def updateResRestraints(self):
  
        textMatrix = []
        objectList = []
        colorMatrix = []

        if self.molPartner: 
            
            # Update the floatEntry dialog with the new air upper distance limit
            self.airUpperDistanceLimit.set(self.molPartner.airUpperDistanceLimit) 
            
            hChains = self.molPartner.sortedChains()

            for hChain in hChains:
                hResidues = [(r.residue.seqCode,r) for r in hChain.residues]
                hResidues.sort()           
         
                for seqCode, hResidue in hResidues:

                    residue     = hResidue.residue
                    flexibility = hResidue.flexibility
                    interaction = hResidue.interaction

                    colours = [None, None, None, None, None]         
                    isActive  = 'No'
                    isPassive = 'No'

                    ccpCode = residue.ccpCode

                    if interaction == 'active':
                        isActive = 'Yes'
                        colours[2] = activeColour
                    elif interaction == 'passive':
                        isPassive= 'Yes'
                        colours[3] = passiveColour

                    if flexibility == 'full': colours[4] = flexibleColor
                    elif flexibility == 'semi': colours[4] = semiflexColor

                    chain = hChain.chain

                    datum = ['%s:%s' % (chain.molSystem.code,chain.code),
                             '%d%s' % (seqCode,ccpCode),
                             isActive, isPassive, flexibility]

                    objectList.append(hResidue)
                    textMatrix.append(datum)               
                    colorMatrix.append(colours)

        self.residueMatrix.update(colorMatrix=colorMatrix,
                                  objectList=objectList,
                                  textMatrix=textMatrix) 
        self.updateSemiFlexModePulldown()     #Update the semiflexibility mode ('manual' or 'automatic')

    def updateRestraintTerms(self):

        if self.run and self.constraintSet:
            if self.run.nmrConstraintStore is not self.constraintSet:
                setRunConstraintSet(self.run, self.constraintSet)

        textMatrix = []
        objectList = []

        if self.run:
            i = 0
            for term in self.run.sortedHaddockEnergyTerms():
                i += 1
                if term.constraintList:
                    fileName = None
                    format  = 'CCPN'
                    constraintList = '%s:(%s)' % (term.constraintList.serial,term.constraintList.className[:-14])
                else:
                    constraintList = None
                    format = 'CNS'
                    fileName = term.fileName

                if term.code in energyTermTypes:    
                    datum = [i,
                             '%s;%d' % (term.code,term.termId),
                             format,
                             constraintList,
                             fileName]

                    objectList.append(term)
                    textMatrix.append(datum)               

        self.restraintMatrix.update(objectList=objectList,textMatrix=textMatrix) 

    def addEnergyTerm(self):
  
        """Add new Haddock energy terms. If the term is a RDC energy term than only 5 entries are allowed.
           For every entry add a corresponding energy protocol.
        """
        if self.run:
            termType = self.termPullldown.getSelected()

            if termType == 'RDC' or termType == 'DANI':
                termlist = [ i.termId for i in self.run.sortedHaddockEnergyTerms() if i.code == termType ]
                termlist.sort()
                if len(termlist) < 5:
                    if len(termlist): termId = termlist[-1] + 1
                    else: termId = 1
                    energyTerm = self.run.newHaddockEnergyTerm(code=termType,termId=termId)
                    if termType == 'RDC': addRdcParam(self.run,termId)
                    if termType == 'DANI': addDaniParam(self.run,termId)
                    self.updateAllAfter()
                else: 
                    showWarning('Warning',('Only 5 %s parameter sets allowed' % termType), parent=self)
                    return
            else:
                termId = 1
                while self.run.findFirstHaddockEnergyTerm(code=termType,termId=termId): termId +=1

                energyTerm = self.run.newHaddockEnergyTerm(code=termType,termId=termId)

                self.updateAllAfter()                
                
        else: showWarning('Warning','No active run. Set this first.',parent=self)

    def removeEnergyTerm(self):
  
        """Remove selected Haddock energy terms. If the energy term is an RDC term the corresponding energy protocol
           is also removed.
        """
        if self.energyTerm:
            if self.energyTerm.code == 'RDC':
                rdcprotocol = self.run.findFirstHaddockEnergyTerm(code='rdcProtocolStore',termId=self.energyTerm.termId)
                if rdcprotocol: rdcprotocol.delete()
            elif self.energyTerm.code == 'DANI':
                daniprotocol = self.run.findFirstHaddockEnergyTerm(code='daniProtocolStore',termId=self.energyTerm.termId)
                if daniprotocol: daniprotocol.delete()
            self.energyTerm.delete()
            self.energyTerm = None
            
            self.updateAllAfter()
            
    def chooseRun(self, index, name):
  
        runs = self.hProject.sortedRuns()

        if runs: 
            run = runs[index]
            if run is not self.run:
                self.run = run
                self.updateRuns()

    def updateRunPulldown(self):
  
        names = []
        index = -1

        if self.hProject:
            runs = self.hProject.sortedRuns()

            if runs:
                if self.run not in runs: self.run = runs[-1]

                names = ['Run %d' % r.serial for r in runs]
                index = runs.index(self.run) 

            else: self.run = None

        self.runPulldown.setup(names, index)
        self.runPulldown2.setup(names, index)

    def updateRuns(self):

        table = self.runMatrix
        if not self.hProject:
                  table.update(objectList=[], textMatrix=[], colorMatrix=[])
                  return

        textMatrix  = []; objectList  = []; colorMatrix = []

        run = self.run
        if run:
            
            #Restraint definitions
            textMatrix.append(['Restraint definition', None, None])
            objectList.append(None)

            textMatrix.append(['Use Centre of Mass Restraints', run.centerOfMassRestraints, 'Whether to define center of mass restraints to\nenforce contact between the molecules'])
            objectList.append([run,'centerOfMassRestraints', None, self.toggleBoolean, None])

            textMatrix.append(['Centre of Mass Constant', run.centerOfMassConstant, 'Force constant for center of mass restraints'])
            objectList.append([run,'centerOfMassConstant', self.floatEntry, self.getFloat, self.setNonNegFloat])

            textMatrix.append(['Use Surface Contact?', run.surfaceContactRestraints, 'Whether to define surface contact restraints to\nenforce contact between the molecules'])
            objectList.append([run,'surfaceContactRestraints', None, self.toggleBoolean, None])

            textMatrix.append(['Surface Contact Constant', run.surfaceContactConstant, 'Force constant for surface contact restraints'])
            objectList.append([run,'surfaceContactConstant', self.floatEntry, self.getFloat, self.setNonNegFloat])

            textMatrix.append(['Random Ambig Restraints?', run.randomAmbigRestraints, 'Whether to define randomly ambiguous interaction\nrestraints from accessible residues'])
            objectList.append([run,'randomAmbigRestraints', None, self.toggleBoolean, None])      

            textMatrix.append(['Random Exlude AIR', run.randomExcludeAir, 'Whether to randomly exclude a fraction of the\nambiguous restraints (AIRs)'])
            objectList.append([run,'randomExcludeAir', None, self.toggleBoolean, None])

            textMatrix.append(['Num Random Exclude Parts', run.randomExclParts, 'Number of partitons for random exclusion'])
            objectList.append([run,'randomExclParts', self.intEntry, self.getInt, self.setNonNegInt])

            textMatrix.append(['Use H-bond Restraints?', run.useHBondRestraints, 'Whether to use the hydrogen bond restraints'])
            objectList.append([run,'useHBondRestraints', None, self.toggleBoolean, None])

            textMatrix.append(['Use DNA Restraints?', run.useDnaRestraints, 'Whether to use DNA/RNA restraints'])
            objectList.append([run,'useDnaRestraints', None, self.toggleBoolean, None])

            #Do or do not use defined symmetry operators
            symRows = 1
            textMatrix.append(['Symmetry Operations', None, None]) 
            objectList.append(None)

            symmetryOps = run.sortedSymmetryRestraints()
            symmetrySets = self.ccpnProject.sortedMolSystemSymmetrySets()

            if symmetrySets:
                for symmetrySet in symmetrySets:
                     for symmOp in symmetrySet.sortedSymmetries():
                        symRows += 1

                        chainRanges = []
                        for segment in symmOp.segments:
                            seqId = segment.firstSeqId
                            chainRanges.append('%s:%d-%d' % (segment.chainCode,seqId,seqId+symmOp.segmentLength))

                        symmText = ','.join(chainRanges)
                        description = 'Whether to use symmetry operation of type "%s"\nover range%s' % (symmOp.symmetryCode,symmText)

                        if symmOp in symmetryOps: inUse = True
                        else: inUse = False

                        textMatrix.append(['%s:%s:%s' % (symmetrySet.molSystem.code,symmOp.symmetryCode,symmText), inUse, description])
                        objectList.append([symmOp,None,None,self.useSymmetry,None])

            else:
                symRows += 1
                textMatrix.append(['None present', None, None])
                objectList.append(None)

            #Docking protocol
            textMatrix.append(['Docking protocol', None, None]) 
            objectList.append(None)
    
            textMatrix.append(['Num It 0 Structures', run.numIt0Structures, 'Number of structures for rigid body docking'])
            objectList.append([run,'numIt0Structures', self.intEntry, self.getInt, self.setNonNegInt])

            textMatrix.append(['Num It 1 Structures', run.numIt1Structures, 'Number of structures for semi-flexible SA'])
            objectList.append([run,'numIt1Structures', self.intEntry, self.getInt, self.setNonNegInt])

            textMatrix.append(['Num Water Structures', run.numWrefStructures, 'Number of structures for water refinement'])
            objectList.append([run,'numWrefStructures', self.intEntry, self.getInt, self.setNonNegInt])    

            textMatrix.append(['Solvent',run.solvent,'Which solvent to use for final explicit\nsolvent refinement'])
            objectList.append([run,'solvent',self.pulldownMenu,self.getSolvent,self.setPulldownMenu])

            textMatrix.append(['Calculate Desolvation', run.calcDesolvation, 'Whether to calculate explicit desolvation energy '])
            objectList.append([run,'calcDesolvation', None, self.toggleBoolean, None])

            textMatrix.append(['Random Start Orient?', run.radomizeStartOriention, 'Whether to randomise the starting orientation'])
            objectList.append([run,'radomizeStartOriention', None, self.toggleBoolean, None])

            textMatrix.append(['Do Initial Rigid?', run.initialRigidBodyMinim, 'Whether to do the initial rigid body minimisation'])
            objectList.append([run,'initialRigidBodyMinim', None, self.toggleBoolean, None])

            textMatrix.append(['Use Rigid Translations?', run.doRigidTranslations, 'Whether to allow translation in rigid body minimisation'])
            objectList.append([run,'doRigidTranslations', None, self.toggleBoolean, None])

            textMatrix.append(['Num trails', run.nTrails, 'Number of trials for rigid body minimisation'])
            objectList.append([run,'nTrails', self.intEntry, self.getInt, self.setNonNegInt])    

            textMatrix.append(['Random Seed', run.randomSeed, 'The seed number for the random number generator'])
            objectList.append([run,'randomSeed', self.intEntry, self.getInt, self.setInt])

            textMatrix.append(['Rotate It0 180', run.rotate180It0, 'Whether to sample 180 degrees rotated solutions\nduring rigid body EM'])
            objectList.append([run,'rotate180It0', None, self.toggleBoolean, None])

            textMatrix.append(['Rotate It1 180', run.rotate180It1, 'Whether to sample 180 degrees rotated solutions\nduring semi-flexible SA'])
            objectList.append([run,'rotate180It1', None, self.toggleBoolean, None])

            textMatrix.append(['Remove Non-polar H?', run.removeNonPolarH, 'Whether to remove non-polar hydrogens\n(speeds up the calculation)'])
            objectList.append([run,'removeNonPolarH', None, self.toggleBoolean, None])

            textMatrix.append(['Skip Structures', run.skipStructures, 'Number of structures to skip in the selection\nof structures in it0'])
            objectList.append([run,'skipStructures', self.intEntry, self.getInt, self.setNonNegInt])

            #Energy and interaction parameters
            textMatrix.append(['Energy and interaction parameters', None, None])
            objectList.append(None)

            textMatrix.append(['Use Dihedral Energy?', run.doIncludeDihEnergy, 'Whether to include dihedral angle energy terms'])
            objectList.append([run,'doIncludeDihEnergy', None, self.toggleBoolean, None])

            textMatrix.append(['Use Rigid Electrostatics?', run.doRigidBodyElectrostatics, 'Whether to include electrostatic during\nrigid body docking (it0)'])
            objectList.append([run,'doRigidBodyElectrostatics', None, self.toggleBoolean, None])

            textMatrix.append(['Use SA Electrostatics?', run.doSAElectrostatics, 'Whether to include electrostatic during\nsemi-flexible SA (it1)'])
            objectList.append([run,'doSAElectrostatics', None, self.toggleBoolean, None])

            textMatrix.append(['Epsilon','%1.1f' % run.epsilon, 'The epsilon constant for the electrostatic energy term'])
            objectList.append([run,'epsilon', self.floatEntry, self.getFloat, self.setFloat])

            textMatrix.append(['Dielectric Type',run.dielectricType,'Whether to use constant (cdie) or distance-dependent\n(rdie) dielectric'])
            objectList.append([run,'dielectricType',self.pulldownMenu,self.getDielectricType,self.setPulldownMenu])

            textMatrix.append(['Non-bonded Type',run.nonBondedType,'specify the type of non-bonded interaction'])
            objectList.append([run,'nonBondedType',self.pulldownMenu,self.getNonBondedType,self.setPulldownMenu])

            textMatrix.append(['Rigid IM Scale','%1.1f' % run.rigidbodyIMinteractScaling, 'Scaling factor of intermolecular interactions\nfor rigid body EM'])
            objectList.append([run,'rigidbodyIMinteractScaling', self.floatEntry, self.getFloat, self.setFloat]) 

            #Solvated Docking
            textMatrix.append(['Solvated Docking', None, None])
            objectList.append(None)
    
            textMatrix.append(['Do Water Docking', run.doWaterDock, 'Whether to perform solvated docking'])
            objectList.append([run,'doWaterDock', None, self.toggleBoolean, None])
    
            textMatrix.append(['Use DB solvation', run.useDbSolvateMethod, 'Wheter to use database driven solvated docking'])
            objectList.append([run,'useDbSolvateMethod', None, self.toggleBoolean, None])

            textMatrix.append(['Do Water analysis', run.doWaterAnalysis, 'Whether to perform some water analysis'])
            objectList.append([run,'doWaterAnalysis', None, self.toggleBoolean, None])

            textMatrix.append(['Rigid body translation', run.doRigidBodyWaterTrans, 'Allows translation of water molecules during rigid-body docking'])
            objectList.append([run,'doRigidBodyWaterTrans', None, self.toggleBoolean, None])

            textMatrix.append(['Initial restraint cutoff','%1.1f' % run.waterInitRestCutoff, 'Initial cutoff for restraints solvating method'])
            objectList.append([run,'waterInitRestCutoff', self.floatEntry, self.getFloat, self.setFloat]) 

            textMatrix.append(['Restraint cutoff','%1.1f' % run.waterRestCutoff, 'Cutoff for restraints solvating method'])
            objectList.append([run,'waterRestCutoff', self.floatEntry, self.getFloat, self.setFloat]) 

            textMatrix.append(['Force constant','%1.1f' % run.waterRestScale, 'Force constant for restraints solvating method'])
            objectList.append([run,'waterRestScale', self.floatEntry, self.getFloat, self.setFloat]) 

            textMatrix.append(['Water to keep','%1.2f' % run.waterToKeep, 'Fraction of water to keep'])
            objectList.append([run,'waterToKeep', self.floatEntry, self.getFloat, self.setFloat]) 

            textMatrix.append(['Water to add','%1.2f' % run.waterToAddRandom, 'Random fraction to be added to the fraction of water to keep'])
            objectList.append([run,'waterToAddRandom', self.floatEntry, self.getFloat, self.setFloat]) 

            textMatrix.append(['Water-surface cutoff','%1.1f' % run.waterSurfaceCutoff, 'Water-protein surface-cutoff'])
            objectList.append([run,'waterSurfaceCutoff', self.floatEntry, self.getFloat, self.setFloat]) 

            #Annealing Protocols and Haddock specific energy terms
            termId = len([ec.code for ec in run.sortedHaddockEnergyTerms()])+1
            for annealProtocolStore in ['dockingProtocolStore','distRestraintEnergyStore','dihRestraintEnergyStore','semiflexInterMolScalingStore']:
                protocol = run.findFirstHaddockEnergyTerm(code=annealProtocolStore)
                if protocol:
                    textMatrix.append([eval(annealProtocolStore)['details'], None, None])
                    objectList.append(None)
            
                    for term in protocol.sortedEnergyTermParameters():
                        textMatrix.append([term.code,'%1.3f' % term.value, 'Energy constant for term %s' % term.code])
                        objectList.append([term,'value', self.floatEntry, self.getFloat, self.setFloat])
            
                else:
                    protocol = eval(annealProtocolStore)
                    energyTermStore = run.newHaddockEnergyTerm(code=annealProtocolStore,termId=termId)
            
                    textMatrix.append([protocol['details'], None, None])
                    objectList.append(None)
            
                    terms = protocol['terms'].keys()
                    terms.sort()
                    for term in terms:
                        energyTerm = energyTermStore.newEnergyTermParameter(code=term,value=protocol['terms'][term]) 
                        textMatrix.append([term,'%1.3f' % energyTerm.value, 'Energy constant for term %s' % term])
                        objectList.append([energyTerm,'value', self.floatEntry, self.getFloat, self.setFloat])
            
                    termId += 1
    
            #Automatic distance restraint weighting    
            textMatrix.append(['Use automated distance restraints weighting', None, None])
            objectList.append(None)
            
            textMatrix.append(['Automatic distance restraints', run.doAirScaling, 'Whether to do automatic distance restraint weighting'])
            objectList.append([run,'doAirScaling', None, self.toggleBoolean, None])

            textMatrix.append(['Num distance restraints', run.numUnambRestautoAir, 'Define the number of distance restraints for automated weighting'])
            objectList.append([run,'numUnambRestautoAir', self.intEntry, self.getInt, self.setNonNegInt])    

            textMatrix.append(['Num AIR restraints', run.numAmbRestautoAir, 'Define the number of AIR restraints for automated weighting'])
            objectList.append([run,'numAmbRestautoAir', self.intEntry, self.getInt, self.setNonNegInt])    
    
            protocol = run.findFirstHaddockEnergyTerm(code='autoDistanceRestraintWeightStore')
            if protocol:
                for term in protocol.sortedEnergyTermParameters():
                    textMatrix.append([term.code,'%1.1f' % term.value, 'Energy constant for term %s' % term.code])
                    objectList.append([term,'value', self.floatEntry, self.getFloat, self.setFloat])
        
            else:
                energyTermStore = run.newHaddockEnergyTerm(code='autoDistanceRestraintWeightStore',termId=termId)
                terms = autoDistanceRestraintWeightStore['terms'].keys()
                terms.sort()
                for term in terms:
                    energyTerm = energyTermStore.newEnergyTermParameter(code=term,value=autoDistanceRestraintWeightStore['terms'][term]) 
                    textMatrix.append([term,'%1.1f' % energyTerm.value, 'Energy constant for term %s' % term])
                    objectList.append([energyTerm,'value', self.floatEntry, self.getFloat, self.setFloat])

            #Scoring Weights
            textMatrix.append(['Scoring weights', None, None])
            objectList.append(None)
            stage_name = ['it0','it1','w']

            scoringWeights = [(sw.term, sw.stage, sw) for sw in run.scoringWeights]
            if scoringWeights:
                scoringWeights.sort()
                for term, stage, scoringWeight in scoringWeights:
                    textMatrix.append(['%s stage %d' % (term,stage), scoringWeight.value, 'Weighting factor for scoring %s terms in stage %s' % (term,stage_name[stage])])
                    objectList.append([scoringWeight,'value', self.floatEntry, self.getFloat, self.setPosFloat])

            else:
                for term in DEFAULT_SCORE:
                    for stage in range(len(DEFAULT_SCORE[term])):
                        scoringWeight = run.newScoringWeight(term=term,stage=stage,value=DEFAULT_SCORE[term][stage])
                        textMatrix.append(['%s stage %d' % (term,stage),DEFAULT_SCORE[term][stage], 'Weighting factor for scoring term %s in stage %s' % (term,stage_name[stage])])
                        objectList.append([scoringWeight,'value', self.floatEntry, self.getFloat, self.setPosFloat])


            #Analysis and clustering
            textMatrix.append(['Analysis and clustering', None, None])
            objectList.append(None)

            textMatrix.append(['Num Analysis Structures?', run.numAnalysisStructures, 'Number of structures to be analysed'])
            objectList.append([run,'numAnalysisStructures', self.intEntry, self.getInt, self.setNonNegInt])    

            textMatrix.append(['Analysis Cluster RMSD','%1.2f' % run.analysisClustRmsd, 'RMSD cutoff for clustering'])
            objectList.append([run,'analysisClustRmsd', self.floatEntry, self.getFloat, self.setPosFloat])

            textMatrix.append(['Analysis Cluster Size', run.analysisClustSize, 'Minimum cluster size'])
            objectList.append([run,'analysisClustSize', self.intEntry, self.getInt, self.setNonNegInt])    

            textMatrix.append(['Analysis H-Bond Dist','%1.2f' % run.analysisDistHBond, 'Cutoff distance (proton-acceptor) to define\nan hydrogen bond'])
            objectList.append([run,'analysisDistHBond', self.floatEntry, self.getFloat, self.setPosFloat])

            textMatrix.append(['Analysis Nonbonded Dist','%1.2f' % run.analysisDistNonbond, 'Cutoff distance (carbon-carbon) to define\nan hydrophobic contact'])
            objectList.append([run,'analysisDistNonbond', self.floatEntry, self.getFloat, self.setPosFloat])    

            #Haddock program settings for local use
            textMatrix.append(['Haddock program settings', None, None])
            objectList.append(None)
            
            textMatrix.append(['Haddock program directory', run.haddockDir, 'Defines the path to the Haddock software directory'])
            objectList.append([run,'haddockDir', self.stringEntry, self.getString, self.setString])
            
            textMatrix.append(['CNS executable', run.cnsExecutable, 'Defines the path to the CNS executable'])
            objectList.append([run,'cnsExecutable', self.stringEntry, self.getString, self.setString])
            
            textMatrix.append(['CPU number', run.cpuNumber, "Defines the number of cpu's used in a cluster \nfor submission of Haddock jobs"])
            objectList.append([run,'cpuNumber', self.intEntry, self.getInt, self.setNonNegInt])
            
            textMatrix.append(['Queue command', run.queueCommand, "Defines the queue type for Haddock job submission.\n 'csh' is default for running on a single machine in C-shell"])
            objectList.append([run,'queueCommand', self.stringEntry, self.getString, self.setString])
            
            # Finally update the ScrolledMatrix table
            n = len(textMatrix)
            m = 3
            colors = [None] * m

            colorMatrix = [colors for i in range(n)]
            colorMatrix[0]  = [headingColor] * m            # Restraint definition
            colorMatrix[10] = [headingColor] * m            # Symmetry Operations
            colorMatrix[10+symRows] = [headingColor] * m    # Docking protocol
            colorMatrix[25+symRows] = [headingColor] * m    # Energy and Interaction Parameters
            colorMatrix[33+symRows] = [headingColor] * m    # Solvated docking
            colorMatrix[44+symRows] = [headingColor] * m    # Docking molecular simulation protocol
            colorMatrix[61+symRows] = [headingColor] * m    # Distance restraints energy constants
            colorMatrix[80+symRows] = [headingColor] * m    # Dihedral restraint energy constants
            colorMatrix[85+symRows] = [headingColor] * m    # Scaling of intermolecular interactions
            colorMatrix[92+symRows] = [headingColor] * m    # Automatic distance restraint weighting
            colorMatrix[112+symRows] = [headingColor] * m    # Scoring weights
            colorMatrix[146+symRows] = [headingColor] * m    # Analysis and clustering
            colorMatrix[152+symRows] = [headingColor] * m    # Haddock program settings

        table.update(colorMatrix=colorMatrix, objectList=objectList, textMatrix=textMatrix)
        self.updateRunPulldown()

    def useSymmetry(self, widget, symmetryOp, attrName):
  
        if self.run:
            if symmetryOp in self.run.symmetryRestraints: self.run.removeSymmetryRestraint(symmetryOp)
            else: self.run.addSymmetryRestraint(symmetryOp)

            self.updateRuns()

    def getDielectricType(self, widget, obj, attrName):
   
        value = getattr(obj, attrName)
        names = self.hProject.metaclass.container.getElement('DielectricType').enumeration
        index = -1

        if value in names: index = names.index(value)
        else: index = 0

        widget.setup(names,index)

    def getNonBondedType(self, widget, obj, attrName):
  
        value = getattr(obj, attrName) 
        names = self.hProject.metaclass.container.getElement('NonBondedType').enumeration
        index = -1

        if value in names: index = names.index(value)
        else: index = 0

        widget.setup(names,index)
        
    def getSolvent(self, widget, obj, attrName):
  
        value = getattr(obj, attrName) 
        names = self.hProject.metaclass.container.getElement('Solvent').enumeration
        index = -1

        if value in names: index = names.index(value)
        else: index = 0

        widget.setup(names,index)

    def addNewRun(self):
  
        """Add new run within the 'Runs' tab"""

        self.run = self.hProject.newRun()
        self.updateAllAfter()
        
    def removeRun(self):
  
        """Remove a run within the 'Runs' tab."""

        if self.run:
            self.run.delete()
            self.run = None
            self.updateAllAfter()
            
    def copyRun(self):
      
        """Copy a run within the 'Runs' tab."""

        runA = self.run
        if runA:
            self.run = copyRun(runA)
            self.updateAllAfter()
   
    def addMolPartner(self):
  
        """Add new HADDOCK molecular partner to project. Cannot define more than 6 partners.
           Automaticly add new partner code and set number of models to maximum number in
           ensemble.
        """
        if self.hProject:
            if not self.ccpnProject.molSystems:
                showWarning('Warning','No molecular systems present in CCPN project',parent=self)
                return

            if len(self.hProject.haddockPartners) < 6: 
                code = 'A'
                while self.hProject.findFirstHaddockPartner(code=code): code = chr(ord(code)+1)

                molSystem = self.ccpnProject.findFirstMolSystem()
                partner = self.hProject.newHaddockPartner(code=code, molSystem=molSystem)
                setPartnerChains(partner, molSystem.chains)
                
                models = []
                ensembles = self.ccpnProject.sortedStructureEnsembles()
                if molSystem: 
                    ensembles = [e for e in ensembles if e.molSystem is molSystem]
                    for ensemble in ensembles:    
                        models += [model for model in ensemble.sortedModels()]    
                    partner.structureEnsemble = ensembles[-1]
                    if not partner.structureEnsemble.getDetails():
                        partner.structureEnsemble.setDetails(str(len(models)))    

                self.updateAllAfter()
            else: showWarning('Warning','Cannot set more than six molecular partners',parent=self)

    def removeMolPartner(self):
  
        if self.molPartner:
            self.molPartner.delete()
            self.updateAllAfter()
            
    def editMolPartners(self, widget, obj, attrName):

        """Basicly set the GUI to tab 1, the 'Molecular Partners' tab."""

        self.tabbedFrame.select(2)

    def editResRestraints(self, widget, obj, attrName):

        """Basicly set the GUI to tab 2, the 'Residues' tab. Update the 'General'
           tab afterwards."""
  
        self.molPartner = obj
        self.updateMolPartnerPulldown() # i.e. select in target tab
        self.updateSemiFlexModePulldown()
        self.tabbedFrame.select(3)
        
    def editRuns(self, widget, obj, attrName):
  
        """Basicly set the GUI to tab 3, the 'Runs' tab. Update the 'General'
           tab afterwards."""

        self.run = obj
        self.updateRunPulldown() # i.e. select in target tab
        self.tabbedFrame.select(4)
        
    def editRestraintTerms(self, widget, obj, attrName):
  
        self.run = obj
        self.updateRunPulldown() # i.e. select in target tab
        self.tabbedFrame.select(5)

  #
  # Specific widget calls
  #
  
    def getConstraintLists(self):
  
        cLists = []

        if self.run and self.run.nmrConstraintStore:
            if self.energyTerm:
                for constraintList in self.run.nmrConstraintStore.sortedConstraintLists():
                    if constraintList.className == energyTermTypes[self.energyTerm.code]: cLists.append(constraintList)

            else: cLists = self.run.nmrConstraintStore.sortedConstraintlists()
            
        return cLists
        
    def getRestrList(self, energyTerm):
  
        names = []
        index = -1

        constraintLists = self.getConstraintLists()

        if constraintLists:
            names = ['%s:(%s)' % (cl.serial,cl.className[:-14]) for cl in constraintLists]
            #if energyTerm.constraintList not in constraintLists: energyTerm.constraintList = None
            #index = constraintLists.index(energyTerm.constraintList)
 
        self.restrListPulldown.setup(names, index)
        
    def getCnsFile(self, energyTerm):
      
        if self.run and self.energyTerm:
            file_types = [ FileType(".tbl", ["*.tbl"]),FileType("All", ["*"]) ]

            fileName = energyTerm.fileName
            if fileName and not os.path.exists(fileName): fileName = None

            popup = FileSelectPopup(self, file_types, file=fileName, dismiss_text='Cancel')
            fileName = popup.getFile()

            if fileName:
                energyTerm.fileName = fileName
                energyTerm.constraintList = None
            
            self.updateAllAfter()
            
    def setRestrList(self, obj, name=None):
 
        index =  self.restrListPulldown.getSelectedIndex()

        if self.run and self.energyTerm:
            constraintList = self.getConstraintLists()[index]
            self.energyTerm.constraintList = constraintList

            if constraintList is not None: self.energyTerm.fileName = None

            self.updateAllAfter()
            
    def setMolSystem(self, partner, name=None):
   
        """Get all molsystems as stored in the project as list. Automaticly set chains and models
           to all available
        """
        index =  self.molSysPulldown.getSelectedIndex()
        molSystems = self.ccpnProject.sortedMolSystems()
   
        if molSystems:
            molSystem = molSystems[index]
            self.molPartner.molSystem = molSystem
            chains = molSystem.sortedChains()
            if chains: setPartnerChains(self.molPartner,chains)

            ensembles = self.ccpnProject.sortedStructureEnsembles()
            if molSystem: 
                ensembles = [e for e in ensembles if e.molSystem is molSystem]    
                self.molPartner.structureEnsemble = ensembles[-1]

        self.updateAllAfter()
        
    def setChains(self, obj):

        """Get the list of chains for the selected molsystem"""
        if self.molPartner and self.molPartner.molSystem:
            if obj is not None:
                chains = self.molPartner.molSystem.sortedChains()
                values = self.chainSelect.get()
                chains = [chains[i] for i in range(len(values)) if values[i]]
                setPartnerChains(self.molPartner,chains)
            else:
                pass    
        self.molPartnerMatrix.keyPressEscape()  
        self.updateAllAfter()
        
    def setModel(self, obj, name=None):
    
        """Set the selected model in the structureEnsemble.details"""

        index =  self.modelPulldown.getSelectedIndex()
        models = self.getModels()
        
        if models:
            if index <= len(models):
                model = models[index]
                self.molPartner.structureEnsemble = model.parent
                self.molPartner.structureEnsemble.setDetails(str(models[index].serial))

        self.updateAllAfter()
    
    def toggleIsDna(self,partner):
        
        """Turn the 'isDna' tag for a partner on/off set the forcefield to TOPALLHDG for proteins and DNA for DNA/RNA"""
        
        partner.isDna = not partner.isDna
    
        if partner.isDna == True: partner.forceFieldCode = 'DNA'
        else: partner.forceFieldCode = 'TOPALLHDG'    
        
        self.updateAllAfter()
        
    def toggleAutoHis(self, partner):

        """Set auto histidine protonation state on/off"""

        partner.autoHistidinePstate = not partner.autoHistidinePstate
        self.updateAllAfter()
        
    def getModels(self):
  
        """Return list of structure ensemble models for selected partner"""

        models = []
        if self.molPartner:
            molSystem = self.molPartner.molSystem
            ensembles = self.ccpnProject.sortedStructureEnsembles()
            if molSystem: ensembles = [e for e in ensembles if e.molSystem is molSystem]
            for ensemble in ensembles:    
                models += [model for model in ensemble.sortedModels()]
        
        return models
        
    def getMolSystem(self, partner):

        """Select molecular system from list of molsystems stored in the project"""

        names = []; index = -1
        molSystem = partner.molSystem
        molSystems = self.ccpnProject.sortedMolSystems()
        if molSystems:
            names = [ms.code for ms in molSystems]
            if molSystem not in molSystems: molSystem = molSystems[0]
            index = molSystems.index(molSystem)

        self.molSysPulldown.setup(names, index)
        
    def getChains(self, partner):
  
        names  = []
        values = []
        
        molSystem = partner.molSystem
        if molSystem:
            for chain in molSystem.sortedChains():
                names.append(chain.code)
                if not partner.chains: values.append(True)
                elif partner.findFirstChain(chain=chain): values.append(True)
                else: values.append(False)
                self.chainSelect.set(values=values,options=names)   
        else:
            showWarning('Warning','Set Mol System or ensemble first',parent=self)
            self.molPartnerMatrix.keyPressEscape()        
            
    def getModel(self, partner):

        names = []; index = -1

        models = self.getModels()
        if models: names = [m.serial for m in models]
        
        storedModelIndex = self.molPartner.structureEnsemble.getDetails()
        if storedModelIndex: index = int(storedModelIndex)
    
        self.modelPulldown.setup(names, index-1)

  #
  # Generic widget calls
  #
    def selectRowObj(self, obj, row, col, table):
  
        self.table = table
        self.rowObj = obj

    def getValue(self, rowObj):
 
        # get correct object, widget and get/set functions
        obj, attrName, widget, getter, setter = rowObj

        # Bit of a hack because widgets are normally set at construction per column
        # Now setting it per row, and hence on-the-fly
        self.table.editWidget = widget
    
        # get current value & setup widget
        getter(widget, obj, attrName)
        
    def setValue(self, event, null=None): # null is for pulldown menu callbacks passing name,index

        # get correct widget and get/set functions
        obj, attrName, widget, getter, setter = self.rowObj
    
        # set and check the appropriate parameter value from current edit widget

        # no setter for boolean toogles - the getter will have done all the toggling already
        # no setter for file selects - the file popup gets and sets the value and cannot be interrupted
        if setter: setter(widget, obj, attrName)

        self.table.refreshFunc() # Update relevant table (ScrolledMatrix)
    
    def getString(self, widget, obj, attrName):

        value = getattr(obj, attrName) 
        widget.set(value)
        
    def setString(self, widget, obj, attrName):

        value = widget.get().strip()

        if value: setattr(obj, attrName, value)
        
    def setPulldownMenu(self, widget, obj, attrName):

        value = widget.getSelected()
        setattr(obj, attrName, value)
        
    def toggleBoolean(self, widget, obj, attrName):

        value = getattr(obj, attrName)
        setattr(obj, attrName, not value)

        self.table.refreshFunc()
        
    def getDirectory(self, widget, obj, attrName):
    
        """Is called by updateGeneral. No widget because table placement is irrelevent and setting process not interruptable"""

        value = getattr(obj, attrName) 
        if not os.path.exists(value): value = None

        popup = FileSelectPopup(self, directory=value, show_file=False)
        value = popup.getDirectory()
        if value:
            setattr(obj, attrName, value)
            self.table.refreshFunc()
            
    def getFile(self, widget, obj, attrName):
        # No widget because table placement is irrelevent and setting process not interruptable

        file_types = [ FileType("All", ["*"]), ]

        value = getattr(obj, attrName) 
        if not os.path.exists(value): value = None

        popup = FileSelectPopup(self, file_types, file=value, dismiss_text='Cancel')
        value = popup.getFile()
        if value:
            setattr(obj, attrName, value)
            self.table.refreshFunc()
    
    def getInt(self, widget, obj, attrName):
  
        value = getattr(obj, attrName) 
        widget.set(value)
        
    def setInt(self, widget, obj, attrName):

        value = widget.get() or 0
        setattr(obj, attrName, value)
        
    def setNonNegInt(self, widget, obj, attrName):

        value = widget.get() or 0
        if value >= 0: setattr(obj, attrName, value)
        
    def setPosInt(self, widget, obj, attrName):

        value = widget.get() or 0
        if value > 0: setattr(obj, attrName, value)
    
    def getFloat(self, widget, obj, attrName):
  
        value = getattr(obj, attrName) 
        widget.set(value)
    
    def setFloat(self, widget, obj, attrName):

        value = widget.get() or 0.0
        setattr(obj, attrName, value)
    
    def setFraction(self, widget, obj, attrName):
   
        value = min(1.0,max(0.0, widget.get() or 0.0))
        setattr(obj, attrName, value)
    
    def setPosFloat(self, widget, obj, attrName):

        value = widget.get() or 0.0
        if value > 0.0: setattr(obj, attrName, value)
    
    def setNonNegFloat(self, widget, obj, attrName):

        value = widget.get() or 0.0
        if value >= 0.0: setattr(obj, attrName, value)
