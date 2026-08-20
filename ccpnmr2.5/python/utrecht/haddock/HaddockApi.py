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
import os

from memops.general.Io import loadProject

from HaddockExportClassic    import    exportClassic
from HaddockExportParam        import    exportParam
from HaddockImportRunCns    import    runCnsImporter
from HaddockServerUpload    import     ServerUpload
from HaddockBasic            import  getStructureFromFile, addRdcParam, addDaniParam
from HaddockLocal            import  *

class HaddockApi(object):
    
    """Description: API for operation of the Haddock package from within a Python script.
       Input      : A valid ccpnProject instance
       Output      : None
       Arguments  : Set debug to True for additional info
    """
    def __init__(self, ccpnProject=None, debug=False):
        
        if ccpnProject: 
            self.ccpnProject = ccpnProject
        else:
            if self.debug: print("ERROR: compulsary ccpn project argument missing.")
            return
        
        self.debug         = debug
        self.projects     = []
    
        self.__initStoredProjects()
    
    def __initStoredProjects(self):
        
        """Description: Instantiates 'HaddockProject' instances from all stored projects
                         in the ccpnProject
           Input      : self.ccpnProject instance as supplied to 'HaddockApi'
           Output      : None, filles the self.projects list with HaddockApi.HaddockProject 
                        instances.
        """
        for hProject in self.ccpnProject.sortedHaddockProjects():
            self.projects.append(HaddockProject(hProject=hProject,
                                                ccpnProject=self.ccpnProject,
                                                debug=self.debug))
    
    def removeHaddockProject(self,hProject=None,name=None):
        
        """Description: Removes provided Haddock project
           Input      : HaddockApi.HaddockProject instance or the name of the project
           Output     : None
        """
        if hProject:
            if hProject in self.projects:
                if self.debug: print("NOTE: Remove Haddock project: %s" % hProject.name)
                self.projects.remove(hProject)
                hProject.hProject.delete()
            else:
                if self.debug: print("ERROR: Project : %s not in list, could not remove" % hProject.name)    
        elif name:
            found = False
            for hProject in self.projects:
                if name == hProject.name:
                    if self.debug: print("NOTE: Remove Haddock project: %s" % hProject.name)
                    self.projects.remove(hProject)
                    hProject.hProject.delete()
                    found = True
            if not found:
                if self.debug: print("ERROR: No project with name %s" % name)    
                    
    def newHaddockProject(self,name='Default',workingDir='.'):

        """Description: Create new default Haddock project
           Input      : None
           Output      : New HaddockApi.HaddockProject project instance
           Arguments  : Project 'name' (Default by default) project working directory
                        'workingDir' (set to current directory by default)
        """
        for project in self.projects:
            if project.name == name:
                if self.debug: print("ERROR: Haddock project with name %s already made" % name)
                return
        
        if self.debug: print("NOTE: Created new Haddock project with name: %s" % name)
        hProject = HaddockProject(self.ccpnProject.newHaddockProject(name=name,workingDir=workingDir),
                                  ccpnProject=self.ccpnProject,
                                  debug=self.debug)
        
        self.projects.append(hProject)
        return hProject
    
    def getHaddockProject(self,name=None):
        
        """Description: Get Haddock project by name
           Input      : Project name
           Output      : HaddockApi.HaddockProject project instance
        """
        request = [project for project in self.projects if project.name == name]
        
        if len(request): 
            if self.debug: print("NOTE: getHaddockProject, return project with name %s" % name)
            return request[0]
        else:
            if self.debug: print("ERROR: getHaddockProject, no project with name %s" % name)

class HaddockProject(object):
    
    """Description: HaddockProject is a convenient storage class for Haddock projects.
                    It takes care of setting default, processing data and ensuring 
                    validitie while working with Haddock projects
       Input      : A valid ccpn Haddock project instance (hProject) and a valid ccpn
                       project instance (ccpnProject)
       Output      : None
       Arguments  : Set debug to True for additional info
    """
    def __init__(self, hProject=None, ccpnProject=None, debug=False):
        
        self.debug = debug
        self.hProject = hProject
        self.ccpnProject = ccpnProject
        
        self.runs = []
        self.partners = []
        self.airUpperDistanceLimit = 2.000
        
        self.__initStoredRuns()
        self.__initStoredPartners()
    
    def __getattr__(self,name):
        
        """Description: Ensures compatibility with the standard ccpn way of 
                         getting data.
           Input      : Ccpn argument name
           Output     : Ccpn argument output
        """
        return getattr(self.hProject,name)
    
    def __initStoredRuns(self):
        
        """Description: Retrieve all stored Haddock runs for the given project and create
                        a Haddock Api 'HaddockRun' instance of each of them.
                        Gets called by default when initiating a HaddockProject instance
           Input      : self.hProject, self.ccpnProject
           Output      : A HaddockApi.HaddockProject.HaddockRun instance, appended to self.runs
        """
        for run in self.hProject.sortedRuns():
            self.runs.append(HaddockRun(run=run,
                                        ccpnProject=self.ccpnProject,
                                        airUpperDistanceLimit=self.airUpperDistanceLimit,
                                        debug=self.debug))
        
        if not len(self.runs): self.newHaddockRun()                            
    
    def __initStoredPartners(self):
        
        """Description: Retrieve all stored Haddock partners for the given project and create
                        a Haddock Api 'HaddockPartner' instance of each of them.
                        Gets called by default when initiating a HaddockProject instance
           Input      : self.hProject, self.ccpnProject
           Output      : A HaddockApi.HaddockProject.HaddockPartner instance, appended to self.runs
        """
        for partner in self.hProject.sortedHaddockPartners():
            self.partners.append(HaddockPartner(hPartner=partner,
                                                  ccpnProject=self.ccpnProject,
                                                  debug=self.debug))
            
    def newHaddockPartner(self,molSystem=None,pdb=None):
        
        """Description: Create a new Haddock partner in the current project. 
           Input      : Create new Haddock partner from ccpn.molSystem instance using
                        the 'molSystem' argument or create new partner from a PDB file 
                        using the 'pdb' argument (define full path to pdb file). In case
                        of a PDB file the new molecular system is either appended to
                        a ccpn.molSystem carrying the same name or a new ccpn.molSystem
                        with the name of the PDB file.
           Output      : New HaddockApi.HaddockProject.HaddockPartner instance
        """
        if molSystem and not pdb:
            partner = HaddockPartner(ccpnProject=self.ccpnProject,
                                     debug=self.debug)
            partner.initDefaultPartner(molSystem=molSystem,hProject=self.hProject)
                                    
            if partner.isValidPartner:
                self.partners.append(partner)
                if self.debug: print("NOTE: added new Haddock Partner with code: %s" % partner.code)
                return partner
        elif pdb:
            newMolsysName = (os.path.basename(pdb)).split('.')[0]
            molSystem = self.ccpnProject.findFirstMolSystem(code=newMolsysName)

            if molSystem:
                i = 1
                name = newMolsysName
                while self.ccpnProject.findFirstMolSystem(code=name):
                    name = '%s_%d' % (newMolsysName, i)
                    i += 1

                molSystem = self.ccpnProject.newMolSystem(code=name,name=name)
            else:
                molSystem = self.ccpnProject.newMolSystem(code=newMolsysName,name=newMolsysName)
            
            if self.debug: print("NOTE: created new molSystem with name: %s" % newMolsysName)    
            structure = getStructureFromFile(molSystem, pdb, fileType='rough', doWarnings=True)    
            
            if structure:
                partner = HaddockPartner(ccpnProject=self.ccpnProject,
                                          debug=self.debug)
                partner.initDefaultPartner(molSystem=molSystem,hProject=self.hProject)
                                    
                if partner.isValidPartner:
                    self.partners.append(partner)
                    if self.debug: print("NOTE: added new Haddock Partner with code: %s" % partner.code)
                    return partner    
        else:
            if self.debug: print('ERROR: Provide ccpn molSystem instance or PDB file to create Haddock partner')
    
    def removeHaddockPartner(self,partner=None,code=None):
        
        """Description: Remove given Haddock Partner
           Input      : HaddockApi.HaddockProject.HaddockPartner instance or partner code
           Output      : None
        """
        if partner:
            if partner in self.partners:
                if self.debug: print("NOTE: removeHaddockPartner, remove Haddock partner: %s" % partner.code)
                partner.hPartner.delete()
                self.partners.remove(partner)
            else:
                if self.debug: print("ERROR: removeHaddockPartner, no partner with code: %s" % partner.code)
        elif code:
            code = code.upper()
            found = False
            for partner in self.partners:
                if partner.code == code:
                    if self.debug: print("NOTE: removeHaddockPartner, remove Haddock partner: %s" % partner.code)
                    partner.hPartner.delete()
                    self.partners.remove(partner)
                    found = True
            if not found:
                if self.debug: print("ERROR: removeHaddockPartner, no partner with code: %s" % code)    
    
    def getHaddockPartner(self,code=None):
        
        """Description: Get Haddock partner by code
           Input      : Partner code
           Output      : HaddockApi.HaddockProject.HaddockPartner instance
        """
        if code: code = code.upper()
        request = [partner for partner in self.partners if partner.code == code]
        if len(request): 
            if self.debug: print("NOTE: getHaddockPartner, return Haddock Partner with code %s" % code)
            return request[0]
        else:
            if self.debug: print("ERROR: getHaddockPartner, no partner with code %s" % code)        
    
    def newHaddockRun(self):
        
        """Description: Create a new default Haddock run in current project
           Input      : None
           Output      : HaddockApi.HaddockProject.HaddockRun instance
        """
        run = HaddockRun(run=self.hProject.newRun(),
                         ccpnProject=self.ccpnProject,
                         airUpperDistanceLimit=self.airUpperDistanceLimit,
                         debug=self.debug)
        
        self.runs.append(run)
        if self.debug: print("NOTE: Create new Haddock run: %i" % run.serial)
                        
        return run
    
    def removeHaddockRun(self,run=None):
        
        """Description: Remove the given Haddock run from the current project
           Input      : HaddockApi.HaddockProject.HaddockRun instance
           Output      : None
        """
        if run in self.runs:
            if self.debug: print("NOTE: Remove current Haddock run: %i" % run.serial)
            self.runs.remove(run)
            run.run.delete()
    
    def copyHaddockRun(self,run=None, nmrConstraintStore=None):

        """Description: Copy a give run to a new run in the same project
           Input      : A HaddockApi.HaddockProject.HaddockRun instance.
           Output      : New HaddockApi.HaddockProject.HaddockRun instance
           Arguments  : 'nmrConstraintStore' to include constraints in the 
                        new run (default is None).
        """
        if run:
            if not run.haddockProject == self.hProject:
                if self.debug: print("WARNING: Cannot copy Haddock runs between projects")
            else:
                runB = self.hProject.newRun(nmrConstraintStore=nmrConstraintStore)
                if self.debug: print("NOTE: Copy run %i to new run %i" % (run.serial,runB.serial))

                for metaAttr in run.metaclass.attributes:
                    if metaAttr.changeability == 'changeable':
                        attrName = metaAttr.name
                        value = getattr(run,attrName)
                        setattr(runB,attrName,value)

                for symmetry in run.symmetryRestraints:
                    runB.addSymmetryRestraint(symmetry)

                for scoringWeight in run.scoringWeights:
                    runB.newScoringWeight(stage=scoringWeight.stage,
                                          term=scoringWeight.term,
                                          value=scoringWeight.value)

                for term in run.haddockEnergyTerms:
                    termB = runB.newHaddockEnergyTerm(code=term.code,
                                                    termId=term.termId,
                                                    name=term.name,
                                                    fileName=term.fileName,
                                                    details=term.details)

                    if run.nmrConstraintStore is nmrConstraintStore:
                        termB.constraintList = term.constraintList

                    for param in term.energyTermParameters:
                        termB.newEnergyTermParameter(code=param.code, value=param.value)

                copyrun = HaddockRun(run=runB,
                                     ccpnProject=self.ccpnProject,
                                     airUpperDistanceLimit=self.airUpperDistanceLimit,
                                     debug=self.debug)

                self.runs.append(copyrun)        
                return copyrun
    
    def getHaddockRun(self,serial=None):

        """Description: Get Haddock Run by serial
           Input      : Run number
           Output      : HaddockApi.HaddockProject.HaddockRun instance
        """
        request = [run for run in self.runs if run.serial == serial]
        if len(request): 
            if self.debug: print("NOTE: getHaddockRun, return Haddock Run %i" % serial)
            return request[0]
        else:
            if self.debug: print("ERROR: getHaddockPartner, no run with %s" % serial)
    
    def importRunCns(self,runCns=None):
        
        """Description: Import some settings from a Haddock run.cns file into a new Haddock run.
           Input      : Full path to run.cns file
           Output      : New HaddockApi.HaddockProject.HaddockRun instance
        """
        if not runCns or not os.path.isfile(runCns):
            if self.debug: print("ERROR: importRunCns, run.cns file not found or not defined")
        else:
            if self.debug: print("NOTE: import run.cns file: %s" % runCns)
            runcsnimport = runCnsImporter(hProject=self.hProject,runCns=runCns)
            run = HaddockRun(run=runcnsimport.newRun,
                             ccpnProject=self.ccpnProject,
                             airUpperDistanceLimit=self.airUpperDistanceLimit,
                             debug=self.debug)
            self.runs.append(run)
            return run

class HaddockPartner(object):
    
    """Description: HaddockPartner is a convenient storage class for Haddock partners.
                    It takes care of setting default, processing data and ensuring 
                    validitie while working with Haddock partners
       Input      : A valid ccpn Haddock partner instance (hPartner) and a valid ccpn
                       project instance (ccpnProject)
       Output      : None
       Arguments  : Set debug to True for additional info
    """
    def __init__(self, hPartner=None, ccpnProject=None, debug=False):

        self.hPartner = hPartner
        if hPartner: self.hProject = hPartner.haddockProject
        self.ccpnProject = ccpnProject
        self.debug = debug
    
        self.isValidPartner = False
    
    def __getattr__(self,name):
        
        """Description: Ensures compatibility with the standard ccpn way of 
                         getting data.
           Input      : Ccpn argument name
           Output     : Ccpn argument output
        """
        return getattr(self.hPartner,name)

    def initDefaultPartner(self,molSystem=None,hProject=None):
        
        """Description: Creates a new deafult partner. This function gets called by the 'newHaddockParter'
                        function in the HaddockProject class.
           Input      : ccpn.MolSystem instance representing the partner and HaddockApi.HaddockProject 
                         instance to include the partner in.
           Output      : None, self.hPartner is set
        """
        if hProject:
            self.hProject = hProject
            if len(self.hProject.haddockPartners) < 6: 
                code = 'A'
                while self.hProject.findFirstHaddockPartner(code=code): code = chr(ord(code)+1)

                partner = self.hProject.newHaddockPartner(code=code, molSystem=molSystem)
                self.hPartner = partner
                if self.debug: print("NOTE: Made new haddock partner:%s, code:%s" % (self.hPartner.molSystem.name,code))

                self.setAutoHistidinePstate()            # Set pState to True by default
            
                ensembles = self.getHaddockEnsemble()
                for ensemble in ensembles:                        # Set the default Haddock ensemble to all
                    self.setHaddockEnsemble(ensemble=ensemble)    # ensembles in the MolSystem
                
                self.isValidPartner = True
            else:
                if self.debug: print('WARNING: Cannot set more than six molecular partners')
        else:
            if self.debug: print('ERROR: Compulsary Haddock project not provided')                                        

    def getPartnerResidues(self):
        
        """Description: Returns a list of all residues belonging to the Haddock partner
           Input      : None, uses self.hPartner
           Output      : List of ccp.molecule.MolSystem.Residue instances
        """
        residues = []
        
        hChains = self.hPartner.sortedChains()
        for hChain in hChains:
            hResidues = [r for r in hChain.residues]
            residues += hResidues
        
        residues.sort()
        return residues    

    def setResidueFlexibility(self,residues=None,state='none'):
        
        """Description: Set the flexibility state for the selected residue of the Haddock Partner
                        Flexibility is used during the (semi)-flexible refinement stage of Haddock.
                        If the partners 'semiFlexMode' parameter is set to 'automatic' (default) than
                        manual set flexibility of the residues will be ignored.
           Input      : List of ccp.molecule.MolSystem.Residue instances. Flexibility state (none,semi,full)
           Output      : None, sets the flexibility state in the model
        """
        if residues:
            if not type(residues) == type([]): residues = [residues]
            
            if state in ['none','semi','full']:
                for residue in residues:
                    residue.flexibility = state
                    if self.debug: print("NOTE: Set flexibility of residue: %s-%s to '%s'" % 
                                        (residue.residue.ccpCode, residue.residue.seqCode, state))
            else:
                if self.debug: print("ERROR: setResidueFlexibility. State %s not allowed (none,semi,full)" % state)        
        else:
            if self.debug: print("ERROR: setResidueFlexibility. No residues provides")            
    
    def setSemiFlexMode(self,mode='automatic'):
        
        """Description: Set the semi-flexibility mode to either manual or automatic
           Input      : String, 'manual' or 'automatic'.
           Output      : None, sets the semiFlexMode in the model
        """
        mode = mode.lower()
        if mode in ['manual','automatic']:
            self.hPartner.semiFlexMode = mode
            if self.debug: print("NOTE: set semiFlexMode to %s" % repr(mode))
        
    def setResidueAirState(self,residues=None,state='none'):            

        """Description: Set the AIR state for the selected residue of the Haddock Partner
           Input      : List of ccp.molecule.MolSystem.Residue instances. AIR state (none,active,passive)
           Output      : None, sets the AIR state in the model
        """
        if residues:
            if not type(residues) == type([]): residues = [residues]
            
            if state in ['none','active','passive']:
                for residue in residues:
                    residue.interaction = state
                    if self.debug: print("NOTE: Set AIR state of residue: %s-%s to '%s'" % 
                                        (residue.residue.ccpCode, residue.residue.seqCode, state))
            else:
                if self.debug: print("ERROR: setResidueAirState. State %s not allowed (none,active,passive)" % state)        
        else:
            if self.debug: print("ERROR: setResidueAirState. No residues provides")

    def setAutoHistidinePstate(self,pState=True):
        
        """Description: Turn the 'autoHistidinePstate' option for the partner to on/off
           Input      : pState True or False
           Output      : None, writes to model
        """
        self.hPartner.autoHistidinePstate = pState
        if self.debug: print("NOTE: AutoHistidinePstate set to :%s" % repr(pState))        

    def setForceField(self,molType=None):
        
        """Description: Defines the forfield type used in Haddock calculation based on the
                        molecule type of the Haddock partner. If partner is RNA or DNA
                        that 'forceFieldCode' is set to RNA or DNA else the code is set to
                        'TOPALLHDG.
           Input      : ccpn.molSystem.molType
           Output      : None, writes to model
        """
        if molType in ['DNA','RNA']:
            self.hPartner.isDna = True
            self.hPartner.forceFieldCode = molType
        else:
            self.hPartner.isDna = False
            self.hPartner.forceFieldCode = 'TOPALLHDG'

        if self.debug: print("NOTE: Haddock Partner isDNA:%s, forcefield code set to:%s" % 
                            (repr(self.hPartner.isDna),self.hPartner.forceFieldCode))        
    
    
                            
    def getHaddockEnsemble(self):
        
        """Description: Returns a list of structures ensembles belonging to the Haddock partner
           Input      : None, uses self.hPartner
           Output       : List of ccpn.molecule.molStructure.StructureEnsemble instance
        """
        molSystem = self.hPartner.molSystem
        ensembles = self.ccpnProject.sortedStructureEnsembles()
        if molSystem: return [e for e in ensembles if e.molSystem is molSystem]

    def setHaddockEnsemble(self,ensemble=None):
        
        """Description: Define ensembles belonging to the Haddock partner
           Input      : (List of) ccpn.molecule.molStructure.StructureEnsemble instance
           Output     : None, adds to model
        """
        if ensemble:
            self.hPartner.structureEnsemble = ensemble
            if not self.hPartner.chains:
                molSystem = ensemble.molSystem
                self.hPartner.ccpMolSysCode = molSystem.code

                chains = [c.chain for c in ensemble.coordChains]
                self.setHaddockPartnerChains(chains=chains)
            if self.debug: print("NOTE: Set Haddock ensemble for partner: %s to %s" % 
                                (self.hPartner.code,ensemble.ensembleId))
        else:
            self.setHaddockPartnerChains(chains=[])                            

    def getHaddockPartnerChains(self):
        
        """Description: Return a list of chains belonging to the Haddock partner
           Input      : None, uses self.hPartner
           Output     : List of ccpn.molecule.MolSystem.Chain instances
        """
        return [chain for chain in self.hPartner.molSystem.sortedChains()]            

    def setHaddockPartnerChains(self,chains=None):
        
        """Description: Define the chains belonging to a Haddock partner
           Input      : ccpn.molSystem.chains instances. You get get the list using 
                        'getHaddockPartnerChains' function
           Output      : None, writes to model
        """    
        if not type(chains) == type([]): chain = [chains]
        
        for hChain in self.hPartner.chains:
            if hChain.chain not in chains: hChain.delete()

        molType = None
        for chain in chains:
            if molType is None: molType = chain.molecule.molType
            elif molType != chain.molecule.molType:
                if self.debug: print('CCPN-HADDOCK setPartnerChains failed: Chains not of same type')
                return

        for chain in chains:
            hChain = self.hPartner.findFirstChain(chain=chain)
            if not hChain: self.hPartner.newChain(chain=chain)

            molType = chain.molecule.molType
        
        self.setForceField(molType=molType)    

        # Curate all haddock residue numbers and CCPN residue links
        hSeqId = 1
        for hChain in self.hPartner.sortedChains():
            chain     = hChain.chain
            residues  = chain.sortedResidues()
            hResidues = hChain.sortedResidues()

            resDict = {}
            for hResidue in hResidues:
                residue = hResidue.residue
                resDict[residue] = hResidue

            for residue in residues:
                hResidue = resDict.get(residue)

                if hResidue is None: hResidue = hChain.newResidue(haddockSeqId=hSeqId,residue=residue)
                else: hResidue.haddockSeqId = hSeqId

                hSeqId += 1    

class HaddockRun(object):
    
    """Description: HaddockRun is a convenient storage class for Haddock runs.
                    It takes care of setting default, processing data and ensuring 
                    validitie while working with Haddock runs
       Input      : A valid ccpn Haddock run instance (run) and a valid ccpn
                       project instance (ccpnProject)
       Output      : None
       Arguments  : Set debug to True for additional info
    """
    def __init__(self, run=None, ccpnProject=None, debug=False, airUpperDistanceLimit=None):

        self.debug = debug
        self.run = run
        self.hProject = run.haddockProject
        self.ccpnProject = ccpnProject
        self.airUpperDistanceLimit = airUpperDistanceLimit
    
        self.__initProtocolStores()
    
    def __getattr__(self,name):
        
        """Description: Ensures compatibility with the standard ccpn way of 
                         getting data. Can also be used to set most of the run parameters
           Input      : Ccpn argument name
           Output     : Ccpn argument output
        """
        return getattr(self.run,name)
    
    def __initProtocolStores(self):
        
        """Description: Initiates the energyProtocolStores that are not in the model by default
           Input      : None
           Output      : None
        """
        termId = len([ec.code for ec in self.run.sortedHaddockEnergyTerms()])+1
        for annealProtocolStore in ['dockingProtocolStore','distRestraintEnergyStore','dihRestraintEnergyStore',
                                    'semiflexInterMolScalingStore','autoDistanceRestraintWeightStore']:
            protocol = self.run.findFirstHaddockEnergyTerm(code=annealProtocolStore)
            if not protocol:
                protocol = eval(annealProtocolStore)
                energyTermStore = self.run.newHaddockEnergyTerm(code=annealProtocolStore,termId=termId)
                
                terms = protocol['terms'].keys()
                terms.sort()
                for term in terms:
                    energyTerm = energyTermStore.newEnergyTermParameter(code=term,value=protocol['terms'][term])
        
                termId += 1
        
        stage_name = ['it0','it1','w']
        scoringWeights = [(sw.term, sw.stage, sw) for sw in self.run.scoringWeights]
        if not scoringWeights:
            for term in DEFAULT_SCORE:
                for stage in range(len(DEFAULT_SCORE[term])):
                    scoringWeight = self.run.newScoringWeight(term=term,stage=stage,value=DEFAULT_SCORE[term][stage])
                
    def exportClassicProject(self):
        
        """Description: Exports current project and run as a 'classical' HADDOCK style project. This means
                           a root directory bearing the projects name containing all PDB structure files that
                        need to be docked, a new.html file, restraint files, a ensemble.list file if multiple 
                        ensembles for a given structure are used. If DNA then set 'useDNARestraints' to True. 
                        In the GUI you have to set this manually.
           Input:        None, uses self.run, self.ccpnProject
           Output:        The various files as described above.                
        """
        if self.debug:
            print("NOTE: Export Haddock project:%s, run:%i as Classic project" % 
                 (self.hProject.name,self.run.serial))
            print("      Storage location: %s" % self.hProject.workingDir)
        
        for partner in self.hProject.sortedHaddockPartners():
            if partner.isDna: self.run.useDnaRestraints = True
            
        classic = exportClassic(hProject=self.hProject,
                                latestRun=self.run,
                                ccpnProject=self.ccpnProject)

    def exportHaddockParameterFile(self):

        """Description:    Export a 'new' type HADDOCK project. The project is exported as a self containing 
                        parameter file. It can be uploaded to the HADDOCK server. If DNA then set 'useDNA
                        -Restraints' to True. In the GUI you have to set this manually. 
           Input:        None, uses self.run, self.ccpnProject
           Output:        Haddock project parameter file
        """
        if self.debug:
            print("NOTE: Export Haddock project:%s, run:%i as Haddock parameter file" % 
                 (self.hProject.name,self.run.serial))
            print("      Storage location: %s" % self.hProject.workingDir)
        
        for partner in self.hProject.sortedHaddockPartners():
            if partner.isDna: self.run.useDnaRestraints = True
        
        paramfile = exportParam(hProject=self.hProject,
                                latestRun=self.run,
                                ccpnProject=self.ccpnProject)
        paramfile.writeToFile()

    def runOnHaddockServer(self,username=None,password=None):
        
        """Description: Class for automatic upload of a HADDOCK webserver parameter file.
           Input      : Valid Hadddock web server username and password
           Output      : Server exeptance or rejection messages. Further communication via
                        users e-mail adress.
        """
        if self.debug: print("NOTE: Prepaire Haddock project:%s, run:%i for server upload" % 
                            (self.hProject.name,self.run.serial))

        if not username or not password:
            if self.debug: print("ERROR: No username or password set for server upload")                    

        paramfile = exportParam(hProject=self.hProject,
                                latestRun=self.run,
                                ccpnProject=self.ccpnProject)    
        
        server = ServerUpload(paramfile.filestring,self.hProject.name,self.run.serial,username,password)                            

    def newRestraintSet(self,termType=None,fileName=None,constraintList=None):
        
        """Description: Add restraint sets to the current run. Supported types: UNAMBIG, RDC, HBOND, 
                        DIHEDRAL, AMBIG, DANI.
           Input      : Type of restraint (termType) and path to restraint file on disk (fileName) or
                        Ccpn constraintList (constraintList).
           Output      : New energyTerm
        """
        termType = termType.upper()
        allowedTerms = ['UNAMBIG','RDC','HBOND','DIHEDRAL','AMBIG','DANI']
        
        if not fileName and not constraintList:
            if self.debug: print("ERROR: newRestraintSet, no filName and no constraintList defined")
            return    
        
        if termType in allowedTerms:
            if termType in ['RDC','DANI']:
                termlist = [ i.termId for i in self.run.sortedHaddockEnergyTerms() if i.code == termType ]
                termlist.sort()
                if len(termlist) < 5:
                    if len(termlist): termId = termlist[-1] + 1
                    else: termId = 1
                    energyTerm = self.run.newHaddockEnergyTerm(code=termType,termId=termId)
                    if termType == 'RDC': addRdcParam(self.run,termId)
                    if termType == 'DANI': addDaniParam(self.run,termId)
                else: 
                    if self.debug: print('WARNING: newRestraintSet, Only 5 %s parameter sets allowed' % termType)
                    return
            else:
                termId = 1
                while self.run.findFirstHaddockEnergyTerm(code=termType,termId=termId): termId +=1
                energyTerm = self.run.newHaddockEnergyTerm(code=termType,termId=termId)
            
            if fileName:
                energyTerm.fileName = fileName
                energyTerm.constraintList = None
            elif constraintList:
                energyTerm.fileName = None
                energyTerm.constraintList = constraintList
            else: pass
            
            return energyTerm    
        else:
            if self.debug: print("ERROR: newRestraintSet, termType %s not allowed, choose from %s" % 
                                 termType," ".join(allowedTypes))
        
        
