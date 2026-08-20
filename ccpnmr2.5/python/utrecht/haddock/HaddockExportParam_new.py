#!/usr/bin/env/ python

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

import     sys, os
from    HaddockBasic    import getPdbString, getAirSegments, getFlexibleResidues, makeBackup
from     HaddockLocal    import rdcProtocolStore, daniProtocolStore

class exportParam:
    
    """Description:    Export a 'new' type HADDOCK project. The project is exported as a self containing parameter file.
                       It can be uploaded to the HADDOCK server.
       Input:        Haddock project instance, CCPN project instance.
       Output:        Haddock project parameter file
    """
    
    def __init__(self,haddockproject,ccpnproject,semiFlexMode):
        
        self.haddockproject = haddockproject
        self.ccpnproject = ccpnproject
        self.workingDir = self.haddockproject.workingDir
        self.partners = self.haddockproject.sortedHaddockPartners()
        self.latestRun = self.haddockproject.sortedRuns()[-1]
        self.HaddockRunParameters = {}
        
        self.parameters = {'analysisClustRmsd':'clust_rmsd','analysisClustSize':'clust_size','analysisDistHBond':'dist_hb',
                           'analysisDistNonbond':'dist_nb','calcDesolvation':'calcdesolv','centerOfMassConstant':'kcont',
                           'centerOfMassRestraints':'cmrest','dielectricType':'dielec','doIncludeDihEnergy':'dihedflag',
                           'doRigidBodyElectrostatics':'elecflag_0','doRigidTranslations':'rigidtrans',
                           'doSAElectrostatics':'elecflag_1','doWaterDock':'waterdock','epsilon':'epsilon',
                           'initialRigidBodyMinim':'rigidmini','nTrails':'ntrails',
                           'nonBondedType':'par_nonbonded','numAnalysisStructures':'anastruc_1','numIt0Structures':'structures_0',
                           'numIt1Structures':'structures_1','numWrefStructures':'waterrefine','radomizeStartOriention':'randorien',
                           'randomAmbigRestraints':'ranair','randomExclParts':'ncvpart','randomExcludeAir':'noecv',
                           'randomSeed':'iniseed','removeNonPolarH':'delenph','rigidbodyIMinteractScaling':'inter_rigid',
                           'rotate180It0':'rotate180_0','rotate180It1':'rotate180_1','skipStructures':'skip_struc',
                           'solvent':'solvent','surfaceContactConstant':'ksurf','surfaceContactRestraints':'surfrest',
                           'useDnaRestraints':'dnarest_on','useDbSolvateMethod':'solvate_method',
                           'waterInitRestCutoff':'water_restraint_initial','waterRestCutoff':'water_restraint_cutoff',
                           'waterRestScale':'water_restraint_scale','waterToKeep':'water_tokeep','waterToAddRandom':'water_randfrac',        
                           'waterSurfaceCutoff':'water_surfcutoff','doWaterAnalysis':'water_analysis','doRigidBodyWaterTrans':'transwater',
                           'numInitWaterShells':'waterensemble'}        
        
        filename = os.path.join(self.workingDir,self.haddockproject.name+'-run'+str(self.latestRun.serial)+'.web')
        print("** Export Haddock project parameter file %s **" % filename)
        makeBackup(filename)
        
        self.writePartnerData(semiFlexMode)
        self.writeRunParameters(self.HaddockRunParameters)
        self.writeScoringWeights(self.HaddockRunParameters)
        self.writeSymmetryData(self.HaddockRunParameters)
        self.writeEnergyProtocolStores(self.HaddockRunParameters)
        
        self.writeToFile(filename)
    
        print("** Export complete **")
        

    def writePartnerData(self,semiFlexMode):
    
        """Write partner specific data"""
        
        ambigstring = None
        constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'AMBIG']
        if constraint:
            for term in constraint: 
                if term.fileName and os.path.isfile(term.fileName):
                    ambigstring = ""
                    ambigfile = file(term.fileName,'r')
                    for line in ambigfile.readlines(): ambigstring += line
                    ambigfile.close()
        
        partnercount = 1
        for partner in self.haddockproject.sortedHaddockPartners():
            self.HaddockRunParameters['p%i' % partnercount] = {}
            self.writeStructureData(partner, self.HaddockRunParameters['p%i' % partnercount])
            if not ambigstring: self.writeActivePassive(partner, self.HaddockRunParameters['p%i' % partnercount])
            self.writeFlexible(partner,semiFlexMode, self.HaddockRunParameters['p%i' % partnercount])
            self.writeHistidinePstate(partner, self.HaddockRunParameters['p%i' % partnercount])
            self.writeForceField(partner, self.HaddockRunParameters['p%i' % partnercount])
            partnercount += 1
    
        if ambigstring: self.HaddockRunParameters['ambigdata'] = ambigstring
    
    def writeStructureData(self,partner,parameters):
        
        """Write structure data for partner"""
    
        parameters['pdb'] = {}
        parameters['pdb']['mode'] = 'submit'
        parameters['pdb']['chain'] = partner.code
        
        molSystem = partner.molSystem
        eId       = partner.structureEnsemble.ensembleId
        ensembles = self.ccpnproject.sortedStructureEnsembles()
        if molSystem: ensembles = [e for e in ensembles if e.molSystem is molSystem]
        if not eId <= len(ensembles): eId = len(ensembles)
        chains = [hc.chain.pdbOneLetterCode.strip() for hc in partner.chains]
        
        pdbstring = ''
        ens = range(0,eId)
        if len(ens) > 1:
            modelcount = 1
            for ensemble in ens:
                if ensemble == ens[-1]:
                    pdbstring += ('MODEL %i' % modelcount)    
                    pdbstring += getPdbString(ensembles[ensemble],chains,chainRename=partner.code,chain=True)
                    pdbstring += 'ENDMDL\nEND\n'
                else:
                    pdbstring += ('MODEL %i\n' % modelcount)    
                    pdbstring += getPdbString(ensembles[ensemble],chains,chainRename=partner.code,chain=True)
                    pdbstring += 'ENDMDL\n'
                    modelcount += 1
        else: pdbstring = getPdbString(ensembles[0],chains,chainRename=partner.code)    
        
        parameters['pdb']['pdbdata'] = pdbstring
    
    def writeFlexible(self,partner,semiFlexMode,parameters):
        
        """Write flexible and semi-flexible residues for the HADDOCK partner."""
        
        flex = getFlexibleResidues(partner)
    
        parameters['semiflex'] = {}
        parameters['semiflex']['mode'] = semiFlexMode[partner.code]
        parameters['semiflex']['segments'] = flex['semi']
        
        parameters['fullyflex'] = {}
        parameters['fullyflex']['segments'] = flex['full']
    
    def writeActivePassive(self,partner,parameters):
        
        """Write active and passive residues for the HADDOCK partner"""
        
        restraints = getAirSegments(partner)
        
        parameters['r'] = {}
        parameters['r']['activereslist'] = restraints['active']
        parameters['r']['passivereslist'] = restraints['passive']
    
    def writeHistidinePstate(self,partner,parameters):
    
        parameters['auto_his'] = partner.autoHistidinePstate
        if partner.autoHistidinePstate == False:
            pstates = []
            molSystem = partner.molSystem
            ensembles = self.ccpnproject.sortedStructureEnsembles()
            if molSystem: ensembles = [e for e in ensembles if e.molSystem is molSystem]
            for chain in ensembles[0].sortedCoordChains():
                for residue in chain.sortedResidues():
                    if residue.residue.ccpCode == 'His':
                        if residue.residue.descriptor == "prot:HD1;deprot:HE2": 
                            pstates.append({residue.residue.seqCode:'HISD'})
                        elif residue.residue.descriptor == "deprot:HD1;prot:HE2":     
                            pstates.append({residue.residue.seqCode:'HISE'})
                        elif residue.residue.descriptor == "prot:HD1;prot:HE2":     
                            pstates.append({residue.residue.seqCode:'HIS+'})
                        else: pass
            
            parameters['his'] = pstates
    
    def writeForceField(self,partner,parameters):
    
        parameters['forcefield'] = partner.forceFieldCode
        if partner.forceFieldCode == 'TOPALLHDG': parameters['moleculetype'] = 'Protein'
        else: parameters['moleculetype'] = partner.forceFieldCode
        
        parameters['dna'] = partner.isDna
        parameters['segid'] = partner.code
    
    def writeScoringWeights(self,parameters):
    
        """Write the scoring weights for the different docking stages"""
        
        scoringWeights = [(sw.term, sw.stage, sw.value) for sw in self.latestRun.scoringWeights]
        scoringWeights.sort()
        scoreTerm = ''
        for term, stage, scoringWeight in scoringWeights:
            if term == scoreTerm:
                parameters['w_%s' % term].append(scoringWeight)
            else:
                scoreTerm = term
                parameters['w_%s' % term] = [scoringWeight]    

    def writeSymmetryData(self,parameters):
        
        """Write symmetry data, type of symmetry, residues and chains involved"""                
    
        symdict = {'ncs':[],'C2':[],'C3':[],'C4':[],'C5':[],'C6':[]}; symw = False
    
        for symmetry in self.latestRun.sortedSymmetryRestraints():
            symdict[symmetry.symmetryCode].append((1,symmetry.segmentLength,'A'))
        
        for symmetry in symdict:
            if len(symdict[symmetry]) > 0:
                if symmetry == 'ncs': 
                    parameters['ncs'] = {}
                    parameters['ncs']['on'] = True
                    parameters['ncs']['constant'] = self.latestRun.get('ncsRestraintConstant')
                    parameters['ncs']['segments'] = []
                    for symrange in symdict[symmetry]: parameters['ncs']['segments'].append(symrange)    
                elif symw == False: 
                    parameters['symmetry'] = {}
                    parameters['symmetry']['on'] = True
                    parameters['symmetry']['constant'] = self.latestRun.get('symmetryRestraintConstant')
                    for symrange in symdict[symmetry]: parameters['ncs']['segments'].append(symrange)
                    symw = True
                else: pass
                

    def writeEnergyProtocolStores(self,parameters):
        
        """Write all protocols that define energy constanst in different stages"""
        
        protocolStores = ['distRestraintEnergyStore','dihRestraintEnergyStore',
                           'semiflexInterMolScalingStore','autoDistanceRestraintWeightStore']
        
        protocolStoreDict = {}
        
        for protocolStore in protocolStores:
            protocol = self.latestRun.findFirstHaddockEnergyTerm(code=protocolStore)
            if protocol:
                for term in protocol.sortedEnergyTermParameters(): protocolStoreDict[term.code] = term.value                

        'Hydrogen bond restraints'
        if self.latestRun.get('useHBondRestraints') == True:
            constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'HBOND']
            if constraint:
                hbondstring = None
                for term in constraint:
                    if term.fileName and os.path.isfile(term.fileName):
                        hbondstring = ""
                        hbondfile = file(term.fileName,'r')
                        for line in hbondfile.readlines(): hbondstring += line
                        hbondfile.close()
            if hbondstring:        
                parameters['hbonds_on'] = self.latestRun.get('useHBondRestraints')
                parameters['hbonddata'] = hbondstring
        
        'Distance restraints energy constants'
        constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'UNAMBIG']
        if constraint:
            unambigstring = None
            for term in constraint:
                if term.fileName and os.path.isfile(term.fileName):
                    unambigstring = ""
                    unambigfile = file(term.fileName,'r')
                    for line in unambigfile.readlines(): unambigstring += line
                    unambigfile.close()
            if unambigstring:        
                parameters['unambigdata'] = unambigstring
        
        for stage in ['unamb','amb','hbond']:
            parameters[stage] = {}
            parameters[stage]['firstit'] = protocolStoreDict[stage+'_firstit']
            parameters[stage]['lastit'] = protocolStoreDict[stage+'_lastit']
            parameters[stage]['stages'] = {}
            parameters[stage]['stages']['hot'] = protocolStoreDict[stage+'_hot']
            parameters[stage]['stages']['cool1'] = protocolStoreDict[stage+'_cool1']
            parameters[stage]['stages']['cool2'] = protocolStoreDict[stage+'_cool2']
            parameters[stage]['stages']['cool3'] = protocolStoreDict[stage+'_cool3']
    
        'Automated distance restraints weighting'
        parameters['air_scaling'] = self.latestRun.get('doAirScaling')
        parameters['tot_unamb'] = self.latestRun.get('numUnambRestautoAir')
        parameters['tot_amb'] = self.latestRun.get('numAmbRestautoAir')
        for stage in ['mrswi','rswi','masy','asy']:
            parameters['stage'] = {}
            parameters['stage']['hot'] = protocolStoreDict[stage+'_hot']
            parameters['stage']['cool1'] = protocolStoreDict[stage+'_cool1']
            parameters['stage']['cool2'] = protocolStoreDict[stage+'_cool2']
            parameters['stage']['cool3'] = protocolStoreDict[stage+'_cool3']
        
        'Dihedral restraint energy constants'
        constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'DIHEDRAL']
        if constraint:
            dihedralstring = None
            for term in constraint:
                if term.fileName and os.path.isfile(term.fileName):
                    dihedralstring = ""
                    dihedralfile = file(term.fileName,'r')
                    for line in dihedralfile.readlines(): dihedralstring += line
                    dihedralfile.close()
            if dihedralstring:        
                parameters['dihedrals_on'] = True
                parameters['dihedraldata'] = dihedralstring
                parameters['stages'] = {}
                parameters['stages']['hot'] = protocolStoreDict['dihedrals_hot']
                parameters['stages']['cool1'] = protocolStoreDict['dihedrals_cool1']
                parameters['stages']['cool2'] = protocolStoreDict['dihedrals_cool2']
                parameters['stages']['cool3'] = protocolStoreDict['dihedrals_cool3']
            
        'RDC constraints'
        rdcs = 1
        constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'RDC']
        for term in constraint:
            energyTermStore = self.latestRun.findFirstHaddockEnergyTerm(code='rdcProtocolStore',termId=term.termId)
            if energyTermStore and term.fileName and os.path.exists(term.fileName):
                storedict = {}
                for energyTerm in energyTermStore.sortedEnergyTermParameters(): storedict[energyTerm.code] = energyTerm.value
                self.writeRdcProtocol(parameters,storedict,term.termId,rdcfile=term.fileName)
                rdcs += 1
            else:
                print("RDC Warning: RDC energyTerm %i defined but no RDC CNS file associated" % term.termId)
        
        while rdcs < 6:
            self.writeRdcProtocol(parameters,rdcProtocolStore['terms'],rdcs,disable=True)
            rdcs += 1
        
        'DANI constraints'
        danis = 1
        constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'DANI']
        for term in constraint:
            energyTermStore = self.latestRun.findFirstHaddockEnergyTerm(code='daniProtocolStore',termId=term.termId)
            if energyTermStore and term.fileName and os.path.exists(term.fileName):
                storedict = {}
                for energyTerm in energyTermStore.sortedEnergyTermParameters(): storedict[energyTerm.code] = energyTerm.value
                self.writeDaniProtocol(parameters,storedict,term.termId,danifile=term.fileName)
                danis += 1
            else:
                print("DANI Warning: DANI energyTerm %i defined but no DANI CNS file associated" % term.termId)

        while danis < 6:
            self.writeDaniProtocol(parameters,daniProtocolStore['terms'],danis,disable=True)
            danis += 1        
        
        'Scaling of intermolecular interactions for semi-flexible SA'
        parameters['scaling_init'] = {}
        parameters['scaling_init']['rigid'] = protocolStoreDict['init_rigid']
        parameters['scaling_init']['cool2'] = protocolStoreDict['init_cool2']
        parameters['scaling_init']['cool3'] = protocolStoreDict['init_cool3']
    
        parameters['scaling_fin'] = {}
        parameters['scaling_fin']['rigid'] = protocolStoreDict['fin_rigid']
        parameters['scaling_fin']['cool2'] = protocolStoreDict['fin_cool2']
        parameters['scaling_fin']['cool3'] = protocolStoreDict['fin_cool3']
    
        'Docking annealing protocol'
        protocol = self.latestRun.findFirstHaddockEnergyTerm(code='dockingProtocolStore')
        for term in protocol.sortedEnergyTermParameters(): 
            parameters['%s' % term.code] = term.value    

    def writeDaniProtocol(self,parameters,storedict,danis,danifile=None,disable=False):        

            parameters['dan%i' % danis] = {}
            if disable == True: parameters['dan%i' % danis]['choice'] = 'NO'
            else: parameters['dan%i' % danis]['choice'] = 'DANI'
            parameters['dan%i' % danis]['constants'] = {}
            parameters['dan%i' % danis]['constants']['firstit'] = storedict['firstIt']
            parameters['dan%i' % danis]['constants']['lastit'] = storedict['lastIt']
            parameters['dan%i' % danis]['constants']['stages'] = {}
            parameters['dan%i' % danis]['constants']['stages']['hot'] = storedict['hot']
            for n in ['1','2','3']: 
                parameters['dan%i' % danis]['constants']['stages']['cool%s' % n] = storedict['cool%s' % n]
            parameters['dan%i' % danis]['tc'] = storedict['tc']
            parameters['dan%i' % danis]['anis'] = storedict['anis']
            parameters['dan%i' % danis]['r'] = storedict['r']
            parameters['dan%i' % danis]['wh'] = storedict['wh']
            parameters['dan%i' % danis]['wn'] = storedict['wn']
            
            if danifile:
                danistring = ""
                openfile = file(danifile,'r')
                for line in openfile.readlines(): danistring += line
                openfile.close()
            
                parameters['dan%i' % danis]['danidata'] = danistring

    def writeRdcProtocol(self,parameters,storedict,rdcs,rdcfile=None,disable=False):    

            parameters['rdc%i' % rdcs] = {}
            if disable == True:
                parameters['rdc%i' % rdcs]['choice'] = 'NO'
            else:    
                parameters['rdc%i' % rdcs]['choice'] = ['NO','SANI','VANGLE'][int(storedict['rdc_choice'])]
            parameters['rdc%i' % rdcs]['r'] = storedict['rdc_r']
            parameters['rdc%i' % rdcs]['d'] = storedict['rdc_d']
            parameters['rdc%i' % rdcs]['constants'] = {}
            parameters['rdc%i' % rdcs]['constants']['firstit'] = storedict['rdc_firstIt']
            parameters['rdc%i' % rdcs]['constants']['lastit'] = storedict['rdc_lastIt']
            parameters['rdc%i' % rdcs]['constants']['stages'] = {}
            parameters['rdc%i' % rdcs]['constants']['stages']['hot'] = storedict['rdc_hot']
            for n in ['1','2','3']: 
                parameters['rdc%i' % rdcs]['constants']['stages']['cool%s' % n] = storedict['rdc_cool%s' % n]
            parameters['rdc%i' % rdcs]['constants']['stages']['ini_bor'] = {}
            parameters['rdc%i' % rdcs]['constants']['stages']['ini_bor']['hot'] = storedict['ini_bor_hot']
            for n in ['1','2','3']: 
                parameters['rdc%i' % rdcs]['constants']['stages']['ini_bor']['cool%s' % n] = storedict['ini_bor_cool%s' % n]
            parameters['rdc%i' % rdcs]['constants']['stages']['fin_bor'] = {}
            parameters['rdc%i' % rdcs]['constants']['stages']['fin_bor']['hot'] = storedict['fin_bor_hot']
            for n in ['1','2','3']: 
                parameters['rdc%i' % rdcs]['constants']['stages']['fin_bor']['cool%s' % n] = storedict['fin_bor_cool%s' % n]
            parameters['rdc%i' % rdcs]['constants']['stages']['ini_cen'] = {}
            parameters['rdc%i' % rdcs]['constants']['stages']['ini_cen']['hot'] = storedict['ini_cen_hot']
            for n in ['1','2','3']: 
                parameters['rdc%i' % rdcs]['constants']['stages']['ini_cen']['cool%s' % n] = storedict['ini_cen_cool%s' % n]
            parameters['rdc%i' % rdcs]['constants']['stages']['fin_cen'] = {}
            parameters['rdc%i' % rdcs]['constants']['stages']['fin_cen']['hot'] = storedict['fin_cen_hot']
            for n in ['1','2','3']: 
                parameters['rdc%i' % rdcs]['constants']['stages']['fin_cen']['cool%s' % n] = storedict['fin_cen_cool%s' % n]

            if rdcfile:
                rdcstring = ""
                openfile = file(rdcfile,'r')
                for line in openfile.readlines(): rdcstring += line
                openfile.close()

                parameters['rdc%i' % rdcs]['rdcdata'] = rdcstring

    def writeRunParameters(self,parameters):
        
        """Write all general docking parameters"""
        
        parameters['runname'] = 'run%i' % self.latestRun.serial
        
        for parameter in self.parameters:
            parameters['%s' % self.parameters[parameter]] = self.latestRun.get(parameter)
    
    def writeToFile(self,filename):
        
        params = file(filename,'w')
        params.write("HaddockRunParameters = %s" % repr(self.HaddockRunParameters))
        params.close()        
