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

import     sys, os
from    HaddockBasic    import getPdbString, getAirSegments, getFlexibleResidues, makeBackup
from     HaddockLocal    import rdcProtocolStore, daniProtocolStore

class exportParam:
    
    """Description:    Export a 'new' type HADDOCK project. The project is exported as a self containing parameter file.
                       It can be uploaded to the HADDOCK server.
       Input:        Haddock project instance, CCPN project instance.
       Output:        Haddock project parameter file
    """
    
    def __init__(self,hProject=None,latestRun=None,ccpnProject=None):
        
        self.haddockproject = hProject
        self.latestRun = latestRun
        self.ccpnproject = ccpnProject
        self.workingDir = self.haddockproject.workingDir
        self.partners = self.haddockproject.sortedHaddockPartners()
        self.filestring = ""
        self.identlevel = 0
        
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
        
        
        print("** Export Haddock project parameter file **")
        
        if len(self.partners):
            self.writeHeader()
            self.writePartnerData()
            self.writeRunParameters()
            self.writeScoringWeights()
            self.writeSymmetryData()
            self.writeEnergyProtocolStores()
            self.writeFooter()
        else: 
            print("ERROR, No partners defined")
            return
        
        print("** Export complete **")
    
    def writeToFile(self):
        
        """Write the parameter file string to file"""
        
        if len(self.filestring):    
            self.filename = os.path.join(self.workingDir,self.haddockproject.name+'-run'+str(self.latestRun.serial)+'.web')
            print("Export parameter to file %s" % self.filename)
            makeBackup(self.filename)
        
            self.file = file(self.filename,'w')
            self.file.write(self.filestring)
            self.file.close()
        
    def writeHeader(self):
    
        """Write header line"""
    
        self.filestring += ("HaddockRunParameters (\n")
    
    def writeFooter(self):
        
        """Write footer line"""
        
        self.identlevel -= 1
        self.filestring += ("%s)\n" % (" "*self.identlevel))
    
    def writePartnerData(self):
    
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
        self.identlevel += 1
        for partner in self.haddockproject.sortedHaddockPartners():
            self.filestring += ("%sp%i = HaddockPartnerParameters (\n" % (" "*self.identlevel,partnercount))
            self.writeStructureData(partner)
            if ambigstring: self.writeActivePassive(partner,deactivate=True)
            else: self.writeActivePassive(partner)
            self.writeFlexible(partner)
            self.writeHistidinePstate(partner)
            self.writeForceField(partner)
            self.identlevel -= 1
            self.filestring += ("%s),\n" % (" "*self.identlevel))
            partnercount += 1
    
        if ambigstring: self.filestring += ("%stbldata = %s,\n" % (" "*self.identlevel,repr(ambigstring)))
    
    def writeStructureData(self,partner):
        
        """Write structure data for partner"""
    
        self.identlevel += 1
        self.filestring += ("%spdb = PDBData (\n" % (" "*self.identlevel))
        self.identlevel += 1
        self.filestring += ("%smode = 'submit',\n" % (" "*self.identlevel))
        self.filestring += ("%schain = %s,\n" % (" "*self.identlevel,repr(partner.code)))
        
        molSystem = partner.molSystem
        modelId   = int(partner.structureEnsemble.getDetails())
        
        ensembles = self.ccpnproject.sortedStructureEnsembles()
        if molSystem: ensembles = [e for e in ensembles if e.molSystem is molSystem]
        models = []
        for e in ensembles: models += [model for model in e.sortedModels()]
        
        chains = [hc.chain.pdbOneLetterCode.strip() for hc in partner.chains]
    
        pdbstring = ''
        if modelId > 1:
            modelcount = 1
            for model in models[0:modelId]:
                if modelcount == modelId:
                    pdbstring += ('MODEL %i\n' % modelcount)    
                    pdbstring += getPdbString(model,
                                              [ch.chain for ch in partner.sortedChains()],
                                              chainRename=partner.code,
                                              fileEnd='ENDMDL\nEND\n') 
                else:
                    pdbstring += ('MODEL %i\n' % modelcount)    
                    pdbstring += getPdbString(model,
                                              [ch.chain for ch in partner.sortedChains()],
                                              chainRename=partner.code,
                                              fileEnd='ENDMDL')
                    modelcount += 1
        else: 
            pdbstring = getPdbString(models[0],[ch.chain for ch in partner.sortedChains()],chainRename=partner.code)
        
        self.filestring += ("%spdbdata = %s,\n" % (" "*self.identlevel,repr(pdbstring)))
        self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
    
    def writeFlexible(self,partner):
        
        """Write flexible and semi-flexible residues for the HADDOCK partner."""
        
        flex = getFlexibleResidues(partner)
        
        self.filestring += ("%ssemiflex = SemiflexSegmentList (\n" % (" "*self.identlevel))
        self.identlevel += 1
        self.filestring += ("%ssegments = RangeArray (\n" % (" "*self.identlevel))
        if partner.semiFlexMode == 'manual':
            self.identlevel += 1
            if not len(flex['semi']) == 0:
                for semiflex in flex['semi']:
                    self.filestring += ("%s%s,\n" % (" "*self.identlevel,repr(semiflex)))
            self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
        self.filestring += ("%smode = %s,\n" % (" "*self.identlevel,repr(partner.semiFlexMode)))
        self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
        self.filestring += ("%sfullyflex = SegmentList (\n" % (" "*self.identlevel))
        self.identlevel += 1
        self.filestring += ("%ssegments = RangeArray (\n" % (" "*self.identlevel))
        self.identlevel += 1
        if not len(flex['full']) == 0:
            for fullyflex in flex['full']:
                self.filestring += ("%s%s,\n" % (" "*self.identlevel,repr(fullyflex)))
        self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
        self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
    
    def writeActivePassive(self,partner,deactivate=False):
        
        """Write active and passive residues for the HADDOCK partner"""
        
        restraints = getAirSegments(partner)
        if deactivate == True: restraints['active'] = []; restraints['passive'] = []
        
        self.filestring += ("%sr = RestraintsInterface (\n" % (" "*self.identlevel))
        self.identlevel += 1
        self.filestring += ("%sactivereslist = IntegerArray (\n" % (" "*self.identlevel))
        self.identlevel += 1
        for activeres in restraints['active']:
            self.filestring += ("%s%i,\n" % (" "*self.identlevel,activeres))
        self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
        self.filestring += ("%spassivereslist = IntegerArray (\n" % (" "*self.identlevel))
        self.identlevel += 1
        for passiveres in restraints['passive']:
            self.filestring += ("%s%i,\n" % (" "*self.identlevel,passiveres))
        self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
        self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
    
    def writeHistidinePstate(self,partner):
    
        self.filestring += ("%sauto_his = %s,\n" % (" "*self.identlevel,partner.autoHistidinePstate))
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
            
            if len(pstates):
                self.filestring += ("%shis = HistidineStateArray (\n" % (" "*self.identlevel)) 
                self.identlevel += 1
                for pstate in pstates:
                    self.filestring +=("%s," % repr(pstate))
                self.identlevel -= 1    
                self.filestring += ("%s),\n" % (" "*self.identlevel))
        else:
            self.filestring += ("%shis = HistidineStateArray (\n" % (" "*self.identlevel)) 
            self.filestring += ("%s),\n" % (" "*self.identlevel))        
    
    def writeForceField(self,partner):
    
        self.filestring += ("%sforcefield = %s,\n" % (" "*self.identlevel,repr(partner.forceFieldCode)))
        if partner.forceFieldCode == 'TOPALLHDG':
            self.filestring += ("%smoleculetype = 'Protein',\n" % (" "*self.identlevel))
        else:
            self.filestring += ("%smoleculetype = %s,\n" % (" "*self.identlevel,repr(partner.forceFieldCode)))    
        self.filestring += ("%sdna = %s,\n" % (" "*self.identlevel,partner.isDna))
        self.filestring += ("%ssegid = %s,\n" % (" "*self.identlevel,repr(partner.code)))
    
    def writeScoringWeights(self):
    
        """Write the scoring weights for the different docking stages"""
        
        scoringWeights = [(sw.term, sw.stage, sw.value) for sw in self.latestRun.scoringWeights]
        scoringWeights.sort()
        scoreTerm = ''
        for term, stage, scoringWeight in scoringWeights:
            if term == scoreTerm:
                self.filestring += ("%s%1.3f,\n" % (" "*self.identlevel,scoringWeight))
            elif scoreTerm == '':
                self.filestring += ("%sw_%s = FloatArray (\n" % (" "*self.identlevel,term))
                self.identlevel += 1
                self.filestring += ("%s%1.3f,\n" % (" "*self.identlevel,scoringWeight))
                scoreTerm = term
            else:
                self.identlevel -= 1
                self.filestring += ("%s),\n" % " "*self.identlevel)
                self.filestring += ("%sw_%s = FloatArray (\n" % (" "*self.identlevel,term))
                self.identlevel += 1
                self.filestring += ("%s%1.3f,\n" % (" "*self.identlevel,scoringWeight))
                scoreTerm = term
        
        self.identlevel -= 1
        self.filestring += ("),\n")        

    def writeSymmetryData(self):
        
        """Write symmetry data, type of symmetry, residues and chains involved"""                
    
        symdict = {'ncs':[],'C2':[],'C3':[],'C4':[],'C5':[],'C6':[]}; symw = False
    
        for symmetry in self.latestRun.sortedSymmetryRestraints():
            symdict[symmetry.symmetryCode].append((1,symmetry.segmentLength,'A'))
        
        for symmetry in symdict:
            if len(symdict[symmetry]) > 0:
                if symmetry == 'ncs': 
                    self.filestring += ("%sncs = SymmetrySpecification (\n" % (" "*self.identlevel))
                    self.identlevel += 1
                    self.filestring += ("%son = True,\n" % (" "*self.identlevel))
                    self.filestring += ("%sconstant = %1.1f,\n" % (" "*self.identlevel,self.latestRun.get('ncsRestraintConstant')))
                    self.identlevel -= 1
                    self.filestring += ("%s),\n" % " "*self.identlevel)
                elif symw == False: 
                    self.filestring += ("%ssymmetry = SymmetrySpecification (\n" % (" "*self.identlevel))
                    self.identlevel += 1
                    self.filestring += ("%son = True,\n" % (" "*self.identlevel))
                    self.filestring += ("%sconstant = %1.1f,\n" % (" "*self.identlevel,self.latestRun.get('symmetryRestraintConstant')))
                    self.identlevel -= 1
                    self.filestring += ("%s),\n" % " "*self.identlevel)
                    symw = True
                else: pass
                
                for symrange in symdict[symmetry]:            
                    self.filestring += ("%s%ssegments = LabeldRangePairArray (\n" % (" "*self.identlevel,symmetry))
                    self.identlevel += 1
                    self.filestring += ("%sLabeledRangePair (\n" % (" "*self.identlevel))
                    self.identlevel += 1
                    self.filestring += ("%sr = LabeledRangeArray (\n" % (" "*self.identlevel))
                    self.identlevel += 1
                    self.filestring += ("%sLabeledRange (\n" % (" "*self.identlevel))
                    self.identlevel += 1
                    self.filestring += ("%sstart = %i,\n" % (" "*self.identlevel,symrange[0]))
                    self.filestring += ("%send = %i,\n" % (" "*self.identlevel,symrange[1]))
                    self.filestring += ("%schain = %s,\n" % (" "*self.identlevel,symrange[2]))    
                    self.identlevel -= 1
                    self.filestring += ("%s),\n" % (" "*self.identlevel))
                    self.identlevel -= 1
                    self.filestring += ("%s),\n" % (" "*self.identlevel))
                    self.identlevel -= 1
                    self.filestring += ("%s),\n" % (" "*self.identlevel))
                    self.identlevel -= 1
                    self.filestring += ("%s),\n" % (" "*self.identlevel))

    def writeEnergyProtocolStores(self):
        
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
                self.filestring += ("%shbonds_on = %s,\n" % (" "*self.identlevel,repr(self.latestRun.get('useHBondRestraints'))))
                self.filestring += ("%shbonddata = %s,\n" % (" "*self.identlevel,repr(hbondstring)))
        
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
                self.filestring += ("%sunambigtbldata = %s,\n" % (" "*self.identlevel,repr(unambigstring)))
        
        for stage in ['unamb','amb','hbond']:
            self.filestring += ("%s%s =  ExtStageConstants (\n" % (" "*self.identlevel,stage))
            self.identlevel += 1
            self.filestring += ("%sfirstit = %i,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_firstit']))
            self.filestring += ("%slastit = %i,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_lastit']))
            self.filestring += ("%sstages = StageConstants (\n" % (" "*self.identlevel))
            self.identlevel += 1
            self.filestring += ("%shot = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_hot']))
            self.filestring += ("%scool1 = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_cool1']))
            self.filestring += ("%scool2 = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_cool2']))
            self.filestring += ("%scool3 = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_cool3']))
            self.identlevel -= 1
            self.filestring += ("%s),\n" % (" "*self.identlevel))
            self.identlevel -= 1
            self.filestring += ("%s),\n" % (" "*self.identlevel))
            
        'Automated distance restraints weighting'
        self.filestring += ("%sair_scaling = %s,\n" % (" "*self.identlevel,repr(self.latestRun.get('doAirScaling'))))
        self.filestring += ("%stot_unamb = %i,\n" % (" "*self.identlevel,self.latestRun.get('numUnambRestautoAir')))
        self.filestring += ("%stot_amb = %i,\n" % (" "*self.identlevel,self.latestRun.get('numAmbRestautoAir')))
        for stage in ['mrswi','rswi','masy','asy']:
            self.filestring += ("%s%s = StageConstants (\n" % (" "*self.identlevel,stage))
            self.identlevel += 1
            self.filestring += ("%shot = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_hot']))
            self.filestring += ("%scool1 = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_cool1']))
            self.filestring += ("%scool2 = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_cool2']))
            self.filestring += ("%scool3 = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict[stage+'_cool3']))
            self.identlevel -= 1
            self.filestring += ("%s),\n" % (" "*self.identlevel))
        
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
                self.filestring += ("%sdihedrals_on = True,\n" % (" "*self.identlevel))
                self.filestring += ("%sdihedraldata = %s,\n" % (" "*self.identlevel,repr(dihedralstring)))
                self.filestring += ("%sstages = StageConstants (\n" % (" "*self.identlevel))
                self.identlevel += 1
                self.filestring += ("%shot = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict['dihedrals_hot']))
                self.filestring += ("%scool1 = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict['dihedrals_cool1']))
                self.filestring += ("%scool2 = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict['dihedrals_cool2']))
                self.filestring += ("%scool3 = %1.1f,\n" % (" "*self.identlevel,protocolStoreDict['dihedrals_cool3']))
                self.identlevel -= 1
                self.filestring += ("%s),\n" % (" "*self.identlevel))
        
        'RDC constraints'
        rdcs = 1
        constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'RDC']
        for term in constraint:
            energyTermStore = self.latestRun.findFirstHaddockEnergyTerm(code='rdcProtocolStore',termId=term.termId)
            if energyTermStore and term.fileName and os.path.exists(term.fileName):
                storedict = {}
                for energyTerm in energyTermStore.sortedEnergyTermParameters(): storedict[energyTerm.code] = energyTerm.value
                self.writeRdcProtocol(storedict,term.termId,rdcfile=term.fileName)
                rdcs += 1
            else:
                print("RDC Warning: RDC energyTerm %i defined but no RDC CNS file associated" % term.termId)
        
        while rdcs < 6:
            self.writeRdcProtocol(rdcProtocolStore['terms'],rdcs,disable=True)
            rdcs += 1
        
        'DANI constraints'
        danis = 1
        constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'DANI']
        for term in constraint:
            energyTermStore = self.latestRun.findFirstHaddockEnergyTerm(code='daniProtocolStore',termId=term.termId)
            if energyTermStore and term.fileName and os.path.exists(term.fileName):
                storedict = {}
                for energyTerm in energyTermStore.sortedEnergyTermParameters(): storedict[energyTerm.code] = energyTerm.value
                self.writeDaniProtocol(storedict,term.termId,danifile=term.fileName)
                danis += 1
            else:
                print("DANI Warning: DANI energyTerm %i defined but no DANI CNS file associated" % term.termId)
        
        while danis < 6:
            self.writeDaniProtocol(daniProtocolStore['terms'],danis,disable=True)
            danis += 1        
        
        'Scaling of intermolecular interactions for semi-flexible SA'
        self.filestring += ("%sscaling_init = ScalingConstants (\n" % (" "*self.identlevel))
        self.identlevel += 1
        self.filestring += ("%srigid = %1.3f,\n" % (" "*self.identlevel,protocolStoreDict['init_rigid']))
        self.filestring += ("%scool2 = %1.3f,\n" % (" "*self.identlevel,protocolStoreDict['init_cool2']))
        self.filestring += ("%scool3 = %1.3f,\n" % (" "*self.identlevel,protocolStoreDict['init_cool3']))
        self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
        
        self.filestring += ("%sscaling_fin = ScalingConstants (\n" % (" "*self.identlevel))
        self.identlevel += 1
        self.filestring += ("%srigid = %1.3f,\n" % (" "*self.identlevel,protocolStoreDict['fin_rigid']))
        self.filestring += ("%scool2 = %1.3f,\n" % (" "*self.identlevel,protocolStoreDict['fin_cool2']))
        self.filestring += ("%scool3 = %1.3f,\n" % (" "*self.identlevel,protocolStoreDict['fin_cool3']))
        self.identlevel -= 1
        self.filestring += ("%s),\n" % (" "*self.identlevel))
        
        'Docking annealing protocol'
        protocol = self.latestRun.findFirstHaddockEnergyTerm(code='dockingProtocolStore')
        for term in protocol.sortedEnergyTermParameters(): 
            self.filestring += ("%s%s = %1.1f,\n" % (" "*self.identlevel,term.code,term.value))    
    
    def writeDaniProtocol(self,storedict,danis,danifile=None,disable=False):        

            self.filestring += ("%sdan%i = DANIParameters (\n" % (" "*self.identlevel,danis))
            self.identlevel += 1
            if disable == True: self.filestring += ("%schoice = 'NO',\n" % (" "*self.identlevel))
            else: self.filestring +=("%schoice = 'DANI',\n" % (" "*self.identlevel))
            self.filestring += ('%sconstants = ExtStageConstants (\n' % (" "*self.identlevel))
            self.identlevel += 1
            self.filestring += ('%sfirstit = %i,\n' % (" "*self.identlevel,storedict['firstIt']))
            self.filestring += ('%slastit = %i,\n' % (" "*self.identlevel,storedict['lastIt']))
            self.filestring += ('%sstages = StageConstants (\n' % (" "*self.identlevel))
            self.identlevel += 1
            self.filestring += ('%shot = %1.1f,\n' % (" "*self.identlevel,storedict['hot']))
            for n in ['1','2','3']: self.filestring += ('%scool%s = %1.1f,\n' % (" "*self.identlevel,n,storedict['cool%s' % n]))
            self.identlevel -= 1         
            self.filestring += ("%s),\n" % (" "*self.identlevel))
            self.identlevel -= 1 
            self.filestring += ("%s),\n" % (" "*self.identlevel))
            self.filestring += ('%stc = %1.1f,\n' % (" "*self.identlevel,storedict['tc']))
            self.filestring += ('%sanis = %1.1f,\n' % (" "*self.identlevel,storedict['anis']))
            self.filestring += ('%sr = %1.1f,\n' % (" "*self.identlevel,storedict['r']))
            self.filestring += ('%swh = %1.1f,\n' % (" "*self.identlevel,storedict['wh']))
            self.filestring += ('%swn = %1.1f,\n' % (" "*self.identlevel,storedict['wn']))
            
            if danifile:
                danistring = ""
                openfile = file(danifile,'r')
                for line in openfile.readlines(): danistring += line
                opnefile.close()
                
                self.filestring += ('%sdanidata = %s,\n' % (" "*self.identlevel,repr(danistring)))    
            
            self.identlevel -= 1 
            self.filestring += ("%s),\n" % (" "*self.identlevel))

    def writeRdcProtocol(self,storedict,rdcs,rdcfile=None,disable=False):    

            self.filestring += ("%srdc%i = RDCParameters (\n" % (" "*self.identlevel,rdcs))
            self.identlevel += 1
            if disable == True:
                self.filestring += ('%schoice = %s,\n' % (" "*self.identlevel,repr('NO')))    
            else:    
                self.filestring += ('%schoice = %s,\n' % (" "*self.identlevel,repr(['NO','SANI','VANGLE'][int(storedict['rdc_choice'])])))
            self.filestring += ('%sr = %1.1f,\n' % (" "*self.identlevel,storedict['rdc_r']))
            self.filestring += ('%sd = %1.1f,\n' % (" "*self.identlevel,storedict['rdc_d']))
            self.filestring += ('%sconstants = ExtStageConstants (\n' % (" "*self.identlevel))
            self.identlevel += 1
            self.filestring += ('%sfirstit = %i,\n' % (" "*self.identlevel,storedict['rdc_firstIt']))
            self.filestring += ('%slastit = %i,\n' % (" "*self.identlevel,storedict['rdc_lastIt']))
            self.filestring += ('%sstages = StageConstants (\n' % (" "*self.identlevel))
            self.identlevel += 1
            self.filestring += ('%shot = %1.1f,\n' % (" "*self.identlevel,storedict['rdc_hot']))
            for n in ['1','2','3']: self.filestring += ('%scool%s = %1.1f,\n' % (" "*self.identlevel,n,storedict['rdc_cool%s' % n]))
            self.identlevel -= 1         
            self.filestring += ("%s),\n" % (" "*self.identlevel))
            self.identlevel -= 1 
            self.filestring += ("%s),\n" % (" "*self.identlevel))
            self.filestring += ('%sini_bor = StageConstants (\n' % (" "*self.identlevel))
            self.identlevel += 1
            self.filestring += ('%shot = %1.1f,\n' % (" "*self.identlevel,storedict['ini_bor_hot']))
            for n in ['1','2','3']: self.filestring += ('%scool%s = %1.1f,\n' % (" "*self.identlevel,n,storedict['ini_bor_cool%s' % n]))
            self.identlevel -= 1 
            self.filestring += ("%s),\n" % (" "*self.identlevel))
            self.filestring += ('%sfin_bor = StageConstants (\n' % (" "*self.identlevel))
            self.identlevel += 1
            self.filestring += ('%shot = %1.1f,\n' % (" "*self.identlevel,storedict['fin_bor_hot']))
            for n in ['1','2','3']: self.filestring += ('%scool%s = %1.1f,\n' % (" "*self.identlevel,n,storedict['fin_bor_cool%s' % n]))
            self.identlevel -= 1 
            self.filestring += ("%s),\n" % (" "*self.identlevel))        
            self.filestring += ('%sini_cen = StageConstants (\n' % (" "*self.identlevel))
            self.identlevel += 1
            self.filestring += ('%shot = %1.1f,\n' % (" "*self.identlevel,storedict['ini_cen_hot']))
            for n in ['1','2','3']: self.filestring += ('%scool%s = %1.1f,\n' % (" "*self.identlevel,n,storedict['ini_cen_cool%s' % n]))
            self.identlevel -= 1 
            self.filestring += ("%s),\n" % (" "*self.identlevel))
            self.filestring += ('%sfin_cen = StageConstants (\n' % (" "*self.identlevel))
            self.identlevel += 1
            self.filestring += ('%shot = %1.1f,\n' % (" "*self.identlevel,storedict['fin_cen_hot']))
            for n in ['1','2','3']: self.filestring += ('%scool%s = %1.1f,\n' % (" "*self.identlevel,n,storedict['fin_cen_cool%s' % n]))
            self.identlevel -= 1 
            self.filestring += ("%s),\n" % (" "*self.identlevel))

            if rdcfile:
                rdcstring = ""
                openfile = file(rdcfile,'r')
                for line in openfile.readlines(): rdcstring += line
                openfile.close()

                self.filestring += ('%srdcdata = %s,\n' % (" "*self.identlevel,repr(rdcstring)))
            
            self.identlevel -= 1 
            self.filestring += ("%s),\n" % (" "*self.identlevel))
        
    def writeRunParameters(self):
        
        """Write all general docking parameters"""
        
        self.filestring += ("%srunname = 'run%i',\n" % (" "*self.identlevel,self.latestRun.serial))
        
        for parameter in self.parameters:
            if parameter == 'useDbSolvateMethod':
                if self.latestRun.get(parameter) == True:
                    self.filestring += ("%s%s = %s,\n" % (" "*self.identlevel,self.parameters[parameter],repr('db')))    
                else:    
                    self.filestring += ("%s%s = %s,\n" % (" "*self.identlevel,self.parameters[parameter],repr('pnt')))    
            else:
                self.filestring += ("%s%s = %s,\n" % (" "*self.identlevel,self.parameters[parameter],repr(self.latestRun.get(parameter))))    
