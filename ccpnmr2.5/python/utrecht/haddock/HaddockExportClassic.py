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

import sys
from     os.path             import join, isdir
from     os                     import makedirs
from    HaddockBasic        import getPdbString, getAirSegments, getFlexibleResidues, makeBackup
from     HaddockDnaRnaRest     import dnaRnaRestraints
from     HaddockLocal        import *

class exportClassic:
    
    """Description: Exports current project and run as a 'classical' HADDOCK style project. This means
                       a root directory bearing the projects name containing all PDB structure files that
                    need to be docked, a new.html file, restraint files, a ensemble.list file if multiple 
                    ensembles for a given structure are used.
       Input:        Haddock project instance, CCPN project instance.
       Output:        The various files as described above.                
    """
    
    def __init__(self,hProject=None,latestRun=None,ccpnProject=None):
        
        self.haddockproject = hProject
        self.latestRun = latestRun
        self.ccpnproject = ccpnProject
        
        self.workingDir = self.haddockproject.workingDir
        self.partners = self.haddockproject.sortedHaddockPartners()
        self.fileNames = {}
        
        if len(self.partners):
            self.__setupDirectoryStructure()
            self.__writeNewHtmlFile()
            self.__writeRunCnsFile()
            print("** Export complete **")
        else:
            print("-->ERROR: Export classic project. No partners defined")    
    
    def __setupDirectoryStructure(self):
        
        """Make a project root directory. If project is allready defined it will not be overwritten."""
        
        print("** Export 'Classic' Haddock project. Name %s, run %i **\n" % (self.haddockproject.name,self.latestRun.serial))
        
        self.projectRoot = join(self.workingDir,self.haddockproject.name)
        if isdir(self.projectRoot):
            print("Project %s allready has a root HADDOCK root directory within the set working directory" % self.haddockproject.name)
        else:
            print("Make project root directory %s within the set working directory" % self.haddockproject.name)
            makedirs(self.projectRoot)
    
    def __writeNewHtmlFile(self):
        
        """Make a new.html file based upon a HADDOCK project object"""
        
        print("Generate 'new.html' file")
        makeBackup(join(self.projectRoot,'new.html'))
        
        file = open(join(self.projectRoot,'new.html'),'w')
        file.write('<html>\n')
        file.write('<head>\n')
        file.write('<title>HADDOCK - start</title>\n')
        file.write('</head>\n')
        file.write('<body bgcolor=#ffffff>\n')
        file.write('<h2>Parameters for the start:</h2> \n')
        file.write('<BR>\n')
        file.write('<h4><!-- HADDOCK -->\n')
        
        ambigfile = False; allowedConstraintStores = ['UNAMBIG','AMBIG','DIHEDRAL','HBOND','RDC','DANI']
        danicount = 1; rdccount = 1
        for constraint in [ i for i in self.latestRun.sortedHaddockEnergyTerms() if i.code in allowedConstraintStores ]:
            if constraint.code == 'AMBIG': ambigfile = True
            if constraint.fileName: 
                if constraint.code == 'DANI':
                    file.write('%s%d_TBL=%s<BR>\n' % (constraint.code,danicount,constraint.fileName))
                    danicount += 1
                elif constraint.code == 'RDC':
                    file.write('%s%d_TBL=%s<BR>\n' % (constraint.code,rdccount,constraint.fileName))
                    rdccount += 1
                else: file.write('%s_TBL=%s<BR>\n' % (constraint.code,constraint.fileName))        
                    
        if ambigfile == False:
            self.__writeAmbigCnsTable()
            file.write('AMBIG_TBL=%s<BR>\n' % join(self.projectRoot,'ambig.tbl'))            
        
        file.write('HADDOCK_DIR=%s<BR>\n' % self.latestRun.haddockDir)
        file.write('N_COMP=%i<BR>\n' % len(self.partners))
        
        i = 1
        for partner in self.partners:
            pdbensemblelist = self.__writePdbFiles(partner)
            file.write('PDB_FILE%d=%s.pdb<BR>\n' % (i, join(self.projectRoot,self.fileNames[partner])))
            if len(pdbensemblelist) > 1:
                filelist = self.__writePdbFileList(partner,pdbensemblelist)
                file.write('PDB_LIST%d=%s<BR>\n' % (i, filelist))
            i += 1
        
        file.write('PROJECT_DIR=%s<BR>\n' % self.projectRoot)
        
        i = 1
        for partner in self.partners:
            file.write('PROT_SEGID_%d=%s<BR>\n' % (i, partner.code))
            i += 1
        
        file.write('RUN_NUMBER=%d<BR>\n' % (self.latestRun.serial))
        file.write('submit_save=Save updated parameters<BR>\n')
        file.write('</h4><!-- HADDOCK -->\n')
        file.write('</body>\n')
        file.write('</html> \n')
        
        file.close()
    
    def __writePdbFiles(self,partner):
        
        """Write PDB files into a given directory (using generic file names) given a list of MolStructrue objects.
           If multiple ensembles, each is exported as a seperate PDB file. The PDB name constitutes the MolSystem
           name, the selected chains in the given Haddock Partner and the ensemble number. If no ensembles the ensemble
           number is always 1, if no chains that the Haddock partner ID is used.
        """
        pdbmodellist = []
        molSystem = partner.molSystem
        modelId   = int(partner.structureEnsemble.getDetails())
        chains       = [hc.chain.pdbOneLetterCode.strip() for hc in partner.chains]
        
        if len("".join(chains)): chainstring = "".join(chains)
        else: chainstring = partner.code
        
        ensembles = self.ccpnproject.sortedStructureEnsembles()
        if molSystem: ensembles = [e for e in ensembles if e.molSystem is molSystem]
        models = []
        for e in ensembles: models += [model for model in e.sortedModels()]
            
        for model in models[0:modelId]:
            fileName = join(self.projectRoot,'%s-%s_%d.pdb' % (molSystem.code,chainstring,model.serial))
            if not self.fileNames.has_key(partner): self.fileNames[partner] = ('%s-%s_%d' % (molSystem.code,chainstring,model.serial))
            print("Export structure %s" % fileName)
            fileObj = open(fileName, 'w')
            fileObj.write(getPdbString(model,[ch.chain for ch in partner.sortedChains()],chainRename=partner.code,blankchain=True))
            fileObj.close()
            pdbmodellist.append(fileName)
        
        return pdbmodellist
    
    def __writePdbFileList(self,partner,pdbensemblelist):
        
        """Write a file containing PDB file names for the ensembles of a HADDOCK run. """

        listFileName = join(self.projectRoot,'PdbEnsemble%s.list' % (partner.code))
        print("Generate ensemble file list %s" % listFileName)
        makeBackup(listFileName)
        
        file = open(listFileName, 'w')
        for pdb in pdbensemblelist: file.write('"%s"\n' % pdb)
        file.close()
        
        return listFileName
    
    def __writeAmbigCnsTable(self):
        
        """Write a CNS style .tbl file containing the residue-residue ambiguous interaction
           restraints. The 'getAirSegments' function retrieves the active and passive
           residues from the haddockpartners. 'upperDistanceLimit' set the ambiguous
           interation restraints upper distance limit. Restraints are always writted to a file
           called 'ambig.tbl'. If only one molecular partner is defined the generated ambig.tbl
           file will be empty.
        """
        airs = {}
        for partner in self.haddockproject.sortedHaddockPartners(): airs[partner] = getAirSegments(partner)
        print("Export ambiguous interaction restraints as CNS style 'ambig.tbl' file")
        makeBackup(join(self.projectRoot,'ambig.tbl'))
        
        file = open(join(self.projectRoot,'ambig.tbl'),'w')
        for partner in airs:
            file.write("!\n")
            file.write("! HADDOCK AIR restraints for molecule %s\n" % partner.code)
            file.write("!\n")
            file.write("!\n")
            active = []
            for p in airs:
                if not p == partner: active += [(r,p.code) for r in airs[p]['active']]
                if not p == partner: active += [(r,p.code) for r in airs[p]['passive']]
            for act in airs[partner]['active']:
                for resid in range(len(active)):
                    if resid == 0:
                        file.write("assign ( resid %i and segid %s)\n" % (act,partner.code))
                        file.write("       (\n")
                        file.write("        ( resid %i and segid %s)\n" % (active[resid][0],active[resid][1]))
                    elif not resid == len(active)-1:
                        file.write("     or\n")
                        file.write("        ( resid %i and segid %s)\n" % (active[resid][0],active[resid][1]))
                    else:
                        file.write("     or\n")
                        file.write("       ( resid %i and segid %s)\n" % (active[resid][0],active[resid][1]))
                        file.write("       )  %1.1f 2.0 0.0\n" % partner.airUpperDistanceLimit)
                        file.write("!\n")
        file.close()
    
    def __writeRunCnsFile(self):
        
        """Write a HADDOCK version 2.1 compliant run.cns file"""
        
        print("Generate 'run.cns' file")
        makeBackup(join(self.projectRoot,'run.cns'))

        self.run = open(join(self.projectRoot,'run.cns'),'w')
        self.run.write(self.__runCnsHeader())
        
        self.run.write("\n{======== number of molecules for docking ==================}\n")
        self.run.write("{* number of components *}\n")
        self.run.write("{===>} ncomponents=%i;\n" % len(self.partners)) 
        
        self.run.write("\n{======================= filenames =========================}\n")
        self.run.write("{*  the name of your current project *}\n")
        self.run.write("{*  this will be used as name for the generated structures *}\n")
        self.run.write('{===>} fileroot="%s";\n' % self.haddockproject.name)

        self.run.write("\n{* RUN directory *}\n")
        self.run.write("{*  the absolute path of your current run, e.g. /home/haddock/run1*}\n")
        self.run.write('{===>} run_dir="%s";\n' % join(self.projectRoot,"run"+str(self.latestRun.serial)))
        
        for partner in self.partners:
            self.run.write('\n{* PDB file of molecule (protein) %s *}\n' % partner.code)
            self.run.write('{===>} prot_coor_%s="%s.pdb";\n' % (partner.code,self.fileNames[partner])) 
            self.run.write('{* PSF file of molecule (protein) %s *}\n' % partner.code)
            self.run.write('{===>} prot_psf_%s="%s.psf";\n' % (partner.code,self.fileNames[partner])) 
            self.run.write('{* segid of molecule (protein) %s *}\n' % partner.code)
            self.run.write('{===>} prot_segid_%s="%s";\n' % (partner.code,partner.code)) 
            self.run.write('{* fileroot of molecule (protein) %s *}\n' % partner.code)
            self.run.write('{===>} prot_root_%s="%s";\n' % (partner.code,self.fileNames[partner])) 
            self.run.write('{* Is molecule %s DNA? *}\n' % partner.code)
            self.run.write('{+ choice: true false +}\n')
            self.run.write('{===>} dna_%s=%s;\n' % (partner.code,str(partner.isDna).lower()))
    
        self.run.write('\n{ Atomname nomenclature }\n')
        self.run.write('{ set true if you have IUPAC (e.g. LEU HB2 and HB3 and not HB2 and HB1) data (e.g. from XEASY) }\n')
        self.run.write('{ choice: true false }\n')
        self.run.write('xplortodiana=false;\n') 

        self.run.write('\n{* Remove non-polar hydrogens? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} delenph=%s;\n' % str(self.latestRun.get('removeNonPolarH')).lower()) 

        self.run.write('\n{* HADDOCK directory *}\n')
        self.run.write('{*  the absolute path of the HADDOCK program files *}\n')
        self.run.write('{===>} haddock_dir="%s";\n' % self.latestRun.haddockDir) 

        self.run.write('\n{* Logfile directory *}\n')
        self.run.write('{* specify a directory for the large CNS log files *}\n')
        self.run.write('{===>} temptrash_dir="%s";\n' % join(self.projectRoot,"run"+str(self.latestRun.serial)))
        
        self.__histidinePatches()
        self.__setFlexInterface()
        self.__setSymmetry()
        self.__setDistanceRestraints()
        self.__setDnaRnaRestraints()
        self.__setDihedrals()
        self.__setKarplusCoupling()
        self.__setRDC()
        self.__setRelaxationData()
        self.__setTopologyFiles()
        self.__setEnergyParams()
        self.__setIterations()
        self.__setDockingProtocol()
        self.__setSolvatedDocking()
        self.__setWaterRefinement()
        self.__setScoringProtocol()
        self.__setAnalysisProtocol()
        self.__setClusterParams()
        
        self.run.write(self.__runCnsFooter())
        self.run.close()
    
    def __histidinePatches(self):
        
        HD1pstates = []; HE2pstates = []; length = 10; pcount = 0
        
        for partner in self.partners:
            
            if partner.autoHistidinePstate == True: 
                for nr in range(1,11): HD1pstates.append((partner.code,nr,0)); HE2pstates.append((partner.code,nr,0))
            else:
                count = 1
                molSystem = partner.molSystem
                ensembles = self.ccpnproject.sortedStructureEnsembles()
                if molSystem: ensembles = [e for e in ensembles if e.molSystem is molSystem]
                for chain in ensembles[0].sortedCoordChains():
                    for residue in chain.sortedResidues():
                        if residue.residue.ccpCode == 'His':
                            if residue.residue.descriptor == "prot:HD1;deprot:HE2": 
                                HD1pstates.append((partner.code,count,residue.residue.seqCode))
                                count += 1
                            elif residue.residue.descriptor == "deprot:HD1;prot:HE2":     
                                HE2pstates.append((partner.code,count,residue.residue.seqCode))
                                count += 1
                            elif residue.residue.descriptor == "prot:HD1;prot:HE2":     
                                HD1pstates.append((partner.code,count,residue.residue.seqCode))
                                HE2pstates.append((partner.code,count,residue.residue.seqCode))
                                count += 1
                            else: pass
            
            count = len(HD1pstates)
            while len(HD1pstates) < length: HD1pstates.append((partner.code,count+1,0)); count += 1
            count = len(HE2pstates)
            while len(HE2pstates) < length: HE2pstates.append((partner.code,count+1,0)); count += 1                    
        
            length += 10
            pcount += 1
    
        self.run.write('\n{==================== histidine patches =====================}\n')
        self.run.write('\n{* Patch to change doubly protonated HIS to singly protonated histidine (HD1) *}\n')
        self.run.write('{* just give the residue number of the histidines for the HISD patch, set them to zero if you do not want them *}\n')
        self.run.write('{+ table: rows=6 "molecule (Protein) A" "molecule (Protein) B"  "molecule (Protein) C"  "molecule (Protein) D"  "molecule (Protein) E"  "molecule (Protein) F" cols=10 "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" +}\n')

        self.run.write('\nnumhisd=%i;\n\n' % (len(HD1pstates)/pcount))

        for pstate in HD1pstates:
            self.run.write('{===>} %s_hisd_resid_%i=%i;\n' % (pstate[0],pstate[1],pstate[2]))

        self.run.write('\n{* Patch to change doubly protonated HIS to singly protonated histidine (HE2) *}\n')
        self.run.write('{* just give the residue number of the histidines for the HISE patch, set them to zero if you do not want them *}\n')
        self.run.write('{+ table: rows=6 "molecule (Protein) A" "molecule (Protein) B"  "molecule (Protein) C"  "molecule (Protein) D"  "molecule (Protein) E"  "molecule (Protein) F" cols=10 "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" +}\n')
        
        self.run.write('\nnumhise=%i;\n\n' % (len(HE2pstates)/pcount))
        
        for pstate in HE2pstates:
            self.run.write('{===>} %s_hisd_resid_%i=%i;\n' % (pstate[0],pstate[1],pstate[2]))

    def __setFlexInterface(self):
        
        self.run.write('\n{========= Definition of semi-flexible interface ============}\n')
        self.run.write('{* Define the interface of each molecule.*}\n')
        self.run.write('{* Side-chains and backbone of these residues will be allowed to move during semi-flexible refinement*}\n')
        
        for partner in self.partners:
            self.run.write('\n{* number of semi-flexible segments for molecule (protein) %s (-1 for automated mode) *}\n' % partner.code)
            self.run.write('{* Note that current max is 10 (edit the run.cns to add more segments *}\n')
            
            flex = getFlexibleResidues(partner)
            
            if partner.semiFlexMode == 'manual':
                self.run.write('\n{===>} nseg_%s=%i;\n' % (partner.code,len(flex['semi'])))    
            else:
                self.run.write('\n{===>} nseg_%s=-1;\n' % partner.code)
            
            self.run.write('\n{* Residues of molecule (protein) %s at interface *}\n' % partner.code)
            self.run.write('{+ table: rows=10 "segment 1" "segment 2" "segment 3" "segment 4" "segment 5" "segment 6" "segment 7" "segment 8" "segment 9" "segment 10" cols=2 "Start residue" "End residue" +}\n')
            self.run.write("\n")
            
            count = 1
            for semiflex in flex['semi']:
                self.run.write('{===>} %s_start_seg_%i="%i";\n' % (partner.code,count,semiflex[0])) 
                self.run.write('{===>} %s_end_seg_%i="%i";\n' % (partner.code,count,semiflex[1]))
                count += 1
        
        self.run.write('\n{=========== Definition of fully flexible segments ==========}\n')
        self.run.write('{* Define the fully flexible segment of each molecule.*}\n')
        self.run.write('{* These segments will be allowed to move at all stages of it1 *}\n')
        
        for partner in self.partners:
            self.run.write('\n{* Number of fully flexible segments for molecule (protein) %s *}\n' % partner.code)
            self.run.write('{* Note that current max is 5 (edit the run.cns to add more segments *}\n')
            
            flex = getFlexibleResidues(partner)
            
            self.run.write('\n{===>} nfle_%s=%i;\n' % (partner.code,len(flex['full'])))    
            
            self.run.write('\n{* Residues of molecule (protein) %s at interface *}\n' % partner.code)
            self.run.write('{+ table: rows=10 "segment 1" "segment 2" "segment 3" "segment 4" "segment 5" "segment 6" "segment 7" "segment 8" "segment 9" "segment 10" cols=2 "Start residue" "End residue" +}\n')
            self.run.write("\n")
            
            count = 1
            for fullyflex in flex['full']:
                self.run.write('{===>} %s_start_seg_%i="%i";\n' % (partner.code,count,fullyflex[0])) 
                self.run.write('{===>} %s_end_seg_%i="%i";\n' % (partner.code,count,fullyflex[1]))
                count += 1
    
    def __setSymmetry(self):
        
        symdict = {'ncs':[],'C2':[],'C3':[],'C5':[]}; symw = False
    
        for symmetry in self.latestRun.sortedSymmetryRestraints():
            symdict[symmetry.symmetryCode].append((1,symmetry.segmentLength,'A'))
        
        self.run.write('\n{====================== NCS restraints  =====================}\n')
        self.run.write('{* Do you want to use NCS restraints? *}\n')
        self.run.write('{+ choice: true false +}\n')
        
        if len(symdict['ncs']) > 0: self.run.write('{===>} ncs_on=true;\n')
        else: self.run.write('{===>} ncs_on=false;\n')
        
        self.run.write('\n{* Force constant for NCS restraints *}\n')
        self.run.write('{===>} kncs=%1.1f;\n' % self.latestRun.get('ncsRestraintConstant'))
        
        self.run.write('\n{* Number of NCS pairs *}\n')
        self.run.write('{===>} numncs=%i;\n' % len(symdict['ncs']))
        
        self.run.write('\n{* Define the segments pairs for NCS restraints *}\n')
        self.run.write('{+ table: rows=5 "pair 1" "pair 2" "pair 3" "pair 4" "pair 5" " cols=6 "Start res seg1" "End res seg1" "Segid seg1" "Start res seg2" "End res seg2" "Segid seg2" +}\n')
    
        symcount = 1
        for ncs in symdict['ncs']:
            self.run.write('{===>} ncs_sta1_%i="%i";\n' % (symcount,ncs[0]))
            self.run.write('{===>} ncs_end1_%i="%i";\n' % (symcount,ncs[1]))
            self.run.write('{===>} ncs_seg1_%i="%s";\n' % (symcount,ncs[2]))
            self.run.write('{===>} ncs_sta2_%i="%i";\n' % (symcount,ncs[0]))
            self.run.write('{===>} ncs_end2_%i="%i";\n' % (symcount,ncs[1]))
            self.run.write('{===>} ncs_seg2_%i="%s";\n' % (symcount,ncs[2]))
            symcount += 1
        
        self.run.write('\n{==================== Symmetry restraints  ==================}\n')
        self.run.write('{* Do you want to use symmetry restraints ? *}\n')
        self.run.write('{+ choice: true false +}\n')
        
        if len(symdict['C2']) > 0 or len(symdict['C3']) > 0 or len(symdict['C5']) > 0:
              self.run.write('{===>} sym_on=true;\n')
        else: self.run.write('{===>} sym_on=false;\n')    

        self.run.write('\n{* Force constant for symmetry restraints ? *}\n')
        self.run.write('{===>} ksym=%1.1f;\n' % self.latestRun.get('symmetryRestraintConstant'))
    
        self.run.write('\n{* Number of C2 symmetry pairs *}\n')
        self.run.write('{===>} numc2sym=%i;\n' % len(symdict['C2']))
        
        self.run.write('\n{* Define the segment pairs C2 symmetry restraints *}\n')
        self.run.write('{+ table: rows=10 "pair 1" "pair 2" "pair 3" "pair 4" "pair 5" "pair 6" "pair 7" "pair 8" "pair 9" "pair 10" cols=6 "Start res seg1" "End res seg1" "Segid seg1" "Start res seg2" "End res seg2" "Segid seg2" +}\n')
        
        symcount = 1
        for c2 in symdict['C2']:
            self.run.write('{===>} c2sym_sta1_%i="%i";\n' % (symcount,c2[0]))
            self.run.write('{===>} c2sym_end1_%i="%i";\n' % (symcount,c2[1]))
            self.run.write('{===>} c2sym_seg1_%i="%s";\n' % (symcount,c2[2]))
            self.run.write('{===>} c2sym_sta2_%i="%i";\n' % (symcount,c2[0]))
            self.run.write('{===>} c2sym_end2_%i="%i";\n' % (symcount,c2[1]))
            self.run.write('{===>} c2sym_seg2_%i="%s";\n' % (symcount,c2[2]))
            symcount += 1
        
        self.run.write('\n{* Number of C3 symmetry pairs *}\n')
        self.run.write('{===>} numc3sym=%i;\n' % len(symdict['C3']))

        self.run.write('\n{* Define the segment pairs C3 symmetry restraints *}\n')
        self.run.write('{+ table: rows=2 "triple 1" "triple 2" cols=9 "Start res seg1" "End res seg1" "Segid seg1" "Start res seg2" "End res seg2" "Segid seg2" "Start res seg3" "End res seg3" "Segid seg3" +}\n')

        symcount = 1
        for c3 in symdict['C3']:
            self.run.write('{===>} c3sym_sta1_%i="%i";\n' % (symcount,c3[0]))
            self.run.write('{===>} c3sym_end1_%i="%i";\n' % (symcount,c3[1]))
            self.run.write('{===>} c3sym_seg1_%i="%s";\n' % (symcount,c3[2]))
            self.run.write('{===>} c3sym_sta2_%i="%i";\n' % (symcount,c3[0]))
            self.run.write('{===>} c3sym_end2_%i="%i";\n' % (symcount,c3[1]))
            self.run.write('{===>} c3sym_seg2_%i="%s";\n' % (symcount,c3[2]))
            self.run.write('{===>} c3sym_sta3_%i="%i";\n' % (symcount,c3[0]))
            self.run.write('{===>} c3sym_end3_%i="%i";\n' % (symcount,c3[1]))
            self.run.write('{===>} c3sym_seg3_%i="%s";\n' % (symcount,c3[2]))
            symcount += 1

        self.run.write('\n{ Number of C4 symmetry quadruples }\n')
        self.run.write('{ Not yet used }\n')
        self.run.write('numc4sym=0;\n')

        self.run.write('\n{* Number of C5 symmetry pairs *}\n')
        self.run.write('{===>} numc5sym=%i;\n' % len(symdict['C5']))

        self.run.write('\n{* Define the segment pairs C5 symmetry restraints *}\n')
        self.run.write('{+ table: rows=5 "Segment1" "Segment2" "Segment3" "Segment4" "Segment5" cols=3 "Start residue" "End residue" "Segid" +}\n')

        symcount = 1
        for c5 in symdict['C5']:
            self.run.write('{===>} c5sym_sta1_%i="%i";\n' % (symcount,c5[0]))
            self.run.write('{===>} c5sym_end1_%i="%i";\n' % (symcount,c5[1]))
            self.run.write('{===>} c5sym_seg1_%i="%s";\n' % (symcount,c5[2]))
            self.run.write('{===>} c5sym_sta2_%i="%i";\n' % (symcount,c5[0]))
            self.run.write('{===>} c5sym_end2_%i="%i";\n' % (symcount,c5[1]))
            self.run.write('{===>} c5sym_seg2_%i="%s";\n' % (symcount,c5[2]))
            self.run.write('{===>} c5sym_sta3_%i="%i";\n' % (symcount,c5[0]))
            self.run.write('{===>} c5sym_end3_%i="%i";\n' % (symcount,c5[1]))
            self.run.write('{===>} c5sym_seg3_%i="%s";\n' % (symcount,c5[2]))
            self.run.write('{===>} c5sym_sta4_%i="%i";\n' % (symcount,c5[0]))
            self.run.write('{===>} c5sym_end4_%i="%i";\n' % (symcount,c5[1]))
            self.run.write('{===>} c5sym_seg4_%i="%s";\n' % (symcount,c5[2]))
            self.run.write('{===>} c5sym_sta5_%i="%i";\n' % (symcount,c5[0]))
            self.run.write('{===>} c5sym_end5_%i="%i";\n' % (symcount,c5[1]))
            self.run.write('{===>} c5sym_seg5_%i="%s";\n' % (symcount,c5[2]))
            symcount += 1    

    def __setDistanceRestraints(self):
        
        self.run.write('\n{=========================== Distance restraints  ========================}\n')
        self.run.write('{* Turn on/off and energy constants for distance restraints *}\n')
        self.run.write('{+ table: rows=3 "distances" "AIR (ambig)" "hbonds" cols=6 "firstIteration" "lastIteration" "hot" "cool1" "cool2" "cool3"+}\n\n')

        distRestraintEnergyStore = self.latestRun.findFirstHaddockEnergyTerm(code='distRestraintEnergyStore')
        for term in distRestraintEnergyStore.sortedEnergyTermParameters():
            self.run.write('{===>} %s=%s;\n' % (term.code,term.value))

        self.run.write('\n{* Do you want to randomly exclude a fraction of the ambiguous restraints (AIRs)? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} noecv=%s;\n' % str(self.latestRun.get('randomExcludeAir')).lower())

        self.run.write('\n{* Number of partitions for random exclusion (%excluded=100/number of partitions)? *}\n')
        self.run.write('{===>} ncvpart=%i;\n' % self.latestRun.get('randomExclParts'))

        self.run.write('\n{* Do you want to use hydrogen bond restraints? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} hbonds_on=%s;\n' % str(self.latestRun.get('useHBondRestraints')).lower()) 

        self.run.write('\n{* Do you want to define randomly ambiguous interaction restraints from accessible residues? *}\n')
        self.run.write('{* Only residues in the defined flexible segments will be considered *}\n')
        self.run.write('{* Note that this option is exclusive with any other distance restraints and only for it0    *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} ranair=%s;\n' % str(self.latestRun.get('randomAmbigRestraints')).lower())

        self.run.write('\n{* Do you want to define center of mass restraints to enforce contact between the molecules? *}\n')
        self.run.write('{* Note that these are only active during it0 and it1 *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} cmrest=%s;\n' % str(self.latestRun.get('centerOfMassRestraints')).lower())

        self.run.write('\n{* Force constant for center of mass restraints *}\n')
        self.run.write('{===>} kcont=%1.1f;\n' % self.latestRun.get('centerOfMassConstant'))

        self.run.write('\n{* Do you want to define surface contact restraints to enforce contact between the molecules? *}\n')
        self.run.write('{* Note that these are only active during it0 and it1 *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} surfrest=%s;\n' % str(self.latestRun.get('surfaceContactRestraints')).lower())

        self.run.write('\n{* Force constant for surface contact restraints *}\n')
        self.run.write('{===>} ksurf=%1.1f;\n' % self.latestRun.get('surfaceContactConstant'))
        
        self.run.write('\n{ Use automated distance restraints weighting }\n')
        self.run.write('{ choice: true false }\n')
        self.run.write('air_scaling=%s;\n' % str(self.latestRun.get('doAirScaling')).lower())  

        self.run.write('\n{ Define the number of distance restraints for automated weighting }\n')
        self.run.write('tot_unamb=%i;\n' % self.latestRun.get('numUnambRestautoAir')) 
        self.run.write('{ Define the number of AIR restraints for automated weighting }\n')
        self.run.write('tot_amb=%i;\n' % self.latestRun.get('numAmbRestautoAir'))

        self.run.write('\n{ potential shape }\n')
        autoDistanceRestraintWeightStore = self.latestRun.findFirstHaddockEnergyTerm(code='autoDistanceRestraintWeightStore')
        for term in autoDistanceRestraintWeightStore.sortedEnergyTermParameters():
            self.run.write('%s=%s;\n' % (term.code,term.value))
                
    def __setDnaRnaRestraints(self):
        
        self.run.write('\n{======================DNA-RNA restraints ============================}\n')
        self.run.write('{* Use DNA/RNA restraints (dna-rna_restraints.def in data/sequence)? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} dnarest_on=%s;\n' % (repr(self.latestRun.get('useDnaRestraints'))).lower())        
        
        if self.latestRun.get('useDnaRestraints') == True:
            for partner in self.partners:
                if partner.isDna == True: 
                    print ("Generate 'dna-rna_restraints.def' file")
                    dnarestraints = dnaRnaRestraints(self.ccpnproject,partner,self.projectRoot)
                    dnarestraints.writeToFile()
        
    def __setDihedrals(self):
        
        self.run.write('\n{=========================== dihedrals ==============================}\n')
        self.run.write('{* energy constants *}\n')
        self.run.write('{+ table: rows=1 "dihedrals" cols=5 "use?" "hot" "cool1" "cool2" "cool3" +}\n')

        self.run.write('\n{+ choice: true false +}\n')
        
        constraint = [c.code for c in self.latestRun.sortedHaddockEnergyTerms()]
        if 'DIHEDRAL' in constraint:
            self.run.write('{===>} dihedrals_on=true;\n')
        else:
            self.run.write('{===>} dihedrals_on=false;\n')    

        dihRestraintEnergyStore = self.latestRun.findFirstHaddockEnergyTerm(code='dihRestraintEnergyStore')
        for term in dihRestraintEnergyStore.sortedEnergyTermParameters():
            self.run.write('{===>} %s=%s;\n' % (term.code,term.value))
    
    def __setKarplusCoupling(self):
        
        self.run.write("""
{=========================== Karplus coupling restraints ====================}

{* Karplus coefficients: edit manually the run.cns file to specify them if needed   *}
{* The jcoupling restraint files should be present in the data/jcouplings directory *}
{* and named c1.tbl, c2.tbl, ... *}

 c1_on=false; 
 c1_karplusa=6.98; 
 c1_karplusb=-1.38; 
 c1_karplusc=1.72; 
 c1_karplusd=-60.0; 
 c1_hot=0.0; 
 c1_cool1=0.2; 
 c1_cool2=1.0; 
 c1_cool3=1.0; 

 c2_on=false; 
 c2_karplusa=6.98; 
 c2_karplusb=-1.38; 
 c2_karplusc=1.72; 
 c2_karplusd=-120.0; 
 c2_hot=0.0; 
 c2_cool1=0.2; 
 c2_cool2=1.0; 
 c2_cool3=1.0; 

 c3_on=false; 
 c3_karplusa=6.98; 
 c3_karplusb=-1.38; 
 c3_karplusc=1.72; 
 c3_karplusd=-120.0; 
 c3_hot=0.0; 
 c3_cool1=0.2; 
 c3_cool2=1.0; 
 c3_cool3=1.0; 

 c4_on=false; 
 c4_karplusa=6.98; 
 c4_karplusb=-1.38; 
 c4_karplusc=1.72; 
 c4_karplusd=-120.0; 
 c4_hot=0.0; 
 c4_cool1=0.2; 
 c4_cool2=1.0; 
 c4_cool3=1.0; 

 c5_on=false; 
 c5_karplusa=6.98; 
 c5_karplusb=-1.38; 
 c5_karplusc=1.72; 
 c5_karplusd=-120.0; 
 c5_hot=0.0; 
 c5_cool1=0.2; 
 c5_cool2=1.0; 
 c5_cool3=1.0;
        """)
    
    def __setRDC(self):
        
        self.run.write('\n{=========================== residual dipolar couplings ======================}\n')
        self.run.write('\n{* Parameters *}\n')
        self.run.write('{+ table: rows=5 "class1" "class2" "class3" "class4" "class5"\n')
        self.run.write('          cols=25 "type" "firstIt" "lastIt" "Ksani<br>(hot)" "Ksani<br>(cool1)" "Ksani<br>(cool2)" "Ksani<br>(cool3)" "R" "D"\n') 
        self.run.write(' "Kvean<br>(ini_bor_hot)" "Kvean<br>(fin_bor_hot)"\n')
        self.run.write(' "Kvean<br>(ini_bor_cool1)" "Kvean<br>(fin_bor_cool1)"\n') 
        self.run.write(' "Kvean<br>(ini_bor_cool2)" "Kvean<br>(fin_bor_cool2)"\n')
        self.run.write(' "Kvean<br>(ini_bor_cool3)" "Kvean<br>(fin_bor_cool3)"\n') 
        self.run.write(' "Kvean<br>(ini_cen_hot)" "Kvean<br>(fin_cen_hot)"\n') 
        self.run.write(' "Kvean<br>(ini_cen_cool1)" "Kvean<br>(fin_cen_cool1)"\n') 
        self.run.write(' "Kvean<br>(ini_cen_cool2)" "Kvean<br>(fin_cen_cool2)"\n') 
        self.run.write(' "Kvean<br>(ini_cen_cool3)" "Kvean<br>(fin_cen_cool3)"+}\n')
        
        constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'RDC']
        
        protocol_count = 1
        if len(constraint):
            for term in constraint:
                if term.fileName:
                    protocol = self.latestRun.findFirstHaddockEnergyTerm(code='rdcProtocolStore',termId=term.termId)
                    if protocol:    
                        self.run.write('\n{+ choice: "NO" "SANI" "VANGLE" +}\n')
                        for energyTerm in protocol.sortedEnergyTermParameters():
                            if energyTerm.code[0:3] == 'rdc':                            
                                if energyTerm.code == 'rdc_choice':
                                    self.run.write('{===>} rdc%i_%s=%s;\n' % (protocol_count,energyTerm.code[4:len(energyTerm.code)],['NO','SANI','VANGLE'][int(energyTerm.value)]))
                                else:
                                    self.run.write('{===>} rdc%i_%s=%s;\n' % (protocol_count,energyTerm.code[4:len(energyTerm.code)],energyTerm.value))
                            else:    
                                self.run.write('{===>} %s_%i=%1.1f;\n' % (energyTerm.code,protocol_count,energyTerm.value))
                    protocol_count += 1
                else:
                    print("RDC Warning: RDC energyTerm %i defined but no RDC CNS file associated" % term.termId)
        else:
            self.run.write('\n{+ choice: "NO" "SANI" "VANGLE" +}\n')
            for termId in range(1,6):    
                self.run.write('\n{===>} rdc%i_choice="NO";\n' % termId)
                for term in rdcProtocolStore['terms']: 
                    if term[0:3] == 'rdc':
                        self.run.write('{===>} rdc%i_%s=%1.1f;\n' % (termId,term[4:len(term)],rdcProtocolStore['terms'][term]))
                    else:
                        self.run.write('{===>} %s_%i=%1.1f;\n' % (term,termId,rdcProtocolStore['terms'][term]))        
            
    def __setRelaxationData(self):
        
        constraint = [c for c in self.latestRun.sortedHaddockEnergyTerms() if c.code == 'DANI']
        
        if len(constraint):
            self.run.write('\n{=========================== relaxation data ======================}\n')
            self.run.write('\n{* Parameters *}\n')
            self.run.write('{+ table: rows=5 "class1" "class2" "class3" "class4" "class5"\n')
            self.run.write('          cols=12 "type" "firstIt" "lastIt" "Kdani(hot)" "Kdani(cool1)" "Kdani(cool2)" "Kdani(cool3)" "Correlation time" "R" "D" "H frequency" "N frequency" +}\n')
            self.run.write('{+ choice: "NO" "DANI" +}\n')
            for term in constraint:
                if term.fileName:
                    protocol = self.latestRun.findFirstHaddockEnergyTerm(code='daniProtocolStore',termId=term.termId)
                    if protocol:
                        self.run.write('{===>} dan%i_choice="DANI";\n' % term.termId)    
                        for energyTerm in protocol.sortedEnergyTermParameters():
                            self.run.write('{===>} dan%i_%s=%1.1f;\n' % (term.termId,energyTerm.code,energyTerm.value))
                else:
                    print("DANI Warning: DANI energyTerm %i defined but no DANI CNS file associated" % term.termId)
                    self.run.write('{===>} dan%i_choice="NO";\n' % term.termId)    
        else:
            self.run.write('{=========================== relaxation data ======================}\n')
            self.run.write('\n{* Parameters *}\n')
            self.run.write('{+ table: rows=5 "class1" "class2" "class3" "class4" "class5"\n')
            self.run.write('          cols=12 "type" "firstIt" "lastIt" "Kdani(hot)" "Kdani(cool1)" "Kdani(cool2)" "Kdani(cool3)" "Correlation time" "R" "D" "H frequency" "N frequency" +}\n')
            self.run.write('{+ choice: "NO" "DANI" +}\n')
            for termId in range(1,6):    
                self.run.write('\n{===>} dan%i_choice="NO";\n' % termId)
                for term in daniProtocolStore['terms']: self.run.write('{===>} dan%i_%s=%1.1f;\n' % (termId,term,daniProtocolStore['terms'][term]))    
    
    def __setTopologyFiles(self):
        
        self.run.write('\n{===================== topology and parameter files ======================}\n')
        
        for partner in self.partners:
            if partner.isDna == True:
                self.run.write('\n{* topology file for molecule (protein) %s *}\n' % partner.code)
                self.run.write('{===>} prot_top_%s="dna-rna-allatom.top";\n' % partner.code)
                self.run.write('{* linkage file for molecule (protein) %s *}\n' % partner.code)
                self.run.write('{===>} prot_link_%s="dna-rna.link";\n' % partner.code)
                self.run.write('{* energy parameter file for molecule (protein) %s *}\n' % partner.code)
                self.run.write('{===>} prot_par_%s="dna-rna-allatom.param";\n' % partner.code)
            else:        
                self.run.write('\n{* topology file for molecule (protein) %s *}\n' % partner.code)
                self.run.write('{===>} prot_top_%s="topallhdg5.3.pro";\n' % partner.code)
                self.run.write('{* linkage file for molecule (protein) %s *}\n' % partner.code)
                self.run.write('{===>} prot_link_%s="topallhdg5.3.pep";\n' % partner.code)
                self.run.write('{* energy parameter file for molecule (protein) %s *}\n' % partner.code)
                self.run.write('{===>} prot_par_%s="parallhdg5.3.pro";\n' % partner.code)
        
        self.run.write('\n{* type of non-bonded parameters *}\n')
        self.run.write('{* specify the type of non-bonded interaction *}\n')
        self.run.write('{+ choice: "PROLSQ" "PARMALLH6" "PARALLHDG" "OPLSX" +}\n')
        self.run.write('{===>} par_nonbonded="%s";\n' % self.latestRun.get('nonBondedType'))

    def __setEnergyParams(self):
        
        self.run.write('\n{===================== energy and interaction parameters ==================}\n')

        self.run.write('\n{ Do you want to include dihedral angle energy terms? }\n')
        self.run.write('{ choice: true false }\n')
        self.run.write('dihedflag=%s;\n' % str(self.latestRun.get('doIncludeDihEnergy')).lower()) 

        self.run.write('\n{* Do you want to include the electrostatic energy term for docking? *}\n')
        self.run.write('{* Note that it will be automatically included in the solvent refinement *}\n')

        self.run.write('\n{* Include electrostatic during rigid body docking (it0)? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} elecflag_0=%s; \n' % str(self.latestRun.get('doRigidBodyElectrostatics')).lower()) 
        self.run.write('{* Include electrostatic during semi-flexible SA (it1)? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} elecflag_1=%s; \n' % str(self.latestRun.get('doSAElectrostatics')).lower()) 

        self.run.write('\n{* Give the epsilon constant for the electrostatic energy term? *}\n')
        self.run.write('{* Note that for explicit solvent refinement cdie with epsilon=1 is used *}\n')
        self.run.write('{===>} epsilon=%1.1f;\n' % self.latestRun.get('epsilon')) 

        self.run.write('\n{* Use constant (cdie) or distance-dependent (rdie) dielectric? *}\n')
        self.run.write('{+ choice: cdie rdie +}\n')
        self.run.write('{===>} dielec=%s;\n'  % self.latestRun.get('dielectricType'))

        self.run.write('\n{* Scaling of intermolecular interactions for rigid body EM*}\n')
        self.run.write('{===>} inter_rigid=%1.1f;\n'  % self.latestRun.get('rigidbodyIMinteractScaling'))

        self.run.write('\n{* Scaling of intermolecular interactions for semi-flexible SA*}\n')
        self.run.write('{+ table: rows=3 "Rigid body dynamic " "SA with flexible side-chains (cool2)" "SA with flexible backbone and side-chains (cool3)" cols=2 "Init value" "Final value" +}\n')

        semiflexInterMolScalingStore = self.latestRun.findFirstHaddockEnergyTerm(code='semiflexInterMolScalingStore')
        for term in semiflexInterMolScalingStore.sortedEnergyTermParameters():
            self.run.write('{===>} %s=%s;\n' % (term.code,term.value))
            
    def __setIterations(self):
        
        self.run.write('\n{===================== Number of structures to dock =======================}\n')
        self.run.write('\n{* Setting for the rigid-body (it0) and semi-flexible refiment (it1) *}\n')

        self.run.write('\n{* number of structures for rigid body docking *}\n')
        self.run.write('{===>} structures_0=%i;\n' % self.latestRun.get('numIt0Structures'))
        self.run.write('       keepstruct_0=&structures_0;\n')
        self.run.write('{* number of structures for refinement *}\n')
        self.run.write('{===>} structures_1=%i;\n' % self.latestRun.get('numIt1Structures')) 
        self.run.write('       keepstruct_1=&structures_1;\n')
        self.run.write('       keepstruct_2=&structures_1;\n')
        self.run.write('{* number of structures to be analysed*}\n')
        self.run.write('{===>} anastruc_1=%i;\n' % self.latestRun.get('numAnalysisStructures')) 
        self.run.write('       anastruc_0=&anastruc_1;\n')
        self.run.write('       anastruc_2=&anastruc_1;\n')

        self.run.write('\n{* Sampling of symmetry related solutions                       *}\n')

        self.run.write('\n{* Sample 180 degrees rotated solutions during rigid body EM?   *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} rotate180_0=%s;\n' % str(self.latestRun.get('rotate180It0')).lower())

        self.run.write('\n{* Sample 180 degrees rotated solutions during semi-flexible SA?*}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} rotate180_1=%s;\n' % str(self.latestRun.get('rotate180It1')).lower())    
    
    def __setDockingProtocol(self):
        
        self.run.write('\n{=========================== DOCKING protocol =============================}\n')
        self.run.write('\n{* Randomize starting orientations? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} randorien=%s;\n' % str(self.latestRun.get('radomizeStartOriention')).lower())

        self.run.write('\n{* Perform initial rigid body minimisation? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} rigidmini=%s;\n' % str(self.latestRun.get('initialRigidBodyMinim')).lower())

        self.run.write('\n{* Allow translation in rigid body minimisation? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} rigidtrans=%s;\n' % str(self.latestRun.get('doRigidTranslations')).lower())

        self.run.write('\n{* Number of trials for rigid body minimisation? *}\n')
        self.run.write('{===>} ntrials=%i;\n' % self.latestRun.get('nTrails'))

        self.run.write('\n{* initial seed for random number generator *}\n')
        self.run.write('{* change to get different initial velocities *}\n')
        self.run.write('{===>} iniseed=%i;\n' % self.latestRun.get('randomSeed')) 

        dockingProtocolStore = self.latestRun.findFirstHaddockEnergyTerm(code='dockingProtocolStore')
        self.dockingProtocolDict = {}
        for term in dockingProtocolStore.sortedEnergyTermParameters():
            self.dockingProtocolDict[term.code] = term.value

        self.run.write('\n{* temperature for rigid body high temperature TAD *}\n')
        self.run.write('{===>} tadhigh_t=%i;\n' % int(self.dockingProtocolDict['tadhigh_t'])) 

        self.run.write('\n{* initial temperature for rigid body first TAD cooling step *}\n')
        self.run.write('{===>} tadinit1_t=%i;\n' % int(self.dockingProtocolDict['tadinit1_t'])) 

        self.run.write('\n{* final temperature after first cooling step *}\n')
        self.run.write('{===>} tadfinal1_t=%i;\n' % int(self.dockingProtocolDict['tadfinal1_t'])) 

        self.run.write('\n{* initial temperature for second TAD cooling step with flexible side-chain at the inferface *}\n')
        self.run.write('{===>} tadinit2_t=%i;\n' % int(self.dockingProtocolDict['tadinit2_t'])) 

        self.run.write('\n{* finale temperature after second cooling step *}\n')
        self.run.write('{===>} tadfinal2_t=%i;\n' % int(self.dockingProtocolDict['tadfinal2_t'])) 

        self.run.write('\n{* initial temperature for third TAD cooling step with fully flexible interface *}\n')
        self.run.write('{===>} tadinit3_t=%i;\n' % int(self.dockingProtocolDict['tadinit3_t'])) 

        self.run.write('\n{* finale temperature after third cooling step *}\n')
        self.run.write('{===>} tadfinal3_t=%i;\n' % int(self.dockingProtocolDict['tadfinal3_t'])) 

        self.run.write('\n{* time step *}\n')
        self.run.write('{===>} timestep=%1.3f;\n' % self.dockingProtocolDict['timestep']) 
        
        self.run.write('\n{* factor for timestep in TAD *}\n')
        self.run.write('{===>} tadfactor=%i;\n' % int(self.dockingProtocolDict['tadfactor'])) 

        self.run.write('\n{* number of MD steps for rigid body high temperature TAD *}\n')
        self.run.write('{===>} initiosteps=%i;\n' % int(self.dockingProtocolDict['initiosteps'])) 

        self.run.write('\n{* number of MD steps during first rigid body cooling stage *}\n')
        self.run.write('{===>} cool1_steps=%i;\n' % int(self.dockingProtocolDict['cool1_steps'])) 

        self.run.write('\n{* number of MD steps during second cooling stage with flexible side-chains at interface *}\n')
        self.run.write('{===>} cool2_steps=%i;\n' % int(self.dockingProtocolDict['cool2_steps'])) 

        self.run.write('\n{* number of MD steps during third cooling stage with fully flexible interface *}\n')
        self.run.write('{===>} cool3_steps=%i;\n' % int(self.dockingProtocolDict['cool3_steps']))
    
    def __setSolvatedDocking(self):
        
        self.run.write('\n{======================= Solvated rigid body docking=======================}\n')
        self.run.write('\n{* perform solvated docking ? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} waterdock=%s;\n' % str(self.latestRun.get('doWaterDock')).lower())
        
        self.run.write('\n{* which method to use for solvating? *}\n')
        self.run.write('{* db: database-based (recommended), restraints: for restrained solvating to amino-acid most often forming\n')
        self.run.write('water mediated contacts and blank (""): for uniform waterlayer *}\n')
        self.run.write('{+ choice: "db" "restraints" "" +}\n')
        
        if self.latestRun.get('useDbSolvateMethod') == True:
            self.run.write('{===>} solvate_method="db";\n')
        else:
            self.run.write('{===>} solvate_method="";\n')    

        self.run.write('\n{* initial cutoff for restraints solvating method *}\n')
        self.run.write('{* all waters further away from a highly occuring water solvated residue will be removed in the generation\n') 
        self.run.write('of the initial solvation shell *}\n')
        self.run.write('{===>} water_restraint_initial=%1.1f;\n' % self.latestRun.get('waterInitRestCutoff'))

        self.run.write('\n{* cutoff for restraints solvating method *}\n')
        self.run.write('{* upper distance limit for defining distance restraints between water and amino-acids often found to be\n') 
        self.run.write('involved in water-mediated contacts *}\n')
        self.run.write('{===>} water_restraint_cutoff=%1.1f;\n' % self.latestRun.get('waterRestCutoff'))

        self.run.write('\n{* force constant for restraints solvating method *}\n')
        self.run.write('{===>} water_restraint_scale=%1.1f;\n' % self.latestRun.get('waterRestScale'))

        self.run.write('\n{* fraction of water to keep *}\n')
        self.run.write('{* this is the fraction of all interface water after the initial rigid body docking that will be kept (note)\n')
        self.run.write('that more waters might be removed if the interaction energy is unfavorable  *}\n')
        self.run.write('{===>} water_tokeep=%1.2f;\n' % self.latestRun.get('waterToKeep'))

        self.run.write('\n{* random fraction to be added to the fraction of water to keep *}\n')
        self.run.write('{===>} water_randfrac=%1.1f;\n' % self.latestRun.get('waterToAddRandom'))

        self.run.write('\n{* water-protein surface-cutoff *}\n')
        self.run.write('{* waters further away than this cutoff distance from any component of the complex will be removed *}\n')
        self.run.write('{===>} water_surfcutoff=%1.1f;\n' % self.latestRun.get('waterSurfaceCutoff'))

        self.run.write('\n{* do some water analysis *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} water_analysis=%s;\n' % str(self.latestRun.get('doWaterAnalysis')).lower())

        self.run.write('\n{* allows translation of water molecules during rigid-body docking, true or false: *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} transwater=%s;\n' % str(self.latestRun.get('doRigidBodyWaterTrans')).lower())

        self.run.write('\n{* number of different initial solvation shells to generate *}\n')
        self.run.write('{===>} waterensemble=%i;\n' % self.latestRun.get('numInitWaterShells'))
        
    def __setWaterRefinement(self):
        
        self.run.write('\n{==================== final explicit solvent refinement  ==================}\n')
        self.run.write('\n{* Do you want to refine your docking models in explicit solvent? *}\n')
        self.run.write('{+ choice: "yes" "no" +}\n')
    
        if self.latestRun.get('numWrefStructures') > 0:
            self.run.write('{===>} firstwater="yes";\n')
        else:
            self.run.write('{===>} firstwater="no";\n')    

        self.run.write('\n{* Which solvent do you want to use? *}\n')
        self.run.write('{+ choice: "water" "dmso" +}\n')
        self.run.write('{===>} solvent="%s";\n' % self.latestRun.get('solvent')) 

        self.run.write('\n{* number of structures for the explicit solvent refinement *}\n')
        self.run.write('{* the n best structures will be refined                    *}\n')
        self.run.write('{===>} waterrefine=%i;\n' % self.latestRun.get('numWrefStructures')) 
        self.run.write('       structures_2=&waterrefine;\n')

        self.run.write('\n{* number of steps for heating phase (100, 200, 300K)?      *}\n')
        self.run.write('{===>} waterheatsteps=%i;\n' % int(self.dockingProtocolDict['waterheatsteps']))

        self.run.write('\n{* number of steps for 300K sampling phase?                 *}\n')
        self.run.write('{===>} watersteps=%i;\n' % int(self.dockingProtocolDict['watersteps']))

        self.run.write('\n{* number of steps for cooling phase (300, 200, 100K)?      *}\n')
        self.run.write('{===>} watercoolsteps=%i;\n' % int(self.dockingProtocolDict['watercoolsteps']))

        self.run.write('\n{* write additional PDB files including solvent ?           *}\n')
        self.run.write('{+ choice: true false +}\n')
        
        if self.latestRun.get('doWaterDock') == True: 
            self.run.write('{===>} keepwater=true;\n')
        else:
            self.run.write('{===>} keepwater=false;\n')     

        self.run.write('\n{ calculate explicit desolvation energy (note this will double the cpu requirements) }\n')
        self.run.write('{ choice: true false }\n')
        self.run.write('{===>} calcdesolv=%s;\n' % str(self.latestRun.get('calcDesolvation')).lower())    
    
    def __setScoringProtocol(self):
        
        """Write the scoring weights for the different docking stages"""
        
        self.run.write('\n{================================ Scoring =================================}\n')
        self.run.write('\n{* Settings for the scoring of the docking solutions *}\n')

        self.run.write('\n{* Define the weights for the various terms for the sorting of structures (scoring) *}\n')
        self.run.write('{+ table: rows=11 "Evdw" "Eelec" "Edist" "Esani" "Edani" "Evean" "Ecdih" "Esym" "BSA" "dEint" "Edesolv"\n')
        self.run.write('         cols=3 "Rigid body EM" "semi-flexible SA" "Water refinement" +}\n')
        
        scoringWeights = [(sw.term, sw.stage, sw.value) for sw in self.latestRun.scoringWeights]
        scoringWeights.sort()
        scoreTerm = ''
        for term, stage, scoreWeight in scoringWeights:
            if term == scoreTerm:
                self.run.write("{===>} w_%s_%i=%1.1f;\n" % (term,stage,scoreWeight))
            elif scoreTerm == '':
                self.run.write("{===>} w_%s_%i=%1.1f;\n" % (term,stage,scoreWeight))
                scoreTerm = term
            else:
                self.run.write("\n")
                self.run.write("{===>} w_%s_%i=%1.1f;\n" % (term,stage,scoreWeight))
                scoreTerm = term
                
        self.run.write('\n{ Perform smoothed-scoring selection for rigid-body docking solutions ? }\n')
        self.run.write('{ choice: true false }\n')
        self.run.write('{ currently not used }\n')
        self.run.write('{===>} smoothing=false;\n')

        self.run.write('\n{* It is possible to skip structures in the selection of structure in it0 *}\n')
        self.run.write('{* Give for this the number of structures to skip: *}\n')
        self.run.write('{===>} skip_struc=%i;\n' % self.latestRun.get('skipStructures'))
    
    def __setAnalysisProtocol(self):
        
        self.run.write('\n{======================= analysis and clustering ==========================}\n')
        self.run.write('\n{* Cutoff distance (proton-acceptor) to define an hydrogen bond? *}\n')
        self.run.write('{===>} dist_hb=%1.1f;\n' % self.latestRun.get('analysisDistHBond'))

        self.run.write('\n{* Cutoff distance (carbon-carbon) to define an hydrophobic contact? *}\n')
        self.run.write('{===>} dist_nb=%1.1f;\n' % self.latestRun.get('analysisDistNonbond'))

        self.run.write('\n{* RMSD cutoff for clustering? *}\n')
        self.run.write('{===>} clust_rmsd=%1.1f;\n' % self.latestRun.get('analysisClustRmsd'))

        self.run.write('\n{* Minimum cluster size? *}\n')
        self.run.write('{===>} clust_size=%i;\n' % self.latestRun.get('analysisClustSize'))
        
        self.run.write('\n{========================= Structure quality analysis =====================}\n')
        self.run.write('\n{* specify location of the executables, e.g.  procheck.scr script. Make sure that your PRODIR and PROSA_BASE system variables are set correctly. Leave fields empty if you do not want to perform these checks *}\n')

        self.run.write('\n{* Procheck executable: *}\n')
        self.run.write('{===>} procheck_exe="";\n')
        self.run.write('{* Procheck_comp executable: *}\n')
        self.run.write('{===>} procheckcomp_exe="";\n')
        self.run.write('{* Whatif executable: *}\n')
        self.run.write('{===>} whatif_exe="";\n')
        self.run.write('{* Prosa executable: *}\n')
        self.run.write('{===>} prosa_exe="";\n')
        self.run.write('{* Number of PDB files to analyse: *}\n')
        self.run.write('{===>} how_many_pdb="20";\n')

        self.run.write('\n{======================= final clean-up ===================================}\n')
        self.run.write('\n{* Clean up the run directory after completion (only files for struct #1 are kept) ? *}\n')
        self.run.write('{+ choice: true false +}\n')
        self.run.write('{===>} cleanup=true;\n')    
    
    def __setClusterParams(self):
        
        self.run.write('\n{============================ parallel jobs ===============================}\n')
        self.run.write('\n{* How many nodes do you want to use in parallel? *}\n')
        self.run.write('{* leave unused fields blank, make sure that the queues are actually running *}\n')
        self.run.write('{+ table: rows=10 "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"\n')
        self.run.write(' cols=3 "queue command" "cns executable" "number of jobs" +}\n')

        self.run.write('\n{===>} queue_1="%s";\n' % self.latestRun.queueCommand) 
        self.run.write('{===>} cns_exe_1="%s";\n' % self.latestRun.cnsExecutable)
        self.run.write('{===>} cpunumber_1=%i;\n' % self.latestRun.cpuNumber) 
        
    def __runCnsHeader(self):
        
        return"""!$Revision: 1.4 $
!$Date: 2010-12-02 10:17:31 $
!$RCSfile: HaddockExportClassic.py,v $


module(
iteration;
filenames;
data;
iterations;
saprotocol;
refine;
toppar; 
analysis;
)

{+ File: run.cns +}
{+ Description: this file contains all necessary information to run HADDOCK. +}

{+ Authors: Alexandre Bonvin<br>
Version: 2.0, January 2007 <br><br>
Initially adapted from ARIA of Nilges and Linge +}

! Please cite the following references when using this protocol: 
{+ reference: Cyril Dominguez, Rolf Boelens and Alexandre M.J.J. Bonvin (2003).  HADDOCK: a protein-protein docking approach 
based on biochemical and/or biophysical information. <i>J. Am. Chem. Soc.</i> <b>125</b>, 1731-1737.
<p>
<b>When using <i>residual dipolar couplings</i> in HADDOCK cite in addition:</b><p>
<LI>A.D.J. van Dijk, D. Fushman and A.M.J.J. Bonvin (2005). Various strategies of using residual dipolar 
couplings in NMR-driven protein docking: Application to Lys48-linked di-ubiquitin and validation against 
15N-relaxation data. <EM>Proteins: Struc. Funct. & Bioinformatics</EM>, <STRONG>60</STRONG>, 367-381.</li>
<p>
<b>When using <i>diffusion anisotropy data</i> in HADDOCK cite in addition:</b><p>
<li>A.D.J. van Dijk, R. Kaptein, R. Boelens and A.M.J.J. Bonvin (2006). Combining NMR relaxation with 
chemical shift perturbation data to drive protein-protein docking. <EM>J. Biomol. NMR</EM>, 
<STRONG>34</STRONG>, 237-244.</li>
<p>
<b>When using <i>solvated docking</i> in HADDOCK cite in addition:</b><p>
<li>A.D.J. van Dijk and A.M.J.J. Bonvin (2006). Solvated docking: introducing water into the modelling 
of biomolecular complexes. <EM>Bioinformatics</EM>,  <STRONG>22</STRONG> 2340-2347.
<p>
<b>When performing <i>flexible protein-DNA docking</i> using HADDOCK cite in addition:</b><p>
<li>M. van Dijk, A.D.J. van Dijk, V. Hsu, R. Boelens and  A.M.J.J. Bonvin (2006).
Information-driven Protein-DNA Docking using HADDOCK: it is a matter of flexibility.
<EM>Nucl. Acids Res.</EM>, <STRONG>34</STRONG> 3317-3325.</li>
+}

{- Guidelines for using this file:
   - all strings must be quoted by double-quotes
   - logical variables (true/false) are not quoted
   - do not remove any evaluate statements from the file
   - pathnames should not exceed 80 characters -}
{- begin block parameter definition -} define(
"""        

    def __runCnsFooter(self):
        
        return """
{===========================================================================}
{        things below this line do not normally need to be changed          }
{===========================================================================}

) {- end block parameter definition -}

!for global parameters (local variables (suffix ) => global variables):
evaluate (&saprotocol.randorien=&randorien)
evaluate (&saprotocol.rigidmini=&rigidmini)
evaluate (&saprotocol.rigidtrans=&rigidtrans)
evaluate (&saprotocol.ntrials=&ntrials)
evaluate (&saprotocol.iniseed=&iniseed)
evaluate (&saprotocol.tadhigh_t=&tadhigh_t)
evaluate (&saprotocol.t1_init=&tadinit1_t)
evaluate (&saprotocol.t2_init=&tadinit2_t)
evaluate (&saprotocol.t3_init=&tadinit3_t)
evaluate (&saprotocol.t1_final=&tadfinal1_t)
evaluate (&saprotocol.t2_final=&tadfinal2_t)
evaluate (&saprotocol.t3_final=&tadfinal3_t)
evaluate (&saprotocol.inter_rigid=&inter_rigid)
evaluate (&saprotocol.inter_init_rigid=&init_rigid)
evaluate (&saprotocol.inter_fin_rigid=&fin_rigid)
evaluate (&saprotocol.inter_init_cool2=&init_cool2)
evaluate (&saprotocol.inter_fin_cool2=&fin_cool2)
evaluate (&saprotocol.inter_init_cool3=&init_cool3)
evaluate (&saprotocol.inter_fin_cool3=&fin_cool3)
evaluate (&saprotocol.tempstep=50)
evaluate (&saprotocol.timestep=&timestep)
evaluate (&saprotocol.tadfactor=&tadfactor)
evaluate (&saprotocol.initiosteps=&initiosteps)
evaluate (&saprotocol.cool1_steps=&cool1_steps)
evaluate (&saprotocol.cool2_steps=&cool2_steps)
evaluate (&saprotocol.cool3_steps=&cool3_steps)
evaluate (&saprotocol.fbeta=100)
evaluate (&saprotocol.mass=100)

evaluate (&filenames.fileroot=&fileroot)
evaluate (&filenames.template=&fileroot + "_1.pdb")

evaluate (&iterations.ini_count    =1) 
evaluate (&iterations.structures   =&structures_$iteration)
evaluate (&iterations.keepstruct   =&keepstruct_$iteration)
evaluate (&iterations.filesort     =&filesort_$iteration)
evaluate (&iterations.w_vdw        =&w_vdw_$iteration)
evaluate (&iterations.w_elec       =&w_elec_$iteration)
evaluate (&iterations.w_dist       =&w_dist_$iteration)
evaluate (&iterations.w_sani       =&w_sani_$iteration)
evaluate (&iterations.w_dani       =&w_dani_$iteration)
evaluate (&iterations.w_vean       =&w_vean_$iteration)
evaluate (&iterations.w_cdih       =&w_cdih_$iteration)
evaluate (&iterations.w_sym        =&w_sym_$iteration)
evaluate (&iterations.w_bsa        =&w_bsa_$iteration)
evaluate (&iterations.w_deint      =&w_deint_$iteration)
evaluate (&iterations.w_desolv     =&w_desolv_$iteration)
evaluate (&iterations.anastruc     =&anastruc_$iteration)
evaluate (&iterations.rotate180    =&rotate180_$iteration)


!topology and parameters, sequence file, template file, interface definition:
evaluate (&toppar.prot_top_1=&prot_top_A )
evaluate (&toppar.prot_top_2=&prot_top_B )
evaluate (&toppar.prot_top_3=&prot_top_C )
evaluate (&toppar.prot_top_4=&prot_top_D )
evaluate (&toppar.prot_top_5=&prot_top_E )
evaluate (&toppar.prot_top_6=&prot_top_F )

evaluate (&toppar.prot_link_1=&prot_link_A )
evaluate (&toppar.prot_link_2=&prot_link_B )
evaluate (&toppar.prot_link_3=&prot_link_C )
evaluate (&toppar.prot_link_4=&prot_link_D )
evaluate (&toppar.prot_link_5=&prot_link_E )
evaluate (&toppar.prot_link_6=&prot_link_F )

evaluate (&toppar.prot_par_1=&prot_par_A )
evaluate (&toppar.prot_par_2=&prot_par_B )
evaluate (&toppar.prot_par_3=&prot_par_C )
evaluate (&toppar.prot_par_4=&prot_par_D )
evaluate (&toppar.prot_par_5=&prot_par_E )
evaluate (&toppar.prot_par_6=&prot_par_F )

evaluate (&toppar.par_nonbonded=&par_nonbonded)

evaluate (&toppar.prot_coor_1=&prot_coor_A)
evaluate (&toppar.prot_coor_2=&prot_coor_B)
evaluate (&toppar.prot_coor_3=&prot_coor_C)
evaluate (&toppar.prot_coor_4=&prot_coor_D)
evaluate (&toppar.prot_coor_5=&prot_coor_E)
evaluate (&toppar.prot_coor_6=&prot_coor_F)

evaluate (&toppar.prot_root_1=&prot_root_A)
evaluate (&toppar.prot_root_2=&prot_root_B)
evaluate (&toppar.prot_root_3=&prot_root_C)
evaluate (&toppar.prot_root_4=&prot_root_D)
evaluate (&toppar.prot_root_5=&prot_root_E)
evaluate (&toppar.prot_root_6=&prot_root_F)

evaluate (&toppar.dna_1=&dna_A)
evaluate (&toppar.dna_2=&dna_B)
evaluate (&toppar.dna_3=&dna_C)
evaluate (&toppar.dna_4=&dna_D)
evaluate (&toppar.dna_5=&dna_E)
evaluate (&toppar.dna_6=&dna_F)

evaluate (&toppar.prot_segid_1=&prot_segid_A)
evaluate (&toppar.prot_segid_2=&prot_segid_B)
evaluate (&toppar.prot_segid_3=&prot_segid_C)
evaluate (&toppar.prot_segid_4=&prot_segid_D)
evaluate (&toppar.prot_segid_5=&prot_segid_E)
evaluate (&toppar.prot_segid_6=&prot_segid_F)

evaluate (&data.ncomponents=&ncomponents) 

evaluate (&toppar.nseg_1=&nseg_A)
evaluate (&toppar.nseg_2=&nseg_B)
evaluate (&toppar.nseg_3=&nseg_C)
evaluate (&toppar.nseg_4=&nseg_D)
evaluate (&toppar.nseg_5=&nseg_E)
evaluate (&toppar.nseg_6=&nseg_F)
evaluate (&toppar.start_seg_1_1=&A_start_seg_1)
evaluate (&toppar.start_seg_1_2=&A_start_seg_2)
evaluate (&toppar.start_seg_1_3=&A_start_seg_3)
evaluate (&toppar.start_seg_1_4=&A_start_seg_4)
evaluate (&toppar.start_seg_1_5=&A_start_seg_5)
evaluate (&toppar.start_seg_1_6=&A_start_seg_6)
evaluate (&toppar.start_seg_1_7=&A_start_seg_7)
evaluate (&toppar.start_seg_1_8=&A_start_seg_8)
evaluate (&toppar.start_seg_1_9=&A_start_seg_9)
evaluate (&toppar.start_seg_1_10=&A_start_seg_10)
evaluate (&toppar.end_seg_1_1=&A_end_seg_1)
evaluate (&toppar.end_seg_1_2=&A_end_seg_2)
evaluate (&toppar.end_seg_1_3=&A_end_seg_3)
evaluate (&toppar.end_seg_1_4=&A_end_seg_4)
evaluate (&toppar.end_seg_1_5=&A_end_seg_5)
evaluate (&toppar.end_seg_1_6=&A_end_seg_6)
evaluate (&toppar.end_seg_1_7=&A_end_seg_7)
evaluate (&toppar.end_seg_1_8=&A_end_seg_8)
evaluate (&toppar.end_seg_1_9=&A_end_seg_9)
evaluate (&toppar.end_seg_1_10=&A_end_seg_10)
evaluate (&toppar.start_seg_2_1=&B_start_seg_1)
evaluate (&toppar.start_seg_2_2=&B_start_seg_2)
evaluate (&toppar.start_seg_2_3=&B_start_seg_3)
evaluate (&toppar.start_seg_2_4=&B_start_seg_4)
evaluate (&toppar.start_seg_2_5=&B_start_seg_5)
evaluate (&toppar.start_seg_2_6=&B_start_seg_6)
evaluate (&toppar.start_seg_2_7=&B_start_seg_7)
evaluate (&toppar.start_seg_2_8=&B_start_seg_8)
evaluate (&toppar.start_seg_2_9=&B_start_seg_9)
evaluate (&toppar.start_seg_2_10=&B_start_seg_10)
evaluate (&toppar.end_seg_2_1=&B_end_seg_1)
evaluate (&toppar.end_seg_2_2=&B_end_seg_2)
evaluate (&toppar.end_seg_2_3=&B_end_seg_3)
evaluate (&toppar.end_seg_2_4=&B_end_seg_4)
evaluate (&toppar.end_seg_2_5=&B_end_seg_5)
evaluate (&toppar.end_seg_2_6=&B_end_seg_6)
evaluate (&toppar.end_seg_2_7=&B_end_seg_7)
evaluate (&toppar.end_seg_2_8=&B_end_seg_8)
evaluate (&toppar.end_seg_2_9=&B_end_seg_9)
evaluate (&toppar.end_seg_2_10=&B_end_seg_10)
evaluate (&toppar.start_seg_3_1=&C_start_seg_1)
evaluate (&toppar.start_seg_3_2=&C_start_seg_2)
evaluate (&toppar.start_seg_3_3=&C_start_seg_3)
evaluate (&toppar.start_seg_3_4=&C_start_seg_4)
evaluate (&toppar.start_seg_3_5=&C_start_seg_5)
evaluate (&toppar.start_seg_3_6=&C_start_seg_6)
evaluate (&toppar.start_seg_3_7=&C_start_seg_7)
evaluate (&toppar.start_seg_3_8=&C_start_seg_8)
evaluate (&toppar.start_seg_3_9=&C_start_seg_9)
evaluate (&toppar.start_seg_3_10=&C_start_seg_10)
evaluate (&toppar.end_seg_3_1=&C_end_seg_1)
evaluate (&toppar.end_seg_3_2=&C_end_seg_2)
evaluate (&toppar.end_seg_3_3=&C_end_seg_3)
evaluate (&toppar.end_seg_3_4=&C_end_seg_4)
evaluate (&toppar.end_seg_3_5=&C_end_seg_5)
evaluate (&toppar.end_seg_3_6=&C_end_seg_6)
evaluate (&toppar.end_seg_3_7=&C_end_seg_7)
evaluate (&toppar.end_seg_3_8=&C_end_seg_8)
evaluate (&toppar.end_seg_3_9=&C_end_seg_9)
evaluate (&toppar.end_seg_3_10=&C_end_seg_10)
evaluate (&toppar.start_seg_4_1=&D_start_seg_1)
evaluate (&toppar.start_seg_4_2=&D_start_seg_2)
evaluate (&toppar.start_seg_4_3=&D_start_seg_3)
evaluate (&toppar.start_seg_4_4=&D_start_seg_4)
evaluate (&toppar.start_seg_4_5=&D_start_seg_5)
evaluate (&toppar.start_seg_4_6=&D_start_seg_6)
evaluate (&toppar.start_seg_4_7=&D_start_seg_7)
evaluate (&toppar.start_seg_4_8=&D_start_seg_8)
evaluate (&toppar.start_seg_4_9=&D_start_seg_9)
evaluate (&toppar.start_seg_4_10=&D_start_seg_10)
evaluate (&toppar.end_seg_4_1=&D_end_seg_1)
evaluate (&toppar.end_seg_4_2=&D_end_seg_2)
evaluate (&toppar.end_seg_4_3=&D_end_seg_3)
evaluate (&toppar.end_seg_4_4=&D_end_seg_4)
evaluate (&toppar.end_seg_4_5=&D_end_seg_5)
evaluate (&toppar.end_seg_4_6=&D_end_seg_6)
evaluate (&toppar.end_seg_4_7=&D_end_seg_7)
evaluate (&toppar.end_seg_4_8=&D_end_seg_8)
evaluate (&toppar.end_seg_4_9=&D_end_seg_9)
evaluate (&toppar.end_seg_4_10=&D_end_seg_10)
evaluate (&toppar.start_seg_5_1=&E_start_seg_1)
evaluate (&toppar.start_seg_5_2=&E_start_seg_2)
evaluate (&toppar.start_seg_5_3=&E_start_seg_3)
evaluate (&toppar.start_seg_5_4=&E_start_seg_4)
evaluate (&toppar.start_seg_5_5=&E_start_seg_5)
evaluate (&toppar.start_seg_5_6=&E_start_seg_6)
evaluate (&toppar.start_seg_5_7=&E_start_seg_7)
evaluate (&toppar.start_seg_5_8=&E_start_seg_8)
evaluate (&toppar.start_seg_5_9=&E_start_seg_9)
evaluate (&toppar.start_seg_5_10=&E_start_seg_10)
evaluate (&toppar.end_seg_5_1=&E_end_seg_1)
evaluate (&toppar.end_seg_5_2=&E_end_seg_2)
evaluate (&toppar.end_seg_5_3=&E_end_seg_3)
evaluate (&toppar.end_seg_5_4=&E_end_seg_4)
evaluate (&toppar.end_seg_5_5=&E_end_seg_5)
evaluate (&toppar.end_seg_5_6=&E_end_seg_6)
evaluate (&toppar.end_seg_5_7=&E_end_seg_7)
evaluate (&toppar.end_seg_5_8=&E_end_seg_8)
evaluate (&toppar.end_seg_5_9=&E_end_seg_9)
evaluate (&toppar.end_seg_5_10=&E_end_seg_10)
evaluate (&toppar.start_seg_6_1=&F_start_seg_1)
evaluate (&toppar.start_seg_6_2=&F_start_seg_2)
evaluate (&toppar.start_seg_6_3=&F_start_seg_3)
evaluate (&toppar.start_seg_6_4=&F_start_seg_4)
evaluate (&toppar.start_seg_6_5=&F_start_seg_5)
evaluate (&toppar.start_seg_6_6=&F_start_seg_6)
evaluate (&toppar.start_seg_6_7=&F_start_seg_7)
evaluate (&toppar.start_seg_6_8=&F_start_seg_8)
evaluate (&toppar.start_seg_6_9=&F_start_seg_9)
evaluate (&toppar.start_seg_6_10=&F_start_seg_10)
evaluate (&toppar.end_seg_6_1=&F_end_seg_1)
evaluate (&toppar.end_seg_6_2=&F_end_seg_2)
evaluate (&toppar.end_seg_6_3=&F_end_seg_3)
evaluate (&toppar.end_seg_6_4=&F_end_seg_4)
evaluate (&toppar.end_seg_6_5=&F_end_seg_5)
evaluate (&toppar.end_seg_6_6=&F_end_seg_6)
evaluate (&toppar.end_seg_6_7=&F_end_seg_7)
evaluate (&toppar.end_seg_6_8=&F_end_seg_8)
evaluate (&toppar.end_seg_6_9=&F_end_seg_9)
evaluate (&toppar.end_seg_6_10=&F_end_seg_10)

evaluate (&toppar.nfle_1=&nfle_A)
evaluate (&toppar.start_fle_1_1=&A_start_fle_1)
evaluate (&toppar.start_fle_1_2=&A_start_fle_2)
evaluate (&toppar.start_fle_1_3=&A_start_fle_3)
evaluate (&toppar.start_fle_1_4=&A_start_fle_4)
evaluate (&toppar.start_fle_1_5=&A_start_fle_5)
evaluate (&toppar.end_fle_1_1=&A_end_fle_1)
evaluate (&toppar.end_fle_1_2=&A_end_fle_2)
evaluate (&toppar.end_fle_1_3=&A_end_fle_3)
evaluate (&toppar.end_fle_1_4=&A_end_fle_4)
evaluate (&toppar.end_fle_1_5=&A_end_fle_5)
evaluate (&toppar.nfle_2=&nfle_B)
evaluate (&toppar.start_fle_2_1=&B_start_fle_1)
evaluate (&toppar.start_fle_2_2=&B_start_fle_2)
evaluate (&toppar.start_fle_2_3=&B_start_fle_3)
evaluate (&toppar.start_fle_2_4=&B_start_fle_4)
evaluate (&toppar.start_fle_2_5=&B_start_fle_5)
evaluate (&toppar.end_fle_2_1=&B_end_fle_1)
evaluate (&toppar.end_fle_2_2=&B_end_fle_2)
evaluate (&toppar.end_fle_2_3=&B_end_fle_3)
evaluate (&toppar.end_fle_2_4=&B_end_fle_4)
evaluate (&toppar.end_fle_2_5=&B_end_fle_5)
evaluate (&toppar.nfle_3=&nfle_C)
evaluate (&toppar.start_fle_3_1=&C_start_fle_1)
evaluate (&toppar.start_fle_3_2=&C_start_fle_2)
evaluate (&toppar.start_fle_3_3=&C_start_fle_3)
evaluate (&toppar.start_fle_3_4=&C_start_fle_4)
evaluate (&toppar.start_fle_3_5=&C_start_fle_5)
evaluate (&toppar.end_fle_3_1=&C_end_fle_1)
evaluate (&toppar.end_fle_3_2=&C_end_fle_2)
evaluate (&toppar.end_fle_3_3=&C_end_fle_3)
evaluate (&toppar.end_fle_3_4=&C_end_fle_4)
evaluate (&toppar.end_fle_3_5=&C_end_fle_5)
evaluate (&toppar.nfle_4=&nfle_D)
evaluate (&toppar.start_fle_4_1=&D_start_fle_1)
evaluate (&toppar.start_fle_4_2=&D_start_fle_2)
evaluate (&toppar.start_fle_4_3=&D_start_fle_3)
evaluate (&toppar.start_fle_4_4=&D_start_fle_4)
evaluate (&toppar.start_fle_4_5=&D_start_fle_5)
evaluate (&toppar.end_fle_4_1=&D_end_fle_1)
evaluate (&toppar.end_fle_4_2=&D_end_fle_2)
evaluate (&toppar.end_fle_4_3=&D_end_fle_3)
evaluate (&toppar.end_fle_4_4=&D_end_fle_4)
evaluate (&toppar.end_fle_4_5=&D_end_fle_5)
evaluate (&toppar.nfle_5=&nfle_E)
evaluate (&toppar.start_fle_5_1=&E_start_fle_1)
evaluate (&toppar.start_fle_5_2=&E_start_fle_2)
evaluate (&toppar.start_fle_5_3=&E_start_fle_3)
evaluate (&toppar.start_fle_5_4=&E_start_fle_4)
evaluate (&toppar.start_fle_5_5=&E_start_fle_5)
evaluate (&toppar.end_fle_5_1=&E_end_fle_1)
evaluate (&toppar.end_fle_5_2=&E_end_fle_2)
evaluate (&toppar.end_fle_5_3=&E_end_fle_3)
evaluate (&toppar.end_fle_5_4=&E_end_fle_4)
evaluate (&toppar.end_fle_5_5=&E_end_fle_5)
evaluate (&toppar.nfle_6=&nfle_F)
evaluate (&toppar.start_fle_6_1=&F_start_fle_1)
evaluate (&toppar.start_fle_6_2=&F_start_fle_2)
evaluate (&toppar.start_fle_6_3=&F_start_fle_3)
evaluate (&toppar.start_fle_6_4=&F_start_fle_4)
evaluate (&toppar.start_fle_6_5=&F_start_fle_5)
evaluate (&toppar.end_fle_6_1=&F_end_fle_1)
evaluate (&toppar.end_fle_6_2=&F_end_fle_2)
evaluate (&toppar.end_fle_6_3=&F_end_fle_3)
evaluate (&toppar.end_fle_6_4=&F_end_fle_4)
evaluate (&toppar.end_fle_6_5=&F_end_fle_5)

evaluate (&data.numncs=&numncs)
evaluate (&toppar.ncs_sta1_1=&ncs_sta1_1)
evaluate (&toppar.ncs_sta1_2=&ncs_sta1_2)
evaluate (&toppar.ncs_sta1_3=&ncs_sta1_3)
evaluate (&toppar.ncs_sta1_4=&ncs_sta1_4)
evaluate (&toppar.ncs_sta1_5=&ncs_sta1_5)
evaluate (&toppar.ncs_end1_1=&ncs_end1_1)
evaluate (&toppar.ncs_end1_2=&ncs_end1_2)
evaluate (&toppar.ncs_end1_3=&ncs_end1_3)
evaluate (&toppar.ncs_end1_4=&ncs_end1_4)
evaluate (&toppar.ncs_end1_5=&ncs_end1_5)
evaluate (&toppar.ncs_seg1_1=&ncs_seg1_1)
evaluate (&toppar.ncs_seg1_2=&ncs_seg1_2)
evaluate (&toppar.ncs_seg1_3=&ncs_seg1_3)
evaluate (&toppar.ncs_seg1_4=&ncs_seg1_4)
evaluate (&toppar.ncs_seg1_5=&ncs_seg1_5)
evaluate (&toppar.ncs_sta2_1=&ncs_sta2_1)
evaluate (&toppar.ncs_sta2_2=&ncs_sta2_2)
evaluate (&toppar.ncs_sta2_3=&ncs_sta2_3)
evaluate (&toppar.ncs_sta2_4=&ncs_sta2_4)
evaluate (&toppar.ncs_sta2_5=&ncs_sta2_5)
evaluate (&toppar.ncs_end2_1=&ncs_end2_1)
evaluate (&toppar.ncs_end2_2=&ncs_end2_2)
evaluate (&toppar.ncs_end2_3=&ncs_end2_3)
evaluate (&toppar.ncs_end2_4=&ncs_end2_4)
evaluate (&toppar.ncs_end2_5=&ncs_end2_5)
evaluate (&toppar.ncs_seg2_1=&ncs_seg2_1)
evaluate (&toppar.ncs_seg2_2=&ncs_seg2_2)
evaluate (&toppar.ncs_seg2_3=&ncs_seg2_3)
evaluate (&toppar.ncs_seg2_4=&ncs_seg2_4)
evaluate (&toppar.ncs_seg2_5=&ncs_seg2_5)

evaluate (&data.numc2sym=&numc2sym)
evaluate (&toppar.c2sym_sta1_1=&c2sym_sta1_1)
evaluate (&toppar.c2sym_sta1_2=&c2sym_sta1_2)
evaluate (&toppar.c2sym_sta1_3=&c2sym_sta1_3)
evaluate (&toppar.c2sym_sta1_4=&c2sym_sta1_4)
evaluate (&toppar.c2sym_sta1_5=&c2sym_sta1_5)
evaluate (&toppar.c2sym_sta1_6=&c2sym_sta1_6)
evaluate (&toppar.c2sym_sta1_7=&c2sym_sta1_7)
evaluate (&toppar.c2sym_sta1_8=&c2sym_sta1_8)
evaluate (&toppar.c2sym_sta1_9=&c2sym_sta1_8)
evaluate (&toppar.c2sym_sta1_10=&c2sym_sta1_10)
evaluate (&toppar.c2sym_end1_1=&c2sym_end1_1)
evaluate (&toppar.c2sym_end1_2=&c2sym_end1_2)
evaluate (&toppar.c2sym_end1_3=&c2sym_end1_3)
evaluate (&toppar.c2sym_end1_4=&c2sym_end1_4)
evaluate (&toppar.c2sym_end1_5=&c2sym_end1_5)
evaluate (&toppar.c2sym_end1_6=&c2sym_end1_6)
evaluate (&toppar.c2sym_end1_7=&c2sym_end1_7)
evaluate (&toppar.c2sym_end1_8=&c2sym_end1_8)
evaluate (&toppar.c2sym_end1_9=&c2sym_end1_9)
evaluate (&toppar.c2sym_end1_10=&c2sym_end1_10)
evaluate (&toppar.c2sym_seg1_1=&c2sym_seg1_1)
evaluate (&toppar.c2sym_seg1_2=&c2sym_seg1_2)
evaluate (&toppar.c2sym_seg1_3=&c2sym_seg1_3)
evaluate (&toppar.c2sym_seg1_4=&c2sym_seg1_4)
evaluate (&toppar.c2sym_seg1_5=&c2sym_seg1_5)
evaluate (&toppar.c2sym_seg1_6=&c2sym_seg1_6)
evaluate (&toppar.c2sym_seg1_7=&c2sym_seg1_7)
evaluate (&toppar.c2sym_seg1_8=&c2sym_seg1_8)
evaluate (&toppar.c2sym_seg1_9=&c2sym_seg1_9)
evaluate (&toppar.c2sym_seg1_10=&c2sym_seg1_10)
evaluate (&toppar.c2sym_sta2_1=&c2sym_sta2_1)
evaluate (&toppar.c2sym_sta2_2=&c2sym_sta2_2)
evaluate (&toppar.c2sym_sta2_3=&c2sym_sta2_3)
evaluate (&toppar.c2sym_sta2_4=&c2sym_sta2_4)
evaluate (&toppar.c2sym_sta2_5=&c2sym_sta2_5)
evaluate (&toppar.c2sym_sta2_6=&c2sym_sta2_6)
evaluate (&toppar.c2sym_sta2_7=&c2sym_sta2_7)
evaluate (&toppar.c2sym_sta2_8=&c2sym_sta2_8)
evaluate (&toppar.c2sym_sta2_9=&c2sym_sta2_9)
evaluate (&toppar.c2sym_sta2_10=&c2sym_sta2_10)
evaluate (&toppar.c2sym_end2_1=&c2sym_end2_1)
evaluate (&toppar.c2sym_end2_2=&c2sym_end2_2)
evaluate (&toppar.c2sym_end2_3=&c2sym_end2_3)
evaluate (&toppar.c2sym_end2_4=&c2sym_end2_4)
evaluate (&toppar.c2sym_end2_5=&c2sym_end2_5)
evaluate (&toppar.c2sym_end2_6=&c2sym_end2_6)
evaluate (&toppar.c2sym_end2_7=&c2sym_end2_7)
evaluate (&toppar.c2sym_end2_8=&c2sym_end2_8)
evaluate (&toppar.c2sym_end2_9=&c2sym_end2_9)
evaluate (&toppar.c2sym_end2_10=&c2sym_end2_10)
evaluate (&toppar.c2sym_seg2_1=&c2sym_seg2_1)
evaluate (&toppar.c2sym_seg2_2=&c2sym_seg2_2)
evaluate (&toppar.c2sym_seg2_3=&c2sym_seg2_3)
evaluate (&toppar.c2sym_seg2_4=&c2sym_seg2_4)
evaluate (&toppar.c2sym_seg2_5=&c2sym_seg2_5)
evaluate (&toppar.c2sym_seg2_6=&c2sym_seg2_6)
evaluate (&toppar.c2sym_seg2_7=&c2sym_seg2_7)
evaluate (&toppar.c2sym_seg2_8=&c2sym_seg2_8)
evaluate (&toppar.c2sym_seg2_9=&c2sym_seg2_9)
evaluate (&toppar.c2sym_seg2_10=&c2sym_seg2_10)

evaluate (&data.numc3sym=&numc3sym)
evaluate (&toppar.c3sym_sta1_1=&c3sym_sta1_1)
evaluate (&toppar.c3sym_sta1_2=&c3sym_sta1_2)
evaluate (&toppar.c3sym_sta2_1=&c3sym_sta2_1)
evaluate (&toppar.c3sym_sta2_2=&c3sym_sta2_2)
evaluate (&toppar.c3sym_sta3_1=&c3sym_sta3_1)
evaluate (&toppar.c3sym_sta3_2=&c3sym_sta3_2)
evaluate (&toppar.c3sym_end1_1=&c3sym_end1_1)
evaluate (&toppar.c3sym_end1_2=&c3sym_end1_2)
evaluate (&toppar.c3sym_end2_1=&c3sym_end2_1)
evaluate (&toppar.c3sym_end2_2=&c3sym_end2_2)
evaluate (&toppar.c3sym_end3_1=&c3sym_end3_1)
evaluate (&toppar.c3sym_end3_2=&c3sym_end3_2)
evaluate (&toppar.c3sym_seg1_1=&c3sym_seg1_1)
evaluate (&toppar.c3sym_seg1_2=&c3sym_seg1_2)
evaluate (&toppar.c3sym_seg2_1=&c3sym_seg2_1)
evaluate (&toppar.c3sym_seg2_2=&c3sym_seg2_2)
evaluate (&toppar.c3sym_seg3_1=&c3sym_seg3_1)
evaluate (&toppar.c3sym_seg3_2=&c3sym_seg3_2)

evaluate (&data.numc4sym=&numc4sym)

evaluate (&data.numc5sym=&numc5sym)
evaluate (&toppar.c5sym_sta1_1=&c5sym_sta1_1)
evaluate (&toppar.c5sym_sta2_1=&c5sym_sta2_1)
evaluate (&toppar.c5sym_sta3_1=&c5sym_sta3_1)
evaluate (&toppar.c5sym_sta4_1=&c5sym_sta4_1)
evaluate (&toppar.c5sym_sta5_1=&c5sym_sta5_1)
evaluate (&toppar.c5sym_end1_1=&c5sym_end1_1)
evaluate (&toppar.c5sym_end2_1=&c5sym_end2_1)
evaluate (&toppar.c5sym_end3_1=&c5sym_end3_1)
evaluate (&toppar.c5sym_end4_1=&c5sym_end4_1)
evaluate (&toppar.c5sym_end5_1=&c5sym_end5_1)
evaluate (&toppar.c5sym_seg1_1=&c5sym_seg1_1)
evaluate (&toppar.c5sym_seg2_1=&c5sym_seg2_1)
evaluate (&toppar.c5sym_seg3_1=&c5sym_seg3_1)
evaluate (&toppar.c5sym_seg4_1=&c5sym_seg4_1)
evaluate (&toppar.c5sym_seg5_1=&c5sym_seg5_1)

evaluate (&toppar.xplortodiana=&xplortodiana)
evaluate (&toppar.delenph=&delenph)

evaluate ($nhisd=&numhisd)
evaluate ($nhise=&numhise)
evaluate ($ncc=1)
while ($ncc <= $nhisd) loop hisd
  evaluate (&toppar.hisd_resid_1_$ncc=&A_hisd_resid_$ncc)
  evaluate (&toppar.hisd_resid_2_$ncc=&B_hisd_resid_$ncc)
  evaluate (&toppar.hisd_resid_3_$ncc=&C_hisd_resid_$ncc)
  evaluate (&toppar.hisd_resid_4_$ncc=&D_hisd_resid_$ncc)
  evaluate (&toppar.hisd_resid_5_$ncc=&E_hisd_resid_$ncc)
  evaluate (&toppar.hisd_resid_6_$ncc=&F_hisd_resid_$ncc)
  evaluate ($ncc = $ncc + 1)
end loop hisd

evaluate ($ncc=1)
while ($ncc <= $nhise) loop hise
  evaluate (&toppar.hise_resid_1_$ncc=&A_hise_resid_$ncc)
  evaluate (&toppar.hise_resid_2_$ncc=&B_hise_resid_$ncc)
  evaluate (&toppar.hise_resid_3_$ncc=&C_hise_resid_$ncc)
  evaluate (&toppar.hise_resid_4_$ncc=&D_hise_resid_$ncc)
  evaluate (&toppar.hise_resid_5_$ncc=&E_hise_resid_$ncc)
  evaluate (&toppar.hise_resid_6_$ncc=&F_hise_resid_$ncc)
  evaluate ($ncc = $ncc + 1)
end loop hise

!Dihedrals, Jcouplings, Residual dipolar couplints, Hbonds, Analysis:
evaluate (&Data.flags.dihed = &dihedflag)
evaluate (&Data.flags.elec0 = &elecflag_0)
evaluate (&Data.flags.elec1 = &elecflag_1)
evaluate (&Data.epsilon = &epsilon)
evaluate (&Data.dielec  = &dielec)
evaluate (&Data.dnarest = &dnarest_on)

evaluate (&Data.flags.noe  =  true)
evaluate (&Data.flags.cdih =  &dihedrals_on)
evaluate (&Data.cdih.on =  &dihedrals_on)

evaluate (&Data.flags.coup =  false)
evaluate (&Data.flags.vean =  false)
if (&rdc1_choice = "VANGLE") then
  evaluate (&Data.flags.vean =  true)
end if
if (&rdc2_choice = "VANGLE") then
  evaluate (&Data.flags.vean =  true)
end if
if (&rdc3_choice = "VANGLE") then
  evaluate (&Data.flags.vean =  true)
end if
if (&rdc4_choice = "VANGLE") then
  evaluate (&Data.flags.vean =  true)
end if
if (&rdc5_choice = "VANGLE") then
  evaluate (&Data.flags.vean =  true)
end if
evaluate (&Data.flags.sani =  false)
if (&rdc1_choice = "SANI") then
  evaluate (&Data.flags.sani =  true)
end if
if (&rdc2_choice = "SANI") then
  evaluate (&Data.flags.sani =  true)
end if
if (&rdc3_choice = "SANI") then
  evaluate (&Data.flags.sani =  true)
end if
if (&rdc4_choice = "SANI") then
  evaluate (&Data.flags.sani =  true)
end if
if (&rdc5_choice = "SANI") then
  evaluate (&Data.flags.sani =  true)
end if
evaluate (&Data.flags.dani =  false)
if (&dan1_choice = "DANI") then
  evaluate (&Data.flags.dani =  true)
end if
if (&dan2_choice = "dani") then
  evaluate (&Data.flags.dani =  true)
end if
if (&dan3_choice = "dani") then
  evaluate (&Data.flags.dani =  true)
end if
if (&dan4_choice = "dani") then
  evaluate (&Data.flags.dani =  true)
end if
if (&dan5_choice = "dani") then
  evaluate (&Data.flags.dani =  true)
end if


evaluate (&Data.flags.plan =  false)
evaluate (&Data.flags.ncs  =  &ncs_on)
evaluate (&Data.flags.sym  =  &sym_on)

evaluate (&data.scaling=&air_scaling)
evaluate (&data.totnoe_unamb=&tot_unamb)
evaluate (&data.unamb_firstit=&unamb_firstit)
evaluate (&data.unamb_lastit=&unamb_lastit)
evaluate (&data.unamb_hot=&unamb_hot)
evaluate (&data.unamb_cool1=&unamb_cool1)
evaluate (&data.unamb_cool2=&unamb_cool2)
evaluate (&data.unamb_cool3=&unamb_cool3)
evaluate (&data.noecv=&noecv)
evaluate (&data.ncvpart=&ncvpart)
evaluate (&data.ranair=&ranair)
if (&data.ranair eq true) then
  evaluate (&data.noecv = false)
end if
evaluate (&data.cmrest=&cmrest)
evaluate (&data.kcont=&kcont)
evaluate (&data.surfrest=&surfrest)
evaluate (&data.ksurf=&ksurf)

evaluate (&data.totnoe_amb=&tot_amb)
evaluate (&data.amb_firstit=&amb_firstit)
evaluate (&data.amb_lastit=&amb_lastit)
evaluate (&data.amb_hot=&amb_hot)
evaluate (&data.amb_cool1=&amb_cool1)
evaluate (&data.amb_cool2=&amb_cool2)
evaluate (&data.amb_cool3=&amb_cool3)

evaluate (&data.kncs=&kncs)
evaluate (&data.ksym=&ksym)

evaluate (&data.hbond_firstit=&hbond_firstit)
evaluate (&data.hbond_lastit=&hbond_lastit)
evaluate (&data.hbond_hot=&hbond_hot)
evaluate (&data.hbond_cool1=&hbond_cool1)
evaluate (&data.hbond_cool2=&hbond_cool2)
evaluate (&data.hbond_cool3=&hbond_cool3)

evaluate (&data.mrswi_hot=&mrswi_hot)
evaluate (&data.mrswi_cool1=&mrswi_cool1)
evaluate (&data.mrswi_cool2=&mrswi_cool2)
evaluate (&data.mrswi_cool3=&mrswi_cool3)

evaluate (&data.rswi_hot=&rswi_hot)
evaluate (&data.rswi_cool1=&rswi_cool1)
evaluate (&data.rswi_cool2=&rswi_cool2)
evaluate (&data.rswi_cool3=&rswi_cool3)

evaluate (&data.masy_hot=&masy_hot)
evaluate (&data.masy_cool1=&masy_cool1)
evaluate (&data.masy_cool2=&masy_cool2)
evaluate (&data.masy_cool3=&masy_cool3)

evaluate (&data.asy_hot=&asy_hot)
evaluate (&data.asy_cool1=&asy_cool1)
evaluate (&data.asy_cool2=&asy_cool2)
evaluate (&data.asy_cool3=&asy_cool3)

evaluate (&data.dihedrals.on=&dihedrals_on)
evaluate (&data.dihedrals_hot=&dihedrals_hot)
evaluate (&data.dihedrals_cool1=&dihedrals_cool1)
evaluate (&data.dihedrals_cool2=&dihedrals_cool2)
evaluate (&data.dihedrals_cool3=&dihedrals_cool3)
evaluate (&data.hbonds_on=&hbonds_on)

evaluate (&data.c1_on=&c1_on)
evaluate (&data.c1_karplusa=&c1_karplusa)
evaluate (&data.c1_karplusb=&c1_karplusb)
evaluate (&data.c1_karplusc=&c1_karplusc)
evaluate (&data.c1_karplusd=&c1_karplusd)
evaluate (&data.c1_hot=&c1_hot)
evaluate (&data.c1_cool1=&c1_cool1)
evaluate (&data.c1_cool2=&c1_cool2)
evaluate (&data.c1_cool3=&c1_cool3)
evaluate (&data.c2_on=&c2_on)
evaluate (&data.c2_karplusa=&c2_karplusa)
evaluate (&data.c2_karplusb=&c2_karplusb)
evaluate (&data.c2_karplusc=&c2_karplusc)
evaluate (&data.c2_karplusd=&c2_karplusd)
evaluate (&data.c2_hot=&c2_hot)
evaluate (&data.c2_cool1=&c2_cool1)
evaluate (&data.c2_cool2=&c2_cool2)
evaluate (&data.c2_cool3=&c2_cool3)
evaluate (&data.c3_on=&c3_on)
evaluate (&data.c3_karplusa=&c3_karplusa)
evaluate (&data.c3_karplusb=&c3_karplusb)
evaluate (&data.c3_karplusc=&c3_karplusc)
evaluate (&data.c3_karplusd=&c3_karplusd)
evaluate (&data.c3_hot=&c3_hot)
evaluate (&data.c3_cool1=&c3_cool1)
evaluate (&data.c3_cool2=&c3_cool2)
evaluate (&data.c3_cool3=&c3_cool3)
evaluate (&data.c4_on=&c4_on)
evaluate (&data.c4_karplusa=&c4_karplusa)
evaluate (&data.c4_karplusb=&c4_karplusb)
evaluate (&data.c4_karplusc=&c4_karplusc)
evaluate (&data.c4_karplusd=&c4_karplusd)
evaluate (&data.c4_hot=&c4_hot)
evaluate (&data.c4_cool1=&c4_cool1)
evaluate (&data.c4_cool2=&c4_cool2)
evaluate (&data.c4_cool3=&c4_cool3)
evaluate (&data.c5_on=&c5_on)
evaluate (&data.c5_karplusa=&c5_karplusa)
evaluate (&data.c5_karplusb=&c5_karplusb)
evaluate (&data.c5_karplusc=&c5_karplusc)
evaluate (&data.c5_karplusd=&c5_karplusd)
evaluate (&data.c5_hot=&c5_hot)
evaluate (&data.c5_cool1=&c5_cool1)
evaluate (&data.c5_cool2=&c5_cool2)
evaluate (&data.c5_cool3=&c5_cool3)

evaluate (&data.rdc1_choice=&rdc1_choice)
evaluate (&data.rdc1_firstIt=&rdc1_firstIt)
evaluate (&data.rdc1_lastIt=&rdc1_lastIt)
evaluate (&data.rdc1_hot=&rdc1_hot)
evaluate (&data.rdc1_cool1=&rdc1_cool1)
evaluate (&data.rdc1_cool2=&rdc1_cool2)
evaluate (&data.rdc1_cool3=&rdc1_cool3)
evaluate (&data.rdc1_r=&rdc1_r)
evaluate (&data.rdc1_d=&rdc1_d)

evaluate (&data.rdc2_choice=&rdc2_choice)
evaluate (&data.rdc2_firstIt=&rdc2_firstIt)
evaluate (&data.rdc2_lastIt=&rdc2_lastIt)
evaluate (&data.rdc2_hot=&rdc2_hot)
evaluate (&data.rdc2_cool1=&rdc2_cool1)
evaluate (&data.rdc2_cool2=&rdc2_cool2)
evaluate (&data.rdc2_cool3=&rdc2_cool3)
evaluate (&data.rdc2_r=&rdc2_r)
evaluate (&data.rdc2_d=&rdc2_d)

evaluate (&data.rdc3_choice=&rdc3_choice)
evaluate (&data.rdc3_firstIt=&rdc3_firstIt)
evaluate (&data.rdc3_lastIt=&rdc3_lastIt)
evaluate (&data.rdc3_hot=&rdc3_hot)
evaluate (&data.rdc3_cool1=&rdc3_cool1)
evaluate (&data.rdc3_cool2=&rdc3_cool2)
evaluate (&data.rdc3_cool3=&rdc3_cool3)
evaluate (&data.rdc3_r=&rdc3_r)
evaluate (&data.rdc3_d=&rdc3_d)

evaluate (&data.rdc4_choice=&rdc4_choice)
evaluate (&data.rdc4_firstIt=&rdc4_firstIt)
evaluate (&data.rdc4_lastIt=&rdc4_lastIt)
evaluate (&data.rdc4_hot=&rdc4_hot)
evaluate (&data.rdc4_cool1=&rdc4_cool1)
evaluate (&data.rdc4_cool2=&rdc4_cool2)
evaluate (&data.rdc4_cool3=&rdc4_cool3)
evaluate (&data.rdc4_r=&rdc4_r)
evaluate (&data.rdc4_d=&rdc4_d)

evaluate (&data.rdc5_choice=&rdc5_choice)
evaluate (&data.rdc5_firstIt=&rdc5_firstIt)
evaluate (&data.rdc5_lastIt=&rdc5_lastIt)
evaluate (&data.rdc5_hot=&rdc5_hot)
evaluate (&data.rdc5_cool1=&rdc5_cool1)
evaluate (&data.rdc5_cool2=&rdc5_cool2)
evaluate (&data.rdc5_cool3=&rdc5_cool3)
evaluate (&data.rdc5_r=&rdc5_r)
evaluate (&data.rdc5_d=&rdc5_d)

evaluate (&data.dan1_choice=&dan1_choice)
evaluate (&data.dan1_firstIt=&dan1_firstIt)
evaluate (&data.dan1_lastIt=&dan1_lastIt)
evaluate (&data.dan1_hot=&dan1_hot)
evaluate (&data.dan1_cool1=&dan1_cool1)
evaluate (&data.dan1_cool2=&dan1_cool2)
evaluate (&data.dan1_cool3=&dan1_cool3)
evaluate (&data.dan1_tc=&dan1_tc)
evaluate (&data.dan1_anis=&dan1_anis)
evaluate (&data.dan1_r=&dan1_r)
evaluate (&data.dan1_wh=&dan1_wh)
evaluate (&data.dan1_wn=&dan1_wn)

evaluate (&data.dan2_choice=&dan2_choice)
evaluate (&data.dan2_firstIt=&dan2_firstIt)
evaluate (&data.dan2_lastIt=&dan2_lastIt)
evaluate (&data.dan2_hot=&dan2_hot)
evaluate (&data.dan2_cool1=&dan2_cool1)
evaluate (&data.dan2_cool2=&dan2_cool2)
evaluate (&data.dan2_cool3=&dan2_cool3)
evaluate (&data.dan2_tc=&dan2_tc)
evaluate (&data.dan2_anis=&dan2_anis)
evaluate (&data.dan2_r=&dan2_r)
evaluate (&data.dan2_wh=&dan2_wh)
evaluate (&data.dan2_wn=&dan2_wn)

evaluate (&data.dan3_choice=&dan3_choice)
evaluate (&data.dan3_firstIt=&dan3_firstIt)
evaluate (&data.dan3_lastIt=&dan3_lastIt)
evaluate (&data.dan3_hot=&dan3_hot)
evaluate (&data.dan3_cool1=&dan3_cool1)
evaluate (&data.dan3_cool2=&dan3_cool2)
evaluate (&data.dan3_cool3=&dan3_cool3)
evaluate (&data.dan3_tc=&dan3_tc)
evaluate (&data.dan3_anis=&dan3_anis)
evaluate (&data.dan3_r=&dan3_r)
evaluate (&data.dan3_wh=&dan3_wh)
evaluate (&data.dan3_wn=&dan3_wn)

evaluate (&data.dan4_choice=&dan4_choice)
evaluate (&data.dan4_firstIt=&dan4_firstIt)
evaluate (&data.dan4_lastIt=&dan4_lastIt)
evaluate (&data.dan4_hot=&dan4_hot)
evaluate (&data.dan4_cool1=&dan4_cool1)
evaluate (&data.dan4_cool2=&dan4_cool2)
evaluate (&data.dan4_cool3=&dan4_cool3)
evaluate (&data.dan4_tc=&dan4_tc)
evaluate (&data.dan4_anis=&dan4_anis)
evaluate (&data.dan4_r=&dan4_r)
evaluate (&data.dan4_wh=&dan4_wh)
evaluate (&data.dan4_wn=&dan4_wn)

evaluate (&data.dan5_choice=&dan5_choice)
evaluate (&data.dan5_firstIt=&dan5_firstIt)
evaluate (&data.dan5_lastIt=&dan5_lastIt)
evaluate (&data.dan5_hot=&dan5_hot)
evaluate (&data.dan5_cool1=&dan5_cool1)
evaluate (&data.dan5_cool2=&dan5_cool2)
evaluate (&data.dan5_cool3=&dan5_cool3)
evaluate (&data.dan5_tc=&dan5_tc)
evaluate (&data.dan5_anis=&dan5_anis)
evaluate (&data.dan5_r=&dan5_r)
evaluate (&data.dan5_wh=&dan5_wh)
evaluate (&data.dan5_wn=&dan5_wn)


!VEAN statement:
evaluate (&data.ini_bor_hot_1=&ini_bor_hot_1)
evaluate (&data.ini_bor_cool1_1=&ini_bor_cool1_1)
evaluate (&data.ini_bor_cool2_1=&ini_bor_cool2_1)
evaluate (&data.ini_bor_cool3_1=&ini_bor_cool3_1)
evaluate (&data.ini_cen_hot_1=&ini_cen_hot_1)
evaluate (&data.ini_cen_cool1_1=&ini_cen_cool1_1)
evaluate (&data.ini_cen_cool2_1=&ini_cen_cool2_1)
evaluate (&data.ini_cen_cool3_1=&ini_cen_cool3_1)
evaluate (&data.fin_bor_hot_1=&fin_bor_hot_1)
evaluate (&data.fin_bor_cool1_1=&fin_bor_cool1_1)
evaluate (&data.fin_bor_cool2_1=&fin_bor_cool2_1)
evaluate (&data.fin_bor_cool3_1=&fin_bor_cool3_1)
evaluate (&data.fin_cen_hot_1=&fin_cen_hot_1)
evaluate (&data.fin_cen_cool1_1=&fin_cen_cool1_1)
evaluate (&data.fin_cen_cool2_1=&fin_cen_cool2_1)
evaluate (&data.fin_cen_cool3_1=&fin_cen_cool3_1)

evaluate (&data.ini_bor_hot_2=&ini_bor_hot_2)
evaluate (&data.ini_bor_cool1_2=&ini_bor_cool1_2)
evaluate (&data.ini_bor_cool2_2=&ini_bor_cool2_2)
evaluate (&data.ini_bor_cool3_2=&ini_bor_cool3_2)
evaluate (&data.ini_cen_hot_2=&ini_cen_hot_2)
evaluate (&data.ini_cen_cool1_2=&ini_cen_cool1_2)
evaluate (&data.ini_cen_cool2_2=&ini_cen_cool2_2)
evaluate (&data.ini_cen_cool3_2=&ini_cen_cool3_2)
evaluate (&data.fin_bor_hot_2=&fin_bor_hot_2)
evaluate (&data.fin_bor_cool1_2=&fin_bor_cool1_2)
evaluate (&data.fin_bor_cool2_2=&fin_bor_cool2_2)
evaluate (&data.fin_bor_cool3_2=&fin_bor_cool3_2)
evaluate (&data.fin_cen_hot_2=&fin_cen_hot_2)
evaluate (&data.fin_cen_cool1_2=&fin_cen_cool1_2)
evaluate (&data.fin_cen_cool2_2=&fin_cen_cool2_2)
evaluate (&data.fin_cen_cool3_2=&fin_cen_cool3_2)


evaluate (&data.ini_bor_hot_3=&ini_bor_hot_3)
evaluate (&data.ini_bor_cool1_3=&ini_bor_cool1_3)
evaluate (&data.ini_bor_cool2_3=&ini_bor_cool2_3)
evaluate (&data.ini_bor_cool3_3=&ini_bor_cool3_3)
evaluate (&data.ini_cen_hot_3=&ini_cen_hot_3)
evaluate (&data.ini_cen_cool1_3=&ini_cen_cool1_3)
evaluate (&data.ini_cen_cool2_3=&ini_cen_cool2_3)
evaluate (&data.ini_cen_cool3_3=&ini_cen_cool3_3)
evaluate (&data.fin_bor_hot_3=&fin_bor_hot_3)
evaluate (&data.fin_bor_cool1_3=&fin_bor_cool1_3)
evaluate (&data.fin_bor_cool2_3=&fin_bor_cool2_3)
evaluate (&data.fin_bor_cool3_3=&fin_bor_cool3_3)
evaluate (&data.fin_cen_hot_3=&fin_cen_hot_3)
evaluate (&data.fin_cen_cool1_3=&fin_cen_cool1_3)
evaluate (&data.fin_cen_cool2_3=&fin_cen_cool2_3)
evaluate (&data.fin_cen_cool3_3=&fin_cen_cool3_3)

evaluate (&data.ini_bor_hot_4=&ini_bor_hot_4)
evaluate (&data.ini_bor_cool1_4=&ini_bor_cool1_4)
evaluate (&data.ini_bor_cool2_4=&ini_bor_cool2_4)
evaluate (&data.ini_bor_cool3_4=&ini_bor_cool3_4)
evaluate (&data.ini_cen_hot_4=&ini_cen_hot_4)
evaluate (&data.ini_cen_cool1_4=&ini_cen_cool1_4)
evaluate (&data.ini_cen_cool2_4=&ini_cen_cool2_4)
evaluate (&data.ini_cen_cool3_4=&ini_cen_cool3_4)
evaluate (&data.fin_bor_hot_4=&fin_bor_hot_4)
evaluate (&data.fin_bor_cool1_4=&fin_bor_cool1_4)
evaluate (&data.fin_bor_cool2_4=&fin_bor_cool2_4)
evaluate (&data.fin_bor_cool3_4=&fin_bor_cool3_4)
evaluate (&data.fin_cen_hot_4=&fin_cen_hot_4)
evaluate (&data.fin_cen_cool1_4=&fin_cen_cool1_4)
evaluate (&data.fin_cen_cool2_4=&fin_cen_cool2_4)
evaluate (&data.fin_cen_cool3_4=&fin_cen_cool3_4)

evaluate (&data.ini_bor_hot_5=&ini_bor_hot_5)
evaluate (&data.ini_bor_cool1_5=&ini_bor_cool1_5)
evaluate (&data.ini_bor_cool2_5=&ini_bor_cool2_5)
evaluate (&data.ini_bor_cool3_5=&ini_bor_cool3_5)
evaluate (&data.ini_cen_hot_5=&ini_cen_hot_5)
evaluate (&data.ini_cen_cool1_5=&ini_cen_cool1_5)
evaluate (&data.ini_cen_cool2_5=&ini_cen_cool2_5)
evaluate (&data.ini_cen_cool3_5=&ini_cen_cool3_5)
evaluate (&data.fin_bor_hot_5=&fin_bor_hot_5)
evaluate (&data.fin_bor_cool1_5=&fin_bor_cool1_5)
evaluate (&data.fin_bor_cool2_5=&fin_bor_cool2_5)
evaluate (&data.fin_bor_cool3_5=&fin_bor_cool3_5)
evaluate (&data.fin_cen_hot_5=&fin_cen_hot_5)
evaluate (&data.fin_cen_cool1_5=&fin_cen_cool1_5)
evaluate (&data.fin_cen_cool2_5=&fin_cen_cool2_5)
evaluate (&data.fin_cen_cool3_5=&fin_cen_cool3_5)

!intermolecular contacts analysis
evaluate (&data.hb_dist=&dist_hb)
evaluate (&data.nb_dist=&dist_nb)

!water refinement
evaluate (&refine.firstwater=&firstwater)
evaluate (&refine.keepwater=&keepwater)
evaluate (&refine.waterrefine=&waterrefine)
evaluate (&refine.solvent=&solvent)
evaluate (&refine.pmrefine_on=&pmrefine_on)
evaluate (&refine.calcdesolv=&calcdesolv)
evaluate (&refine.heatsteps=&waterheatsteps)
evaluate (&refine.steps=&watersteps)
evaluate (&refine.coolsteps=&watercoolsteps)


!for the non-bonded parameters (the section was taken out of 
!parallhdg5.0.pro and parallhdg5.1.pro, so be careful!):
if (&toppar.par_nonbonded eq "PROLSQ") then
    evaluate (&toppar.repel_radius = 1.0)
    evaluate (&toppar.repel_rcons = 20)
    evaluate (&toppar.repel_rexpo  = 4)
    evaluate (&toppar.repel_irexp  = 1)
elseif (&toppar.par_nonbonded eq "PARMALLH6") then
    evaluate (&toppar.repel_radius = 0.8)
    evaluate (&toppar.repel_rcons = 5.0)
    evaluate (&toppar.repel_rexpo  = 2)
    evaluate (&toppar.repel_irexp  = 2)
elseif (&toppar.par_nonbonded eq "OPLSX") then
    evaluate (&toppar.repel_radius = 0.0)
else        {...now the standard PARALLHDG parameters}
    evaluate (&toppar.repel_radius = 0.78)
    evaluate (&toppar.repel_rcons = 5.0)
    evaluate (&toppar.repel_rexpo  = 2)
    evaluate (&toppar.repel_irexp  = 2)
end if

!Procheck analysis:
evaluate (&analysis.procheckdir=&procheckdir)

! Water in rigid body docking
evaluate (&data.waterdock=&waterdock)
evaluate (&data.water_tokeep=&water_tokeep)
evaluate (&data.water_randfrac=&water_randfrac)
evaluate (&data.solvate_method=&solvate_method)
evaluate (&data.water_surfcutoff=&water_surfcutoff)
evaluate (&data.water_analysis=&water_analysis)
evaluate (&data.transwater=&transwater)
evaluate (&data.water_restraint_initial=&water_restraint_initial)
evaluate (&data.water_restraint_cutoff=&water_restraint_cutoff)
evaluate (&data.water_restraint_scale=&water_restraint_scale)
evaluate (&data.waterensemble=&waterensemble)

if (&data.waterdock eq true) then
  evaluate (&iterations.rotate180   = false)
  evaluate (&SaProtocol.initiosteps = 0)
  evaluate (&SaProtocol.cool1_steps = 0)
  display SOLVATED DOCKING TURNED ON: initiosteps and cool1_steps set to 0, rotate180 set to false
end if
"""
