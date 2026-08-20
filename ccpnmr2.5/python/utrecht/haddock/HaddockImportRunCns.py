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

import os, sys, re
from HaddockLocal import *

class runCnsImporter:
    
    """Description: Import run specific parameters and their settings from a Haddock run.cns file.
       Input:        A valid Haddock run.cns parameter file and the CCPN Haddock project under witch
                    to create the new run.
       Output:        A new haddock run in the current project.                
    """
    
    def __init__(self,hProject=None,runCns=None):
        
        self.hProject = hProject
        self.newRun    = None
        self.rawparams = {}
        self.formattedparams = {}
        self.boolTrue = ['true','"true"','yes','"db"']
        self.boolFalse = ['false','"false"','no','""']
        self.energyStores = ['dockingProtocolStore','distRestraintEnergyStore','dihRestraintEnergyStore','scoringWeights',
                             'semiflexInterMolScalingStore','autoDistanceRestraintWeightStore','rdcProtocolStore']
        self.parameters = {'delenph': 'removeNonPolarH', 'clust_rmsd': 'analysisClustRmsd', 'noecv': 'randomExcludeAir', 
                           'structures_1': 'numIt1Structures', 'structures_0': 'numIt0Structures', 'skip_struc': 'skipStructures', 
                           'water_surfcutoff': 'waterSurfaceCutoff', 'solvent': 'solvent', 'dielec': 'dielectricType', 
                           'rigidtrans': 'doRigidTranslations', 'ncvpart': 'randomExclParts','air_scaling':'doAirScaling', 
                           'water_restraint_scale': 'waterRestScale', 'waterrefine': 'numWrefStructures', 
                           'ksurf': 'surfaceContactConstant', 'anastruc_1': 'numAnalysisStructures','tot_unamb':'numUnambRestautoAir', 
                           'cmrest': 'centerOfMassRestraints', 'surfrest': 'surfaceContactRestraints','tot_amb':'numAmbRestautoAir', 
                           'solvate_method': 'useDbSolvateMethod', 'dihedflag': 'doIncludeDihEnergy', 
                           'dist_nb': 'analysisDistNonbond', 'transwater': 'doRigidBodyWaterTrans', 
                           'water_restraint_initial': 'waterInitRestCutoff', 'dnarest_on': 'useDnaRestraints', 
                           'elecflag_0': 'doRigidBodyElectrostatics', 'elecflag_1': 'doSAElectrostatics', 
                           'par_nonbonded': 'nonBondedType', 'rigidmini': 'initialRigidBodyMinim', 
                           'water_analysis': 'doWaterAnalysis', 'epsilon': 'epsilon', 'ntrials': 'nTrails', 
                           'waterensemble': 'numInitWaterShells', 'water_restraint_cutoff': 'waterRestCutoff', 
                           'randorien': 'radomizeStartOriention', 'inter_rigid': 'rigidbodyIMinteractScaling', 
                           'ranair': 'randomAmbigRestraints', 'calcdesolv': 'calcDesolvation', 'dist_hb': 'analysisDistHBond', 
                           'water_randfrac': 'waterToAddRandom', 'water_tokeep': 'waterToKeep', 'rotate180_0': 'rotate180It0', 
                           'rotate180_1': 'rotate180It1', 'kcont': 'centerOfMassConstant', 'clust_size': 'analysisClustSize', 
                           'iniseed': 'randomSeed', 'waterdock': 'doWaterDock','hbonds_on':'useHBondRestraints'}
        
        print("** Import Haddock run.cns parameter file: %s **" % runCns)                
        self.readRunCns(runCns)
        if not len(self.rawparams) == 0:
            self.formatSingleParams()
            self.formatEnergyStores()
            if not len(self.formattedparams) == 0:
                self.createNewRun()
                print("** Run.cns import finished **")
            else:
                print(" - Error while importing run.cns nothing imported")    
        else:
            print(" - Error while importing run.cns nothing imported")
    
    def createNewRun(self):
        
        run = self.hProject.newRun()
        self.newRun = run
        termId = len([ec.code for ec in run.sortedHaddockEnergyTerms()])+1
        print(" - Create new run with number: %i" % run.serial)
        
        for param in self.formattedparams:
            if not param in self.energyStores:
                run.set(param,self.formattedparams[param])
            elif param == 'scoringWeights':
                scoringWeights = [(sw.term, sw.stage, sw) for sw in run.scoringWeights]
                if scoringWeights:
                    scoringWeights.sort()
                    for term, stage, scoringWeight in scoringWeights:
                        scoringWeight = self.formattedParams['scoringWeights'][term][stage]
                else:
                    for term in self.formattedparams['scoringWeights']:
                        for stage in range(len(self.formattedparams['scoringWeights'][term])):
                            scoringWeight = run.newScoringWeight(term=term,stage=stage,value=self.formattedparams['scoringWeights'][term][stage])    
            elif param == 'rdcProtocolStore':
                for termId in self.formattedparams[param]:
                    energyTermStore = run.newHaddockEnergyTerm(code='rdcProtocolStore',termId=termId)
                    terms = self.formattedparams[param][termId].keys()
                    for term in terms:
                        energyTerm = energyTermStore.newEnergyTermParameter(code=term,value=self.formattedparams[param][termId][term])
            else:
                protocol = run.findFirstHaddockEnergyTerm(code=param)
                if protocol:
                    for term in protocol.sortedEnergyTermParameters():
                        if term.code in self.formattedparams[param]: term.value = self.formattedparams[param][term.code]
                else:
                    energyTermStore = run.newHaddockEnergyTerm(code=param,termId=termId)
                    terms = self.formattedparams[param].keys()
                    terms.sort()
                    for term in terms:
                        energyTerm = energyTermStore.newEnergyTermParameter(code=term,value=self.formattedparams[param][term])
                    termId += 1
    
    def processLine(self,line):
        
        splitline = line[1].split('=')
        code = splitline[0].strip()
        value = splitline[1][0:-1].strip()
        
        convert = False
        if convert == False:
            try: 
                value = int(value)
                convert = True
            except: pass
        
        if convert == False:
            try:
                value = float(value)
                convert = True
            except: pass
        
        if convert == False:
            if value in self.boolTrue: 
                value = True
                convert = True
            elif value in self.boolFalse: 
                value = False
                convert = True
            else: pass    
        
        if convert == False:
            if value[0] == '"' and value[-1] == '"': 
                value = value[1:-1]
                convert = True
    
        self.rawparams[code] = value                    
    
    def readRunCns(self,runcns):
        
        param = re.compile('{===>}')
        runfile = open(runcns,'r')
        for line in runfile.readlines():
            line = line.strip()
            if param.match(line): self.processLine(param.split(line))

    def formatSingleParams(self):
        
        for param in self.parameters:
            if self.rawparams.has_key(param): self.formattedparams[self.parameters[param]] = self.rawparams[param]
            else: print(" - Parameter %s not found in run.cns. Use model default" % param)
                    
    def formatEnergyStores(self):
        
        for energyStore in self.energyStores:
            if energyStore == 'rdcProtocolStore':    
                formatted = {}; count = 1
                while not count == 6:
                    if self.rawparams.has_key('rdc%i_choice' % count):
                        if self.rawparams['rdc%i_choice' % count] == 'SANI' or self.rawparams['rdc%i_choice' % count] == 'VANGLE':
                            rdc = {}; choise = {'SANI':1.0,'VANGLE':2.0}
                            rdc['rdc_choice'] = choise[self.rawparams['rdc%i_choice' % count]]
                            rdc['rdc_firstIt'] = float(self.rawparams['rdc%i_firstIt' % count])
                            rdc['rdc_lastIt'] = float(self.rawparams['rdc%i_lastIt' % count]) 
                            rdc['rdc_hot'] = float(self.rawparams['rdc%i_hot' % count]) 
                            rdc['rdc_cool1'] = float(self.rawparams['rdc%i_cool1' % count]) 
                            rdc['rdc_cool2'] = float(self.rawparams['rdc%i_cool2' % count]) 
                            rdc['rdc_cool3'] = float(self.rawparams['rdc%i_cool3' % count])
                            rdc['rdc_r'] = float(self.rawparams['rdc%i_r' % count]) 
                            rdc['rdc_d'] = float(self.rawparams['rdc%i_d' % count]) 
                            rdc['ini_bor_hot'] = float(self.rawparams['ini_bor_hot_%i' % count])
                            rdc['fin_bor_hot'] = float(self.rawparams['fin_bor_hot_%i' % count])
                            rdc['ini_bor_cool1'] = float(self.rawparams['ini_bor_cool1_%i' % count])
                            rdc['fin_bor_cool1'] = float(self.rawparams['fin_bor_cool1_%i' % count])
                            rdc['ini_bor_cool2'] = float(self.rawparams['ini_bor_cool2_%i' % count])
                            rdc['fin_bor_cool2'] = float(self.rawparams['fin_bor_cool2_%i' % count])
                            rdc['ini_bor_cool3'] = float(self.rawparams['ini_bor_cool3_%i' % count])
                            rdc['fin_bor_cool3'] = float(self.rawparams['fin_bor_cool3_%i' % count])
                            rdc['ini_cen_hot'] = float(self.rawparams['ini_cen_hot_%i' % count])
                            rdc['fin_cen_hot'] = float(self.rawparams['fin_cen_hot_%i' % count])
                            rdc['ini_cen_cool1'] = float(self.rawparams['ini_cen_cool1_%i' % count])
                            rdc['fin_cen_cool1'] = float(self.rawparams['fin_cen_cool1_%i' % count])
                            rdc['ini_cen_cool2'] = float(self.rawparams['ini_cen_cool2_%i' % count])
                            rdc['fin_cen_cool2'] = float(self.rawparams['fin_cen_cool2_%i' % count])
                            rdc['ini_cen_cool3'] = float(self.rawparams['ini_cen_cool3_%i' % count])
                            rdc['fin_cen_cool3'] = float(self.rawparams['fin_cen_cool3_%i' % count])
                            formatted[count] = rdc
                    count += 1

                self.formattedparams['rdcProtocolStore'] = formatted
            elif energyStore == 'scoringWeights':
                formatted = {}
                for param in DEFAULT_SCORE:
                    stagelist = []
                    for stage in [0,1,2]: 
                        key = 'w_'+param+'_'+str(stage)
                        if self.rawparams.has_key(key): stagelist.append(float(self.rawparams[key]))
                        else:
                            stagelist.append(DEFAULT_SCORE[param][stage])
                            print(" - Parameter %s not found on run.cns. Use default value %1.3f" % (param,DEFAULT_SCORE[param][stage]))
                    formatted[param] = stagelist
                    
                self.formattedparams['scoringWeights'] = formatted
            else:    
                formatted = {}
                terms = eval(energyStore)['terms']
                for param in terms: formatted[param] = terms[param]
                for param in formatted:
                    if self.rawparams.has_key(param): formatted[param] = float(self.rawparams[param])
                    else: print(" - Parameter %s not found on run.cns. Use default value %1.3f" % (param,formatted[param]))

                self.formattedparams['distRestraintEnergyStore'] = formatted        
            