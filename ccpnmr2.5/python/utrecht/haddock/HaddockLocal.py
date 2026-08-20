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

"""
The 'HaddockLocal' file contains various environmental variables that do
not have a storage location in the CCPN data model. They are however
necessary for the propper functioning of HADDOCK when run in a local
server environment.

The various 'EneryStore' and 'ProtocolStore' dictionaries are necessary 
for the JIT initiation of the parameters when not yet stored in the model.
Onces stored they can be accessed as every other parameter in the model. 
"""

# Default protein-protein scoring weights
DEFAULT_SCORE = {'vdw':(0.01,1.0,1.0), 'elec':(1.0,1.0,0.2),
                 'dist':(0.01,0.1,0.1), 'sani':(0.1,0.1,0.1),
                 'dani':(0.01,0.1,0.1), 'vean':(0.1,0.1,0.1), 
                 'cdih':(0.0,0.0,0.0), 'sym':(0.1,0.1,0.1), 
                 'bsa':(-0.01,-0.01,0.0), 'deint':(0.0,0.0,0.0), 
                 'desolv':(1.0,1.0,1.0)}

# Distance restraints energy constants
distRestraintEnergyStore = {'details':'Distance restraints energy constants',
                            'terms':{    
                            'unamb_firstit':0.0,
                            'unamb_lastit':2.0, 
                            'unamb_hot':10.0, 
                            'unamb_cool1':10.0, 
                            'unamb_cool2':50.0, 
                            'unamb_cool3':50.0, 
                            'amb_firstit':0.0, 
                            'amb_lastit':2.0, 
                            'amb_hot':10.0, 
                            'amb_cool1':10.0, 
                            'amb_cool2':50.0, 
                            'amb_cool3':50.0, 
                            'hbond_firstit':1.0,
                            'hbond_lastit':2.0, 
                            'hbond_hot':10.0, 
                            'hbond_cool1':10.0, 
                            'hbond_cool2':50.0,
                            'hbond_cool3':50.0}}

# Dihedral restraint energy constants 
dihRestraintEnergyStore = {'details':'Dihedral restraint energy constants',
                           'terms':{
                           'dihedrals_hot':5.0, 
                           'dihedrals_cool1':5.0, 
                           'dihedrals_cool2':50.0, 
                           'dihedrals_cool3':200.0}}

# Scaling of intermolecular interactions for semi-flexible SA
semiflexInterMolScalingStore = {'details':'Scaling of intermolecular interactions for semi-flexible SA',
                                'terms':{
                                'init_rigid':0.001,
                                'fin_rigid':0.001,
                                'init_cool2':0.001,
                                'fin_cool2':1.0,
                                'init_cool3':0.05,
                                'fin_cool3':1.0}}

# Automated distance restraints weighting
autoDistanceRestraintWeightStore = {'details':'Use automated distance restraints weighting',
                                    'terms':{
                                    'mrswi_hot':0.5,
                                    'mrswi_cool1':0.5,
                                    'mrswi_cool2':0.5,
                                    'mrswi_cool3':0.5,
                                    'rswi_hot':0.5,
                                    'rswi_cool1':0.5,
                                    'rswi_cool2':0.5,
                                    'rswi_cool3':0.5,
                                    'masy_hot':-1.0,
                                    'masy_cool1':-1.0,
                                    'masy_cool2':-0.1,
                                    'masy_cool3':-0.1,
                                    'asy_hot':1.0,
                                    'asy_cool1':1.0,
                                    'asy_cool2':0.1,
                                    'asy_cool3':0.1,}}

# RDC parameter settings
rdcProtocolStore = {'terms':{
                    'rdc_choice':1.0, 
                    'rdc_firstIt':2.0, 
                    'rdc_lastIt':2.0, 
                    'rdc_hot':0.001, 
                    'rdc_cool1':0.02,
                    'rdc_cool2':0.2, 
                    'rdc_cool3':0.2, 
                    'rdc_r':0.057, 
                    'rdc_d':-11.49, 
                    'ini_bor_hot':1.0, 
                    'fin_bor_hot':10.0, 
                    'ini_bor_cool1':10.0, 
                    'fin_bor_cool1':40.0, 
                    'ini_bor_cool2':40.0, 
                    'fin_bor_cool2':40.0, 
                    'ini_bor_cool3':40.0, 
                    'fin_bor_cool3':40.0, 
                    'ini_cen_hot':0.25, 
                    'fin_cen_hot':2.5, 
                    'ini_cen_cool1':2.5, 
                    'fin_cen_cool1':10.0, 
                    'ini_cen_cool2':10.0, 
                    'fin_cen_cool2':10.0, 
                    'ini_cen_cool3':10.0, 
                    'fin_cen_cool3':10.0}}

# Relaxation Data
daniProtocolStore = {'terms':{
                    'firstIt':0, 
                    'lastIt':2, 
                    'hot':1,
                    'cool1':5.0,
                    'cool2':10.0,
                    'cool3':10.0, 
                    'tc':9.771,
                    'anis':1.557,
                    'r':0.455,
                    'wh':599.91, 
                    'wn':60.82}}

# Docking annealing protocol
dockingProtocolStore = {'details':'Docking molecular simulation protocol',
                        'terms':{
                        'tadhigh_t':2000.0, 
                        'tadinit1_t':2000.0, 
                        'tadfinal1_t':500.0, 
                        'tadinit2_t':1000.0, 
                        'tadfinal2_t':50.0, 
                        'tadinit3_t':300.0, 
                        'tadfinal3_t':50.0, 
                        'timestep':0.002, 
                        'tadfactor':4.0, 
                        'initiosteps':500.0, 
                        'cool1_steps':500.0, 
                        'cool2_steps':1000.0, 
                        'cool3_steps':1000.0,
                        'waterheatsteps':100.0,
                        'watersteps':750.0,
                        'watercoolsteps':500.0
                        }}                
