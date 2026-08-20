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
A Python example script for the use of the Haddock API.
"""

from HaddockApi    import *

if __name__ == '__main__':

  # Load project and create an instance of HaddockApi (Note that debugging is set to True)
  project = loadProject('/Users/marcvandijk/Documents/PHD/ExtendNMR/demo/george')
  api = HaddockApi(project,debug=True)

  # Upon creation of the HaddockApi instance a hierarchical class tree is generated of all
  # projects and runs saved in the CCPN project.

  # Loop over stored projects
  print("\n** Summary of stored Haddock projects and their content **")
  for pr in api.projects:
      print "project: ", pr.name
 
      # Loop over stored Haddock partners
      for partner in pr.partners:
          print "  partners: ", partner.code
 
      # Loop over stored runs
      for run in pr.runs:
          print "  run: ", run.serial

  # Projects, partners and runs can be retrieved specificly by name
  print("\n** Query the Haddock projects **")

  p = api.getHaddockProject(name='2GZK-demo')
  h = p.getHaddockPartner(code='A')
  r = p.getHaddockRun(serial=1)

  print("--> Fetched project: %s, partner: %s, run: %i" % (p.name,h.code,r.serial))

  # The generated classes still give access to all CCPN methods.
  print("\n--> Retrieve run specific data by model name")
  print("    numIt0Structures: %s" % (repr(r.numIt0Structures)))
  print("    numIt1Structures: %s" % (repr(r.numIt1Structures)))
  print("    numWrefStructures: %s" % (repr(r.numWrefStructures)))
  print("    useDnaRestraints: %s" % (repr(r.useDnaRestraints)))

  # In addition several convenience functions are included
  print("\n--> Retrieve all partner %s residues:" % h.code)
  residues = h.getPartnerResidues()
  for residue in residues:
      print residue.residue.seqId, residue.residue.ccpCode

  print("\n--> Set partner semi-flexibility mode")
  h.semiFlexMode = 'manual'

  print("\n--> Set partner residue semi-flexibility state")
  for res in residues[0:5]:
      h.setResidueFlexibility(residues=res,state='semi')

  print("\n--> Set partner residue AIR state")
  for res in residues[0:5]:
      h.setResidueAirState(residues=res,state='active')

  print("\n--> Set forcefield to be used")
  h.setForceField(molType='DNA')

  print("\n** creating a new project **")
  # Creating a new Haddock project automaticly initiates a default run
  # populated with the most commonly used parameters
  newp = api.newHaddockProject(name='NewProject')

  # Lets get some structure to populate the Haddock partners
  p1 = newp.newHaddockPartner(pdb='/Users/marcvandijk/Documents/PHD/ExtendNMR/1AZP_bound-DNA.pdb')
  p2 = newp.newHaddockPartner(molSystem=project.findFirstMolSystem())

  # Configure the run
  residues = p1.getPartnerResidues()
  for res in residues[1:-2]:
      p1.setResidueAirState(residues=res,state='active')

  p1.setSemiFlexMode(mode='manual')

  residues = p2.getPartnerResidues()
  for res in[residues[n] for n in [12,14,16,17,20,23,32]]:
      p2.setResidueAirState(residues=res,state='active')
  for res in[residues[n] for n in [11,15,18,19,21,22,24]]:
      p2.setResidueAirState(residues=res,state='passive')

  newr = newp.getHaddockRun(serial=1)
  newr.run.numIt0Structures = 100
  newr.run.numIt1Structures = 20
  newr.run.numWrefStructures = 20
 
  print("\n** Export haddock projects **")

  # Exporting is performed at the run level
  # Export all files needed to run Haddock on the local computer infrastructure (Classic mode)
  newr.exportClassicProject()

  # Export a selfcontaining parameter file you can upload to the server. In future releases of
  # Haddock this parameter file might also be used as input for a local run
  newr.exportHaddockParameterFile()



