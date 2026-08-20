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

import     sys
from     os.path         import     join
from     HaddockBasic    import    evalWcPairing

class dnaRnaRestraints:
    
    def __init__(self,ccpnproject,partner,projectRoot,arg=None):
        
        self.defaults         = {'debug':0,'verbose':True,'bpplan':False,'bplan':True,'pickpuc':True,
                                'pickbacdih':True,'c1pick':False,'c1lower':0.05,'c1upper':0.05,'wcpairing':True,
                               'wc_up':0.05,'wc_low':0.05,'wc_uri_up':0.01,'wc_uri_low':0.01}
        self.partner         = partner
        self.residueZones     = []
        self.fileString        = ""
        self.projectRoot     = projectRoot
        
        molSystem = self.partner.molSystem
        ensembles = ccpnproject.sortedStructureEnsembles()
        if molSystem: ensembles = [e for e in ensembles if e.molSystem is molSystem]
        for chain in ensembles[0].sortedCoordChains():
            residues = [r.residue.seqCode for r in chain.sortedResidues()]
            self.residueZones.append((residues[0],residues[-1]))

        self.__writeHeader()
        self.__writeBPplanarity()
        self.__writeBasePlanarity()
        self.__writePucker()
        self.__writeSPBackbone()
        self.__writeC1C1restraint()
        self.__writeWCpairing()
        self.__writeFooter()
        
    def writeToFile(self):
        
        outfile = open(join(self.projectRoot,'dna-rna_restraints.def'),'w')
        outfile.write(self.fileString)
        outfile.close()    

    def __writeBPplanarity(self):

        self.fileString += ("{=========================================== base-pair planarity ===========================================}\n")
        
        self.fileString += ("{* Use planarity restraints for Watson-Crick base pairing *}\n")
        self.fileString += ("{+ choice: true false +}\n\n")

        self.fileString += ("{===>} basepair_planar=%s;\n\n" % str(self.defaults['bpplan']).lower())

    def __writeBasePlanarity(self):    

        self.fileString +=  ("{============================================== base planarity =============================================}\n\n")
        
        self.fileString +=  ("{* Restrain base planarity. This selection must only include nucleotide residues *}\n\n")

        if self.defaults['bplan'] == True:
            zone = ""
            for rzone in self.residueZones:
                if len(zone) == 0: zone += ("(resid %i:%i and segid %s)" % (rzone[0],rzone[1],self.partner.code))
                else: zone += (" or (resid %i:%i and segid %s)" % (rzone[0],rzone[1],self.partner.code))    
            self.fileString +=  ("{===>} bases_planar=(%s);\n\n" % zone)    
        else:
            self.fileString +=  ("{* Base planarity not restraint *}\n\n")

    def __writePucker(self):

        self.fileString += ("{=================================== sugar-pucker dihedral angle restraints ================================}\n\n")
        
        self.fileString += ("{* Pick the dihedral angles of the sugar pucker from the input structure\n")
        self.fileString += ("   and restrain them within the given error range *}\n")
        self.fileString += ("{+ choice: true false +}\n\n")

        self.fileString += ("{===>} dna_pick_pucdih=%s;\n" % str(self.defaults['pickpuc']).lower())

        puckercount = 1
        for rzone in self.residueZones:
            pucker_group = ("resid %i:%i and segid %s" % (rzone[0],rzone[1],self.partner.code))

            self.fileString += ("{* residues with sugar pucker restrained - group %i *}\n" % puckercount)
            self.fileString += ("{===>} pucker_%i=(%s);\n\n" % (puckercount,pucker_group))

            self.fileString += ("{* conformation of group %i *}\n" % puckercount)
            self.fileString += ('{+ choice: "a-form" "b-form" "other" +}\n')
            
            self.fileString += ("{===>} pform_%i=\"other\";\n\n" % puckercount)
                
            self.fileString += ("{* user defined sugar pucker for group %i *}\n" % puckercount)

            self.fileString += ("{* dihedral C1'-C2'-C3'-C4' *}\n")
            self.fileString += ("{===>} dihedral_nu2_%i=-34.9;\n" % puckercount) 
            self.fileString += ("{* dihedral C1'-C2'-C3'-C4' error range *}\n")
            self.fileString += ("{===>} error_nu2_%i=0.0;\n" % puckercount)
            self.fileString += ("{* dihedral C5'-C4'-C3'-C2' *}\n")
            self.fileString += ("{===>} dihedral_nu3_%i=-86.4;\n" % puckercount) 
            self.fileString += ("{* dihedral C5'-C4'-C3'-C2' error range *}\n")
            self.fileString += ("{===>} error_nu3_%i=0.0;\n" % puckercount)
            self.fileString += ("{* dihedral C1'-O4'-C4'-C5' *}\n")
            self.fileString += ("{===>} dihedral_nu4_%i=106.4;\n" % puckercount) 
            self.fileString += ("{* dihedral C1'-O4'-C4'-C5' error range *}\n")
            self.fileString += ("{===>} error_nu4_%i=0.0;\n\n" % puckercount)

            puckercount += 1

    def __writeSPBackbone(self):

        self.fileString += ("{================================ phosphate backbone dihedral angle restraints =============================}\n\n")
        
        self.fileString += ("{* Pick the dihedral angles of the phosphate backbone from the input structure and\n")
        self.fileString += ("   restrain them within the given error range *}\n")
        self.fileString += ("{+ choice: true false +}\n\n")

        self.fileString += ("{===>} dna_pick_bacdih=%s;\n\n" % str(self.defaults['pickbacdih']).lower())

        bacdihcount = 1
        for rzone in self.residueZones:
            bacdih_group = ("resid %i:%i and segid %s" % (rzone[0],rzone[1],self.partner.code))

            self.fileString += ("{* residues with phosphate backbone restrained - group %i *}\n" % bacdihcount)
            self.fileString += ("{===>} dihedral_%i=(%s);\n\n" % (bacdihcount,bacdih_group))

            self.fileString += ("{* conformation of group %i *}\n" % bacdihcount)
            self.fileString += ('{+ choice: "a-form" "b-form" "other" +}\n')
            
            self.fileString += ("{===>} dform_%i=\"other\";\n\n" % bacdihcount)
        
            self.fileString += ("{* user defined posphate backbone for group %i *}\n" % bacdihcount)

            self.fileString += ("{* alpha dihedral O3'-P-O5'-C5' *}\n")
            self.fileString += ("{===>} dihedral_alpha_%i=-10.0;\n" % bacdihcount) 
            self.fileString += ("{* alpha dihedral range *}\n")
            self.fileString += ("{===>} error_alpha_%i=10.0;\n" % bacdihcount) 
            self.fileString += ("{* beta dihedral P-O5'-C5'-C4' *}\n")
            self.fileString += ("{===>} dihedral_beta_%i=136.4;\n" % bacdihcount) 
            self.fileString += ("{* beta dihedral range *}\n")
            self.fileString += ("{===>} error_beta_%i=40.0;\n" % bacdihcount) 
            self.fileString += ("{* gamma dihedral O5'-C5'-C4'-C3' *}\n")
            self.fileString += ("{===>} dihedral_gamma_%i=31.1;\n" % bacdihcount) 
            self.fileString += ("{* gamma dihedral range *}\n")
            self.fileString += ("{===>} error_gamma_%i=20.0;\n" % bacdihcount) 
            self.fileString += ("{* delta dihedral C5'-C4'-C3'-O3' *}\n")
            self.fileString += ("{===>} dihedral_delta_%i=-165.0;\n" % bacdihcount) 
            self.fileString += ("{* delta dihedral range *}\n")
            self.fileString += ("{===>} error_delta_%i=50.0;\n" % bacdihcount) 
            self.fileString += ("{* epsilon dihedral C4'-C3'-O3'-P *}\n")
            self.fileString += ("{===>} dihedral_eps_%i=-165.0;\n" % bacdihcount) 
            self.fileString += ("{* epsilon dihedral range *}\n")
            self.fileString += ("{===>} error_eps_%i=10.0;\n" % bacdihcount) 
            self.fileString += ("{* zeta dihedral C3'-O3'-P-O5' *}\n")
            self.fileString += ("{===>} dihedral_zeta_%i=-150.8;\n" % bacdihcount) 
            self.fileString += ("{* zeta dihedral range *}\n")
            self.fileString += ("{===>} error_zeta_%i=50.0;\n\n" % bacdihcount)

            bacdihcount += 1

    def __writeC1C1restraint(self):

        self.fileString += ("{============================================= C1'-C1' restraints ==========================================}\n\n")
        
        self.fileString += ("{* Have the length of the C1'-C1' virtual bonds measured and restraints. *}\n")
        self.fileString += ("{+ choice: true false +}\n")

        self.fileString += ("{===>} dna_pick_c1=%s;\n\n" % str(self.defaults['c1pick']).lower())

        self.fileString += ("{* Error range used for C1'-C1' virtual bonds  *}\n")
        self.fileString += ("{===>} c1_low=%1.3f;\n" % self.defaults['c1lower'])
        self.fileString += ("{===>} c1_up=%1.3f;\n\n" % self.defaults['c1upper'])    

    def __writeWCpairing(self):

        self.fileString += ("{=========================================== Watson-Crick base pairs =======================================}\n\n")

        self.fileString += ("{* pick Watson-Crick restraint values from structure *}\n")
        self.fileString += ("{+ choice: true false +}\n")

        self.fileString += ("{===>} dna_pick_wc=%s;\n" % str(self.defaults['wcpairing']).lower())
        self.fileString += ("{* error range used for dna_pick_wc defined Watson-Crick restraints *}\n")
        self.fileString += ("{===>} wc_low=%1.3f;\n" % self.defaults['wc_low'])
        self.fileString += ("{===>} wc_up=%1.3f;\n" % self.defaults['wc_up'])
        self.fileString += ("{* for URI, for default much lower range... why?*}\n")
        self.fileString += ("{===>} wc_low_uri=%1.3f;\n" % self.defaults['wc_uri_low'])
        self.fileString += ("{===>} wc_up_uri=%1.3f;\n\n" % self.defaults['wc_uri_up'])

        self.fileString += ("{* residues which form Watson-Crick pairs *}\n\n")
        
        paircount = 1
        pairs = evalWcPairing(self.partner)
        for pair in pairs.pairs:
            self.fileString += ("{* selection for pair %i base A *}\n" % paircount)
            self.fileString += ("{===>} base_a_%i=(resid %i and segid %s);\n" % (paircount,pair[0].seqCode,self.partner.code))
            self.fileString += ("{* selection for pair %i base B *}\n" % paircount)
            self.fileString += ("{===>} base_b_%i=(resid %i and segid %s);\n\n" % (paircount,pair[1].seqCode,self.partner.code))

            paircount += 1

    def __writeHeader(self):

        self.fileString += ("""{+ file: dna-rna_restraints.def       directory: protocols +}
{+ description: Creates restraints to maintain conformation of DNA/RNA +}
{+ comment:This file is to be read by refinement files that modify atom coordinates +}
{+ authors: Axel T. Brunger, and Paul D. Adams, <br>
            modified by Alexandre Bonvin and Marc van Dijk for HADDOCK use
            <br><br>
Additions and changes were made to allow for flexibility during docking<br><br>
Changes include: <br>
<ul>
<li> flags to turn all options on or off
<li> separation of sugar
<li> pucker restraints and phosphate backbone restraints
<li> option to have sugar-phosphate backbone dihedrals measured and restrained within a user defined error range
<li> option to have the length of the Watson-Crick hydrogen bonds measured from the structure measured and restrained within a user defined error range.
</ul>
+}

set message=normal echo=on end

{- begin block parameter definition -} define(\n\n""")

    def __writeFooter(self):

        self.fileString += ("""{=========================================================================================================}
{                        things below this line do not normally need to be changed                        }
{=========================================================================================================}

 ) {- end block parameter definition -}

{- the planarity restraints for Watson-Crick base pairing -}

if (&basepair_planar=true) then
 evaluate ($pair=1)
 evaluate ($done=false)
 while ( $done = false ) loop plan_paired
   if ( &exist_base_a_$pair = true ) then
     if ( &exist_base_b_$pair = true ) then
       show (segid) ( &base_a_$pair and name C1' ) 
       evaluate ($Asegid=$result)
       show (resid) ( &base_a_$pair and name C1' ) 
       evaluate ($Aresid=$result)
       show (segid) ( &base_b_$pair and name C1' ) 
       evaluate ($Bsegid=$result)
       show (resid) ( &base_b_$pair and name C1' ) 
       evaluate ($Bresid=$result)
       evaluate ($plweight = 20)            ! Enforce planarity by increasing plweight value.

       restraints plane

         group
           selection=(((segid $Asegid and resid $Aresid) or (segid $Bsegid and resid $Bresid)) and
                      (resname THY or resname CYT or resname GUA or
                       resname ADE or resname URI) and
                       not (name c#' or name h#' or name h#'' or name o#p or
                            name h7# or name o#' or name p or name h#t or name o#t))
           weight=$plweight
         end
       end
     end if
   else
     evalute ($done = true)
   end if
     evaluate ($pair = $pair + 1)
 end loop plan_paired
else
end if
flag include plan end

{- the planarity restraints single bases -}

 for $id in id ( &bases_planar and tag ) loop plan
   show (segid) (id $id)
   evaluate ($segid=$result)
   show (resid) (id $id)
   evaluate ($resid=decode($result))
   evaluate ($plweight = 20)

   restraints plane

     group
       selection=( segid $segid and resid $resid and
                  (resname THY or resname CYT or resname GUA or
                   resname ADE or resname URI) and
                   not (name c#' or name h#' or name h#'' or name o#p or
                        name h7# or name o#' or name p or name h#t or name o#t))
       weight=$plweight
     end
   end
 end loop plan

{- Dihedral restraints for the sugar pucker -}

if (&dna_pick_pucdih=true) then
  evaluate ($group=1)
  evaluate ($done=false)
  while ( $done = false ) loop dihe
   if ( &exist_pucker_$group = true ) then
     show sum(1) ( &pucker_$group )
     if ( $result > 0 ) then
       evaluate ($min_resid_$group = 99999)
       evaluate ($max_resid_$group = -99999)
       evaluate ($error_nu2=&error_nu2_$group)
       evaluate ($error_nu3=&error_nu3_$group)
       evaluate ($error_nu4=&error_nu4_$group)     
       for $id in id ( &pucker_$group and tag ) loop resid
         show (segid) (id $id)
         evaluate ($segid=$result)
         show (resid) ( id $id )
         evaluate ($resid=decode($result))
     evaluate ($min_resid_$group = max($min_resid_$group,$resid))
     evaluate ($max_resid_$group = max($max_resid_$group,$resid))
         pick dihedral
                   ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name c2' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c4' ) 
          geometry
     evaluatate ($dihedral_nu2=$result)
     pick dihedral
                   ( segid $segid and resid $resid and name c5' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c2' ) 
          geometry
     evaluatate ($dihedral_nu3=$result)
     pick dihedral
                   ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name o4' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c5' ) 
          geometry
     evaluatate ($dihedral_nu4=$result)

      restraints dihedral
           assign  ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name c2' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c4' ) 
                                                       20.0 $dihedral_nu2 $error_nu2 2
           assign  ( segid $segid and resid $resid and name c5' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c2' ) 
                                                       20.0 $dihedral_nu3 $error_nu3 2
           assign  ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name o4' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c5' ) 
                                                       20.0 $dihedral_nu4 $error_nu4 2
           scale=20.0
         end
       end loop resid
     end if
   else
     evaluate ($done=true)
   end if
   evaluate ($group=$group+1)
  end loop dihe

 else

 evaluate ($group=1)
 evaluate ($done=false)
 while ( $done = false ) loop dihe
   if ( &exist_pucker_$group = true ) then
     show sum(1) ( &pucker_$group )
     if ( $result > 0 ) then
       if ( &pform_$group = "a-form" ) then
         evaluate ($dihedral_nu2=37.053)
         evaluate ($dihedral_nu3=-155.59)
         evaluate ($dihedral_nu4=144.26)
       elseif ( &pform_$group = "b-form" ) then
         evaluate ($dihedral_nu2=-34.9)
         evaluate ($dihedral_nu3=-86.4)
         evaluate ($dihedral_nu4=106.4)
       elseif ( &pform_$group = "other" ) then
         evaluate ($dihedral_nu2=&dihedral_nu2_$group)
         evaluate ($dihedral_nu3=&dihedral_nu3_$group)
         evaluate ($dihedral_nu4=&dihedral_nu4_$group)
       end if

       evaluate ($min_resid_$group = 99999)
       evaluate ($max_resid_$group = -99999)

       for $id in id ( &pucker_$group and tag ) loop resid

         show (segid) (id $id)
         evaluate ($segid=$result)
         show (resid) ( id $id )
         evaluate ($resid=decode($result))
     evaluate ($min_resid_$group = max($min_resid_$group,$resid))
     evaluate ($max_resid_$group = max($max_resid_$group,$resid))

         restraints dihedral
           assign  ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name c2' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c4' ) 
                                                       20.0 $dihedral_nu2 0.0 2
           assign  ( segid $segid and resid $resid and name c5' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c3' )
                   ( segid $segid and resid $resid and name c2' ) 
                                                       20.0 $dihedral_nu3 0.0 2
           assign  ( segid $segid and resid $resid and name c1' )
                   ( segid $segid and resid $resid and name o4' )
                   ( segid $segid and resid $resid and name c4' )
                   ( segid $segid and resid $resid and name c5' ) 
                                                       20.0 $dihedral_nu4 0.0 2

           scale=20.0
         end
       end loop resid
     end if
   else
     evaluate ($done=true)
   end if
   evaluate ($group=$group+1)
 end loop dihe
end if
flags include cdih end

{- Dihedral restraints for the phosphate backbone -}

if (&dna_pick_bacdih=true) then
  evaluate ($group=1)
  evaluate ($done=false)
  while ( $done = false ) loop bdihe
   if ( &exist_dihedral_$group = true ) then
     show sum(1) ( &dihedral_$group )
     if ( $result > 0 ) then
       evaluate ($resid=$min_resid_$group)
       evaluate ($nres=$max_resid_$group - $min_resid_$group + 1)
       evaluate ($error_alpha=&error_alpha_$group)
       evaluate ($error_beta=&error_beta_$group)
       evaluate ($error_gamma=&error_gamma_$group)
       evaluate ($error_zeta=&error_zeta_$group)
       evaluate ($error_epsilon=&error_eps_$group)
       evaluate ($error_delta=&error_delta_$group)
       for $id in id ( &dihedral_$group and tag ) loop resid
         show (segid) (id $id)
         evaluate ($segid=$result)
         show (resid) ( id $id )
         evaluate ($resid=decode($result))
         if ($resid > $min_resid_$group) then
           evaluate ($rprec = $resid - 1)
       pick dihedral
                     ( segid $segid and resid $rprec and name O3' )
                     ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' ) 
          geometry
       evaluatate ($dihedral_alpha=$result)
       pick dihedral
                     ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' ) 
          geometry
       evaluatate ($dihedral_beta=$result)

           restraint dihedral
        ! alpha
             assign  ( segid $segid and resid $rprec and name O3' )
                     ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' ) 
                                                       1.0 $dihedral_alpha $error_alpha 2
        ! beta                           
             assign  ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' ) 
                                                       1.0 $dihedral_beta $error_beta 2
             scale 200.0
           end
         end if

     pick dihedral
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' ) 
            geometry
     evaluatate ($dihedral_gamma=$result)
         pick dihedral
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' ) 
          geometry
       evaluatate ($dihedral_delta=$result)

     restraints dihedral
        ! gamma
             assign  ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' ) 
                                                       1.0 $dihedral_gamma $error_gamma 2
        ! delta                           
             assign  ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' ) 
                                                       1.0 $dihedral_delta $error_delta 2            
          scale=200.0
            end

      if ($resid < $max_resid_$group) then
           evaluate ($rfoll = $resid + 1)
       pick dihedral
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P ) 
          geometry
       evaluatate ($dihedral_epsilon=$result)
       pick dihedral
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P )
                     ( segid $segid and resid $rfoll and name O5' ) 
          geometry
       evaluatate ($dihedral_zeta=$result)
           restraint dihedral
             ! epsilon
         assign  ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P ) 
                                                       1.0 $dihedral_epsilon $error_epsilon 2
             ! zeta
         assign  ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P )
                     ( segid $segid and resid $rfoll and name O5' ) 
                                                       1.0 $dihedral_zeta $error_zeta 2
             scale 200.0
           end
         end if
       end loop resid
     end if
   else
     evaluate ($done=true)
   end if
     evaluate ($group=$group+1)
 end loop bdihe

 else

 evaluate ($group=1)
 evaluate ($done=false)
 while ( $done = false ) loop bdihe
 if ( &exist_dihedral_$group = true ) then
     show sum(1) ( &dihedral_$group )
     if ( $result > 0 ) then
       evaluate ($resid=$min_resid_$group)
       evaluate ($nres=$max_resid_$group - $min_resid_$group + 1)
       if ( &dform_$group = "a-form" ) then
         evaluate ($dihedral_alpha=-70)
     evaluate ($error_alpha=50)
         evaluate ($dihedral_beta=180)
     evaluate ($error_beta=50)
         evaluate ($dihedral_gamma=60)
     evaluate ($error_gamma=35)
         evaluate ($dihedral_delta=81)
         evaluate ($error_delta=20)
         evaluate ($dihedral_zeta=-85)
     evaluate ($error_zeta=50)
         evaluate ($dihedral_epsilon=180)
     evaluate ($error_epsilon=35)
       elseif ( &dform_$group = "b-form" ) then
         evaluate ($dihedral_alpha=-63.6)
     evaluate ($error_alpha=6)
         evaluate ($dihedral_beta=176)
     evaluate ($error_beta=7)
         evaluate ($dihedral_gamma=51.4)
     evaluate ($error_gamma=7)
         evaluate ($dihedral_delta=128)
         evaluate ($error_delta=13)
         evaluate ($dihedral_epsilon=-171.7)
     evaluate ($error_epsilon=3.7)
         evaluate ($dihedral_zeta=-103.8)
     evaluate ($error_zeta=10)
       elseif ( &dform_$group = "other" ) then
         evaluate ($dihedral_alpha=&dihedral_alpha_$group)
     evaluate ($error_alpha=&error_alpha_$group)
         evaluate ($dihedral_beta=&dihedral_beta_$group)
     evaluate ($error_beta=&error_beta_$group)
         evaluate ($dihedral_gamma=&dihedral_gamma_$group)
     evaluate ($error_gamma=&error_gamma_$group)
         evaluate ($dihedral_delta=&dihedral_delta_$group)
     evaluate ($error_delta=&error_delta_$group)
     evaluate ($dihedral_zeta=&dihedral_zeta_$group)
     evaluate ($error_zeta=&error_zeta_$group)
         evaluate ($dihedral_epsilon=&dihedral_eps_$group)
     evaluate ($error_epsilon=&error_eps_$group)
       end if

       for $id in id ( &dihedral_$group and tag ) loop resid
         show (segid) (id $id)
         evaluate ($segid=$result)
         show (resid) ( id $id )
         evaluate ($resid=decode($result))
         if ($resid > $min_resid_$group) then
           evaluate ($rprec = $resid - 1)
           restraint dihedral
             ! alpha
         assign  ( segid $segid and resid $rprec and name O3' )
                     ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' ) 
                                                       1.0 $dihedral_alpha $error_alpha 2
             ! beta
         assign  ( segid $segid and resid $resid and name P )
                     ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' ) 
                                                       1.0 $dihedral_beta $error_beta 2
             scale 200.0
           end
         end if

         restraints dihedral
             ! gamma
         assign  ( segid $segid and resid $resid and name O5' )
                     ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' ) 
                                                       1.0 $dihedral_gamma $error_gamma 2
             !delta
         assign  ( segid $segid and resid $resid and name C5' )
                     ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' ) 
                                                       1.0 $dihedral_delta $error_delta 2             
         scale=200.0
           end

     if ($resid < $max_resid_$group) then
           evaluate ($rfoll = $resid + 1)
           restraint dihedral
             ! epsilon
         assign  ( segid $segid and resid $resid and name C4' )
                     ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P ) 
                                                       1.0 $dihedral_epsilon $error_epsilon 2
             ! zeta
         assign  ( segid $segid and resid $resid and name C3' )
                     ( segid $segid and resid $resid and name O3' )
                     ( segid $segid and resid $rfoll and name P )
                     ( segid $segid and resid $rfoll and name O5' ) 
                                                       1.0 $dihedral_zeta $error_zeta 2
             scale 200.0
           end
         end if
       end loop resid
     end if
   else
     evaluate ($done=true)
   end if
   evaluate ($group=$group+1)
  end loop bdihe
 end if
flags include cdih end

{- C1'-C1' virtual bond length restraints -}

noe
   class hres
   averaging hres cent
   potential hres square
   sqconstant hres 1.
   sqexponent hres 2
   scale hres 70.
 end           

if (&dna_pick_c1 = true) then
  evaluate ($pair=1)
  evaluate ($done=false)
  while ( $done = false ) loop noe
   if ( &exist_base_a_$pair = true ) then
     if ( &exist_base_b_$pair = true ) then
       show ( resname ) ( &base_a_$pair and name C1' ) 
       evaluate ($ares=$result)
       show ( resname ) ( &base_b_$pair and name C1' ) 
       evaluate ($bres=$result)
        pick bond
            (&base_a_$pair and name C1') 
            (&base_b_$pair and name C1')
       geometry
        evaluate ($c1c1=$result)
        noe
        assign (&base_a_$pair and name C1') 
               (&base_b_$pair and name C1') $c1c1 &c1_low &c1_up 
        end
     end if
   else
     evaluate ($done=true)
   end if         
     evaluate ($pair=$pair+1)
  end loop noe
 else
end if
flags include noe end

{- Watson-Crick base pairing -}

 noe
   class hres
   averaging hres cent
   potential hres square
   sqconstant hres 1.
   sqexponent hres 2
   scale hres 70.
 end           

 if (&dna_pick_wc = true) then
  evaluate ($pair=1)
  evaluate ($done=false)
  while ( $done = false ) loop noe
   if ( &exist_base_a_$pair = true ) then
     if ( &exist_base_b_$pair = true ) then
       show ( resname ) ( &base_a_$pair and name C1' ) 
       evaluate ($ares=$result)
       show ( resname ) ( &base_b_$pair and name C1' ) 
       evaluate ($bres=$result)
       if ( $ares = THY ) then
        pick bond
                  (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6)
          geometry
        evaluate ($o4n6=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1)
          geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_a_$pair and name h3) 
                  (&base_b_$pair and name n1)
          geometry
        evaluate ($h3n1=$result)
        pick bond
                  (&base_a_$pair and name o2) 
                  (&base_b_$pair and name h2)
          geometry
        evaluate ($o2h2=$result)
        pick bond
                  (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1)
          geometry
        evaluate ($o4n1=$result)
        pick bond
                  (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n1)
          geometry
        evaluate ($o2n1=$result)
       elseif ( $ares = URI ) then
        pick bond
                  (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6)
          geometry
        evaluate ($o4n6=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1)
          geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1)
          geometry
        evaluate ($o4n1=$result)
        pick bond
                  (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n6)
          geometry
        evaluate ($o2n6=$result)
       elseif ( $ares = ADE ) then
        pick bond
                  (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n6)
          geometry
        evaluate ($o4n6=$result)
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1)
          geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_b_$pair and name h3) 
                  (&base_a_$pair and name n1)
          geometry
        evaluate ($h3n1=$result)
        pick bond
                  (&base_b_$pair and name o2) 
                  (&base_a_$pair and name h2)
          geometry
        evaluate ($o2h2=$result)
        pick bond
                  (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n1)
          geometry
        evaluate ($o4n1=$result)
        pick bond
                  (&base_b_$pair and name o2) 
                  (&base_a_$pair and name n1)
          geometry
        evaluate ($o2n1=$result)
       elseif ( $ares = CYT ) then
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1)
          geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name h1)
          geometry
        evaluate ($n3h1=$result)
        pick bond
                  (&base_a_$pair and name n4) 
                  (&base_b_$pair and name o6)
          geometry
        evaluate ($n4o6=$result)
        pick bond
                  (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n2)
          geometry
        evaluate ($o2h2=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name o6)
          geometry
        evaluate ($n3o6=$result)
        pick bond
                  (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n2)
          geometry
        evaluate ($n3n2=$result)
       elseif ( $ares = GUA ) then
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1)
          geometry
        evaluate ($n3n1=$result)
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name h1)
          geometry
        evaluate ($n3h1=$result)
        pick bond
                  (&base_b_$pair and name n4) 
                  (&base_a_$pair and name o6)
          geometry
        evaluate ($n4o6=$result)
        pick bond
                  (&base_b_$pair and name o2) 
                  (&base_a_$pair and name n2)
          geometry
        evaluate ($o2n2=$result)
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name o6)
          geometry
        evaluate ($n3o6=$result)
        pick bond
                  (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n2)
          geometry
        evaluate ($n3n2=$result)

       end if
       noe
         if ( $ares = THY ) then
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6) $o4n6 &wc_low &wc_up 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) $n3n1 &wc_low &wc_up
           assign (&base_a_$pair and name h3) 
                  (&base_b_$pair and name n1) $h3n1 &wc_low &wc_up 
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name h2) $o2h2 &wc_low &wc_up
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1) $o4n1 &wc_low &wc_up
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n1) $o2n1 &wc_low &wc_up
         elseif ( $ares = URI ) then
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6) $o4n6 &wc_low_uri &wc_up_uri 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) $n3n1 &wc_low_uri &wc_up_uri
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1) $o4n1 &wc_low_uri &wc_up_uri 
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n6) $o2n6 &wc_low_uri &wc_up_uri
         elseif ( $ares = ADE ) then
           assign (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n6) $o4n6 &wc_low &wc_up 
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1) $n3n1 &wc_low &wc_up
           assign (&base_b_$pair and name h3) 
                  (&base_a_$pair and name n1) $h3n1 &wc_low &wc_up 
           assign (&base_b_$pair and name o2) 
                  (&base_a_$pair and name h2) $o2h2 &wc_low &wc_up
           assign (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n1) $o4n1 &wc_low &wc_up
           assign (&base_b_$pair and name o2) 
                  (&base_a_$pair and name n1) $o2n1 &wc_low &wc_up
         elseif ( $ares = CYT ) then
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) $n3n1 &wc_low &wc_up 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name h1) $n3h1 &wc_low &wc_up
           assign (&base_a_$pair and name n4)
                  (&base_b_$pair and name o6) $n4o6 &wc_low &wc_up
           assign (&base_a_$pair and name o2)
                  (&base_b_$pair and name n2) $o2n2 &wc_low &wc_up 
           assign (&base_a_$pair and name n3)
                  (&base_b_$pair and name o6) $n3o6 &wc_low &wc_up
           assign (&base_a_$pair and name n3)
                  (&base_b_$pair and name n2) $n3n2 &wc_low &wc_up
         elseif ( $ares = GUA ) then
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1) $n3n1 &wc_low &wc_up 
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name h1) $n3h1 &wc_low &wc_up
           assign (&base_b_$pair and name n4)
                  (&base_a_$pair and name o6) $n4o6 &wc_low &wc_up
           assign (&base_b_$pair and name o2)
                  (&base_a_$pair and name n2) $o2n2 &wc_low &wc_up 
           assign (&base_b_$pair and name n3)
                  (&base_a_$pair and name o6) $n3o6 &wc_low &wc_up
           assign (&base_b_$pair and name n3)
                  (&base_a_$pair and name n2) $n3n2 &wc_low &wc_up
         end if
       end
     end if
   else
     evaluate ($done=true)
   end if         
   evaluate ($pair=$pair+1)
  end loop noe

 else

 evaluate ($pair=1)
 evaluate ($done=false)
 while ( $done = false ) loop noe
   if ( &exist_base_a_$pair = true ) then
     if ( &exist_base_b_$pair = true ) then
       show ( resname ) ( &base_a_$pair and name C1' ) 
       evaluate ($ares=$result)
       show ( resname ) ( &base_b_$pair and name C1' ) 
       evaluate ($bres=$result)
       noe
         if ( $ares = THY ) then
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6) 2.89 0.2 0.2 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) 2.92 0.2 0.2
           assign (&base_a_$pair and name h3) 
                  (&base_b_$pair and name n1) 1.87 0.2 0.2 
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name h2) 2.94 0.2 0.2
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1) 3.69 0.2 0.2
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n1) 3.67 0.2 0.2        
         elseif ( $ares = URI ) then
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n6) 2.95 0.01 0.01 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) 2.82 0.01 0.01
           assign (&base_a_$pair and name o4) 
                  (&base_b_$pair and name n1) 3.63 0.01 0.01 
           assign (&base_a_$pair and name o2) 
                  (&base_b_$pair and name n6) 5.40 0.01 0.01
         elseif ( $ares = ADE ) then
           assign (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n6) 2.89 0.2 0.2 
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1) 2.92 0.2 0.2
           assign (&base_b_$pair and name h3) 
                  (&base_a_$pair and name n1) 1.87 0.2 0.2 
           assign (&base_b_$pair and name o2) 
                  (&base_a_$pair and name h2) 2.94 0.2 0.2
           assign (&base_b_$pair and name o4) 
                  (&base_a_$pair and name n1) 3.69 0.2 0.2
           assign (&base_b_$pair and name o2) 
                  (&base_a_$pair and name n1) 3.67 0.2 0.2
         elseif ( $ares = CYT ) then
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name n1) 2.87 0.2 0.2 
           assign (&base_a_$pair and name n3) 
                  (&base_b_$pair and name h1) 1.86 0.2 0.2
           assign (&base_a_$pair and name n4)
                  (&base_b_$pair and name o6) 2.81 0.2 0.2
           assign (&base_a_$pair and name o2)
                  (&base_b_$pair and name n2) 2.81 0.2 0.2 
           assign (&base_a_$pair and name n3)
                  (&base_b_$pair and name o6) 3.58 0.2 0.2
           assign (&base_a_$pair and name n3)
                  (&base_b_$pair and name n2) 3.63 0.2 0.2
         elseif ( $ares = GUA ) then
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name n1) 2.87 0.2 0.2 
           assign (&base_b_$pair and name n3) 
                  (&base_a_$pair and name h1) 1.86 0.2 0.2
           assign (&base_b_$pair and name n4)
                  (&base_a_$pair and name o6) 2.81 0.2 0.2
           assign (&base_b_$pair and name o2)
                  (&base_a_$pair and name n2) 2.81 0.2 0.2 
           assign (&base_b_$pair and name n3)
                  (&base_a_$pair and name o6) 3.58 0.2 0.2
           assign (&base_b_$pair and name n3)
                  (&base_a_$pair and name n2) 3.63 0.2 0.2
         end if
       end
     end if
   else
     evaluate ($done=true)
   end if         
   evaluate ($pair=$pair+1)
 end loop noe
 end if
 flags include noe end

    set message=off echo=off end""")