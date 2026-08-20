"""
======================COPYRIGHT/LICENSE START==========================

constants.py: Constants used for nmrStar file reading.writing.

Copyright (C) 2006 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)
- PDBe website (http://www.ebi.ac.uk/pdbe/)

- contact Wim Vranken (wim@ebi.ac.uk)
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

# Only contains fully valid mappings - incomplete ones are in presetDict.py
# for bmrb mapping

bmrbCodeToCcpCode = {

  'Aib':      (None,('protein','ABA'),None),
  'AIB':      (None,('protein','ABA'),None),
  'Glp':      (None,('protein','PCA'),None),
  '5HP':      (None,('protein','PCA'),None),
  'HYPR':     (None,('protein','PCA'),None),
  'Nle':      (None,('protein','NLE'),None),
  'Phol':     (None,('other','PHL'),None),
  'MYR':      (None,('other','MYR'),None),
  'fMET':     (None,('other','FME'),None),
  'FMET':     (None,('other','FME'),None),
  'C_heme_c': (None,('other','HEC'),None),
 
  'PFF':      (None,('other','PFF'),None),
  'NAL':      (None,('other','NAL'),None),
 
  # N-terminal acetylated
  'Ac_A':   (('other','ACE'),('protein','ALA'),None),
  'ALA_AC': (('other','ACE'),('protein','ALA'),None),
  
  'Ac_Aib': (('other','ACE'),('protein','ABA'),None),
  
  'Ac_C':   (('other','ACE'),('protein','CYS'),None),
  'CYS_AC': (('other','ACE'),('protein','CYS'),None),

  'Ac_D':   (('other','ACE'),('protein','ASP'),None),
  'ASP_AC': (('other','ACE'),('protein','ASP'),None),
  
  'Ac_E':   (('other','ACE'),('protein','GLU'),None),
  'GLU_AC': (('other','ACE'),('protein','GLU'),None),

  'PHE_AC': (('other','ACE'),('protein','PHE'),None),
  
  'Ac_G':   (('other','ACE'),('protein','GLY'),None),
  'GLY_AC': (('other','ACE'),('protein','GLY'),None),
  'GLY-ACE':(('other','ACE'),('protein','GLY'),None),
  
  'Ac_H':   (('other','ACE'),('protein','HIS'),None),
  '$Ac_H':  (('other','ACE'),('protein','HIS'),None),
  
  'ILE_AC': (('other','ACE'),('protein','ILE'),None),
  'ILE-ACE':(('other','ACE'),('protein','ILE'),None),
  
  'Ac_L':   (('other','ACE'),('protein','LEU'),None),
  'LEU_AC': (('other','ACE'),('protein','LEU'),None),
  'LEU-ACE':(('other','ACE'),('protein','LEU'),None),
  
  'Ac_M':   (('other','ACE'),('protein','MET'),None),
  
  'Ac_R':   (('other','ACE'),('protein','ARG'),None),
  'ARG_AC': (('other','ACE'),('protein','ARG'),None),

  'Ac_S':   (('other','ACE'),('protein','SER'),None),
  '$Ac_S':  (('other','ACE'),('protein','SER'),None),
  
  'Ac_T':   (('other','ACE'),('protein','THR'),None),
  'THR_AC': (('other','ACE'),('protein','THR'),None),
  
  'Ac_Y':   (('other','ACE'),('protein','TYR'),None),
  'TYR_AC': (('other','ACE'),('protein','TYR'),None),

  'TRP_AC': (('other','ACE'),('protein','TRP'),None),

  
  # C-terminal aminated
  'A_NH2':   (None,('protein','ALA'),('other','NH2')),
  '$A_NH2':  (None,('protein','ALA'),('other','NH2')),
  'ALA-NH2': (None,('protein','ALA'),('other','NH2')),

  'C_NH2':   (None,('protein','CYS'),('other','NH2')),
  'CYS-NH2': (None,('protein','CYS'),('other','NH2')),
  'CYS_NH2': (None,('protein','CYS'),('other','NH2')),

  'E_NH2':   (None,('protein','GLU'),('other','NH2')),
  'GLU-NH2': (None,('protein','GLU'),('other','NH2')),
  'L-GLU-NH2':(None,('protein','GLU'),('other','NH2')),

  'F_NH2':   (None,('protein','PHE'),('other','NH2')),
  'PHE_NH2': (None,('protein','PHE'),('other','NH2')),

  'G_NH2':   (None,('protein','GLY'),('other','NH2')),
  'GLY_NH2': (None,('protein','GLY'),('other','NH2')),
  'GLY-NH2': (None,('protein','GLY'),('other','NH2')),

  'H_NH2':   (None,('protein','HIS'),('other','NH2')),
  'HIS-NH2': (None,('protein','HIS'),('other','NH2')),

  'I_NH2':   (None,('protein','ILE'),('other','NH2')),
  'ILE-NH2': (None,('protein','ILE'),('other','NH2')),
  'ILE_NH2': (None,('protein','ILE'),('other','NH2')),

  'K_NH2':   (None,('protein','LYS'),('other','NH2')),
  'LYS-NH2': (None,('protein','LYS'),('other','NH2')),
  'LYS_NH2': (None,('protein','LYS'),('other','NH2')),

  'L_NH2':   (None,('protein','LEU'),('other','NH2')),
  'LEU-NH2': (None,('protein','LEU'),('other','NH2')),
  'LEU_NH2': (None,('protein','LEU'),('other','NH2')),

  'M_NH2':   (None,('protein','MET'),('other','NH2')),
  'MET-NH2': (None,('protein','MET'),('other','NH2')),

  'N_NH2':   (None,('protein','ASN'),('other','NH2')),
  'ASN-NH2': (None,('protein','ASN'),('other','NH2')),

  'P_NH2':   (None,('protein','PRO'),('other','NH2')),
  'PRO-NH2': (None,('protein','PRO'),('other','NH2')),

  'Q_NH2':   (None,('protein','GLN'),('other','NH2')),
  'GLN-NH2': (None,('protein','GLN'),('other','NH2')),

  'R_NH2':   (None,('protein','ARG'),('other','NH2')),
  'AARR':    (None,('protein','ARG'),('other','NH2')),
  'ARG-NH2': (None,('protein','ARG'),('other','NH2')),
  'ARG_NH2': (None,('protein','ARG'),('other','NH2')),

  'S_NH2':   (None,('protein','SER'),('other','NH2')),
  'SER-N':   (None,('protein','SER'),('other','NH2')),
  'SER-NH2': (None,('protein','SER'),('other','NH2')),
  'SER_NH2': (None,('protein','SER'),('other','NH2')),

  'T_NH2':   (None,('protein','THR'),('other','NH2')),
  '$T_NH2':  (None,('protein','THR'),('other','NH2')),
  'THR_NH2': (None,('protein','THR'),('other','NH2')),

  'Y_NH2':   (None,('protein','TYR'),('other','NH2')),
  'TYR-NH2': (None,('protein','TYR'),('other','NH2')),
  'TYR_NH2': (None,('protein','TYR'),('other','NH2')),
  
  'TRP_NH2': (None,('protein','TRP'),('other','NH2')),

}
