
"""
======================COPYRIGHT/LICENSE START==========================

Constants.py: Useful constants for scripts in this directory and its subdirectories

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

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

#
# Translation between 1 letter code3Letter and '3 letter code'
# This '3 letter code' is used by Cns.
#

#nuclAcidCodes = {
bioPolymerCodes = {

     'RNA': [['A','C','G','U'],
             ['ADE','CYT','GUA','URA']],
     'DNA': [['A','DA','C','DC','G','DG','T','DT','U','DU'],
             ['ADE','ADE','CYT','CYT','GUA','GUA','THY','THY','URA','URA']],
     'protein': [['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'],
                 ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']]     
            }

nucleotideList = ['RNA','DNA']

defaultMolCode = ' '

defaultLowerDist = 1.8

defaultSeqInsertCode = ' '
defaultAltLoc = ' '

indent = "  "

atomOrder = ['A','B','G','D','E','Z','H']

chainCodeString = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'

stableBondTypes = ('covalent','disulfide','link')
bondTypes = stableBondTypes + ('hydrogen','salt')
