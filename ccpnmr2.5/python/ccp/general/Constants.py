"""
======================COPYRIGHT/LICENSE START==========================

Constants.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================
"""

# Tolerance for determining that two shifts are identical
shiftIdentityTolerance = 0.001
 
# conversion for nucleus NMR frequencies
# from :
# Recommendations for the presentation of NMR structures of proteins and nucleic acids
# John L. Markley, Ad Bax, Yoji Arata, C.W.Hilbers, Robert Kaptein, Brian D. Sykes, 
# Peter E. Wrightg and Kurt Wuthrich
# Journal of Biomolecular NMR, 12: 1-23, 1998
#
# The value for 19F is from the BioMagResBank web site (TFA in a sphere)
# citing Maurer, T. and Kalbitzer, H. R., 
# "Indirect Referencing of 31P and 19F NMR Spectra," 
# J. Magn. Reson. B 113, 177-178 (1996)

chemShiftRefRatios = {'1H':  1.0,
                      '2H':  0.153506088,
                      '13C': 0.251449530,
                      '15N': 0.101329118,
                      '31P': 0.404808636,
                      '19F': 0.940866982}
shiftRatio_13C_TMS = 0.25145020 # From 

#
# Reference repository name
#

refRepositoryName = 'refData'

#
# HTTP location for chemComp download
#

chemCompServer = "www.ebi.ac.uk"
chemCompWebPath = "/pdbe/docs/NMR/chemCompXml_2.0/files/"

#
# String to indicate void AtomSysName for ChemComps. This is to avoid problems with inheritance in some cases.
#

voidAtomSysName = " "

#
# Standard residue codes
#

standardResidueCcpCodes = {

    'protein':
   
      ['Ala','Arg','Asn','Asp','Cys',
       'Gln','Glu','Gly','His','Ile',
       'Leu','Lys','Met','Phe','Pro',
       'Ser','Thr','Trp','Tyr','Val'],
        
    'DNA':
    
      ['A','C','G','T'],
      
    'RNA':
    
      ['A','C','G','U']
      
                          }

code1LetterToCcpCodeDict = {

  'protein': {
     
     'A': 'Ala',
     'C': 'Cys',
     'D': 'Asp',
     'E': 'Glu',
     'F': 'Phe',
     'G': 'Gly',
     'H': 'His',
     'I': 'Ile',
     'K': 'Lys',
     'L': 'Leu',
     'M': 'Met',
     'N': 'Asn',
     'P': 'Pro',
     'Q': 'Gln',
     'R': 'Arg',
     'S': 'Ser',
     'T': 'Thr',
     'V': 'Val',
     'W': 'Trp',
     'Y': 'Tyr'
     
             },
             
  'DNA': {
     
     'A': 'A',
     'C': 'C',
     'G': 'G',
     'T': 'T',
     'I': 'I',
     'U': 'U',
     
          },
             
  'RNA': {
     
     'A': 'A',
     'C': 'C',
     'G': 'G',
     'T': 'T',
     'I': 'I',
     'U': 'U',
     
          }

    }

ccpCodeToCode1LetterDict = {

  'protein': {}, # set up below
  'DNA': code1LetterToCcpCodeDict['DNA'].copy(),
  'RNA': code1LetterToCcpCodeDict['RNA'].copy()
}

dict1 = code1LetterToCcpCodeDict['protein']
dict2 = ccpCodeToCode1LetterDict['protein']
for key in dict1.keys():
  dict2[dict1[key]] = key

#
# Information for pseudoatom corrections (from Wuthrich 1986)
#

amxPseudo =  {'HB*': (1.0, {'H': 0.6}),
              'HG*': (1.0, {'H': 1.0, 'HA': 0.6}),
              'HD*': (1.0, {'HA': 1.0}),
              'HE*': (1.0, {})}
              
aroPseudo = {'HB*': (1.0, {'H': 0.6}),
             'HD*': (2.0, {'H': 2.0, 'HA': 1.0}),
             'CG':  (2.0, {'H': 2.0, 'HA': 1.0}),
             'HE*': (2.0, {'H': 1.0, 'HA': 0.6}),
             'CZ':  (2.0, {'H': 1.0, 'HA': 0.6}),
             'HD*|HE*': (2.4, {'H': 2.4, 'HA': 2.0})}

pseudoAtomCorrectionsWuthrich = {

'Gly': {'HA*': (1.0,{})},

'Ala': {'HB*': (1.0, {'H': 0.6})},

'Ile': {'HG2*':(1.0, {'H': 1.0, 'HA': 0.6}),
        'HG1*':(1.0, {'H': 1.0, 'HA': 0.6}),
        'HD1*':(1.0, {'HA': 1.0}),
        },

'Thr': {'HG2*':(1.0, {'H': 1.0, 'HA': 0.6})},

'Val': {'HG1*':(1.0, {'H': 1.0, 'HA': 0.6}),
        'HG2*':(1.0, {'H': 1.0, 'HA': 0.6}),
        'HG*': (2.4, {'H': 1.7, 'HA': 0.6})
        },

'Leu': {'HB*': (1.0, {'H': 0.6}),
        'HD1*':(1.0, {'HA': 1.0}),
        'HD2*':(1.0, {'HA': 1.0}),
        'HD*': (2.4, {'HA': 1.7})
        },

'Lys':  amxPseudo,
'Ser':  amxPseudo,
'Asp':  amxPseudo,
'Asn':  amxPseudo,
'Cys':  amxPseudo,
'His':  amxPseudo,
'Trp':  amxPseudo,
'Glu':  amxPseudo,
'Gln':  amxPseudo,
'Met':  amxPseudo,
'Arg':  amxPseudo,
'Pro':  amxPseudo,
'Phe':  aroPseudo,
'Tyr':  aroPseudo
}


#
# Information on PDB secondary structure codes, and associated short codes
# Sometimes stored as application data on Residue level
#

secStrucInfo_kw = 'secStrucInfo'

# THIS SHOULD BECOME OBSOLETE!
ssPdbNamesDict = {             
      'Polyproline'    :    'ppr',
      'Right-handed 310':   '310',
      'Right-handed alpha': 'alf',
      'Right-handed pi':    'pih',
      'sheet'          :    'bsh',
      'turn'           :    'trn',
      None             :    'xxx'
    }


# TODO THIS SHOULD GO ELSEWHERE!!
"""
Not handled on DSSP side:
  * B = residue in isolated beta-bridge
  * S = bend 
Not handled on PDB side:
    'Polyproline'    :    '',  
"""
pdbToDsspNamesDict = {

      'Right-handed 310':   'G',
      'Right-handed alpha': 'H',
      'Right-handed pi':    'I',
      'sheet'          :    'E',
      'turn'           :    'T',
      None             :    'C',  # Is really space in DSSP.

      'Polyproline'    :    'C'  # No type available!
}


#
# Information to set up linkEnds (Wim 15/09/2004)
#

standardBackboneAtoms = {

 ('protein',):  ["H","N","CA","C","O"],
 ('DNA','RNA'): ["P","OP1","OP2","O5'","C5'","C4'","C3'","O3'","H3'"]
 
}

linkCodes = [('prev',2),('next',1)]

linkTorsions = {

  ('protein',):    [('OMEGA',['prev_2','prev_1','N','CA']),
                    ('PHI',  ['prev_1','N','CA','C']),
                    ('PSI',  ['N','CA','C','next_1'])],
             
  ('DNA','RNA'):  [('ZETA',  ['prev_2','prev_1','P',"O5'"]),
                   ('ALPHA', ['prev_1','P',"O5'","C5'"]),
                   ('EPSI',  ["C4'","C3'","O3'",'next_1'])]
  }


linkBonds = {

  ('protein',):    [('singleplanar','unknown',['prev_1','N']),
                    ('singleplanar','unknown',['C','next_1'])],
             
  ('DNA','RNA'):  [('single',None,['prev_1','P']),
                   ('single',None,["O3'",'next_1'])]
  }

linkEndDict = {

  ('protein',):    {'prev': ['N','CA'],
                    'next': ['C','CA']}, # These are 'backwards'
             
  ('DNA','RNA'):  {'prev': ['P',"O5'"],
                   'next': ["O3'","C3'"]} # These are 'backwards'
  
  }
          
linkingEnds = {

  'middle': ['prev','next'],
  'start': ['next'],
  'end': ['prev']

}
