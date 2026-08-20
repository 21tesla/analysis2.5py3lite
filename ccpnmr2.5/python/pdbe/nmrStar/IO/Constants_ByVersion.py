
#
# Constants class, can be specified by CCPN and/or NMR-STAR version
#
# Note that these dictionary are reversed in Util.py for NMR-STAR to CCPN!
#

constants_Ccpn_To_NmrStar = {

  #
  # Warning: this is copied from NmrStarFormat.
  #
  # The main mapping is kept here, but this could be version-specific it can
  # also be versioned specifically (see further down for example)
  #

  'molTypes': {

        'protein': 'polypeptide(L)',
        'DNA':     'polydeoxyribonucleotide',
        'RNA':     'polyribonucleotide',
        'DNA/RNA': 'DNA/RNA hybrid',
        'carbohydrate': 'carbohydrates'

                        }, 

  'boolean': {

          False: 'no',
          True:  'yes'

                        }
} 
    

constants_Ccpn__1_1_a2__To_NmrStar__3_0__ = constants_Ccpn_To_NmrStar
#constants_Ccpn__1_1_a2__To_NmrStar__3_0__['redefine'] = {'redefine': 'redefine'}
constants_Ccpn__1_1_2__To_NmrStar__3_1__  = constants_Ccpn_To_NmrStar
constants_Ccpn__1_1_a2__To_NmrStar__3_1__ = constants_Ccpn_To_NmrStar
constants_Ccpn__1_1_a3__To_NmrStar__3_0__ = constants_Ccpn_To_NmrStar
constants_Ccpn__1_1_a3__To_NmrStar__3_1__ = constants_Ccpn_To_NmrStar
