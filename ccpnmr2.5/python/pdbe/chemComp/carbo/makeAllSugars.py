import sys, shutil, os

from pdbe.chemComp.carbo.makeFullSugar import makeFullSugar
from pdbe.chemComp.carbo.Constants import carboBaseNames

from pdbe.chemComp.modify.addSubstituent import addSubstituentToBaseUnit


if __name__ == '__main__':
   
  if '-create' in sys.argv:
    print "Warning: creating new sugar in edit/ directory!"
    testMode = False
  else:
    print "Creating in test directory!"
    testMode = True
 
  
  if '-replace' in sys.argv:
    print "Warning: replacing existing chemComps!"
    # TODO should really separate base and subst generation, seems to get mixed up...
    replace = True
  else:
    replace = False
 
  #
  # Remove local directory for creating subsituted chemComps
  #
  
  if os.path.exists('chemComp'):
    shutil.rmtree('chemComp')
   
  #
  # Make base sugars
  #
  
  coordSystem = 'euroCarbDb'
    
  for (carboBaseName,baseGlycoCtCode) in carboBaseNames: 
    makeFullSugar(carboBaseName,coordSystem,baseGlycoCtCode,testMode,saveData = True,replace=replace)

  #
  # Make modified sugars
  #
  
  modifiedCarbs = {
  
    'dgal-hex-1-5': ['C1_OMe',
                     'C2_NAc','C2_N',
                     'O3_Me','O4_Me','O6_Me',
                     'O3_SO3','O4_SO3','O6_SO3',
                     'C1_OMe:C2_NAc','C2_NAc:O3_Ac','O6_Ac',
                     'C2_NAc:O4_SO3','C2_NAc:O6_SO3','O1_P:C2_NAc','C2_NAc:O4_SO3:O6_SO3',
                     'O3_P','O4_P','O6_P'],
                     
    'dgal-hex-1-5-6d':    ['C2_NAc','C3_NAc','C4_NAc',
                           'C3_N','C2_NAc:C4_NAc'],
                     
    'dgal-hex-0-0-1aldi': [#'C1_OMe',
                           'C2_NAc',
                          ],
                     
    'lgal-hex-1-5-6d': ['C1_OMe','C1_OMe:C2_NAc','C2_NAc','O2_Me','O2_SO3','O3_Me','O4_Me','O2_SO3:O4_SO3','C2_NAc:O3_Ac','O3_Ac'],

    'dglc-hex-1-5': ['C1_OMe',
                     
                     'C1_OMe:C2_NAc','C2_NAc','C2_N','C2_NAc:O6_SO3','C1_OMe:C2_N','C2_NAc:O1_P','C2_NAc:O3_Me',
                     
                     'O2_P','O3_P','O4_P','O6_P',
                     
                     'O1_P:C2_NAc','C1_OMe:C2_NAc:O3_Me'],

    'dglc-hex-0-0-1aldi': [#'C1_OMe',
                           'C2_NAc',
                           'C2_NAc:O6_SO3',
                          ],

    'dgro-dgal-non-2-6-1a-2keto-3d': [#'C1_OMe',
                     #'C2_NAc',
                     #'O3_Me',
                     'C5_NAc',
                     'C5_NGlycol',
                     'O4_Ac:C5_NAc'#'C1_OMe:C2_NAc'
                     ],

    'dman-hex-1-5': ['C1_OMe',
                     'C2_NAc',
                     'O3_Me',
                     'C1_OMe:C2_NAc',
                     'O6_P'],
                     
    'lman-hex-1-5-6d': ['C1_OMe','C1_OMe:C2_NAc','C2_NAc','O2_Ac','C2_NAc:C3_NAc'],

    'dxyl-pen-1-5': ['C1_OMe', 'O2_P'],

    #'dgal-hex-1-5-6d': ['C1_OMe'],
    #'dman-hex-1-5-6d': ['C1_OMe'],
    
  }
  
  baseUnitCcpCodes = modifiedCarbs.keys()
  baseUnitCcpCodes.sort()
  
  for baseUnitCcpCode in baseUnitCcpCodes:
    
    for substituentInfo in modifiedCarbs[baseUnitCcpCode]:
    
      substituentList = substituentInfo.split(':')
    
      mergeInfoList = []
      
      for substituentElement in substituentList:
      
        (baseBindingAtomName,substituent) = substituentElement.split("_")
        baseBindingIndex = int(baseBindingAtomName[-1])
        
        if baseBindingAtomName[0] == 'C':
          removeBaseAtomNames = ("O%d" % baseBindingIndex, "H%dO" % baseBindingIndex)
        elif baseBindingAtomName[0] == 'O':
          # Have to remove the linking atom here...
          removeBaseAtomNames = ("O%d_1" % baseBindingIndex, "H%dO" % baseBindingIndex)
        

         # If base binding atom is C1, this is a methoxy modification at the reducing end, directly on the C1!!
        if substituent == 'OMe':
          mergeInfoList.append({

            'baseBindingAtomName': baseBindingAtomName,
            'removeBaseAtomNames': removeBaseAtomNames,
            'substituent': 'methoxy',
            # TODO: not sure about this one below...
            'renameSubstituentAtomNames': {'OM': 'O%d' % baseBindingIndex}

            })

        elif substituent == 'NAc':
          mergeInfoList.append({

            'baseBindingAtomName': baseBindingAtomName,
            'removeBaseAtomNames': removeBaseAtomNames,
            'substituent': 'n-acetyl',
            'renameSubstituentAtomNames': {}

            })
            
        elif substituent == 'NGlycol':
          mergeInfoList.append({

            'baseBindingAtomName': baseBindingAtomName,
            'removeBaseAtomNames': removeBaseAtomNames,
            'substituent': 'n-glycolyl',
            'renameSubstituentAtomNames': {}

            })

        elif substituent == 'N':
          mergeInfoList.append({

            'baseBindingAtomName': baseBindingAtomName,
            'removeBaseAtomNames': removeBaseAtomNames,
            'substituent': 'n',
            'renameSubstituentAtomNames': {}

            })
            
        elif substituent == 'Me':
          mergeInfoList.append({

            'baseBindingAtomName': baseBindingAtomName,
            'removeBaseAtomNames': removeBaseAtomNames,
            'substituent': 'methyl',
            'renameSubstituentAtomNames': {}

            })
            
        elif substituent == 'P':
          mergeInfoList.append({

            'baseBindingAtomName': baseBindingAtomName,
            'removeBaseAtomNames': removeBaseAtomNames,
            'substituent': 'phosphate',
            'renameSubstituentAtomNames': {}

            })

        elif substituent == 'Ac':
          mergeInfoList.append({

            'baseBindingAtomName': baseBindingAtomName,
            'removeBaseAtomNames': removeBaseAtomNames,
            'substituent': 'acetyl',
            'renameSubstituentAtomNames': {}

            })

        elif substituent == 'SO3':
          mergeInfoList.append({

            'baseBindingAtomName': baseBindingAtomName,
            'removeBaseAtomNames': removeBaseAtomNames,
            'substituent': 'sulfate',
            'renameSubstituentAtomNames': {}

            })

      addSubstituentToBaseUnit(baseUnitCcpCode,'carbohydrate',testMode,mergeInfoList,coordSystem, saveData = True, replace=replace,namingSystemName = 'EuroCarbDb', resetGlycoCtCode = True)
