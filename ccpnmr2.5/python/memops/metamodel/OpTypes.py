""" Data and functions for checking and manipulating operations and opTypes
    NB must be usable in Python 2.1
"""

from memops.universal import Util as uniUtil
from memops.metamodel import ImpConstants
from memops.general import Constants as genConstants

True = not 0
False = not True


# input for operations
#
# Dictionary key is opType. Special cases where the operation is
# different from the standard for the opType are handled in ModelAdapt
#
# subOps dictionary key is opSubType
#
# opSubType==None is the default form of the operation (with no suffix)
#
# The 'group' classifies the operation with respect to access control
#
# suffix,  opDocTemplate,  parDocumentation
# all default to empty strings.
#
# Operation targets (main dictionary) serve to derive parameters etc. 
# and is also involved in access control.
# 
# Parameter targets (subOps parmeters) must be a (qualified name of a) 
# Model Element
# If it is a Class or DataType this will be the type of the final parameter
# If it is an Attribute or Role the final parameter will be calculated
# from the target.
#
# useCollection determines whether something pertains to the collection of a
# ClassElement, or to a single member of it.
# For operations (main dictionary) it determins whether target name or basename
# is used for the operation name. For parameters it setermines if the type is a
# collection or a member element.
#
#
# MetaParameters default to hicard:1,locard:1,isOrdered:False,isUnique:False
#
# The basic modelGen handles 'result' (return value), 'value' (input value)
# and 'parent' (parent class of given class)
# whose specification is given by operation and context.
#
# parameter names are substituted according to self.varNames
# Also parDocumentation is % formatted with self.varNames and documentations are

# NBNB TBD parameter target strings - where if ever converted to real types?

operationData = {
 'get':{
 'group':'query',
 'targetTag':'ClassElement',
 'useCollection':True,
# 'name':'useInsteadOfOpType',
# 'scope':ImpConstants.classifier_level # default is 'instance_level'
# 'taggedValues':{'key:'value', 'key:'value'}
 },
 'findFirst':{
 'group':'query',
 'targetTag':'ClassElement',
 'useCollection':False,
 },
 'findAll':{
 'group':'query',
 'targetTag':'ClassElement',
 'useCollection':True,
 },
 'sorted':{
 'group':'query',
 'targetTag':'MetaRole',
 'useCollection':True,
 },
 'set':{
 'group':'modify',
 'targetTag':'ClassElement',
 'useCollection':True,
 },
 'add':{
 'group':'modify',
 'targetTag':'ClassElement',
 'useCollection':False,
 },
 'remove':{
 'group':'modify',
 'targetTag':'ClassElement',
 'useCollection':False,
 },
 'new':{
 'group':'create',
 'targetTag':'ChildClass',
 },
 'init':{
 'group':'create',
 'targetTag':'container',
 },
 'clone':{
 'group':'create',
 'targetTag':'container',
 },
 'fullDelete':{
 'group':'delete',
 'targetTag':'container',
 'name':'delete',
 },
 'checkDelete':{
 'group':'delete',
 'targetTag':'container',
 'name':'_checkDelete',
 },
 'singleDelete':{
 'group':'delete',
 'targetTag':'container',
 'name':'_singleDelete',
 },
 'checkValid':{
 'group':'query',
 'targetTag':'container',
 },
 'checkAllValid':{
 'group':'query',
 'targetTag':'container',
 },
 'getAttr':{
 'group':'query',
 'targetTag':'container',
 'useCollection':True,
 'name':'get',
 },
 'setAttr':{
 'group':'modify',
 'targetTag':'container',
 'useCollection':True,
 'name':'set',
 },
 'getByKey':{
 'group':'query',
 'targetTag':'container',
 'scope':ImpConstants.classifier_level,
 },
 'getFullKey':{
 'group':'query',
 'targetTag':'container',
 },
 'getLocalKey':{
 'group':'query',
 'targetTag':'container',
 },
 'otherQuery':{
 'group':'query',
 'targetTag':'masterOp',
 },
 'otherModify':{
 'group':'modify',
 'targetTag':'masterOp',
 },
 'otherCreate':{
 'group':'create',
 'targetTag':'masterOp',
 },
 'otherDelete':{
 'group':'delete',
 'targetTag':'masterOp',
 },
 'other':{
 'group':'other',
 'targetTag':'masterOp',
 },
}

# OperationData subOps and parameters - given separately for readability
#
# Note that parameters named 'parent' or 'self' are given as target
# op.container.parentRole and op.container, respectively.
# Other paremeters without a target are given op.target by default

# NBNB TBD recheck mechanism for setting default types of parameters

# get
operationData['get']['subOps'] = {
 None:{
  # 'suffix':'var1', # suffix defaults to the subOpType
  #                  # or to empty string if that is None
                     # - only given where different from default
  # 'taggedValues':{'key:'value', 'key:'value'}
  # 'opDocTemplate':"some text %s other text",
  'parameters':[
   {'name':'result', 'direction':ImpConstants.return_direction,
    'useCollection':True,
    # 'parDocumentation':"some text, maybe with %(value)s formatting"
    #                     #making use of the self.varNames dictionary
    # 'target':'memops.Implementation.String', # dataType or class, qualName. 
                                               # May be None.
    #'hicard':2, 'locard':1, isOrdered:True, isUnique:True,
    # 'defaultValue':aValue
   },
  ],
 },
}

operationData['findFirst']['subOps'] = {
 None:{
  'parameters':[
   {'name':'conditions', 'direction':ImpConstants.in_direction,
    'target':'memops.Implementation.StringKeyDict',
    'taggedValues':{'isSubdivided':'True'},
   },
   {'name':'result', 'direction':ImpConstants.return_direction,
    'parDocumentation':"the first %(value)s that satisfies the %(conditions)s ",
   },
  ],
 },
}

operationData['findAll']['subOps'] = {
 None:{
  'parameters':[
   {'name':'conditions', 'direction':ImpConstants.in_direction,
    'target':'memops.Implementation.StringKeyDict',
    'taggedValues':{'isSubdivided':'True'},
   },
   {'name':'result', 'direction':ImpConstants.return_direction,
    'parDocumentation':"all %(value)s that satisfy the %(conditions)s ",
    'locard':0, 'useCollection':True,
   },
  ],
 },
}

# sorted NBNB TBD see how the parameter is generated
operationData['sorted']['subOps'] = {
 None:{
  'parameters':[
   {'name':'result', 'direction':ImpConstants.return_direction,
    'isOrdered':True, 'isUnique':True, 
    'locard':0, 'hicard':genConstants.infinity,
   },
  ],
 },
}

# set
operationData['set']['subOps'] = {
 None:{
  'parameters':[
   {'name':'value', 'direction':ImpConstants.in_direction,
    'useCollection':True,
    'parDocumentation':"the %(value)s(s) to set",
   },
  ],
 },
}

# add
operationData['add']['subOps'] = {
 None:{
  'parameters':[
   {'name':'value', 'direction':ImpConstants.in_direction,
    'parDocumentation':"the %(value)s to add",},
  ],
 },
}

# remove
operationData['remove']['subOps'] = {
 None:{
  'parameters':[
   {'name':'value', 'direction':ImpConstants.in_direction,
    'parDocumentation':"the %(value)s to remove",},
  ],
 },
}

# new
operationData['new']['subOps'] = {
 None:{
  'opDocTemplate':'Factory function to create %s',
  'parameters':[
   {'name':'attrlinks', 'direction':ImpConstants.in_direction,
    'target':'memops.Implementation.StringKeyDict',
    'taggedValues':{'isSubdivided':'True'},
    'parDocumentation':"the attribute and link parameters",
   },
   {'name':'result', 'direction':ImpConstants.return_direction,
    'parDocumentation':"the new object",
   }
  ],
 },
}

# init
operationData['init']['subOps'] = {
 None:{
  'opDocTemplate':'Constructor for %s',
  'parameters':[
   {'name':'parent', 'direction':ImpConstants.in_direction,
    'parDocumentation':"the parent object",},
   {'name':'attrlinks', 'direction':ImpConstants.in_direction,
    'target':'memops.Implementation.StringKeyDict',
    'taggedValues':{'isSubdivided':'True'},
    'parDocumentation':"the attribute and link parameters",
   },
   {'name':'result', 'direction':ImpConstants.return_direction,
    'parDocumentation':"the new object",
   }
  ],
 },
}

# clone
operationData['clone']['subOps'] = {
 None:{
  'opDocTemplate':'Clone function for %s',
  'parameters':[
   {'name':'attrDict', 'direction':ImpConstants.in_direction,
    'target':'memops.Implementation.StringKeyDict',
    'taggedValues':{'isSubdivided':'True'},
    'parDocumentation':"the attribute parameters",
   },
   {'name':'result', 'direction':ImpConstants.return_direction,
    'parDocumentation':"the new object clone",
   }
  ],
 },
}

# checkValid
operationData['checkValid']['subOps'] = {
 None:{
  'parameters':[
   {'name':'complete', 'direction':ImpConstants.in_direction,
     'parDocumentation':"whether to do a complete but slow check",
     'target':'memops.Implementation.Boolean',
     'defaultValue':False
    },
   ],
 },
}

# checkAllValid
operationData['checkAllValid']['subOps'] = {
 None:{
  'parameters':[
   {'name':'complete', 'direction':ImpConstants.in_direction,
     'parDocumentation':"whether to do a complete but slow check",
     'target':'memops.Implementation.Boolean',
     'defaultValue':False
    },
   ],
 },
}

# getAttr
operationData['getAttr']['subOps'] = {
 None:{
  'parameters':[
   {'name':'name', 'direction':ImpConstants.in_direction,
    'parDocumentation':"the attribute or link name",
    'target':'memops.Implementation.String',
   },
   {'name':'result', 'direction':ImpConstants.return_direction,
    'parDocumentation':"get attribute or link by name",
    'target':'memops.Implementation.Any',
   },
  ],
 },
}

# setAttr
operationData['setAttr']['subOps'] = {
 None:{
  'parameters':[
   {'name':'name', 'direction':ImpConstants.in_direction,
    'parDocumentation':"the attribute or link name",
    'target':'memops.Implementation.String',
   },
   {'name':'value', 'direction':ImpConstants.in_direction,
    'parDocumentation':"value to set attribute or link to",
    'target':'memops.Implementation.Any',
   },
  ],
 },
}

# getFullKey
operationData['getFullKey']['subOps'] = {
 None:{
  'parameters':[
   {'name':'useGuid', 'direction':ImpConstants.in_direction,
     'parDocumentation':"whether to use guid instead of key for TopObject.",
     'target':'memops.Implementation.Boolean', 
     'defaultValue':False
    },
   {'name':'result', 'direction':ImpConstants.return_direction,
     'parDocumentation':"list containing full object key",
     'target':'memops.Implementation.Any', 'hicard':genConstants.infinity,  
     'locard':0, 'isOrdered':True, 'isUnique':False,
    },
   ],
 },
}

# geLocalKey
operationData['getLocalKey']['subOps'] = {
 None:{
  'parameters':[
   {'name':'result', 'direction':ImpConstants.return_direction,
     'parDocumentation':"Local object key",
     'target':'memops.Implementation.Any', 
    },
   ],
 },
}

# getByKey
operationData['getByKey']['subOps'] = {
 None:{
  'parameters':[
   {'name':'startObj', 'direction':ImpConstants.in_direction,
     'parDocumentation':"TopObject that is the start of the search",
     'target':'memops.Implementation.MemopsObject',
    },
   {'name':'fullKey', 'direction':ImpConstants.in_direction,
     'parDocumentation':"list containing full object key",
     'target':'memops.Implementation.Any', 'hicard':genConstants.infinity, 
     'locard':0, 'isOrdered':True, 'isUnique':False,
    },
   {'name':'result', 'direction':ImpConstants.return_direction,
     'parDocumentation':"object corresponding to key",
    },
   ],
 },
}

# fullDelete
operationData['fullDelete']['subOps'] = {
 None:{'parameters':[],},
}

# singleDelete
operationData['singleDelete']['subOps'] = {
 None:{
  'opDocTemplate':
  '''singleDelete for %s:   deletes objects
*Implementation function* - will CORRUPT DATA if called outside the API delete function.''',
  'parameters':[
   {'name':'objsToBeDeleted', 'direction':ImpConstants.in_direction,
    'parDocumentation':"Set of objects to be deleted",
    'target':'memops.Implementation.MemopsObject',
    'hicard':genConstants.infinity, 'locard':0,
    'isOrdered':False, 'isUnique':True,
   },
  ],
 },
}

# checkDelete
operationData['checkDelete']['subOps'] = {
 None:{
  'opDocTemplate':
  '''checkDelete for %s:
   determines cascading deletes to follow from delete of object.
   *Implementation function* - should be called only by API delete function.''',
  'parameters':[
   {'name':'objsToBeDeleted', 'direction':ImpConstants.in_direction,
    'parDocumentation':"Set of objects to be deleted (input/ouput)",
    'target':'memops.Implementation.MemopsObject',
    'hicard':genConstants.infinity, 'locard':0,
    'isOrdered':False, 'isUnique':True,
   },
   {'name':'objsToBeChecked', 'direction':ImpConstants.in_direction,
    'parDocumentation':
    "List of objects to be checked for deletion (input/ouput)",
    'target':'memops.Implementation.MemopsObject',
    'hicard':genConstants.infinity, 'locard':0,
    'isOrdered':True, 'isUnique':False,
   },
   {'name':'linkCounter', 'direction':ImpConstants.in_direction,
    'parDocumentation':
    "Dictionary to track links with locard checks (input/output)",
    'target':'memops.Implementation.Dict',
   },
  ],
 },
}


def getTarget(metaOp, opData=None):
  """ gets target for MetaOperation based on target class and operation name.
  opData overrides the standard operationData, for use e.g. in ModelAdapt.
  NB used within ObjectDomain to make targets for Operations entered in UML
  
  The  converse operation is ModelAdapt.getOperationName
  
  NB 
  - for target 'container' assumes nothing
  - for target 'masterOp' assumes names of the form somename_somesubtypesuffix
  - for ClassElement and ChildClass targets assumes names of the form 
    prefixElementname_somesubtypesuffix
  """
  from memops.metamodel.MetaModel import MemopsError
 
  # set up
  opType = metaOp.opType
  if opData:
    infoDict = opData.get(opType)
  else:
    infoDict = operationData.get(opType)
  if infoDict is None:
    raise MemopsError("%s: operation has unrecognised opType %s"
                      % (metaOp,`opType`))
  targetTag = infoDict['targetTag']
  prefix = infoDict.get('name')
  if prefix is None:
    prefix = opType
 
  # process, depending on targettag
  if targetTag == 'container':
    result = metaOp.container
 
  else:
 
    # get master operation
    if metaOp.opSubType is None:
      masterOp = metaOp
    else:
      # remove '_suffix' from name
      opName = '_'.join(metaOp.name.split('_')[:-1])
      #get equivalent master operation
      masterOp = metaOp.container.getElement(opName)
      if not masterOp:
        raise MemopsError("%s: no master operation (subType:None) found"
                          % (metaOp,))
 
    if targetTag == 'masterOp':
      result = masterOp
 
 
    elif targetTag == 'ChildClass':
      targetName = masterOp.name[len(prefix):]
      result = None
      for childRole in metaOp.container.getChildRoles():
        for childClass in childRole.valueType.getAllSubtypes():
          if childClass.name == targetName:
            if result is None:
              result = childClass
            else:
              raise MemopsError(
               "%s: Ambiguous child class target for operation %s"
               % (metaOp, masterOp.name)
              )
      if result is None:
        raise MemopsError("%s: No child class target for operation %s"
         % (metaOp, masterOp.name)
        )
      
    else:
      # target is a ClassElement
      targetName = uniUtil.lowerFirst(masterOp.name[len(prefix):])
      if infoDict.get('useCollection'):
        result = metaOp.container.getElement(targetName)
        
      else:
        ll =[x for x in metaOp.container.getAllClassElements()
                 if x.baseName == targetName]
        if ll:
          result = ll[0]
        else:
          result = None
 
      if not result:
        raise MemopsError(
         "%s, opType %s, targetTag %s: no target named %s found"
         % (metaOp, opType, targetTag, targetName)
        )
  #
  return result

  ###########################################################################
