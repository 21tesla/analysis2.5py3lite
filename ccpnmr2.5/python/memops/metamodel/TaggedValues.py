
from memops.universal.Constants import trueString, falseString

# used for packages with no packageGroup defined
defaultPackageGroup = 'core'

allowedTags = {
 'MetaModelElement':{
 },
 'ConstrainedElement':{
 },
 'ComplexDataType':{
  'specificImplementation':('db','file'),
 },
 'AbstractDataType':{
  'repositoryId':None,
  'repositoryTag':None,
 },
 'MetaDataObjType':{
 },
 'ClassElement':{
 },
 'AbstractValue':{'originalGuid':None,
 },
 'HasParameters':{
 },
 'MetaAttribute':{
 },
 'MetaClass':{
 },
 'MetaDataType':{
 },
 'MetaOperation':{'originalGuid':None,
 },
 'MetaPackage':{
  'docDiagramNames':None,
  'packageGroup':('nmr','pp',defaultPackageGroup),
  'repositoryId':None,
  'repositoryTag':None,
  'isDraft':(trueString,), # NBNB TBD do we need this?
  'isReferenceData':(trueString,), # NBNB TBD do we need this?
 },
 'MetaParameter':{
  'isSubdivided':(trueString,)  # slightly obscure. Used to identify
                                # collections of undeclared parameters
                                # in languages that support them.
                                # If the valueType is a stringKeyDict 
                                # the parameters are assumed to be keyword-value
                                # else they are assumed to be normal.
                                # In practice used for Python * and ** params
 },
 'MetaRole':{
 },
 'MetaConstraint':{
 },
 'MetaConstant':{
  'repositoryId':None,
  'repositoryTag':None,
 },
 'MetaException':{
  'repositoryId':None,
  'repositoryTag':None,
 }
} 
