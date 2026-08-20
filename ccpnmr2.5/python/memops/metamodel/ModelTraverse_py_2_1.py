"""  Version for python version >= 2.1

Operations on the model in memory, adapting the basic model for a 
particular language and implementation. 

Includes general model querying functions the depend on context.
Includes default version of standard information, such as operations data,
variable names, etc.

Does no output.
"""

from memops.universal import Util as uniUtil
from memops.metamodel import MetaModel
from memops.metamodel import Util as metaUtil
MemopsError = MetaModel.MemopsError
from memops.metamodel import ImpConstants
from memops.metamodel.ModelPortal import ModelPortal

True = not 0
False = not True
  
mandatoryAttributes = ('modelPortal',)

class ModelTraverse_py_2_1:
  """ Tools for traversing a model - loops and calls empty functions
  Also model query functions that may need overriding elsewhere.
  """
    
  # set some general info
  # special metaClasses and metaPackages
  implPackageName = '%s.%s' % (ImpConstants.modellingPackageName,
                               ImpConstants.implementationPackageName)
  dataRootName = '%s.%s' % (implPackageName, ImpConstants.dataRootName)
  baseClassName = '%s.%s' % (implPackageName, ImpConstants.baseClassName)
  rootClassName = '%s.%s' % (implPackageName, ImpConstants.rootClassName)
  baseDataTypeObjName = '%s.%s' % (implPackageName, 
                                   ImpConstants.baseDataTypeObjName)
  topObjClassName = '%s.%s' % (implPackageName, ImpConstants.topObjClassName)
  dataObjClassName = '%s.%s' % (implPackageName, ImpConstants.dataObjClassName)
  implObjClassName = '%s.%s' % (implPackageName, 
                              ImpConstants.implObjClassName)
    
  ###########################################################################
  
  ###########################################################################
  ###
  ### code executing model processing.
  ### To be overridden in subclasses
  ### Does no active work, except for checking the model flavour
  ###
  ###########################################################################

  ###########################################################################
  
  # Only function in this section to be called from ourside the class
  def processModel(self):
    """ top function that does all actions. Only function that should
    be called from outside the module.
    """
    
    # check model flavours
    for (key, val) in self.modelFlavours.items():
      self.modelPortal.checkModelFlavour(key, val)

    pp = self.modelPortal.topPackage
    # Normally we rely on topObjClassName etc. but in this file only
    # this kind of attribute is needed
    self.implPackage = impl = pp.metaObjFromQualName(self.implPackageName)
    self.topObject = impl.getElement(ImpConstants.topObjClassName)
    self.dataRoot = impl.getElement(ImpConstants.dataRootName)
    self.rootClass = impl.getElement(ImpConstants.rootClassName)
    self.baseClass = impl.getElement(ImpConstants.baseClassName)
    self.dataObject = impl.getElement(ImpConstants.dataObjClassName)
    self.implObject = impl.getElement(ImpConstants.implObjClassName)
    self.baseDataType = impl.getElement(ImpConstants.baseDataTypeObjName)
    self.anyObject = impl.getElement('Any')
    
    # branch packages
    for package in self.modelPortal.branchPackages():
      self.processBranchPackage(package)
    
    # leaf packages.
    # Note that import order means that superclasses, datatypes etc.
    # are processed before they are needed
    for package in self.modelPortal.leafPackagesByImport():
      self.initLeafPackage(package)
      self.processLeafPackage(package)
      self.endLeafPackage(package)
    
  ###########################################################################
    
  ###########################################################################

  def processBranchPackage(self, package):
    """ processing actions for branch package
    """
    pass
    
  ###########################################################################

  ###########################################################################
    
  def initLeafPackage(self, package):
    """ processing actions for start of LeafPackage
    """
    pass
    
  ###########################################################################

  ###########################################################################
  
  def processLeafPackage(self, package):
    """  processing actions for body of leaf package
 
    NB programmers note:
    Processing may rely on the model being valid and all subclasses being known
    at the time of processing the superclass.
    Any overriding that adds classes must maintain this.
 
    Also note that DataTypes are needed in many other elements, and that
    exceptions are needed for operations. The processing order matters.
    """
    
    # data types
    # Note that inheritance order means that 
    # supertypes are processed before subtypes
    for xx in self.modelPortal.dataTypesByInheritance(package):
      self.processDataType(xx)
    
    # constants
    for xx in self.modelPortal.constantsAlphabetic(package):
      self.processConstant(xx)
    
    # exceptions
    for xx in self.modelPortal.exceptionsByInheritance(package):
      self.processException(xx)
    
    # dataObjTypes
    for xx in self.modelPortal.dataObjTypesByInheritance(package):
      self.initDataObjType(xx)
      self.processDataObjType(xx)
      self.endDataObjType(xx)
      
    # classes
    for xx in self.modelPortal.classesByInheritance(package):
      self.initClass(xx)
      self.processClass(xx)
      self.endClass(xx)
    
  ###########################################################################
 
  ###########################################################################
    
  def endLeafPackage(self, package):
    """ processing actions for end of LeafPackage
    """
    pass
    
  ###########################################################################

  ###########################################################################
  
  def processDataType(self, dataType):
    """ processing actions for data type
    """
    pass
    
  ###########################################################################

  ###########################################################################
    
  def processConstant(self, constant):
    """ processing actions for constant
    """
    pass
    
  ###########################################################################

  ###########################################################################
    
  def processException(self, exception):
    """ processing actions for exception
    """
    for xx in exception.parameters:
      self.processParameter(xx, exception)
    
    
    
  ###########################################################################

  ###########################################################################
    
  def initClass(self, clazz):
    """ processing actions for start of class
    """
    pass
    
  ###########################################################################

  ###########################################################################
    
  def processClass(self, clazz):
    """ processing actions for body of class
    """
    
    # attributes
    for xx in clazz.getAllAttributes():
      self.processAttribute(xx, clazz)
    
    # roles
    for xx in clazz.getAllRoles():
      self.processRole(xx, clazz)
    
    # operations
    for xx in clazz.getAllOperations():
      self.initOperation(xx, clazz)
      self.processOperation(xx, clazz)
      self.endOperation(xx, clazz)
    
  ###########################################################################

  ###########################################################################
    
  def endClass(self, clazz):
    """ processing actions for end of class
    """
    pass
    
  ###########################################################################

  ###########################################################################
    
  def initDataObjType(self, clazz):
    """ processing actions for start of class
    """
    pass
    
  ###########################################################################

  ###########################################################################
    
  def processDataObjType(self, clazz):
    """ processing actions for body of class
    """
    
    # attributes
    for xx in clazz.getAllAttributes():
      self.processAttribute(xx, clazz)
    
    # operations
    for xx in clazz.getAllOperations():
      self.initOperation(xx, clazz)
      self.processOperation(xx, clazz)
      self.endOperation(xx, clazz)
    
  ###########################################################################

  ###########################################################################
    
  def endDataObjType(self, clazz):
    """ processing actions for end of class
    """
    pass
    
  ###########################################################################

  ###########################################################################
    
  def processAttribute(self, attribute, inClass):
    """ processing actions for attribute
    """
    
    pass
    
  ###########################################################################

  ###########################################################################
    
  def processRole(self, role, inClass):
    """ processing actions for role
    """
    
    pass
    
  ###########################################################################
 
  ###########################################################################
    
  def initOperation(self, op, inClass):
    """ processing actions for start of operation
    """
    pass
    
  ###########################################################################

  ###########################################################################
    
  def processOperation(self, op, inClass):
    """ processing actions for body of operation
    """
    for xx in op.parameters:
      self.processParameter(xx, op)
    
  ###########################################################################
 
  ###########################################################################
    
  def endOperation(self, op, inClass):
    """ processing actions for end of operation
    """
    pass
    
  ###########################################################################

  ###########################################################################
    
  def processParameter(self, parameter, inElement):
    """ processing actions for operation/exception parameter
    """
    
    pass
    
  ###########################################################################
 
  ###########################################################################
  ###
  ### code for calling from subclasses - not intended for overriding
  ###
  ###########################################################################

  ###########################################################################
  
  def addModelFlavour(self, key, val):
    """ Add model flavour to local requirements, checking for conflicts
    """
    
    if not hasattr(self,'modelFlavours'):
      self.modelFlavours = {}
    
    flav = self.modelFlavours.get(key)
    if flav == val:
      return
    elif flav is None:
      self.modelFlavours[key] = val
    else:
      raise MemopsError(
       "attempt to reset modelFlavour %s:%s to %s)" % (key, flav, val)
      )
    
  ###########################################################################
 
  ###########################################################################
  ###
  ### overridable model queries
  ### for use outside the class
  ###
  ###########################################################################

  ###########################################################################
    
  def getFuncname(self, op, inClass=None):
    """ get function calling name for operation
    
    the name is the operation name of the corresponding operation with
    opSubType==None
    """
    return metaUtil.getFuncname(op, inClass)
    
  ###########################################################################

  ###########################################################################

  def valueVar(self, element, doInstance=False, prefix=''):

    var = self.varNames['value']
    if (not doInstance and element.hicard != 1):
      var = self.varNames['values']

    if prefix:
      var = prefix + uniUtil.upperFirst(var)

    return var

  ###########################################################################

  ###########################################################################
    
  def defineFunc(self, op):
    """write function definition
    """
    
    funcName = self.getFuncname(op)
    params = self.getFuncParams(op)
    docString = self.getDocString(op)
    returnType = self.getReturnType(op)
    throwsTypes = self.getThrowsTypes(op)
    # TBD: exceptions
    self.startFunc(funcName, params, docString, returnType, throwsTypes)
    
  ###########################################################################

  ###########################################################################

  def getImportName(self, oo, pp=None, subDirs=None):
    
    if subDirs is None:
      if self.codeDirName is None:
        subDirs = []
      else:
        subDirs = [self.codeDirName]
        
    if oo is self.modelPortal.topPackage:
      # root package - special case - use memops package instead
      oo = oo.getElement(ImpConstants.modellingPackageName)
      
        
    return metaUtil.getReferenceName(oo, pp, subDirs=subDirs)

  ###########################################################################
  
  ###########################################################################
  
  def getPartitioningLinkages(self, inputRole):
    """ get chains of roles leading upwards from inputRole.container
    and inputRole.valueType to reach an object that must be shared from both
    sides. 
    Return values are (uplinks, partitionRoles)
    
    Examples: if <xyz> indicates a MetaRole names 'xyz', we have  :
    
    For ccp.molecule.Molecule.MolResLinkEnd.molResLink, the return is
    [[<molResidue>,<molecule>],[<molecule>]],[None,None]
    
    while for ccp/nmr.NmrConstraint.AbstractConstraintList.experiments it is
    [[<nmrConstraintStore>],[<nmrProject>]],[<nmrProject>,<nmrConstraintStores>])
    
    Note that the roles in uplinks[0], uplinks[1] and *either*  
    partitionRoles[0] *or* partitionRoles[1] (if either is set) combine to make 
    a chain of links that connect inputRole.container and inputRole.valueType
    
    For MetaRoles contained in classes without a parentRole, the uplinks are
    [<topObject>,<memopsRoot>]
    
    NB - as written this code checks that:
    - no partitioning links between parent and descendant
    - No two partitioning links (inherited or not) between same pair of classes 
    - If two partitioning links span the same up-down path in the parent tree,
      one link cannot be lower in the tree in one branch and higher in the other
    Doing the checks in MetaModel would be better (and faster) but would
    require partially duplicating some quite complicated code.
    The code will be slow, considering that it is called for every writeModifier
    and some getValue for roles.   
    """
    
    # get list of uplinks and list of the classes they point to
    startClasses = (inputRole.container,inputRole.valueType)
    classes = ([],[])
    uplinks = ([],[])
    for ic in (0,1):
      ul = uplinks[ic]
      cl = classes[ic]
      cc = startClasses[ic]
      pr = cc.parentRole
      if pr is None:
        # class has no parent link 
        # the best we can do is to go directly to the TopObject
        # we will miss some things, but it cannot be helped
        ul.append(cc.getElement(self.varNames['topObject']))
        topObjClass = cc.container.topObjectClass
        cl.append(topObjClass)
        pr = topObjClass.parentRole
      while pr:
        ul.append(pr)
        cc = pr.valueType
        cl.append(cc)
        pr = cc.parentRole
    
    # Find index of lowest common partitioning class
    # - loops starts from end of list
    # NB, MemopsRoot is partitioning so we always get one
    sentinel = -min(len(classes[0]), len(classes[1]))
    ii = -1
    while ii >= sentinel:
      cc = classes[0][ii]
      if cc is classes[1][ii]:
        if cc.partitionsChildren:
          indexcc = ii
        ii = ii - 1
      else:
        break

    # get lists of non-shared parents
    separate = (classes[0][:indexcc], classes[1][:indexcc])

    # find partitioning links, if any
    partitionRoles = [None,None]
    indices = [[None, None],[None, None]]
    for ic in (0,1):
      spx = separate[ic]
      spz = separate[1-ic]
      for ii in range(len(spx)):
        cc = spx[ii]
        found = False
        for role in cc.getAllRoles():
          if role.partitionsChildren:
            valueType = role.valueType
            for jj in range(len(spz)):
              cc2 = spz[jj]
              if valueType in cc2.getAllSupertypes():
                if (indices[ic][0] == ii and indices[ic][1] == jj):
                  raise MemopsError(
                   "Partitioning roles %s and %s incompatible for %s"
                   % (partitionRoles[ic], role, inputRole)
                  )
                else:
                  found = True
                  indices[ic][0] = ii
                  indices[ic][1] = jj
                  partitionRoles[ic] = role
                  break
            else:
              for cc3 in classes[ic][ii:]:
                if valueType in cc3.getAllSupertypes():
                  raise MemopsError(
                   "Partitioning role %s connects ancestor and descendant"
                   % (role,)
                  )
              continue
            break
      else:
        continue
      break
     
    # select partitioning link, get both roles, and further check compatibility
    if partitionRoles[0] is None and partitionRoles[1] is None:
      # no link 
      resultIndices = (indexcc, indexcc)
    
    elif partitionRoles[0] is None:
      # one-way link
      resultIndices = indices[1]
      resultIndices.reverse()
    
    elif partitionRoles[1] is None:
      # one-way link
      resultIndices = indices[0]
      
    elif indices[0][0] <= indices[1][1] and indices[0][1] <= indices[1][0]:
      # link starting on 0 has priority (might be the same link, of course)
      partitionRoles[1] = partitionRoles[0].otherRole
      resultIndices = indices[0]
    
    elif indices[0][0] >= indices[1][1] and indices[0][1] >= indices[1][0]:
      # link starting on 1 has priority
      partitionRoles[0] = partitionRoles[1].otherRole
      resultIndices = indices[1]
      resultIndices.reverse()
      
    else:
      raise MemopsError("Partitioning roles %s and %s are incompatible"
                        % (partitionRoles[0], partitionRole[1]))
    
    # set uplinks lists
    result = ([],partitionRoles)
    for ic in (0,1):
      
      topObjRole = startClasses[ic].getElement(self.varNames['topObject'])
      topObjClass = startClasses[ic].container.topObjectClass
      ir = resultIndices[ic]
      cl = classes[ic]
      if cl[ir] is cl[-1]:
        # lowest parent is last element - Memopsroot.
        result[0].append([startClasses[ic].getElement(self.varNames['root'])])
      elif cl[ir] is topObjClass:
        # lowest parent is TopObject. Go there directly
        result[0].append([topObjRole])
      else:
        # normal case
        result[0].append(uplinks[ic][:ir+1])
    #
    return result

  ###########################################################################

  ###########################################################################
  
  #  Convenience functions, calling other Language functions only
  # (Language-independent)

  ###########################################################################

  ###########################################################################

  def valueIsNone(self, var):

    return self.comparison(var, 'is', self.noneValue)

  ###########################################################################

  ###########################################################################

  def valueIsNotNone(self, var):
  
    return self.comparison(var, 'is not', self.noneValue)

  ###########################################################################
