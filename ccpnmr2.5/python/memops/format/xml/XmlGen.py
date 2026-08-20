"""
======================COPYRIGHT/LICENSE START==========================

XmlGen.py: Code generation for CCPN framework

Copyright (C) 2005 Wayne Boucher (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../license/GPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
 
You should have received a copy of the GNU General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

Top level classes for generation of I/O maps.

Subclasses in other files will then complete the I/O maps and write them, or 
derived products like XML schemas, to disk. Note that the mapping dictionaries
generated in these classes are identical in structure to those written to disk
(in Python) and could be used directly for XML I/O once language specific 
extensions were added


Structure of I/O maps

A top level dictionary (globalMap) contains:
 
 package specific mapping dictionaries for every package 
  (key is package.shortName)

 A dictionary 'loadMaps' of XML tag : element I/O map

   This has maps for all legal XML tags.
 
   Note that all XML tags are unique across the entire model.
   They are of the following forms:
   - PPPP.Cccccc for AbstractTypes, where the prefix PPPP is the
     package.shortName and Cccccc is the AbstractType.name
   - PPPP.Cccccc.aaaa for ClassElements, where PPPP and Cccccc are as above,
     and aaaa is the ClassElement.name
   - PPPP.exo-Cccccc for exolinks. The tqg describes an out-of-package link *to*
     Class Cccccc in package with shortName PPPP.

 A dictionary 'mapsByGuid' of model element guid : element I/O map
   which is used for model comparison and for crossreferences between I/O maps.
   These maps include all abstract types and all attributes or links,
   but not the exolinks descriptors
  
The maps from different packages refer to each others contents, 
importing package to imported package, and Implementation package to all others.

All dictionaries contain enough information to allow you to write to file
a reference to the mapping dictionary given the globalMap

Every package specific mapping contains:
- 'globalRelease' - the global release version
- '.prefix' - the package prefix (== package.shortName)
- '.qName' - the package qualifiedName
- '.name':'mapping' (fixed value) - describes the kind of dictionary 
   for use in subsequent generation
- '.isReference' - optional Boolean 'are these a reference data?' 
- 'abstractTypes' - the AbstractTypes dictionary
- 'exolinks' - the dictionary of out-of-package link descriptors


the abstractTypes dictionary contains 
- '.prefix' - the package prefix (== package.shortName)
- '.qName' - the package qualifiedName
- '.name':'abstractTypes' (fixed value) - describes the kind of dictionary 
It further contains an I/O map for every AbstractType 
(Class (also abstract), DataObjType (also abstract) or 
DataType (not abstract) ). The dictionary key is the 
AbstractType.name.


the exolinks dictionary contains 
- '.prefix' - the package prefix (== package.shortName)
- '.qName' - the package qualifiedName
- '.name':'exolinks' (fixed value) - describes the kind of dictionary 
It further contains an I/O map for out-of-package links to every 
non-abstract class in the package. The dictionary key is the Class.name.


All I/O maps contain the following:
- 'tag': The XML tag to use when writing. May differ from the corresponding 
  key in the 'loadMaps' which is used for reading.
- 'guid': The guid for the corresponding model element
- 'type' the element type, which is one of the following:
  'class', 'cplx', 'simple', 'exo', 'attr', 'link', 'child', 
  'exotop', 'exolink', 'dobj'
- 'proc': Optional. Allowed values differ, but all maps allow the value 'skip'.
  'Skip is used only for compatibility mappings, and indicates that the element
  is ignored.
  Further content depends on the type - some are a bit arcane, they are chosen
  to make XML I/O fast, not to make the description simple:


- type 'cplx' - for DataObjTypes (and 'class' - for Classes) further has
 'eType' : 'cplx' (fixed value)
 'class' : The actual API class corresponding to the DataObjType (or Class)
           Also serves as a constructor.
           If set to 'skip'the element and its contents will be skipped
 'content' : a dictionary of ClassElement.name : I/O map for classElement
 'headerAttrs' : A list of attribute names for attributes that
                 are written as XML attributes in the Class XML element.
 'simpleAttrs' : A list of ClassElement names for attributes with DataType
                 valueType and intrapackage links that are not included in
                 'headerAttrs'
  'cplxAttrs'  : A list of ClassElement names for child links, attributes of
                 DataObjType valueType, and out-of-package links

  For use in version conversion, two more elements are optional:
  'cnvrt': a function that takes an object as a parameter. It is called on
           the relevant object(s) after reading and link dereferencing, but
           before validity checking. The function converts from the then 
           (possibly intermediate) state of the object to the final state.
  'priority' : An integer that determines the order (lowest value first) 
            in which the cnvrt functions are called. Intended to ensure
            that later conversions can rely on earlier conversions to be 
            finished.

- type 'class' - has all the tags used for type 'cplx', plus the following
   'children': A list of ClassElement names for child links
   'optLinks': A list of role names for intrapackage  roles where 
               role.hicard == 1 and role.otherRole.hicard != 1.
               These may be treated as headerAttrs, simpleAttrs, or skipped,
               depending on the settings of the save function.
   'fromParent': equal to clazz.parentRole.otherRole.name
   'objkey' :  Optional. The name of the single simple attribute that forms
               the local key of the class. Where appropriate it allows 
               faster object creation
   'isTop'   : Optional. True (fixed value).Is this a topObject? 
   'singleKid' : Optional. True (fixed value).Is this an only child? 

type 'simple' - for DataTypes, further has toStr, cnvrt, toEtype
  'toStr' : a function to convert a value to a string, 
   or the string 'text' for string types
  'cnvrt' : a function to convert a value from a string, 
   or the string 'text' for string types 
  'toEtype' : optional. == 'cplx' (fixed value)
            toEtype == 'cplx' is set for DataTypes 
            if attributes of that type must have eType == 'cplx'
             - in practice strings that can contain whitespace or
            that are long or unlimited in length.


type 'exo' - for exolinks
 'eType' : 'cplx' (fixed value)
 'class' : The actual API class corresponding to the Class
 'keyMaps': a list of I/O maps that describes the fullKey of the object -
            used for writing.
            For key elements of DataType type the I/O map of the dataType is 
            used. For key elements of DataObjType type, the mapping of the
            relevant package is used. For key elements of type Class the
            exolinks directory of the relevant package is used.


Finally there are the types for ClassElements - attributes and links


all maps have tag, type and guid
Further (square brackets may be needed for schemas, not used in I/O

type 'attr' - attributes of DataType valueType
 'hicard'   : attr.hicard
 'locard'   : attr.locard. Not needed for I/O but for subtree copying and 
              maybe be for XML schemas.
 'name'     : attr.name. Used for writing XML attributes and for setting the
              result when reading
 'eType'    : Optional. 'cplx' (fixed value). If set attribute is read and
              written with an XMl element containing other elements. Used (e.g.)
              for DataTypes that may contain whitespace (like String).
 'data'     : The I/O map dictionary of the valueType. Must be type 'simple'
 'default'  : attribute default value.
 'implSkip' : Optional. True (fixed value). If set, do not write this attribute
              if the package being written is Implementation. Used for writing
              TopObjects in Implementation package.
 'isDerived': attr.isDerived. Not needed for I/O but may be so for XML schemas.
              Optional. Only set if true.
 'proc'     : Optional. Processing information, 'direct', 'delay', or 'skip'.
              If proc=='direct' the attribute is set into the object bypassing 
              the API. Can only happen when hicard == 1 and there are no 
              constraints.
              If proc == 'delay' the attribute is converted and set into the
              compDataDict dictionary for later conversion - Backwards
              compatibility only
              If proc == 'skip' the element and its contents will be skipped.
              If proc is not set the value is converted if necessary 
              and set through a setter.

type 'link' - intrapackage links
 'hicard'   : link.hicard
 'locard'   : link.locard. Not needed for I/O but for subtree copying and 
              maybe be for XML schemas.
 'name'     : link.name. Used for writing XML attributes and for setting the
              result when reading

 'implSkip' : Optional. True (fixed value). If set, do not write this link
              if the package being written is Implementation. Used for writing
              TopObjects in Implementation package.
 'isDerived': link.isDerived. Not needed for I/O but may be so for XML schemas.
              Optional. Only set if true.
 'proc'     : Optional. Processing information, proc=='delay'  or 'skip'
              If proc == 'delay' the linkIDs are set into the compDataDict
              dictionary for later conversion - Backwards compatibility only.
              If proc == 'skip' the element and its contents will be skipped.
              If proc is not set the value is put into the crossLinkData
              dictionary for later dereferencing.

type 'child' - child links
 'hicard'   : link.hicard
 'locard'   : link.locard. Not needed for I/O but may be so for XML schemas.
 'eType'    : 'cplx' (fixed value). 
 'content'  : The abstractDataTypes dictionary for the package containing the 
              child class.
 'implSkip' : Optional. True (fixed value). If set, do not write this link
              if the package being written is Implementation. Used for writing
              TopObjects in Implementaion package.
 'proc'     : Optional. 'loadDelayed' or 'skip'. 
              If proc == 'loadDelayed'If set, the end of this  link will trigger
              a call to loadDelayedData.  Needed for the special case of reading
              Implementation.xml, where certain objects must be set and linked
              before the rest of the  data can be processed.
              If proc == 'skip' the element and its contents will be skipped.

types 'exolink' (out of package), 'exotop' (out-of-package TopObject), 
      'dobj' (DataObjType
All have the same attributes: 
 'hicard'   : attr.hicard
 'locard'   : attr.locard. Not needed for I/O but for subtree copying and 
              maybe be for XML schemas.
 'name'     : attr.name. Used for setting the result when reading
 'eType'    : 'cplx' (fixed value). 
 'content'  : The abstractDataTypes dictionary for the package containing the 
              DataObjType (for DataObjType attributes) or the exolinks
              dictionary for the package containing the valueType 
              (for exolinks).
 'implSkip' : Optional. True (fixed value). If set, do not write this link
              if the package being written is Implementation. Used for writing
              TopObjects in Implementation package.
 'isDerived': attr.isDerived. Not needed for I/O but may be so for XML schemas.
              Optional. Only set if true.
 'proc'     : Optional. Processing information, proc=='delay'  or 'skip'
              If proc == 'delay' the linkIDs are set into the compDataDict
              dictionary for later conversion - Backwards compatibility only.
              If proc == 'skip' the element and its contents will be skipped.
              If proc is not set the value is put into the crossLinkData
              dictionary for later dereferencing.


"""
from memops.metamodel import MetaModel
MemopsError = MetaModel.MemopsError
from memops.metamodel import ImpConstants

import memops.general.Constants as genConstants

from memops.metamodel.ModelTraverse import ModelTraverse

headerTextTypeNames = ['Token','Word']
implPackageName = '%s.%s' % (ImpConstants.modellingPackageName,
                             ImpConstants.implementationPackageName)
dataObjClassName = '%s.%s' % (implPackageName,
                             ImpConstants.dataObjClassName)


mandatoryAttributes = ()

repositoryTag = '$Name:  $'
repositoryId  = '$Id: XmlGen.py,v 1.10.2.2 2011/05/04 14:03:57 rhfogh Exp $'

# Requires other writers also to be implemented in subclass
class XmlGen(ModelTraverse):

  ###########################################################################

  ###########################################################################

  def __init__(self):
    
    # init handling
    super(XmlGen, self).__init__()
    pp = self.modelPortal.topPackage
    self.implPackage = pp.metaObjFromQualName(self.implPackageName)
    #self.baseClass = self.implPackage.getElement(ImpConstants.baseClassName)
    self.booleanType = self.implPackage.getElement(genConstants.boolean_code)
    self.stringType = self.implPackage.getElement(genConstants.string_code)
    self.intType = self.implPackage.getElement(genConstants.int_code)
    self.anyType = self.implPackage.getElement('Any')
    self.dataObject = self.implPackage.getElement(ImpConstants.dataObjClassName)
    self.topObject = self.implPackage.getElement(ImpConstants.topObjClassName)
    self.dataRoot = self.implPackage.getElement(ImpConstants.dataRootName)
    self.dictType = self.implPackage.getElement('StringKeyDict')
    
    for tag in mandatoryAttributes:
      if not hasattr(self, tag):
        raise MemopsError(" XmlGen lacks mandatory %s attribute" % tag)
        
    # has to be done this way to allow for different initialisation orders
    if not hasattr(self, 'modelFlavours'):
      self.modelFlavours = {}
    

  ###########################################################################

  ###########################################################################
  
  # overrides ModelTraverse
  def initLeafPackage(self, package):
    """ Set up dictionaries and write topObject maps
    """
    
    self.prefix = package.shortName
    
    if package is self.implPackage:
      
      self.globalMap = dd = {}
      dd['mapsByGuid'] = {}
      dd['loadMaps'] = {}
      
      for pp in self.modelPortal.leafPackagesByImport():
    
        prefix = pp.shortName
        qName = pp.qualifiedName()
        self.globalMap[prefix] = mapping = {}
 
        # set top level map
        mapping['abstractTypes'] = abstractTypes = {'.prefix':prefix,
                                                    '.qName':qName,
                                                    '.name':'abstractTypes'}
        mapping['exolinks'] = exolinks = {'.prefix':prefix,
                                          '.qName':qName, 
                                          '.name':'exolinks'}
        mapping['guid'] = pp.guid
        mapping['.prefix'] = prefix
        mapping['.qName'] = qName
        mapping['.name'] = 'mapping'
        mapping['globalRelease'] = genConstants.currentModelVersion

  ###########################################################################

  ###########################################################################
  
  # overrides ModelTraverse
  def processLeafPackage(self, package):
    """ make data dictionaries for package containing actual code
    """
    
    # data types
    # Note that inheritance order means that 
    # supertypes are processed before subtypes
    for xx in self.modelPortal.dataTypesByInheritance(package):
      if not xx.isAbstract:
        self.processDataType(xx)
    
    # dataObjTypes
    for xx in self.modelPortal.dataObjTypesByInheritance(package):
      self.processDataObjType(xx)
      self.endComplexDataType(xx)
      
    # classes
    for xx in self.modelPortal.classesByInheritance(package):
      self.processClass(xx)
      self.endComplexDataType(xx)
    
    # exolink maps
    if package is not self.implPackage:
      for xx in self.modelPortal.classesByInheritance(package):
        if not xx.isAbstract:
          self.processExoLink(xx)
    
  ###########################################################################

  ###########################################################################
    
  # overrides ModelTraverse
  def endLeafPackage(self, package):
    """ processing actions for end of LeafPackage
    """
    
    del self.prefix
    
  ###########################################################################

  ###########################################################################
  
  # overrides ModelTraverse
  def processDataType(self, xx):
    """ toStr, cnvrt, and baseType postponed, being language dependent
    """
    
    guid = xx.guid
    tag = xmlTagFromElem(xx)
    result = {'tag':tag, 'type':'simple', 'guid':guid,}
    
    if xx.typeCodes.get('python') == genConstants.string_code:
      for xx2 in xx.getAllSupertypes():
        if xx2.name in ('Token', 'GuidString'):
          # Token, spaceless max 32 chars 
          # or GuidString, spaceless max 80 chars
          break
        if xx2.name == 'SpacelessString':
          if xx.length or (xx.enumeration and not xx.isOpen):
            # spaceless, limited length or closed enumeration
            break
      else:
        # string, none of the above - use complex element
        result['toEtype'] = 'cplx'
        
    
    self.globalMap[self.prefix]['abstractTypes'][xx.name] = result
    self.globalMap['mapsByGuid'][guid] = result
    self.globalMap['loadMaps'][tag] = result
    
  ###########################################################################

  ###########################################################################
  
  # overrides ModelTraverse
  def processDataObjType(self, xx):
    
    guid = xx.guid
    
    result = {'type':'cplx', 'eType':'cplx', 'guid':guid, 'content':{}, 
              'headerAttrs':[], 'simpleAttrs':[], 'cplxAttrs':[],}
    
    self.globalMap[self.prefix]['abstractTypes'][xx.name] = result
    self.globalMap['mapsByGuid'][guid] = result
    
    if not xx.isAbstract:
      tag = xmlTagFromElem(xx)
      result['tag'] = tag
      self.globalMap['loadMaps'][tag] = result
    
  ###########################################################################

  ###########################################################################
  
  # overrides ModelTraverse
  def endComplexDataType(self, xx):
    
    self.processClassAttrs(xx)
    if isinstance(xx, MetaModel.MetaClass):
      self.processClassRoles(xx)
      
      if self.dataRoot in xx.getAllSupertypes():
        classMap = self.globalMap[self.prefix]['abstractTypes'][xx.name]
        cplxAttrs = (classMap['cplxAttrs'])
        
        # these roles must be stored and loaded before anything else
        ic = ImpConstants
        cplxAttrs.remove(ic.packageLocatorRole)
        cplxAttrs.remove(ic.repositoryRole)
        cplxAttrs[:0] = [ic.repositoryRole, ic.packageLocatorRole]
        
        
  ###########################################################################

  ###########################################################################
  
  # overrides ModelTraverse
  def processClass(self, xx):
    
    guid = xx.guid
    
    result = {'type':'class',  'eType':'cplx', 'guid':guid, 'content':{},
     'headerAttrs':[], 'simpleAttrs':[], 'optLinks':[], 'cplxAttrs':[], 
     'children':[],
    }
    
    pr = xx.parentRole
    if pr is not None:
      
      if self.topObject in xx.getAllSupertypes():
        result['isTop'] = True
        
      if not pr.isAbstract:
        result['fromParent'] = pr.otherRole.name
      
        if pr.otherRole.hicard == 1:
          result['singleKid'] = True
    
    ll = xx.keyNames
    if len(ll) == 1:
      tag = ll[0]
      if xx.getElement(tag).hicard == 1:
        result['objkey'] = tag
    
    self.globalMap[self.prefix]['abstractTypes'][xx.name] = result
    self.globalMap['mapsByGuid'][guid] = result
    
    if not xx.isAbstract:
      tag = xmlTagFromElem(xx)
      result['tag'] = tag
      self.globalMap['loadMaps'][tag] = result
    
  ###########################################################################

  ###########################################################################
  
  def processClassAttrs(self, clazz):
    """ process atttributes for class *or* DataObjType
    """
    
    mapping = self.globalMap[self.prefix]
    classMap = mapping['abstractTypes'][clazz.name]
    classContent = classMap['content']
    
    for elem in clazz.getAllAttributes():
    
      elemMap = self.processAttribute(elem, clazz)
      
      if elemMap is not None:
      
        # connect up
        classContent[elemMap['name']] = elemMap
        
        # add to appropriate list for writing order
        if not elemMap.get('isDerived'):
        
          typ = elemMap['type']
          
          if typ == 'dobj':
            # typ == 'dobj' DataObjType
            classMap['cplxAttrs'].append(elemMap['name'])
          
          else:
            # typ == 'attr', simple attribute
            if elem.hicard != 1 or elemMap.get('eType') == 'cplx':
              classMap['simpleAttrs'].append(elemMap['name'])
            else:
              classMap['headerAttrs'].append(elemMap['name'])
    
  ###########################################################################

  ###########################################################################
  
  def processClassRoles(self, clazz):
    """ process roles for class
    """
    mapping = self.globalMap[clazz.container.shortName]
    classMap = mapping['abstractTypes'][clazz.name]
    classContent = classMap['content']
    
    children = classMap['children']
    
    for elem in clazz.getAllRoles():
    
      elemMap = self.processRole(elem, clazz)
      
      if elemMap is not None:
      
        # connect up
        classContent[elemMap['name']] = elemMap
        
        # add to apprropriate list for writing order
        if not elemMap.get('isDerived'):
        
          typ = elemMap['type']
 
          if typ == 'link':
            # intrapackage link
            otherRole = elem.otherRole
            if elem.hicard == 1 and otherRole and otherRole.hicard != 1:
              classMap['optLinks'].append(elemMap['name'])
            else:
              classMap['simpleAttrs'].append(elemMap['name'])
 
          elif typ == 'child':
            # intrapackage child link
            children.append(elemMap['name'])

          else:
            # exolink
            classMap['cplxAttrs'].append(elemMap['name'])
        
    classMap['cplxAttrs'].extend(children)

  ###########################################################################

  ###########################################################################
  # overrides ModelTraverse
  
  def processAttribute(self, elem, inClass):
    
    globalMap = self.globalMap
    
    result = getAttrMap(elem, inClass, globalMap)
    
    if result is not None:
      # connect up
      globalMap['mapsByGuid'][elem.guid] = result
      globalMap['loadMaps'][result['tag']] = result
    #
    return result

  ###########################################################################

  ###########################################################################
  # overrides ModelTraverse
  
  def processRole(self, elem, inClass):
    
    globalMap = self.globalMap
    
    result = getRoleMap(elem, inClass, globalMap)
    
    if result is not None:
      # connect up
      globalMap['mapsByGuid'][elem.guid] = result
      globalMap['loadMaps'][result['tag']] = result
    #
    return result
    
  ###########################################################################

  ###########################################################################
  
  def processExoLink(self, clazz):
    """
    """
    
    globalMap = self.globalMap
    
    tag = xmlTagFromElem(clazz, var='exo')
    guid = clazz.guid
    
    keyMaps = []
    
    result = {'tag':tag, 
              'name':clazz.name,
              'type':'exo', 
              'eType':'cplx', 
              'guid':guid,
              'keyMaps':keyMaps}
    
    # get parent list
    # set up. First parent list
    parents = [clazz]
    for pp in parents:
      pr = pp.parentRole
      if pr is None:
        break
      else:
        parents.append(pr.valueType)
    
    # remove MemopsRoot
    parents.pop()
    
    topObject = parents.pop()
    keyElements = [topObject.getElement('guid')]
    
    parents.reverse()
    
    for pp in parents:
      for ss in pp.keyNames:
        elem = pp.getElement(ss)
        keyElements.append(elem)
    
    for elem in keyElements:
    
      vt = elem.valueType
 
      packageMap = globalMap[vt.container.shortName]
 
      if isinstance(vt, MetaModel.MetaClass):
        # exolink
        elemMap = packageMap['exolinks']
 
      elif isinstance(vt, MetaModel.MetaDataObjType):
        # MetaDataObjType
        elemMap = packageMap['abstractTypes']
 
      else:
        # simple type
        elemMap = packageMap['abstractTypes'][vt.name]
 
      keyMaps.extend(elemMap for dummy in range(elem.hicard))
    
    self.globalMap[self.prefix]['exolinks'][clazz.name] = result
    self.globalMap['loadMaps'][tag] = result
    
    #
    return result
    
###########################################################################

# generally useful functions . Used e.g. for comaptibility generation

###########################################################################
  
def xmlTag(prefix, typeName, elemName=None, var=None):
  """ make xml tag from component strings
  """
  sep = '.'
  
  # package short5name prefix
  ll = [prefix]
  
  # AbstractDataType field
  if var is None:
    ll.append(typeName)
  else:
    ll.append('%s-%s' % (var, typeName))
  
  # element name
  if elemName is not None:
    ll.append(elemName)
  
  #
  return sep.join(ll)
  
    
###########################################################################

###########################################################################
  
  
def xmlTagFromElem(elem, var=None):
  """ make xml tag for element corresponding to MetaObject
  """
  
  if isinstance(elem, MetaModel.AbstractDataType):
    clazz = elem
    elemName = None
  else:
    # assert isinstance(elem, MetaModel.ClassElement)
    clazz = elem.container
    elemName = elem.name
  
  return xmlTag(clazz.container.shortName, clazz.name, elemName, var)
  
    
###########################################################################

###########################################################################

def getClassElemMap(elem, globalMap):
  """ Get I/O map for class element
  """
  
  result = {}
  inClass = elem.container
 
  if elem.isImplementation or elem.isAbstract:
    return None
 
  # derived elements
  if elem.isDerived:
    if elem.changeability == ImpConstants.frozen:
      return None
    else:
      result['isDerived'] = True
  
  name = elem.name
  result['name'] = name
  result['tag'] = xmlTagFromElem(elem)
  result['guid'] = elem.guid
  result['hicard'] = elem.hicard
  result['locard'] = elem.locard
  
  if inClass.__class__.__name__ == 'MetaClass':
  
    # proc:delay for derived classes (for backwards compatibility)
    if inClass.isDerived:
      #if name not in inClass.keyNames:
      result['proc'] = 'delay'
      result['isDerived'] = True
      
      
    # for topObject shells in Impl package only guid and keys are stored.
    parent = inClass.parentRole
    if (inClass.qualifiedName() == dataObjClassName or 
        parent and parent.valueType.container is not inClass.container):
      if isinstance(elem, MetaModel.MetaRole):
        vpr = elem.valueType.parentRole
        if name in inClass.keyNames:
          # role part of key - do not skip
          pass
        elif vpr and vpr.valueType.container is not elem.valueType.container:
          # link to other TopObject - do not skip
          pass
        else:
          #none of the above. Skip for Topobject in impl package
          result['implSkip'] = True
  
  # 
  return result
    
###########################################################################

###########################################################################

def getAttrMap(elem, inClass, globalMap, valueTypeGuid=None):
  """ Get Attribute-specific parts of ClassElement I/O map
  """
  
  if elem.container is inClass:
    # make new map
  
    result = getClassElemMap(elem, globalMap)
 
    if result is not None:
 
      # set up
      valueType = elem.valueType
      vtPrefix = valueType.container.shortName
      isDataType = valueType.__class__.__name__ == 'MetaDataType'
      hicard = elem.hicard
      
      if valueTypeGuid is None:
        valueTypeGuid = valueType.guid
 
      # set type
      if isDataType:
        result['type'] = 'attr'
        dd = globalMap['mapsByGuid'][valueTypeGuid]
        result['data'] = dd
        toEtype = dd.get('toEtype')
        if toEtype:
          result['eType'] = toEtype
          
      else:
        # complex data type
        result['type'] = 'dobj'
        result['eType'] = 'cplx'
        result['content'] = globalMap[vtPrefix]['abstractTypes']
 
      # default value
      default = elem.defaultValue
      if default:
        if hicard == 1:
          result['default'] = default[0]
        elif elem.isUnique:
          result['default'] = set(default)
        else:
          result['default'] = default
 
      # processing tag
      if hicard == 1 and not elem.isDerived:
        if isDataType:
          if (not valueType.getAllConstraints() and valueType.length is None
              and (valueType.isOpen or not valueType.enumeration)
              and elem.name != ImpConstants.serial_attribute):
            # primitive, unconstrained type. set directly into memory
            result['proc'] = 'direct'
 
  else:
    # map probably exists - get it
    result = globalMap['mapsByGuid'].get(elem.guid)
    
  #
  return result
    
###########################################################################

###########################################################################

def getRoleMap(elem, inClass, globalMap):
  """ Get Role-specific parts of ClassElement I/O map
  """
  
  if elem.container is inClass:
    # make new map
  
    result = getClassElemMap(elem, globalMap)
 
    if result is not None:
 
      # set up
 
      valueType = elem.valueType
      vtPrefix = valueType.container.shortName
      package = elem.container.container
      
      # Correct for derived target class NBNB New Rasmus 11 Apr 2011
      if valueType.isDerived:
        result['isDerived'] = True
        result['proc'] = 'delay'
 
      # processing tag
      hierarchy = elem.hierarchy
      if hierarchy == ImpConstants.parent_hierarchy:
        # parent link
        result = None
 
      elif hierarchy == ImpConstants.child_hierarchy:
        # intra- or inter-package childlink
        result['type'] = 'child'
        result['eType'] = 'cplx'
        result['content'] = globalMap[vtPrefix]['abstractTypes']
        if (not inClass.isAbstract and inClass.parentRole is None
            and elem.name == ImpConstants.packageLocatorRole):
          # MemopsRoot.packageLocator
          result['proc'] = 'loadDelayed'
      
      else:
        # crosslink
      
        if valueType.container is package:
          # intra-package crosslink
          result['type'] = 'link'
 
        elif (valueType.container in package.accessedPackages
              and elem.container.qualifiedName() != dataObjClassName and
              elem.container.container.qualifiedName() != implPackageName):
          # wrong way exolink
          # The exceptions are DataObject.access and MemopsRoot.currentXyz links
          result = None
 
        elif valueType.qualifiedName() == dataObjClassName:
          # AccessObject.dataObjects
          result = None
 
        else:
          # normal exolink
          parent = valueType.parentRole
          if (parent and parent.valueType.container is not valueType.container):
            # is TopObject
            result['type'] = 'exotop'
          else:
            result['type'] = 'exolink'
          result['eType'] = 'cplx'
          result['content'] = globalMap[vtPrefix]['exolinks']
        
        if result is not None:
          # add copyOverride
          copyOverride = True
          otherRole = elem.otherRole
          if otherRole:
            if ( otherRole.hicard == 1
                 or otherRole.hicard == otherRole.locard
                 or otherRole.changeability != ImpConstants.changeable
            ):
              copyOverride = False
          #
          result['copyOverride'] = copyOverride
 
  else:
    # map probably exists - get it
    result = globalMap['mapsByGuid'].get(elem.guid)
        
  #
  return result

###########################################################################

