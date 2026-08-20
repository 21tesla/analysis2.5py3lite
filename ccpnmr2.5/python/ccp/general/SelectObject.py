
"""b
======================COPYRIGHT/LICENSE START==========================

SelectObjectPopup.py: Part of the CcpNmr Analysis program

Copyright (C) 2005 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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

# NBNB comments: rhfogh 15/1/09
# Differences to pickObjectByAttribute:
#   labels do not default to attr 'name'



import Tkinter
from memops.general import Implementation
ApiError = Implementation.ApiError
from memops.gui.Label import Label
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.BasePopup import BasePopup
from memops.metamodel import ImpConstants

class SelectObjectPopup(BasePopup):

  def __init__(self, parent, objects, key=None, labels=None, objectName=None, 
               refClassName=None, **kw):
    
    self.parent = parent
    
    self.resetObjectType(objects, key, labels, objectName, refClassName)
    
    self.doRegisters()
    
    BasePopup.__init__(self, parent, title= "Object Selection",
                       borderwidth=5,  **kw)
  
  def resetObjectType(self, objects, key=None, labels=None, objectName=None, 
                      refClassName=None, **kw):
    if refClassName:
      apiClass = Implementation.getClassFromFullName(refClassName)
      refClass = apiClass._metaclass
    else:
      apiClass = objects[0].__class__
      refClass = objects[0].metaclass
      refClassName = refClass.qualifiedName()
      
    self.refClassName = refClassName
    self.apiClass = apiClass

    for obj in objects:
      if not isinstance(obj, self.apiClass):
        raise ApiError("Object %s is not instance of %s" % (obj, refClassName))
    
    self.objects = objects
    
    self.key    = key
    if labels is None:
      labels = {}
    self.labels = labels.copy()
    self.objectName = objectName or refClass.name 
        
    self.headingList = self.getAttributeHeadings(refClass)
    
    # add 'fullKey to labels
    self.labels['fullKey'] = [str(getattr(obj, 'getFullKey')()) 
                              for obj in self.objects]
    
    if refClass.isAbstract:
      # add extra column to define the type
      self.labels['Class'] = [obj.__class__.__name__ for obj in self.objects]
      self.headingList.append('Class')
    
  
  def body(self, guiFrame):
    
    x = self.parent.winfo_rootx() + 50
    y = self.parent.winfo_rooty() + 50
    location = '500x400+%d+%d' % (x, y)
    self.geometry(location)
    
    guiFrame.expandGrid(1,0)
    
    textMatrix = self.getTextMatrix(self.objects)
        
    self.label = Label(guiFrame, text='Select a %s ' % self.objectName, grid=(0,0))
    
    self.scrolledMatrix = ScrolledMatrix(guiFrame, initialRows=10,
                                         headingList=self.headingList,
                                         callback=self.selectCell,
                                         objectList=self.objects,
                                         textMatrix=textMatrix,
                                         grid=(1,0))

    
    self.protocol('WM_DELETE_WINDOW', self.close)
  
  def getAttributeHeadings(self, refClass):
  
    attrNames, keyNames = getDisplayAttrNames(refClass, self.labels)
  
    # set up for key parameter and fullkey
    # fullKey will be used as key if none othe is present,
    # and will be put at the end otherwise
    if self.key is None:
      if keyNames:
        attrNames.append('fullKey')
      else:
        # key in class not available, no key passed in. Use fullKekey = 'fullKey'
        attrNames.insert(0,'fullKey')
    else:
      try:
        attrNames.remove(self.key)
      except ValueError:
        raise ApiError("key %s not in attributes or labels: %s"
                        % (self.key, attrNames))
      attrNames.insert(0,self.key)
      attrNames.append('fullKey')
     
    return attrNames
      
  def getTextMatrix(self,objects):
    
    headings = self.headingList
    nHeadings = len(headings)
    matrix = []
    
    for ii,object in enumerate(objects):
      line = nHeadings * ['None']
      for jj,attrName in enumerate(headings):
      
        if self.labels.get(attrName) is None:
          
          if hasattr(object,attrName):
            val = getattr(object,attrName)
            if isinstance(val, frozenset):
              val = list(val)
              val.sort()
            ss = str(val)
            if len(ss) > 32:
              ss = ss[:32] + '...'
            line[jj] = ss
        
        else:
          line[jj] = self.labels[attrName][ii]
                        
      matrix.append(line)
    
    return matrix
    
  def update(self, objects, key=None, labels=None, objectName=None, 
             refClassName=None, remove=None, reset=True):
    
    if reset:
      self.doUnregisters()
      self.resetObjectType(objects, key, labels, objectName, refClassName)
      self.doRegisters()
    
    if remove:
      ii = self.objects.index(remove)
      del self.objects[ii]
      for ll in self.labels.values():
        del ll[ii]
    
    if not self.objects:
      self.selectCell(None)
    
    if len(objects) == 1:
      return objects[0]

    self.label.set('Select %s ' % self.objectName)
    textMatrix = self.getTextMatrix(self.objects)
    self.scrolledMatrix.update(objectList=self.objects, 
                               headingList=self.headingList, 
                               textMatrix=textMatrix)

  def close(self):
  
    self.doUnregisters()
    self.parent.argumentVar.set(1)
 
    BasePopup.close(self)

  def destroy(self):
  
    self.doUnregisters()
    try:
      self.parent.argumentVar.set(1)
    except AttributeError, e:
      pass
      
    BasePopup.destroy(self)
  
  def updateAfter(self, object):
  
    # only for deleted object at the mo
    self.after_idle(lambda s = self:s.update(s.objects, remove=object,
                                             reset=False))
    
  def selectCell(self, object, row=None, col=None):
  
    self.doUnregisters()
    #self.destroy()
    self.object = object
    self.parent.argumentVar.set(1)
    self.parent.argumentObject = object
    self.close()
    return self.object
    
  def doRegisters(self):

    try:
      Implementation.registerNotify(self.updateAfter, self.refClassName,
                                    'delete')
    except AttributeError:
      # refClassName may be an abstract class, 
      # in which case registering will fail
      return
  
  def doUnregisters(self):

    try:
      Implementation.unregisterNotify(self.updateAfter, self.refClassName,
                                    'delete')
    except AttributeError:
      # refClassName may be an abstract class, 
      # in which case registering will fail
      return

def getDisplayAttrNames(refClass, labels=None, minimal=False):
  """ get list of attributes to display for object selection.
  Heuristic, designed to provide informative display.
  returns all attributs, and, separately, key attributes if all keyNames are
  attributes."""
  
  ignoreClasses = ('DataObject', 'MemopsObject', 'ComplexDataType')

  if labels is None:
    labels = {}
  
  attributes = refClass.getAllAttributes()
  
  # filter out undesirables
  refPackage = refClass.container
  ll = []
  for attr in attributes:
    if attr.isDerived:
      continue
    elif attr.container.name in ignoreClasses:
      continue
    elif attr.container.container is not refPackage:
      continue
    else:
      ll.append(attr)
  attributes = ll
  
  keyNames = refClass.keyNames
  for tag in keyNames:
    attr = refClass.getElement(tag)
    if attr not in attributes:
      # we have a role in the key - key can not be used
      keyNames = []
      break
  
  # Now divide attributes into classes
  # for minimal annotations we only want the frozen ones
  lldict = {}
  lfrozen = lldict['frozen'] = []
  lsingle = lldict['single'] = []
  lmultiple = lldict['multiple'] = []
  for attr in attributes:
    if attr.hicard == attr.locard and attr.changeability == ImpConstants.frozen:
      lfrozen.append(attr.name)
    elif minimal:
      continue
    elif attr.hicard == 1:
      lsingle.append(attr.name)
    else:
      lmultiple.append(attr.name)
  lfrozen.sort()
  lsingle.sort()
  lmultiple.sort()
  if 'details' in lsingle:
    lsingle.remove('details')
    lsingle.append('details')
 
  # build up main result
  attrNames = []
  
  # first add keyNames
  for tag in keyNames:
    if tag not in attrNames:
      attrNames.append(tag)
  
  # then add external labels
  attrNames.extend(list(sorted(labels.keys())))
  
  # then name, if any
  if 'name' not in attrNames and refClass.getElement('name') is not None:
    attrNames.append('name')
  
  # then rest of attributes in order froaen, single, multiple
  for ss in ('frozen', 'single', 'multiple'):
    for tag in lldict[ss]:
      if tag not in attrNames:
        attrNames.append(tag)
  
  #
  return (attrNames, keyNames)
    
  
