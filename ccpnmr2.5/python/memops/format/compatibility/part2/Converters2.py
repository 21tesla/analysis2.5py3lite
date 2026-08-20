""" Code split off from compatibility.converters to clear namespace for
duplicate module problem

NBNB this was crucial to avoid e.g. memops.api to end up duplicated under
both 'memops.api' and 'memops.format.compatibility.memops.api'
in sys.modules. 

SPlitting the code away from the code in Converters and removing top level 
imports in this file proved necessary to make this work.

The effect is some kind of interaction between the import and our
trick for razing sys.modules.
"""

# needed  to make it work in Python 2.5
from __future__ import absolute_import

def getOldData(oldPath, oldTag, currentVersionStr, oldName=None):
  import os
  from memops.universal import Io as uniIo
  from memops.metamodel.Util import topologicalSortSubgraph

  # find old implementation top directory and change temporarily to it.
  currentTopDir = uniIo.getTopDirectory()
  ss1, lastdir = os.path.split(currentTopDir)
  ss2, newTag =  os.path.split(ss1)
  prevTopDir = os.path.join(ss2, oldTag, lastdir)
  
  uniIo.setTopDirectory(prevTopDir)
  
  # NB these things are put inside the function to avoid problems with
  # the uniIo.setTopDirectory approach
  
  def getCurrentVersionStr():
    """ Encapsulated to get import out of main. Hoping to fix import problem.
    """
    from memops.general.Constants import currentModelVersion
    return str(currentModelVersion)
  
    
  def downlinkTagsByImport(root, sorter=topologicalSortSubgraph):
    """ downlink tags from MemopsRoot object for leafPackages in import order
    This is the post-v2.0 version
    """
    
    leafPackages = []
    packages = [root.metaclass.container.topPackage()]
    for pp in packages:
      childPackages = pp.containedPackages
      if childPackages:
        packages.extend(childPackages)
      else:
        leafPackages.append(pp)
  
    # sort leafPackages by import (imported before importing)
    leafPackages = sorter(leafPackages, 'accessedPackages')
    tags = []
    for pp in leafPackages:
      cc = pp.topObjectClass
      if cc is not None:
        pr = cc.parentRole
        if pr is not None:
          tags.append(pr.otherRole.name)
    #
    return tags 
  
  try:
    # get old version as string
    oldVersionStr = getCurrentVersionStr()
    from memops.general import Io as oldIo
 
    # post-version-2 behaviour - get repository dir and project name
    if oldName is None:
      head,tail = os.path.split(oldPath)
      oldName = os.path.splitext(tail)[0]
      oldPath = os.path.dirname(os.path.dirname(head))
 
    # get root Object and force load all contents
    oldRoot = oldIo.loadProject(oldPath, projectName=oldName)
 
    if oldVersionStr != currentVersionStr:
      # We do have to convert - force load all TopObjects
      tags = downlinkTagsByImport(oldRoot)
      for tag in tags:
        for topObj in getattr(oldRoot, tag):
          if not topObj.isLoaded:
            topObj.load()
  finally:
    # set path back to current
    uniIo.setTopDirectory(currentTopDir)
  
  return (oldRoot, oldVersionStr, tags)
