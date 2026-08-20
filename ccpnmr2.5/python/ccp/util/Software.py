"""
======================COPYRIGHT/LICENSE START==========================

Software.py: Utility functions

Copyright (C) 2009 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

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

def getMethodStore(project):
  """Descrn: Get project.currentMethodStore, creating or setting it if need be.
     Inputs: memops.api.Implementation.MemopsRoot
     Output: ccp.api.general.Method.MethodStore
  """

  if not project.currentMethodStore:
    methodStores = project.sortedMethodStores()
    if methodStores:
      project.currentMethodStore = methodStores[0]
    else:
      # this automatically sets project.currentMethodStore
      project.newMethodStore(name=project.name)

  return project.currentMethodStore

def getSoftware(project, name, version, vendorName='', vendorAddress='', vendorWebAddress='', details=''):
  """Descrn: Get software if it exists, otherwise create a new software.
             Two softwares the same if name and version are the same.
     Inputs: memops.api.Implementation.MemopsRoot, String, String [, String, String, String]
     Output: ccp.api.general.Method.Software
  """

  methodStore = getMethodStore(project)
  software = methodStore.findFirstSoftware(name=name, version=version)
  if not software:
    software = methodStore.newSoftware(name=name, version=version)

  if vendorName:
    software.vendorName = vendorName
  
  if vendorAddress:
    software.vendorAddress = vendorAddress
  
  if vendorWebAddress:
    software.vendorWebAddress = vendorWebAddress

  if details:
    software.details = details

  return software

def getMethod(software, task, procedure=None, parameters=None, details=None):
  """Descrn: Get method if it exists, otherwise create a new method.
             Two methods the same if software, procedure and parameters are the same.
     Inputs: ccp.api.general.Method.Software, String [, String, List((String, Object)), String]
     Output: ccp.api.general.Method.Method
  """

  if parameters is None:
    parameters = []
  else:
    parameters = list(parameters)

  methodStore = getMethodStore(software.root)
  methods = methodStore.findAllMethods(software=software, procedure=procedure)
  for method in methods:
    params = [(p.name,p.value) for p in method.parameters]
    if params == parameters:
      return method

  method = methodStore.newMethod(software=software, procedure=procedure,
                         details=details)
  for (name, value) in parameters:
    method.newParameter(name=name, value=value)

  return method

