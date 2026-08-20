"""
======================COPYRIGHT/LICENSE START==========================

General.py: Data compatibility handling

Copyright (C) 2007 Rasmus Fogh (CCPN project)
 
=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license.
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

======================COPYRIGHT/LICENSE END============================

To obtain more information about this code:

- CCPN website (http://www.ccpn.ac.uk)

- contact Rasmus Fogh (ccpn@bioc.cam.ac.uk)

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following reference:

===========================REFERENCE START=============================
Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and 
automated software development. Bioinformatics 21, 1678-1684.
===========================REFERENCE END===============================
 
"""

from memops.format.compatibility.upgrade.v_2_0_a2 import General as General20a2


# functions to convert fullKeys before application
fullKeyConverters = {}
fullKeyConverters.update(General20a2.fullKeyConverters)

# guids of elements that should be treated as old
elemsTreatedAsOld = set((
 'www.ccpn.ac.uk_Fogh_2008-05-05-17:24:14_00002',
 'www.ccpn.ac.uk_Fogh_2008-05-05-18:37:55_00020',
 'www.ccpn.ac.uk_Fogh_2008-05-05-18:37:55_00022',
 'www.ccpn.ac.uk_Fogh_2008-05-05-17:24:14_00008',
 'www.ccpn.ac.uk_Fogh_2006-09-14-18:48:26_00002',
 'www.ccpn.ac.uk_Fogh_2006-08-16-18:20:11_00004')
).union(General20a2.elemsTreatedAsOld)

# pairs of element guids that should be treated as matching, e.g. when
# a single element must match with several elements in subclasses

# these are links from MemopsRoot to individual TopObjects that all have changed their guids in a model change.
elementPairings = [
 ('www.ccpn.ac.uk_Fogh_2006-09-04-17:21:42_00003', 'ccpn_automatic_memops.Implementation.MemopsRoot.accessControlStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:10:51_00005', 'ccpn_automatic_memops.Implementation.MemopsRoot.affiliationStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-14:16:25_00009', 'ccpn_automatic_memops.Implementation.MemopsRoot.analysisProject'),
 ('www.ccpn.ac.uk_Fogh_2007-06-05-13:42:07_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.annealProtocolStore'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:55:49_00003', 'ccpn_automatic_memops.Implementation.MemopsRoot.annotationStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:10:50_00005', 'ccpn_automatic_memops.Implementation.MemopsRoot.chemCompCharge'),
 ('www.ccpn.ac.uk_Fogh_2006-09-14-18:48:28_00026', 'ccpn_automatic_memops.Implementation.MemopsRoot.chemCompCoord'),
 ('www.ccpn.ac.uk_Fogh_2006-09-14-18:48:28_00024', 'ccpn_automatic_memops.Implementation.MemopsRoot.chemComp'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:10:50_00009', 'ccpn_automatic_memops.Implementation.MemopsRoot.chemElementStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:12:00_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.citationStore'),
 ('www.ccpn.ac.uk_Fogh_2006-09-04-17:21:42_00001', 'ccpn_automatic_memops.Implementation.MemopsRoot.classification'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:39:15_00001', 'ccpn_automatic_memops.Implementation.MemopsRoot.crystallizationStore'),
 ('www.ccpn.ac.uk_Fogh_2007-10-01-14:14:00_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.dangleStore'),
 ('www.ccpn.ac.uk_Fogh_2007-04-03-19:10:06_00001', 'ccpn_automatic_memops.Implementation.MemopsRoot.dataLocationStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:17:05_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.database'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:39:16_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.expBlueprintStore'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:39:16_00004', 'ccpn_automatic_memops.Implementation.MemopsRoot.experimentStore'),
 ('www.ccpn.ac.uk_Fogh_2006-12-31-09:03:01_00008', 'ccpn_automatic_memops.Implementation.MemopsRoot.haddockProject'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:39:16_00006', 'ccpn_automatic_memops.Implementation.MemopsRoot.holderStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:12:01_00001', 'ccpn_automatic_memops.Implementation.MemopsRoot.instrumentStore'),
 ('www.ccpn.ac.uk_Fogh_2007-11-14-17:02:56_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.keywordDefinitionStore'),
 ('www.ccpn.ac.uk_Fogh_2006-12-31-09:03:01_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.labeledMolecule'),
 ('www.ccpn.ac.uk_Fogh_2007-01-24-12:23:19_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.labelingScheme'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:39:17_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.locationStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:12:01_00003', 'ccpn_automatic_memops.Implementation.MemopsRoot.methodStore'),
 ('www.ccpn.ac.uk_Fogh_2008-02-20-18:17:58_00003', 'ccpn_automatic_memops.Implementation.MemopsRoot.molSystemSymmetrySet'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-14:16:24_00026', 'ccpn_automatic_memops.Implementation.MemopsRoot.molSystem'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-14:16:24_00024', 'ccpn_automatic_memops.Implementation.MemopsRoot.molecule'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-14:16:25_00007', 'ccpn_automatic_memops.Implementation.MemopsRoot.nameMappingStore'),
 ('www.ccpn.ac.uk_Fogh_2006-09-14-18:48:28_00011', 'ccpn_automatic_memops.Implementation.MemopsRoot.nmrConstraintStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:10:50_00013', 'ccpn_automatic_memops.Implementation.MemopsRoot.nmrEntryStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:10:51_00001', 'ccpn_automatic_memops.Implementation.MemopsRoot.nmrExpPrototype'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-14:23:33_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.nmrProject'),
 ('www.ccpn.ac.uk_Fogh_2006-09-14-18:48:28_00014', 'ccpn_automatic_memops.Implementation.MemopsRoot.nmrReferenceStore'),
 #('www.ccpn.ac.uk_Fogh_2007-06-05-13:42:07_00005', 'ccpn_automatic_memops.Implementation.MemopsRoot.nmrSimStore'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:39:17_00004', 'ccpn_automatic_memops.Implementation.MemopsRoot.protocolStore'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:39:17_00006', 'ccpn_automatic_memops.Implementation.MemopsRoot.refSampleComponentStore'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:39:17_00008', 'ccpn_automatic_memops.Implementation.MemopsRoot.sampleStore'),
 ('www.ccpn.ac.uk_Fogh_2006-08-17-16:10:50_00011', 'ccpn_automatic_memops.Implementation.MemopsRoot.stereochemistryStore'),
 ('www.ccpn.ac.uk_Fogh_2006-09-14-18:48:28_00021', 'ccpn_automatic_memops.Implementation.MemopsRoot.structureEnsemble'),
 ('www.ccpn.ac.uk_Fogh_2007-11-13-15:55:51_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.validationStore'),
 ('www.ccpn.ac.uk_Fogh_2006-09-05-12:39:18_00002', 'ccpn_automatic_memops.Implementation.MemopsRoot.targetStore'),
] + General20a2.elementPairings

def extraMapChanges(globalMapping):
  """ Extra map changes specific for a given step
  """
  
  # change string attributes to cplx
  for dd in globalMapping['mapsByGuid'].values():
    data = dd.get('data')
    if data:
      if data['toStr'] == 'text':
        dd['eType'] = 'cplx'
  
  # Set links Nmr -> nmrExpPrototype to delayed
  for guid in ('www.ccpn.ac.uk_Fogh_2006-08-16-18:20:06_00008',
   'www.ccpn.ac.uk_Fogh_2006-08-16-18:23:00_00002',
   'www.ccpn.ac.uk_Fogh_2006-08-16-18:20:05_00025'
  ):
    globalMapping['mapsByGuid'][guid]['proc'] = 'delay'
  
  #
  General20a2.extraMapChanges(globalMapping)
