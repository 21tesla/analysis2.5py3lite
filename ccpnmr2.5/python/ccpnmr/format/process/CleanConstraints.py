#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

CleanConstraints.py: Cleans up constraint lists in the data model

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
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

- contact Wim Vranken (wim@ebi.ac.uk)
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

import sys

from ccpnmr.format.general.userInteraction import setupMessageReporter

from ccp.general.Util import getOtherAtom, getResonancesFromPairwiseConstraintItem

import ccp.api.nmr.Nmr as Nmr

#
# Set up a class for this
#

class CleanConstraints:

  def __init__(self,project):
  
    self.project = project
    
  def cleanAllDistanceConstraintLists(self):
    
    for strucGen in self.project.currentNmrProject.findAllStructureGenerations():

      for constraintList in strucGen.nmrConstraintStore.findAllConstraintLists(className = 'DistanceConstraintList'):
    
        self.cleanDistanceConstraintList(constraintList)
        
  def cleanAllDihedralConstraintLists(self):
    
    for strucGen in self.project.currentNmrProject.findAllStructureGenerations():

      for constraintList in strucGen.nmrConstraintStore.findAllConstraintLists(className = 'DihedralConstraintList'):
    
        self.cleanDihedralConstraintList(constraintList)
         
  def cleanAllHBondConstraintLists(self):
    
    for strucGen in self.project.currentNmrProject.findAllStructureGenerations():

      for constraintList in strucGen.nmrConstraintStore.findAllConstraintLists(className = 'HBondConstraintList'):
    
        self.cleanDistanceConstraintList(constraintList)
        
  def cleanHBondConstraintList(self,constraintList):
  
    self.cleanDistanceConstraintList(constraintList)
    
  
  def cleanDihedralConstraintList(self,constraintList):
  
    print "Cleaning dihedral constraint list '%s'." % constraintList.name        

    constraints = []
    constraintRefs = []
    
    #
    # Set up some reference data and do initial checking
    #
    
    for constraint in constraintList.constraints:
      
      constraints.append([])
      constraintRefs.append(constraint)
      
      constraintDeleted = 0
      
      resonances = list(constraint.resonances)
          
      #
      # Check if restraint between covalently linked atoms. Assuming resonances
      # to be listed in order.
      # TODO: need handling of phi, psi, ... before can do this based on reference data
      #       can't find covalent connections in current case!!
      #
      """
      for i in range(0,3):
      
        resonance = resonances[i]
      
        if resonance.resonanceSet and len(resonance.resonanceSet.atomSets) == 1:
            
          nextResonance = resonances[i+1]

          if nextResonance.resonanceSet and len(nextResonance.resonanceSet.atomSets) == 1:
          
            nmrAtoms = resonance.resonanceSet.atomSets[0].atoms
            nextNmrAtoms = nextResonance.resonanceSet.atomSets[0].atoms
            
            if len(nmrAtoms) == 1 and len(nextNmrAtoms) == 1:
              
              nextChemAtom = nextNmrAtoms[0].chemAtom
              chemAtom = nmrAtoms[0].chemAtom
              
              atomsConnected = 0
              
              for chemBond in chemAtom.chemBonds:
              
                if getOtherAtom(chemAtom,chemBond) == nextChemAtom:
                
                  atomsConnected = 1
                  
              if not atomsConnected:
                
                print "  Deleting dihedral constraint %s (no connection via atoms %s,%s)" % (constraint.serial,chemAtom.name,nextChemAtom.name)

                constraint.delete()
                constraintDeleted = 1

      """
      #
      # TODO: sorting all resonances now... not entirely correct but should be OK
      #

      if not constraintDeleted:
      
        resonances.sort()
        constraints[-1].append(resonances)


    #
    # Now check whether there are any duplicates
    #
    
    i = 0
    
    while (i < len(constraints)):

      refConstraint = constraints[i]
      
      j = i + 1
      
      while (j < len(constraints)):
        
        if refConstraint == constraints[j]:
          items = ""
          for (res1,res2,res3,res4) in constraints[j]:
            items += "[%s,%s,%s,%s]," % (res1.name,res2.name,res3.name,res4.name)
          items = items[:-1]
          
          print "  Constraint %s (items %s) is duplicate of %s - deleted" % (constraintRefs[j].serial,items,constraintRefs[i].serial)

          #
          # Check if items are the same (ONLY checking upper, lower!), copy over if not
          #
          
          for otherDihedralItem in constraintRefs[j].items:
          
            itemExists = 0
            
            for dihedralItem in constraintRefs[i].items:
              
              if otherDihedralItem.upperLimit == dihedralItem.upperLimit and otherDihedralItem.lowerLimit == dihedralItem.lowerLimit:
              
                itemExists = 1
                break
                
            if not itemExists:
              keywds = {} 
              for attribute in ['upperLimit','targetValue','lowerLimit','error']:
                keywds[attribute] = getattr(otherDihedralItem,attribute)
              
              Nmr.DihedralConstraintItem(constraintRefs[i],**keywds)
              print "  Copied item with values %s" % (str(keywds))
          
          #
          # Delete after copying information
          #
          
          constraintRefs[j].delete()
          constraints.pop(j)
          constraintRefs.pop(j)
 
        else:
          j += 1
          
      i+= 1

  def cleanDistanceConstraintList(self,constraintList):
  
    print "Cleaning distance constraint list '%s'." % constraintList.name        

    constraints = []
    constraintRefs = []
    
    #
    # Set up some reference data and do initial checking
    #
    
    for constraint in constraintList.constraints:
      
      constraints.append([])
      constraintRefs.append(constraint)
      
      for item in constraint.sortedItems():
        
        resonances = getResonancesFromPairwiseConstraintItem(item)
        
        #
        # Check for items between the same resonance (or resonanceSet)
        # Remove if present.
        #
        
        if resonances[0] == resonances[1]:
          print "  Deleting constraint item (between same resonance)"
          item.delete()
          continue
        
        elif resonances[0].resonanceSet and resonances[0].resonanceSet == resonances[1].resonanceSet:
          print "  Deleting constraint item (between same resonanceSets)"
          # This item can be deleted because ambiguity present on resonanceSet level -
          # need not be repeated on constraint level!
          item.delete()
          continue
          
        #
        # Check if restraint between directly covalently linked atoms
        #
        
        if resonances[0].resonanceSet and len(resonances[0].resonanceSet.atomSets) == 1:
            
          if resonances[1].resonanceSet and len(resonances[1].resonanceSet.atomSets) == 1:
          
            #
            # Need direct link to atomSet for both resonances...
            #
        
            nmrAtoms0 = list(resonances[0].resonanceSet.findFirstAtomSet().atoms)
            
            refChemAtom0 = self.getRefChemAtom(nmrAtoms0)
        
            if refChemAtom0:
            
              nmrAtoms1 = list(resonances[1].resonanceSet.findFirstAtomSet().atoms)
            
              refChemAtom1 = self.getRefChemAtom(nmrAtoms1)

              if refChemAtom1 and nmrAtoms0[0].residue == nmrAtoms1[0].residue:
              
                itemDeleted = 0
                
                if refChemAtom1 == refChemAtom0:
                
                  print "  Deleting constraint item (connected via atom %s)" % (refChemAtom0.name)

                  item.delete()
                  itemDeleted = 1
              
                #else:
                
                #  for chemBond in refChemAtom0.chemBonds:
                
                #    if refChemAtom1 == getOtherAtom(refChemAtom0,chemBond):
                  
                #      print "  Deleting constraint item (between atoms %s, %s directly connected via two covalent bonds)" % (refChemAtom0.name,refChemAtom1.name)

                #      item.delete()
                #      itemDeleted = 1
                #      break
                  
                if itemDeleted:
                
                  continue
        
        #
        # Keep track of valid items
        #

        constraints[-1].append(resonances)
              
      #
      # Check if there are any items left
      #
      
      if not constraint.items:
        print "  Deleting constraint %d (no items left)" % constraint.serial
        constraint.delete()    
        constraints.pop(-1)
        constraintRefs.pop(-1)
        
      else:
      
        constraints[-1].sort()

    #
    # Now check whether there are any duplicates
    #
    
    i = 0
    
    while (i < len(constraints)):

      refConstraint = constraints[i]
      
      j = i + 1
      
      while (j < len(constraints)):
        
        if refConstraint == constraints[j]:
          items = ""
          for (res1,res2) in constraints[j]:
            items += "[%s,%s]," % (res1.name,res2.name)
          items = items[:-1]
          
          #
          # Always keep tightest value if available
          #
          
          for attribute in ['upperLimit','targetValue','lowerLimit']:
            if getattr(constraintRefs[i],attribute) and getattr(constraintRefs[j],attribute):
              setattr(constraintRefs[i],attribute,min(getattr(constraintRefs[i],attribute),getattr(constraintRefs[j],attribute)))

          print "  Constraint %s (items %s) is duplicate of %s - deleted" % (constraintRefs[j].serial,items,constraintRefs[i].serial)
          constraintRefs[j].delete()
          constraints.pop(j)
          constraintRefs.pop(j)
 
        else:
          j += 1
          
      i+= 1


  def getRefChemAtom(self,nmrAtoms):
  
    refChemAtom = nmrAtoms[0].chemAtom
          
    if len(refChemAtom.chemBonds) == 1:

      refChemAtom = getOtherAtom(refChemAtom,refChemAtom.chemBonds[0])
      
    else:
        
      refChemAtom = None
        
    return refChemAtom
  
if __name__ == "__main__":

  inFile = sys.argv[1]  # Try ../../info.reference/ccpNmr/aartUtrecht/1d3z/tlink.xml
  #outFile = sys.argv[2]
  
  # Read xml file
  from memops.format.xml import XmlIO
  proj = XmlIO.loadProjectFile(inFile)
  
  cc = CleanConstraints(proj)
  
  cc.cleanAllDistanceConstraintLists()
  cc.cleanAllDihedralConstraintLists()

  # Write out as xml!
  #PyXmlIO.exportXml(outFile,proj)

  
