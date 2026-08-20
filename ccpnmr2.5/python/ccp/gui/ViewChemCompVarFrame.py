
"""
======================COPYRIGHT/LICENSE START==========================

ViewChemCompVarFrame.py: <write function here>

Copyright (C) 2005 Tim Stevens (University of Cambridge)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in license/LGPL.license
 
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

from ccp.gui.ViewStructureFrame import ViewStructureFrame, symbolColor, symbolMultiplier, symbolSize
from ccp.general.Io import getChemCompCoord

from memops.gui.MessageReporter import showWarning

try:
  from ccp.c import StructAtom
  from ccp.c import StructBond
  from ccp.c import StructStructure
except Exception, e:
  print 'Error, the ViewChemCompVar module will not work, something is wrong with the C code.'
  print 'Exception:', e
  print 'Will continue without viewing functionality'
  StructAtom = StructBond = StructStructure = TkHandler = None
  
PROTEIN_BACKBONE = ('N','C','CA')
CARBOHYDRATE_BACKBONE = ("C1","C2","C3","C4","C5","C6","C7","C8")
NUCLEIC_ACID_BACKBONE = ("OP1","P","O3'","O5'","C3'","C4'","C5'")

SOURCE_NAMES = ('ccpnImport', 'pdb','ideal','euroCarbDb')

class ViewChemCompVarFrame(ViewStructureFrame):

  def __init__(self, parent, chemCompVar=None, **kw):

    self.showAssign  = False
    self.chemCompVar = chemCompVar

    kw['zoom'] = 2.0
    kw['bondWidth'] = 4.0
    kw['radiiScale'] = 0.4

    ViewStructureFrame.__init__(self, parent,  **kw)
  
  def configMenu(self):
 
    items = [{ 'kind': 'command', 'label': 'Reset View', 'command' : self.resetView,
               'tipText': 'Reset size and orientation to original values'},
             { 'kind': 'command', 'label': 'Spin', 'command' : self.spin,
               'tipText': 'Slowly rotate the structure' },
             { 'kind': 'command', 'label': 'Print to File', 'command' : self.printCanvas,
               'tipText': 'Save the graphical display as a PostScript file' }]
      
    self.menu.setMenuItems(items)
      
  def displayStructure(self, model=0):
    
    self.cStructure = None
    self.cStructureRotate = None
    self.cStructureZoom = None
    self.clearConnections()
    self.clearHighlights()
    self.backbone = []
    getAtomDisplayScheme = self.getAtomDisplayScheme
    newAtom = StructAtom.StructAtom
       
    if self.chemCompVar:
    
    
      cAtomDict = {}
      cBondDict = {}
      self.cStructure = StructStructure.StructStructure()
      self.cStructureRotate = self.cStructure.rotate
      self.cStructureZoom = self.cStructure.zoom
      cStructure = self.cStructure
      allAtoms = []
      allAtomsAppend = allAtoms.append
    
      project = self.chemCompVar.root
      molType = self.chemCompVar.chemComp.molType
      ccpCode = self.chemCompVar.chemComp.ccpCode

      chemCompVarCoord = None
      for sourceName in SOURCE_NAMES:
        chemCompCoord = project.findFirstChemCompCoord(sourceName=sourceName,
                                                       molType=molType,
                                                       ccpCode=ccpCode)
        
        if chemCompCoord:
          break
      
      if not chemCompCoord:
        chemCompCoord = project.findFirstChemCompCoord(molType=molType, ccpCode=ccpCode)

      if chemCompCoord:
        linking          = self.chemCompVar.linking
        descriptor       = self.chemCompVar.descriptor
        chemCompVarCoord = chemCompCoord.findFirstChemCompVarCoord(linking=linking,
                                                                   descriptor=descriptor)
        if not chemCompVarCoord:
          chemCompVarCoord = chemCompCoord # just a default
   
      if not chemCompVarCoord:
        msg = 'Chemical component %s available but no atom coordinates, so structure cannot be drawn' % self.chemCompVar.chemComp.ccpCode
        showWarning('Warning', msg, parent=self)
        return
 
      n = 0
      bonds = {}
      for atom1 in self.chemCompVar.chemAtoms:
      
        # no isolated atoms unless a single atom
        # Now allow this to see what is really in the chemCompCoord
        #if len(self.chemCompVar.chemAtoms) > 1:
        #  if not atom1.chemBonds:
        #    continue
      
        if atom1.className == 'LinkAtom':
          continue
          
        chemAtomCoord = chemCompVarCoord.findFirstChemAtomCoord(name=atom1.name)
        if not chemAtomCoord:
          continue
        
        coords = (chemAtomCoord.x,chemAtomCoord.y,chemAtomCoord.z)
        if None in coords: 
          continue
        
        symbol = atom1.elementSymbol
        bonds1 = atom1.chemBonds
        color, label, size = getAtomDisplayScheme(atom1)
        cAtom  = newAtom(size, '', label, coords, color)
        
        self.atomColors[cAtom] = color
        cAtomDict[atom1] = cAtom
        cAtomDict[cAtom] = atom1
        cBondDict[atom1] = []
        cStructure.addAtom(cAtom)
        allAtomsAppend(cAtom)
        n += 1
 
        for bond in bonds1:
          key = list(bond.chemAtoms)
          
          if cAtomDict.get(key[0]) is None:
            continue
          elif cAtomDict.get(key[1]) is None:
            continue
            
          key.sort()
          bonds[tuple(key)] = True
 
      self.cAtomDict = cAtomDict
      self.cBondDict = cBondDict
      
      molType = self.chemCompVar.chemComp.molType
      
      # StructBonds
      for atom1,atom2 in bonds.keys():
        cBond = StructBond.StructBond(cAtomDict[atom1],cAtomDict[atom2])
        cBond.setLineWidth(self.bondWidth)
        cBondDict[atom1].append(cBond)
        cBondDict[atom2].append(cBond)
        cBondDict[cBond] = (atom1, atom2)
        cStructure.addBond(cBond)
        
        if molType == 'protein':
          if atom1.name in PROTEIN_BACKBONE:
            if atom2.name in PROTEIN_BACKBONE:
              self.backbone.append( (atom1,atom2,cBond) )

        elif molType == 'carbohydrate':
          if atom1.name in CARBOHYDRATE_BACKBONE:
            if atom2.name in CARBOHYDRATE_BACKBONE:
              self.backbone.append( (atom1,atom2,cBond) )

        elif molType in ('RNA','DNA'):
          if atom1.name in NUCLEIC_ACID_BACKBONE:
            if atom2.name in NUCLEIC_ACID_BACKBONE:
              self.backbone.append( (atom1,atom2,cBond) )
 
      if n > 1:
        cStructure.moveToCenter()
        self.drawStructure()
    
      self.allAtoms = allAtoms
    
    else:
      self.drawStructure()
  
  def getAtomDisplayScheme(self, atom):
      
    symbol = atom.elementSymbol
    color  = symbolColor.get(symbol, (0.6, 0.6, 0.6))
    label  = atom.name
    size   = self.radiiScale * symbolMultiplier.get(symbol, 1.0)
    
    return color, label, size

  def update(self, chemCompVar):
  
    self.chemCompVar = chemCompVar or None
    self.displayStructure()
 
  def highlightAtom(self, atom, color=(0.0,1.0,0.0), atomSize=0.2 ):

    cAtomDict = self.cAtomDict
    cBondDict = self.cBondDict
    cAtom = cAtomDict.get(atom)
    if cAtom:
      # cAtom.setAnnotation(text)
      cAtom.setAnnotationColor(color)
      self.selectedAtoms.add(cAtom)


  def md(self, atoms, maxCycles, repScale, bondScale):
     
    R2 = self.bondLength*self.bondLength
    cycleEnergy = 0
    for c in range(maxCycles):
      prevCycleEnergy = cycleEnergy
      cycleEnergy = 0
      meanPosition = [0,0,0]
      totalMass    = 0
      for chemAtom in atoms:
        bonded = 0
        for chemAtom2 in atoms:
          if chemAtom2 is chemAtom:
            continue
          if chemAtom.bondedDict.get(chemAtom2):
            bonded = 1
        if not bonded:
          continue
        
        (x,y,z) = chemAtom.drawCoords
        
        minEnergy = None
        minIJ = (1,1,1)
        for i in range(3):
          x1 = x + i - 1
          for j in range(3):
            y1 = y + j - 1
            for k in range(3):
              z1 = z + k - 1
              energy = 0.0
              bonded = 0
              for chemAtom2 in atoms:
                if chemAtom2 is chemAtom:
                  continue

                (x2,y2,z2) = chemAtom2.drawCoords
                if z1 == z2:
                  z2 += 1
                dx = x1-x2
                dy = y1-y2
                dz = z1-z2
                r2 = (dx*dx) + (dy*dy) + (dz*dz)
                if chemAtom.bondedDict.get(chemAtom2):
                  dlength = abs(R2-r2)*bondScale
                  energy += dlength
                  bonded = 1

                else:
                  energy += repScale/float(r2)

              if minEnergy is None:
                minEnergy = energy
                minIJ = (i,j,k)
              
              elif energy <= minEnergy:
                minIJ = (i,j,k)
                minEnergy = energy
                      
        x = x - 1 + minIJ[0]
        y = y - 1 + minIJ[1]
        z = z - 1 + minIJ[2]
        if minEnergy > 0:
          chemAtom.drawCoords = [x,y,z]

        cycleEnergy += minEnergy
        
        massContrib = [chemAtom.mass*coord for coord in chemAtom.drawCoords]
        for i in range(3):
          meanPosition[i] += massContrib[i]
        
        totalMass += chemAtom.mass
      
      for i in range(3):
        meanPosition[i] /= max(1, totalMass)

      for chemAtom in atoms:
        coords = chemAtom.drawCoords
        coords = Geometry.vectorsSubtract(meanPosition, coords)
        coords = Geometry.vectorsAdd(self.center, coords)
        #chemAtom.drawCoords = coords
              
      if c % 3 == 0:
        pass
        #self.drawAtoms(atoms)
        
      # quit early if nothing happens
      if cycleEnergy == prevCycleEnergy:
        break
  
