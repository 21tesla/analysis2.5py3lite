
"""
======================COPYRIGHT/LICENSE START==========================

ViewIsotopomerFrame.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

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
from ccp.gui.ViewChemCompVarFrame import ViewChemCompVarFrame
from ccp.gui.ViewStructureFrame import symbolColor, symbolMultiplier, symbolSize
from ccp.util.LabeledMolecule import getIsotopomerSingleAtomFractions
from ccpnmr.analysis.core.MoleculeBasic import DEFAULT_ISOTOPES

class ViewIsotopomerFrame(ViewChemCompVarFrame):

  def __init__(self, parent, isotopomers, *args, **kw):
  
    self.chemCompVar = None
    self.atomLabelDict = {}
    
    ViewChemCompVarFrame.__init__(self, parent, self.chemCompVar, *args, **kw)

    self.setIsotopomers(isotopomers)

  def setIsotopomers(self, isotopomers=None):
    
    if isotopomers:
      chemComp = isotopomers[0].chemCompLabel.chemComp
      for isotopomer in isotopomers[1:]:
        if isotopomer.chemCompLabel.chemComp is not chemComp:
          raise Exception('Isotopomers from different ChemComps')
 
      chemCompVar = chemComp.findFirstChemCompVar(isDefaultVar=True) \
                     or chemComp.findFirstChemCompVar(linking='none') \
                     or chemComp.findFirstChemCompVar()
 
      self.chemCompVar = chemCompVar
    
    else:
      self.chemCompVar = None
        
    self.isotopomers = isotopomers
    self.displayStructure()
  
  
  def getAtomDisplayScheme(self, atom):
  
    fracDict = getIsotopomerSingleAtomFractions(self.isotopomers,
                                                atom.name, atom.subType)
    
    symbol = atom.elementSymbol
    default = DEFAULT_ISOTOPES.get(symbol)
    
    if default:
      v = fracDict.get(default, 0.0)
      r = 0.3 + (v*0.5)
      g = 0.3 + (v*0.7)
      color = (r, g, 0.3) 
      
    else:
      color = (0.3, 0.3, 0.3)
      
    label  = atom.name
    size   = self.radiiScale * symbolMultiplier.get(symbol, 1.0)
    
    return color, label, size

  
