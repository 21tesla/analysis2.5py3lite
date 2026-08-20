def listChemCompInfo(chemComp):

  linkings = []
  chemAtomsByLinking = {}
  chemBondsByLinking = {}
         
  for ccv in chemComp.sortedChemCompVars():
    chemAtomNames =  [ca.name for ca in ccv.sortedChemAtoms()]
    if ccv.linking not in linkings:
      linkings.append(ccv.linking)
      chemAtomsByLinking[ccv.linking] = chemAtomNames
      chemBondsByLinking[ccv.linking] = list(ccv.sortedChemBonds())
    else:
      caIndex = 0
      while (caIndex < len(chemAtomsByLinking[ccv.linking])):
        caName = chemAtomsByLinking[ccv.linking][caIndex]
        if caName not in chemAtomNames:
          chemAtomsByLinking[ccv.linking].pop(caIndex)
          for chemBond in chemBondsByLinking[ccv.linking]:
            if chemBond.chemAtoms[0].name == caName or  chemBond.chemAtoms[1].name == caName:
              chemBondsByLinking[ccv.linking].pop(chemBondsByLinking[ccv.linking].index(chemBond))
        else:
          caIndex += 1
    
    print chemComp.molType, chemComp.ccpCode

    for linking in linkings:
      print
      print "Linking:",linking
      chemAtomsByLinking[linking].sort()
      print
      print "ChemAtoms:"
      print chemAtomsByLinking[linking]
      print
      print "ChemBonds:"
      for cb in chemBondsByLinking[linking]:
        chemAtoms = cb.sortedChemAtoms()
        print "  (\"%s\",\"%s\"): '%s'," % (chemAtoms[0].name,chemAtoms[1].name,cb.bondType)


if __name__ == '__main__':

  import os, shutil
  
  from ccpnmr.format.general.Conversion import FormatConversion
  
  #
  # Set variables for this project
  #
  
  mol2File = "data/gnp.mol2"
  projectName = 'testImport'
  
  #
  # Remove project if it already exists - otherwise will interfere with creation of data
  # as it will be picked up automatically.
  #
  
  if os.path.exists(projectName):
    shutil.rmtree(projectName)
    
  #
  # Create the FormatConversion object (is a wrapper around the FormatConverter that takes care
  # of CCPN object creation and tracking.
  #

  fc = FormatConversion(identifier=projectName, useGui=True)
  
  #
  # Import the chemComp - returns a list of objects
  #
  
  chemComps = fc.importFile('chemComps','mol2',mol2File)
  
  #
  # If the import worked, list info on the chemComps that were created
  #
  
  for chemComp in chemComps:
    listChemCompInfo(chemComp)

  #
  # Save the project. The XML files that you'll need are in:
  #
  #  - ChemComp info in ccp/molecule/ChemComp/
  #  - ChemCompCoord info in ccp/molecule/ChemCompCoord/
  #
  # You can either copy these into the project where you need them (they'll be picked up
  # automatically), or you can put them in the central data directory of your current
  # CCPN installation:
  #
  #  - ccpnmrXXX/data/ccp/molecule/ChemComp(Coord)
  #
  # in this case they'll be picked up by all projects using this installation.
  #
  
  fc.saveCcpnProject()
  
