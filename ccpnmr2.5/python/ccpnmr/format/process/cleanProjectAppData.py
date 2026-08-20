from ccpnmr.format.general.Constants import allFormatsDict

def findAppData(object,appDataObjectDict):
  
  appData = object.findAllApplicationData()
  
  for appDatum in appData:
    if not appDataObjectDict.has_key(appDatum.application):
      appDataObjectDict[appDatum.application] = []
    appDataObjectDict[appDatum.application].append((object,appDatum))

def cleanAppData(project):
  
  appDataObjectDict = {}
  
  print "Collecting application data... this might take a while."
  
  #findAppData(project,appDataObjectDict)
  
  for nmrExp in project.currentNmrProject.experiments:
    findAppData(nmrExp,appDataObjectDict)
    for dataSource in nmrExp.dataSources:
      findAppData(dataSource,appDataObjectDict)
      for peakList in dataSource.peakLists:
        findAppData(peakList,appDataObjectDict)
        for peak in peakList.peaks:
          findAppData(peak,appDataObjectDict)
        
    for expDim in nmrExp.expDims:
      # Not this one
      for expDimRef in expDim.expDimRefs:
        findAppData(expDimRef,appDataObjectDict)
  
  for resonance in project.currentNmrProject.resonances:
    findAppData(resonance,appDataObjectDict)
  
  for resGroup in project.currentNmrProject.resonanceGroups:
    findAppData(resGroup,appDataObjectDict)
  
  for molSys in project.molSystems:
    findAppData(molSys,appDataObjectDict)
    for strucEnsemble in molSys.structureEnsembles:
      findAppData(strucEnsemble,appDataObjectDict)
      for model in strucEnsemble.models:
        findAppData(model,appDataObjectDict)
      for cChain in strucEnsemble.coordChains:
        findAppData(cChain,appDataObjectDict)
        for cResidue in cChain.residues:
          findAppData(cResidue,appDataObjectDict)
          for cAtom in cResidue.atoms:
            findAppData(cAtom,appDataObjectDict)
            
          
    for chain in molSys.chains:
      findAppData(chain,appDataObjectDict)
      for residue in chain.residues:
        findAppData(residue,appDataObjectDict)
        for atom in residue.atoms:
          findAppData(atom,appDataObjectDict)
  
  print "Finished collecting application data."
  
  for format in ['ccpNmr'] + allFormatsDict.keys():
    if appDataObjectDict.has_key(format):
      for (object,appData) in appDataObjectDict[format]:
        object.removeApplicationData(appData)
      print "  Deleted application data for %s..." % format

if __name__ == '__main__':

  inFile = '/ebi/msd/nmrqual/reference/ccpNmr/magnusLund/to_wim/d28k_test.xml'

  # Read xml file
  from memops.format.xml import XmlIO
  project = XmlIO.loadProjectFile(inFile)
  
  cleanAppData(project)
  
  project.saveModified()
  
