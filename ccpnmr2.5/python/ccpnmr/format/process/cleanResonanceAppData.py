from ccpnmr.format.general.Constants import assign_kw, atomSerial_kw

def cleanResonanceAppData(project,formatName):

  if project.currentNmrProject:
    for resonance in project.currentNmrProject.resonances:
      for keyword in (assign_kw,atomSerial_kw,'shiftList','mappingSingleMatch'):
        appData = resonance.findFirstApplicationData(application = formatName, keyword = keyword)
        if appData:
          resonance.removeApplicationData(appData)

if __name__ == '__main__':

  from memops.general.Io import loadProject
  import sys

  projectDir = sys.argv[1]

  project = loadProject(projectDir)
  
  cleanResonanceAppData(project,'xeasy')
  
  #project.saveModified()
