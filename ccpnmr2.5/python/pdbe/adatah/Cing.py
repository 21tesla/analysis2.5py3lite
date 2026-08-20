import sys, os, re

from pdbe.adatah.Io import getDataFromHttp

def getCingSummaryTextInfo(filePath):

  """

  Note: not getting WHATIF and PROCHECK info, but can get at it if required!!

  """

  fin = open(filePath)
  lines = fin.readlines()
  fin.close()
  
  
  cingSummary = {'distanceRestraints': {}, 'dihedralRestraints': {}}
  
  curDict = None
  curType = None
  
  for line in lines:
    
    cols = line.split()
    
    if line.count(" DistanceRestraintList "):
   
      distanceRestraintListName = cols[2][1:-1]
      
      curType = 'distanceRestraints'
      cingSummary[curType][distanceRestraintListName] = {'violations': {}}
      curDict = cingSummary[curType][distanceRestraintListName]

    elif line.count(" DihedralRestraintList "):
   
      dihedralRestraintListName = cols[2][1:-1]
      
      curType = 'dihedralRestraints'
      cingSummary[curType][dihedralRestraintListName] = {'violations': {}}
      curDict = cingSummary[curType][dihedralRestraintListName]
      
    elif line.count(" RDCRestraintList "):
   
      curType = 'rdcRestraints'
      #cingSummary[curType][dihedralRestraintListName] = {'violations': {}}
      #curDict = cingSummary[curType][dihedralRestraintListName]

    elif line.count("CING ROG analysis "):
      curType = 'cingRog'
      curDict = cingSummary[curType] = {}
      
    elif line.count(" WHAT IF Summary "):
      curType = 'whatIf'
      curDict = cingSummary[curType] = {}

    elif line.count(" Procheck Summary "):
      curType = 'procheck'
      curDict = cingSummary[curType] = {}

    elif line.count(" Wattos Summary "):
      curType = 'wattos'
      curDict = cingSummary[curType] = {}
    
    elif curType == 'distanceRestraints' and cols:
    
      if cols[0] == 'rmsd:':
        if cols[1] == '.':
          average = None
        else:
          average = float(cols[1])
        
        if cols[3] == '.':
          dev = None
        else:
          dev = float(cols[3])
       
        curDict['rmsd'] = (average,dev)

      elif line.count("total all"):
        curDict['numRestraints'] = int(cols[2])
      elif line.count(" A: "):
        if cols[0] == '<-0.1':
          curDict['violations']['lower'] = int(cols[2])
        else:
          if cols[1] == '0.1':
            curDict['violations']['small'] = int(cols[3])
          elif cols[1] == '0.3':
            curDict['violations']['medium'] = int(cols[3])
          elif cols[1] == '0.5':
            curDict['violations']['large'] = int(cols[3])
            
    elif curType == 'dihedralRestraints' and cols:
    
      if cols[0] == 'rmsd:':
        (part1,value2) = line.split("+-")
        (crap,value1) = part1.split()
        
        curDict['rmsd'] = (float(value1.strip()),float(value2.strip()))

      elif line.count(" degree: "):
        if cols[1] == '>1':
          curDict['violations']['small'] = int(cols[3])
        elif cols[1] == '>3':
          curDict['violations']['medium'] = int(cols[3])
        elif cols[1] == '>5':
          curDict['violations']['large'] = int(cols[3])
            
    elif curType == 'cingRog' and cols:
    
      if not line.count('---'):
        if cols[0][-1] == ':':
          cols[0] = cols[0][:-1]
        curDict[cols[0]] = int(line[15:18])
        
    elif curType == 'whatIf' and cols:
    
      if line.count('+/-'):
        index = cols.index('+/-')
        chkType = ' '.join(cols[:index-2])
        curDict[chkType] = (getCingTextFloat(cols[index-1]),getCingTextFloat(cols[index+1]))
        
    elif curType == 'wattos' and cols:
    
      if cols[0] == 'Overall':
        curDict['completeness'] = float(cols[4])

  return cingSummary

def getCingTextFloat(value):

  if value == '.':
    value = None
  else:
    value = float(value)
    
  return value
  
def getCingCcpnTgzFileList(nrgCingCcpnUrl=None):
  
  print "Getting information about current NRG-CING files..."
  
  from pdbe.adatah.Io import getTextFromHttp
  
  dirNamePatt = re.compile("\<a href\=\"[a-z0-9]+\/\"\>([a-z0-9]+)\/\<\/a\>")
  ccpnProjPatt = re.compile("\<a href\=\"[a-z0-9]+\.tgz\"\>([a-z0-9]+)\.tgz\<\/a\>\s+\<\/td\>\<td align\=\"right\"\>(\d+-\w+-\d+ \d+:\d+)\s+\<\/td\>")
  
  if not nrgCingCcpnUrl:
    nrgCingCcpnUrl = "http://nmr.cmbi.ru.nl/NRG-CING/input/"
  
  #
  # Get the list of directories containing the CCPN projects
  #
  
  directoriesHttpPageLines = getTextFromHttp(nrgCingCcpnUrl)

  directoryList = []
  for line in directoriesHttpPageLines:
    dirNameSearch = dirNamePatt.search(line)
    if dirNameSearch:
      directoryList.append(dirNameSearch.group(1))
  
  #
  # Get the actual locations of the CCPN project .tgz files
  #
  
  ccpnProjects = {}  
  for pdbCodeDirectory in directoryList:
    
    nrgCingCcpnProjectsUrl = os.path.join(nrgCingCcpnUrl,pdbCodeDirectory)
    
    directoriesHttpPageLines = getTextFromHttp(nrgCingCcpnProjectsUrl)
    
    for line in directoriesHttpPageLines:
      ccpnProjSearch = ccpnProjPatt.search(line)
      if ccpnProjSearch:
        # Track date and URL for download
        pdbCode = ccpnProjSearch.group(1)
        ccpnProjects[pdbCode] = (ccpnProjSearch.group(2),pdbCodeDirectory,os.path.join(nrgCingCcpnProjectsUrl,"%s.tgz" % pdbCode))

  return ccpnProjects


# Note: needs to move or available to PDBe before this can work there!
try:  
  from sbb.Io import Pickler

  class CingWeeklyDownload(Pickler):


    def __init__(self,pickleDataDir,archiveDir):

      self.pickleDataDir = pickleDataDir
      self.archiveDir = archiveDir

    def updateFiles(self):

      self.getPickleFile("nrgCingPdbCodes")

      newNrgCingPdbCodes = getCingCcpnTgzFileList()

      pdbCodes = newNrgCingPdbCodes.keys()
      pdbCodes.sort()

      for pdbCode in pdbCodes:

        print pdbCode

        downloadTgz = False

        (date,pdbCodeDirectory,downloadUrl) = newNrgCingPdbCodes[pdbCode]

        if not self.nrgCingPdbCodes.has_key(pdbCode):
          downloadTgz = True

        elif self.nrgCingPdbCodes[pdbCode][0] != date:
          # Date change...
          print "  Date change on NRG-CING file for %s..." % (pdbCode)
          downloadTgz = True

        #
        # Create directories and check if file exists, also download if not (could be earlier script problem)
        #

        saveDir = os.path.join(self.archiveDir,pdbCodeDirectory)
        if not os.path.exists(saveDir):
          os.mkdir(saveDir)

        savePath = os.path.join(saveDir,"%s.tgz" % pdbCode)
        if not os.path.exists(savePath):
          downloadTgz = True

        #
        # Now download and update
        #

        if downloadTgz:

          print "  Downloading..."

          data = getDataFromHttp(downloadUrl)

          if os.path.exists(savePath):
            os.remove(savePath)

          fout = open(savePath,'w')
          fout.write(data)
          fout.close()

          # Update pickled information and save... bit slower this but safer
          self.nrgCingPdbCodes[pdbCode] = newNrgCingPdbCodes[pdbCode]
          self.createPickleFile("nrgCingPdbCodes")

except:
  print "  Warning need to install sbb/ directory for NRG-CING downloads!"
  
if __name__ == '__main__':

  cd = getCingSummaryTextInfo("summary.txt")
  
  print cd
