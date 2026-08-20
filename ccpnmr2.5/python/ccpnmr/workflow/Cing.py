"""

Code for plugging CING into workflow

"""

import time, sys, os, webbrowser

from memops.universal.Io import joinPath

from pdbe.adatah.Io import posturl
from pdbe.adatah.Io import getDataFromHttp

from ccpnmr.workflow.Util import WorkFlow

class CingWorkFlow(WorkFlow):
  
  """
  Should this also work from non-CCPN input?
  """
  
  componentList = ['peaks','restraints']
  cingSoftwareName = 'cing'
  
  def cingValidation(self):
  
    """
    Defines full program specific workflow, running code and analysing output
    """  
  
    self.cingRun()
    
    self.analyseCingResults()
  
  def cingRun(self):

    """
    Run CING. Defined in subclasses (web or local)
    """  
  
    pass
  
  def cingAnalyseResults(self):
  
    """
    Analyse results of CING. Defined in here. Might have to be version-specific.
    """  

    pass
    
 
class WebCing(CingWorkFlow):

  """
  Code to run the online version of CING
  """

  # TODO Should become https at one point?
  cmbiSite = "http://nmr.cmbi.ru.nl/"

  cingService = os.path.join(cmbiSite,"icing/serv/iCingServlet")
  cingDownload = os.path.join(cmbiSite,"tmp/cing")
  
  def cingRequest(self,action,uploadFileName=None):
  
    if not action:
      print "Need an action!"
      return
  
    if uploadFileName:

      #
      # First create a file handle
      #
    
      fileHandle = open(uploadFileName)
      fileData = fileHandle.read()
      fileHandle.close()
      
      uploadList = [('UploadFile',"@%s" % uploadFileName,fileData)]

    else:
      uploadList = []
      
    returnPage = posturl(
    
                  self.cingService,
    
                  [('UserId',self.userId),
                   ('AccessKey',self.accessKey),
                   ('Action',action)],
                   
                  uploadList)
    
    return returnPage

  def cingProjectDownload(self,outputDir):
  
    #
    # Put in CCPN project as default. Always usea time-stamped directory name.
    #
    
    if not outputDir:
      outputDir = self.ccpnProject.findFirstRepository(name='userData').url.path
      outputDir = os.path.join(outputDir,'cing',self.getFileTimeString())
  
    #
    # Make sure output dir exists
    #
    
    if not os.path.exists(outputDir):
      os.makedirs(outputDir)
    
    #
    # Get the CING project name - safer in case changed on server side
    #
    
    projectPage = self.cingRequest("ProjectName")

    pyDict = eval(projectPage)
    projectName = pyDict['Result']
    
    #
    # Now download the file - first wait a bit just in case
    #
    
    time.sleep(20)
    
    cingZipFile = "%s_CING_report.zip" % projectName
    cingZipFileHttp = os.path.join(self.cingDownload,self.userId,self.accessKey,cingZipFile)
    
    # Try to re-download if not found.
    outputData = None
    
    while (not outputData or outputData.count("Not Found")):

      time.sleep(20)
      outputData = getDataFromHttp(cingZipFileHttp)
    
    #
    # And write it out...
    #
    
    self.cingZipFileOut = joinPath(outputDir,cingZipFile)

    fout = open(self.cingZipFileOut,'w')
    fout.write(outputData)
    fout.close()

  def cingRun(self,userId,accessKey,downloadDir=None,openBrowser=False,purgeData=True,verbose = True,baseNameSuffix="_testCing"):
    
    #
    # First pack up project
    #
    
    ccpnTgzFile = self.makeCcpnProjectTgz(baseNameSuffix=baseNameSuffix)    
    
    #
    # Initialize
    #
    
    self.userId = userId
    self.accessKey = accessKey
    
    #
    # Upload the file
    #

    initPage = self.cingRequest("Save",uploadFileName = ccpnTgzFile)
    
    if not initPage.count('Success'):
      return None
      
    #
    # Start the CING run
    #
    
    runPage = self.cingRequest("Run")
        
    if not runPage.count('Success'):
      return None    
    
    #
    # Wait for the CING run to finish
    #
    
    statusPage = None
    
    if verbose:
      print "Waiting for termination",
      
    while (not statusPage or 'notDone' in statusPage):
    
      statusPage = self.cingRequest("Status")
      time.sleep(2)
      
      if verbose:
        print ".",
        sys.stdout.flush()
    
    if verbose: 
      print
    
    #
    # Get the log (bit obsolete as also done in file download)
    #
       
    logPage = self.cingRequest("Log")
    
    if verbose:
      pyDict = eval(logPage)
      print pyDict['Result']
        
    #
    # Now get the .tgz of the CING project...
    #
    
    if verbose:
      print "Downloading CING information..."
    
    self.cingProjectDownload(downloadDir)
    
    #
    # Can visualise results if required
    #
    
    if openBrowser:
      url = os.path.join(self.cingDownload,self.userId,self.accessKey)
      browser = webbrowser.get()  # get default browser
      browser.open(url)
  
    #
    # Remove data from server when done (if desired)
    #
    
    if purgeData:
      purgePage = self.cingRequest("Purge")
    
    return True
    
if __name__ == '__main__':

  #
  # Test setup.
  #
  # Make it possible to store userId and accessKey somewhere!!
  #

  if 0:
    userId = "jd3"
    accessKey = "123456"
    ccpnTgzFile = "/Users/wim/tmp/basp.tgz"
  elif 1:
    userId = "wim"
    accessKey = "1test1"
    ccpnTgzFile = "/Users/wim/workspace/stable/all/python/ccpnmr/workflow/local/taf3_tagged.tgz"
  elif 0:
    # This one takes forever...
    userId = "wim"
    accessKey = "2test2"
    ccpnTgzFile = "/Users/wim/tmp/cippalippa.tgz"

  wc = WebCing(ccpnProjectTgz = ccpnTgzFile)
  wc.cingRun(userId,accessKey,downloadDir='./test/',openBrowser=True,purgeData=False)
