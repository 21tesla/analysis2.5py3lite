"""

Code for plugging ARIA into workflow

"""

import os
import re
import time

from ccpnmr.workflow.Util import WorkFlow


class AriaWorkFlow(WorkFlow):
    componentList = ["peaks", "restraints"]
    ariaSoftwareName = "aria"

    def ariaCalculation(self):
        """
        Defines full program specific workflow, running code and analysing output
        """

        self.ariaRun()

        self.ariaAnalyseResults()

    def ariaRun(self):
        """
        Run ARIA. Defined in subclasses (web or local)
        """

        pass

    def ariaAnalyseResults(self):
        """
        Analyse results of ARIA. Defined in here. Might have to be version-specific.
        """

        pass

    def convertNmrSimIntoTags(self):
        """
        Temporary function!
        """

        ariaNmrSimStore = self.ccpnProject.findFirstNmrSimStore(name="ARIA")

        if ariaNmrSimStore:
            runs = ariaNmrSimStore.sortedRuns()
            runs.reverse()

            for run in runs:
                #
                # Use latest run that does not have any output
                #

                if not run.outputConstraintStore:
                    self.molSystem = run.molSystem
                    for chain in self.molSystem.sortedChains():
                        self.tagCcpnObjectForExternalExec(chain, self.ariaSoftwareName)

                    self.nmrConstraintStore = run.inputConstraintStore

                    for constraintList in self.nmrConstraintStore.sortedConstraintLists():
                        self.tagCcpnObjectForExternalExec(constraintList, self.ariaSoftwareName)

                    for peakList in run.inputPeakLists:
                        self.tagCcpnObjectForExternalExec(peakList, self.ariaSoftwareName)

                    # Note: structures are NOT done.

                    break

                """
        print(run.inputConstraintStore)
        print(run.inputMeasurementLists)
        print(run.inputPeakLists)
        print(run.inputStructures)
        print(run.molSystem)
        print(run.outputConstraintStore)
        print(run.outputEnsemble)
        print(run.outputMeasurementLists)
        print(run.outputPeakLists)
        """


# TODO: change all url handling to this code! Much nicer... so posturl stuff can be changed?
import http.cookiejar as cookielib  # py2 cookielib moved to http.cookiejar
import urllib
import urllib.request

from pdbe.adatah.Io import MultipartPostHandler


# Note: might not need this subclass as long as methods don't overlap between local and other
# execs!
class AriaCcpnGrid(AriaWorkFlow):
    """
    Code to run the online CCPNGRID version of ARIA
    """

    ccpnGridSite = "http://webapps.ccpn.ac.uk/"

    ccpnGridLoginUrl = os.path.join(ccpnGridSite, "accounts/login/")
    ccpnGridLocal = "ccpngrid"

    ccpnGridLocalUpload = os.path.join(ccpnGridLocal, "upload")

    ccpnGridStatusUrl = os.path.join(ccpnGridSite, ccpnGridLocal, "status")

    ccpnGridHiddenFieldPatt = re.compile('\\<input.+type\\="hidden"')
    # ccpnGridCheckBoxPatt = re.compile('input.+type\=\"checkbox\"')
    # ccpnGridSelectStartPatt = re.compile('\<select')
    # ccpnGridSelectEndPatt = re.compile('\<\/select')

    ccpnGridNamePatt = re.compile('name\\="([^"]+)"')
    ccpnGridValuePatt = re.compile('value\\="([^"]*)"')

    ccpnGridLinkPatt = re.compile(r"\<a\>([^<]+)\<\/a\>.+\<a\>([^<]+)\<\/a\>")
    ccpnGridStatusPatt = re.compile(r"\<a style\=[^>]+\>([^<]+)\<\/a\>\<a\>([^<]+)\<\/a\>")

    userId = "test"
    password = "test123"

    def ariaRun(self, userId, password, inputFile, uniqueIdentifier=None):

        # Nov 2011: Django 1.3 requires csrf tokens to be passed around

        uploadUrl = self.ccpnGridLogin(userId, password)

        returnPage = self.ccpnGridUploadProject(uploadUrl, inputFile, uniqueIdentifier=uniqueIdentifier)

        return returnPage

    def ccpnGridLogin(self, userId, password):

        #
        # Use default test userId, password...
        #

        if not userId:
            userId = self.userId

        if not password:
            password = self.password

        #
        # Log in to the website
        #

        loginToCcpnGridUrl = os.path.join(self.ccpnGridLoginUrl, "?next=/", self.ccpnGridLocalUpload)

        cookies = cookielib.CookieJar()
        cookieHandler = urllib.request.HTTPCookieProcessor(cookies)
        self.connection = urllib.request.build_opener(cookieHandler, MultipartPostHandler)
        urllib.request.install_opener(self.connection)

        # do this GET only so that can get hold of csrf token
        req = urllib.request.Request(loginToCcpnGridUrl)
        handle = urllib.request.urlopen(req)

        self.csrf_cookie = None
        for cookie in cookies:
            if cookie.name == "csrftoken":
                self.csrf_cookie = cookie
                break
        if not self.csrf_cookie:
            raise OSError("No csrf cookie found")

        # login using the usr, pwd, and csrf token
        dd = {
            "username": userId,
            "password": password,
        }
        login_data = urllib.urlencode(dd)
        """  perhaps do not need to urlencode but play safe
    login_data = urllib.urlencode( dict(
        username=userId, password=password ) )
    """

        req = urllib.request.Request(loginToCcpnGridUrl, login_data)
        # below is the key line; this is how you pass csrf token to POST
        req.add_header("X-CSRFToken", self.csrf_cookie.value)
        returnPage = urllib.request.urlopen(req)

        uploadUrl = returnPage.url

        return uploadUrl

    def ccpnGridUploadProject(self, uploadUrl, inputFile, uniqueIdentifier=None):

        #
        # This unique identifier should be time-stamped (but in a way that CCPN can handle it)
        # Should do this for all runs consistently - also fits into John's stuff I imagine
        #

        if not uniqueIdentifier:
            (path, baseName) = os.path.split(inputFile)
            uniqueIdentifier = baseName.replace(".tgz", "")

        uniqueIdentifier = "%s.%s" % (uniqueIdentifier, self.getFileTimeString())
        self.uniqueIdentifier = uniqueIdentifier

        uploadInfo = {"title": uniqueIdentifier, "file": open(inputFile, "rb"), "tagged": "on"}

        #
        # So what comes back here depends on whether project is tagged or not, will only allow tagged projects!!!
        #

        ##returnPage = self.connection.open(uploadUrl, uploadInfo)
        req = urllib.request.Request(uploadUrl, uploadInfo)
        if self.csrf_cookie:
            req.add_header("X-CSRFToken", self.csrf_cookie.value)
        returnPage = urllib.request.urlopen(req)

        # print returnPage
        # print returnPage.url
        # print returnPage.readlines()

        return returnPage

    def getStatusPageInfo(self, userId="test", password="test123"):

        #
        # Connect to status page - should probably keep this open and re-access!
        #
        # Can do this by setting as self (in dict), then reopening?
        #

        # statusPage = self.connection.open( self.ccpnGridStatusUrl,  {} )
        req = urllib.request.Request(self.ccpnGridStatusUrl)
        if self.csrf_cookie:
            req.add_header("X-CSRFToken", self.csrf_cookie.value)
        statusPage = urllib.request.urlopen(req)

        #
        # So what comes back here depends on whether project is tagged or not, will only allow tagged projects!!!
        #
        # from ccpnmr.workflow.local.ccpnGridStatusPage import linesProjectRunning as lines

        lines = statusPage.readlines()

        infoDict = {}

        fieldType = None
        name = None
        status = None

        lineIndex = 0

        while lineIndex < len(lines):
            dataLine = lines[lineIndex]

            if dataLine.count('<table class="status"'):
                status = None
                while not status:
                    lineIndex += 1
                    dataLine = lines[lineIndex]

                    searchLink = self.ccpnGridLinkPatt.search(dataLine)
                    if searchLink:
                        (statusId, time) = searchLink.groups()
                        continue

                    statusLink = self.ccpnGridStatusPatt.search(dataLine)
                    if statusLink:
                        (status, iteration) = statusLink.groups()

                        infoDict[statusId] = {"time": time, "status": status, "iteration": iteration}
                        break

            elif self.ccpnGridHiddenFieldPatt.search(dataLine):
                name = self.ccpnGridNamePatt.search(dataLine).group(1)
                infoDict[statusId][name] = self.ccpnGridValuePatt.search(dataLine).group(1)

            lineIndex += 1

        return infoDict


if __name__ == "__main__":
    aria = AriaCcpnGrid(identifier="test")

    aria.ariaRun(None, None, "local/taf3_tagged.tgz", uniqueIdentifier="test")

    status = None

    while status not in ("Finished", "Failed"):
        infoDict = aria.getStatusPageInfo()
        print(aria.uniqueIdentifier, infoDict[aria.uniqueIdentifier])

        status = infoDict[aria.uniqueIdentifier]["status"]

        time.sleep(5)
