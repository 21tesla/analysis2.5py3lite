#!/usr/bin/env python

"""
=========================================================================
Package:    - Code for the Graphical User Interface fontend to the 
              Haddock model package in the CCPN data model. 
            - Code for the export of a Haddock compatible project. A 
              Haddock compatible project can either be a parameter file
              ready for submission to the Haddock webserver or a
              directory structure with necessary files for use with a 
              localy installed version of Haddock.

Dependencies: The CCPN Haddock package requires CCPN data model version
              2.0 or higher. The export of a webserver compatible 
              parameter file requires Haddock webserver version 2.1 or 
              higher and a valid user account. The export of a 'classic' 
              Haddock project requires Haddock version 2.0 or higher.

Copyright and License information:
              The Haddock data model as implemented in the CCPN data
              model as well as the use of CCPN GUI code elements is 
              licenced to the CCPN Projects (Copyright (C) 2008) and
              distributed under the terms of the GNU Lesser General
              Public License.
            
              The Haddock project export code as well as the use of 
              Haddock software is covert in the Haddock License
              agreement (Copyright (C) 2008 Haddock Project, Bijvoet
              Center for Biomolecular Research, Utrecht University,
              The Netherlands).

GNU LGPL:        This library is free software; you can redistribute it 
              and/or modify it under the terms of the GNU Lesser General 
              Public License as published by the Free Software 
              Foundation; either version 2.1 of the License, or (at 
              your option) any later version.
 
              A copy of this license can be found in LGPL.license
 
              This library is distributed in the hope that it will be 
              useful, but WITHOUT ANY WARRANTY; without even the implied 
              warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
              PURPOSE. See the GNU Lesser General Public License for 
              more details.
 
              You should have received a copy of the GNU Lesser General 
              Public License along with this library; if not, write to 
              the Free Software Foundation, Inc., 59 Temple Place, Suite 
              330, Boston, MA 02111-1307 USA.

Information:  For further information regarding CCPN, please contact:
              - CCPN website (http://www.ccpn.ac.uk/)
              - email: ccpn@bioc.cam.ac.uk
              
              For further information regarding Haddock, please contact
              Alexandre M.J.J. Bonvin:
              - http://haddock.chem.uu.nl
              - email: a.m.j.j.bonvin@uu.nl    

Citing:          If you are using this software for academic purposes, we 
                suggest quoting the following references:

              For CCPN:    
              Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
              Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. 
              Ionides and Ernest D. Laue (2005). A framework for 
              scientific data modeling and automated software development. 
              Bioinformatics 21, 1678-1684.
            
              For Haddock:
              Cyril Dominguez, Rolf Boelens and Alexandre M.J.J. Bonvin 
              (2003). HADDOCK: a protein-protein docking approach based 
              on biochemical and/or biophysical information. 
              J. Am. Chem. Soc. 125, 1731-1737.
            
              S.J. de Vries, A.D.J. van Dijk, M. Krzeminski, M. van Dijk, 
              A. Thureau, V. Hsu, T. Wassenaar and A.M.J.J. Bonvin (2007) 
              HADDOCK versus HADDOCK: New features and performance of 
              HADDOCK2.0 on the CAPRI targets. 
              Proteins: Struc. Funct. & Bioinformatic 69, 726-733.    
=========================================================================
"""
import httplib, urlparse, urllib, sys, os

from memops.gui.Button          import Button
from memops.gui.LabeledEntry    import LabeledEntry
from memops.gui.LabelFrame      import LabelFrame
from memops.editor.BasePopup    import BasePopup
from memops.editor.Util         import createDismissHelpButtonList
from memops.gui.MessageReporter import showWarning

from HaddockExportParam         import exportParam

class HaddockServerUpload(BasePopup):

    def __init__(self,parent,hProject=None,latestRun=None,ccpnProject=None):

        self.parent         = parent  #GUI parent
        self.hProject        = hProject
        self.latestRun         = latestRun
        self.ccpnProject     = ccpnProject
        
        self.username        = None
        self.password        = None
        
        BasePopup.__init__(self, parent=parent, title='HADDOCK server upload')

    def body(self, guiFrame):

        guiFrame.grid_columnconfigure(1, weight=1)
        guiFrame.grid_rowconfigure(1, weight=1)

        self.usernameEntry = LabeledEntry(guiFrame, 'username')
        self.usernameEntry.grid(row=0)

        self.passwordEntry = LabeledEntry(guiFrame, 'password', show="*")
        self.passwordEntry.grid(row=1)
        
        texts = ['Upload']
        commands = [self.serverUpload]
        buttonList = createDismissHelpButtonList(guiFrame, texts=texts, commands=commands, expands=True)
        buttonList.grid(row=2)

    def serverUpload(self):
        
        """Process username and password. If both valid then start uploading otherwise issue warning.
           Export the project and run as a .web file to disk. Than upload this file
        """
        warning = ''
        
        if self.usernameEntry.getEntry(): self.username = self.usernameEntry.getEntry()
        else: warning += ' No username defined ' 
        
        if self.passwordEntry.getEntry(): self.password = self.passwordEntry.getEntry()
        else: warning += ' No password defined '            
        
        if len(warning): showWarning('Warning', warning, parent=self)
        else:
            exportParams = exportParam(hProject=self.hProject,
                                       latestRun=self.latestRun,
                                       ccpnProject=self.ccpnProject)
            if len(exportParams.filestring): upload = ServerUpload(exportParams.filestring,
                                                                     self.hProject.name,
                                                                     str(self.latestRun.serial),
                                                                     self.username,
                                                                     self.password)
            else: print("ERROR: Export project as parameter file failed")
                
    def open(self):

        BasePopup.open(self)    

    def destroy(self):

        BasePopup.destroy(self)    

class ServerUpload(object):

    """Class for automatic upload of a HADDOCK webserver parameter file.
       Input: filestring - a compatible parameter file as string
              project - a HADDOCK project name
              runname - a run number a string
              username - a valid HADDOCK server user name
              password - a valid HADDOCK server password
              server - URL of HADDOCK server (predefined)
              cgi - cgi script accepting the parameter file (predefined)
       Output: Server exeptance or rejection messages
    """

    def __init__(self,filestring,project,runname,username,password,
                 server='http://haddock.chem.uu.nl/',cgi='cgi-bin/haddockserver4.cgi'):
        
        self.files = [('params',project,filestring)]
        self.url = server+cgi
        
        print("** Upload project file to HADDOCK server **")
        print("Project: %s\nRun: %s" % (project,runname))
        
        self.__upload(runname,username,password)
        
        print("** Upload complete **")

    def __encode_multipart_formdata(self,fields,files):    

        BOUNDARY = '----------1234567890abcdefghij_$'
        CRLF = '\r\n'
        L = []

        for file in files:
            (filekey, filename, filevalue) = file
            L.append('--' + BOUNDARY)
            L.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (filekey, filename))
            L.append('Content-Type: text/plain')
            L.append('')
            for line in filevalue.split('\n'): L.append(line)    

        for (key, value) in fields:
            if value == None: value = ""
            L.append('--' + BOUNDARY)
            L.append('Content-Disposition: form-data; name="%s"' % str(key))
            L.append('')
            L.append(str(value))

        L.append('--' + BOUNDARY + '--')
        L.append('')
        body = CRLF.join(L)
        content_type = 'multipart/form-data; boundary=%s' % BOUNDARY

        return content_type, body

    def __upload(self,runname,username,password):    

        fields = [('runname', runname), ('username', username), ('password', password)]

        content_type, body = self.__encode_multipart_formdata(fields,self.files)
        urlparts = urlparse.urlsplit(self.url)
        host = urlparts[1]
        selector = urlparts[2]    
        http = httplib.HTTPConnection(host)
        headers = {'User-Agent':'anonymous', 'Content-Type':content_type}
        http.request('POST', selector, body, headers)
        res = http.getresponse() 
        response = res.read()
        
        print("HADDOCK server response:")
        print response

