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

import sys

from memops.gui.Frame           import Frame
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.Util            import createDismissHelpButtonList
from memops.editor.BasePopup    import BasePopup
from memops.gui.FileSelectPopup import FileSelectPopup

from HaddockFrame               import HaddockFrame

"""
Description: Graphical user interface around the Utrecht HADDOCK molecular docking 
             package.
Use           : Package can be used as macro for the Analysis software suite by loading
             the function 'HaddockGUI' from the macro import menu or it can be used
             in stand alone mode. 
"""

def HaddockGUI(argServer):

    """Function to load Haddock GUI as CCPN Macro"""

    popup = HaddockPopup(argServer.parent, argServer.getProject())
    popup.open()

def standAloneUse(argv=None):
    
    """Use the Haddock GUI in stand alone mode"""
    
    import Tkinter, os
    from memops.general.Io import loadProject
    
    guiRoot = Tkinter.Tk()
    project = None
    
    if len(argv) > 1:
        project = loadProject(argv[1])
    else:
        popup = FileSelectPopup(guiRoot, show_file=False)
        if popup.getDirectory():
            project = loadProject(popup.getDirectory())
        else:
            print("No project defined means exit")
        popup.destroy()
    
    if project:
        guiRoot.withdraw()
        popup = HaddockPopup(guiRoot, project)    
        guiRoot.mainloop()
    
class HaddockPopup(BasePopup):
    
    def __init__(self, parent, ccpnProject):

        self.parent      = parent
        self.ccpnProject = ccpnProject

        BasePopup.__init__(self, parent=parent, title='HADDOCK Project Setup')

        self.font = 'Helvetica 12'
        self.setFont()
                    
    def body(self, guiFrame):

        # Ensure that the first row and column in popup expand
        self.geometry('820x500')
        guiFrame.grid_rowconfigure(0, weight=1)
        guiFrame.grid_columnconfigure(0, weight=1, minsize=300)

        self.hFrame = HaddockFrame(guiFrame, self.ccpnProject)
        self.hFrame.grid(row=0, column=0, sticky='nsew')

        texts    = ['Save Project']
        commands = [self.save]
        self.bottomButtons = createDismissHelpButtonList(guiFrame, texts=texts, expands=True, commands=commands)
        self.bottomButtons.grid(row=1, column = 0, sticky = 'ew')

    def save(self):
        
        self.hFrame.ccpnProject.saveModified()
        
    def close(self):
  
        BasePopup.destroy(self)
        #sys.exit(0)
 
    def destroy(self):
    
        BasePopup.destroy(self)
    
if __name__ == '__main__':

    standAloneUse(sys.argv)

  
  
  
  
  
