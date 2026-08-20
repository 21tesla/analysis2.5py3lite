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

from memops.gui.IntEntry        import IntEntry
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.Label           import Label
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.MessageReporter import showWarning

from memops.editor.BasePopup    import BasePopup
from memops.editor.Util         import createDismissHelpButtonList

headingColor  = '#80C080'

class HaddockDaniPopup(BasePopup):

    def __init__(self, parent, run):

        self.parent         = parent
        self.latestRun      = run
        self.daniprotocol   = None
        self.waiting        = False

        self.getStoredDaniProtocols()
        BasePopup.__init__(self, parent=parent, title='Relaxation data')

    def body(self, guiFrame):
        
        self.geometry('500x500')
        guiFrame.grid_columnconfigure(0, weight=1)
        guiFrame.grid_rowconfigure(1, weight=1)

        frame = LabelFrame(guiFrame, text='Relaxation protocol settings')
        frame.grid(row=1,column=0,sticky='nsew')
        frame.grid_columnconfigure(0, weight=1)
        frame.grid_rowconfigure(0, weight=1)
        
        self.floatEntry  = FloatEntry(self, returnCallback=self.setValue)
        
        headingList      = ['Parameter','Value','Description']
        justifyList      = ['center','center', 'left']
        editWidgets      = [None, self.floatEntry, None]    
        editGetCallbacks = [None, self.getValue, None]
        editSetCallbacks = [None, self.setValue, None]
        self.daniMatrix = ScrolledMatrix(frame, headingList=headingList,
                                        justifyList=justifyList,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks, 
                                        editWidgets=editWidgets,
                                        multiSelect=False, initialRows=10,
                                        passSelfToCallback=True,
                                        callback=self.selectDani)

        self.daniMatrix.grid(row=0, column=0, sticky='nsew')
        self.daniMatrix.refreshFunc = self.updateDaniSet

        self.updateDaniSet()
        
    def getStoredDaniProtocols(self):
        
        self.daniprotocols = [i for i in self.latestRun.sortedHaddockEnergyTerms() if i.code == 'daniProtocolStore' ]

    def selectDani(self, obj, row, col, table):

        self.table = table
        self.daniprotocol = obj

    def updateDaniSet(self):
        
        textMatrix  = []; objectList  = []; colorMatrix = []
        
        for protocol in self.daniprotocols:
            textMatrix.append(['DANI protocol %i' % protocol.termId, None, None])
            objectList.append(None)
        
            for term in protocol.sortedEnergyTermParameters():    
                textMatrix.append([term.code, term.value, 'Energy constant for term %s' % term.code])
                objectList.append([term,'value', self.floatEntry, self.getFloat, self.setFloat])
    
        n = len(textMatrix)
        colors = [None] * 3
        colorMatrix = [colors for i in range(n)]
        if len(colorMatrix) > 0:
            stages = range(0,n,12)
            for s in stages: colorMatrix[s]  = [headingColor] * 3
        
        self.daniMatrix.update(colorMatrix=colorMatrix, objectList=objectList, textMatrix=textMatrix)

    def getValue(self, rowObj):

        obj, attrName, widget, getter, setter = rowObj
        
        self.table.editWidget = widget
        
        getter(widget, obj, attrName)

    def setValue(self, event, null=None): 
        
        obj, attrName, widget, getter, setter = self.daniprotocol

        if setter: setter(widget, obj, attrName)

        self.table.refreshFunc()

    def getFloat(self, widget, obj, attrName):

        value = getattr(obj, attrName) 
        widget.set(value)

    def setFloat(self, widget, obj, attrName):

        value = widget.get() or 0.0
        setattr(obj, attrName, value)    

    def destroy(self):

        BasePopup.destroy(self)

