
"""
======================COPYRIGHT/LICENSE START==========================

ReferenceFrame.py: Part of the CcpNmr Analysis program

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
from memops.gui.Button import Button
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showYesNo

from ccpnmr.analysis.core.Util import stringFromExperimentSpectrum

class ReferenceFrame(Frame):

  def __init__(self, parent, peakList):
 
    apply(Frame.__init__, (self, parent) )

    #parent.grid_columnconfigure(0, weight=1)
    #parent.grid_rowconfigure(1, weight=1)
    #self.grid_columnconfigure(1, weight=1)
    #self.grid_columnconfigure(2, weight=1)

    self.parent      = parent
    self.peakList    = peakList
    self.specLabel   = Label(self, text='Spectrum: '  )
    self.listLabel   = Label(self, text='Peak List: ' )
    self.dimLabel    = Label(self, text='Dimension' )
    self.adjustLabel = Label(self, text='Adjustment')
    self.specLabel.grid(row = 0, column = 0, columnspan=2, sticky='nsew')
    self.listLabel.grid(row = 0, column = 2, columnspan=1, sticky='nsew')
    self.dimLabel.grid(row = 1, column = 0, columnspan=1, sticky='nsew')
    self.adjustLabel.grid(row = 1, column = 1, columnspan=2, sticky='nsew')
    
    self.goButton     = Button(self, text='Go!', command=self.go,
                               tipText='Commit the re-referencing adjustment')
    self.clearButton  = Button(self, text='Clear', command=self.update,
                               tipText='Clear the input re-referencing values')
    self.cancelButton = Button(self, text='Cancel', command=self.close,
                               tipText='Abort operation without making any changes')
   
    self.update()
   
  def update(self):
  
    spec = self.peakList.dataSource.name
    expt = self.peakList.dataSource.experiment.name
    es = stringFromExperimentSpectrum(expt,spec)
    self.specLabel.set('Spectrum: '+es)
    self.listLabel.set('Peak List: List %d' % self.peakList.serial)

    numDims = self.peakList.dataSource.numDim
    self.numDims = numDims
    self.dimLabel=numDims * ['']
    self.adjustEntry=numDims * ['']
    self.adjusts=numDims * [0]
    
    for i in range(numDims):
      self.dimLabel[i] = Label(self, text = "F%d" % (i+1) )
      self.dimLabel[i].grid(row = i+2, column = 0, columnspan=1, sticky='w')
      self.adjustEntry[i] = FloatEntry(self, text = "%f" % self.adjusts[i],
                                       width=10, tipText='The referencing adjustment to add for this dimension')
      self.adjustEntry[i].grid(row = i+2, column = 1, columnspan=2, sticky='w')
      
    self.goButton     .grid(row = i+3, column = 0, columnspan=1, sticky='nsew')
    self.clearButton  .grid(row = i+3, column = 1, columnspan=1, sticky='nsew')
    self.cancelButton .grid(row = i+3, column = 2, columnspan=1, sticky='nsew')
  
  def go(self):
  
    if showYesNo('Re-reference', "Are you sure?"):
            
      for i in range(self.numDims):
        
        self.adjusts[i] = self.adjustEntry[i].get()
	      
      for peak in self.peakList.peaks:
        
	dims = peak.peakDims
        
        for i, peakDim in enumerate(peak.sortedPeakDims()):
          
	  peakDim.position = peakDim.position + self.adjusts[i]
        
    self.parent.parent.owner.peaksUpdate()
    self.close() 
    
      
