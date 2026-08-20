
"""
======================COPYRIGHT/LICENSE START==========================

ExperimentCreatePopup.py: GUI window to create new experiment

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)
- PDBe website (http://www.ebi.ac.uk/pdbe/)

- contact Wim Vranken (wim@ebi.ac.uk)
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

import string, Tkinter, sys

from memops.universal.Io import joinPath

from ccpnmr.format.general.Io import getHelpUrlDir
from ccpnmr.format.gui.BasePopup import TemporaryBasePopup

from ccp.general.Util import setCurrentStore

from memops.gui.ScrolledMatrix import ScrolledMatrix, FilterMatrix
from memops.gui.HelpPopup import showHelpUrl
from memops.gui.Label import Label
from memops.gui.DataEntry import askString
from memops.gui.MessageReporter import showError

from memops.gui.Util import createDismissHelpButtonList

class ExperimentCreatePopup(TemporaryBasePopup):

  help_url = joinPath(getHelpUrlDir(),'ExperimentCreate.html')

  def __init__(self, parent, project, numDim = None, topText = None, expName = ''):
  
    self.parent = parent
    self.project = project
    self.numDim = numDim
    self.topText = topText
    self.skip = 1
    self.expName = expName
    
    setCurrentStore(project,'ChemElementStore')

    TemporaryBasePopup.__init__(self, parent=parent, title="Project '%s': " % project.name + 'Create Experiment', modal=False, transient=True)

  def body(self, parent):

    row = 0
    columnspan = 3
    
    if self.numDim != None:
      columnspan += 1

    if self.topText:
      label = Label(parent, text= self.topText)
      label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)
    
      row = row + 1
      
    if self.expName:
      expNameText = " for experiment '%s'" % self.expName
    else:
      expNameText = ''

    label = Label(parent, text= 'Select a reference experiment%s from the list below, then press Create' % expNameText)
    label.grid(row=row, column=0, columnspan = columnspan, sticky=Tkinter.EW)

    row = row + 1
      
    colHeadings = self.getHeadings()
      
    parent.grid_rowconfigure(row, weight=1, minsize=300)
    parent.grid_columnconfigure(0, weight=1, minsize=100)
    
    self.table = ScrolledMatrix(parent, headingList=colHeadings,callback=None)
    self.table.grid(row = row, column = 0, columnspan = columnspan, sticky=Tkinter.NSEW)
    
    row = row + 1

    texts = [ 'Create' ]
    commands = [ self.ok ]   # This calls 'ok' in BasePopup, this then calls 'apply' in here

    if self.numDim != None:
      texts.append('Show all experiments')
      commands.append(self.showAll)

    buttons = createDismissHelpButtonList(parent, texts=texts, commands=commands, dismiss_text = 'Skip', help_url=self.help_url)
    buttons.grid(row=row, column=0)
    
    self.update()
        
  def update(self):

    textMatrix = []
    objectList = []
  
    for nep in self.project.sortedNmrExpPrototypes():
    
      #print "%d: %s (%s)   %s" % (nep.serial, nep.name, nep.category, nep.synonym)
    
      for re in nep.refExperiments:
        """
        #print "  %s:" % re.name,
        atomDimList = []
        for refExpDim in re.refExpDims:
          for refExpDimRef in refExpDim.refExpDimRefs:
            expMeasurement = refExpDimRef.expMeasurement
            atoms = []
            addText = ""
            for atomSite in expMeasurement.atomSites:
              atoms.append(atomSite.isotopeCode[-1])
            if expMeasurement.measurementType != "Shift":
              addText = "(%s)" % expMeasurement.measurementType
            atomDimList.append("%s%s" % (string.join(atoms,''),addText))
            string.join(atomDimList,'-'),
        """ 
        valueList = [nep.synonym,re.name,len(re.refExpDims),nep.name,nep.category]
        objectList.append(re)
        textMatrix.append(valueList)
    
    self.table.update(textMatrix = textMatrix, objectList = objectList)

    if self.numDim != None:
      filterClass = FilterMatrix(self.table.parent, self.table, searchString = str(self.numDim), column = 2)
      filterClass.apply()
  
  def showAll(self):
  
    self.numDim = None
    self.update()
  
  def getHeadings(self):
  
    colHeadings = ['Common pathway name','Experiment','Dimensions','Pathway name','Category']
    
    return colHeadings
    
  def apply(self):
  
    self.refExperiment = self.table.currentObject
    
    names = [exp.name for exp in self.project.currentNmrProject.sortedExperiments()]

    while 1:
      self.name = askString('Give the name of the experiment','Experiment name:',initial_value = self.expName, parent = self.parent)
      
      if (self.name in names):
        showError('Repeated experiment name', 'Name already used - choose another.')
      else:
        break

    self.skip = 0
    
    return True

if __name__ == '__main__':

  root = Tkinter.Tk()

  from memops.api import Implementation
  from ccp.api.nmr import Nmr
  project = Implementation.MemopsRoot(name = 'test')
  nmrProject = Nmr.NmrProject(project, name = project.name)
  
  popup = ExperimentCreatePopup(root,project,numDim = 2)
  
  root.mainloop()
