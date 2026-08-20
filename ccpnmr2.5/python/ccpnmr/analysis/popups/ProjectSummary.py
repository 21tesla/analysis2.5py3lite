
"""
======================COPYRIGHT/LICENSE START==========================

ProjectSummary.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2014 Wayne Boucher and Tim Stevens (University of Cambridge)

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

from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.LabelDivider    import LabelDivider
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.WebBrowser      import WebBrowser

from ccp.general.Io import getDataSourceFileName

from ccpnmr.analysis.popups.BasePopup import BasePopup

from ccpnmr.analysis.core.AssignmentBasic import isPeakAssigned, getAtomSetShifts

TABLE_BASE_URL = 'http://www2.ccpn.ac.uk/table'

class Table:

  def __init__(self, colHeadings):

    self.colHeadings = colHeadings
    self.rowObjects = None
    self.rowTexts = []

  def update(self, rowTexts, rowObjects=None):

    self.rowTexts = rowTexts
    self.rowObjects = rowObjects

class ScrolledTable(Table, ScrolledMatrix):

  def __init__(self, guiParent, colHeadings, **kw):

    Table.__init__(self, colHeadings)
    kw['headingList'] = colHeadings
    ScrolledMatrix.__init__(self, guiParent, **kw)
    
  def update(self, rowTexts, rowObjects):

    Table.update(self, rowTexts, rowObjects)
    ScrolledMatrix.update(self, textMatrix=rowTexts, objectList=rowObjects)

class ProjectSummaryPopup(BasePopup):
  """
  **An Executive Summary of the Project**
  
  This popup window provides a executive summary of the project.
  It gives informaton about the project name, the project (non-spectra)
  directories, the spectra, the peak lists, the chains, the restraint
  lists and the structure ensembles.

  The "Update" button will update the information, when the project has
  changed.

  The "Show in Browser" button will open up a browswer window with the
  information in HTML, which can then be printed (which on many computers
  also allows the output to be saved to a file).
  """

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    title = "Project : Summary"

    BasePopup.__init__(self, parent=parent, title=title, **kw)
   
  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)
    
  def body(self, guiFrame):
    
    self.geometry("+250+250")

    guiFrame.grid_columnconfigure(0, weight=1)
    
    row = 0

    frame = Frame(guiFrame, grid=(row,0))
    frame.expandGrid(0,0)
    guiFrame.grid_rowconfigure(row, weight=1)
    row += 1

    self.projectNameLabel = Label(frame, text='Project name:', grid=(0,0), sticky='w')
    utilButtons = UtilityButtonList(frame, helpUrl=self.help_url, grid=(0,1), sticky='e')

    div = LabelDivider(guiFrame, text='Project directories', grid=(row,0))
    row += 1

    tipTexts = ['The record number',
                'The repository name',
                'The full path of the repository',
               ]
    colHeadings = ['#','Name','Path']
                   
    self.projectTable = ScrolledTable(guiFrame, colHeadings, initialRows=3,
                                      grid=(row,0), tipTexts=tipTexts)
    guiFrame.grid_rowconfigure(row, weight=1)
    row += 1
    
    div = LabelDivider(guiFrame, text='Spectra', grid=(row,0))
    row += 1

    tipTexts = ['The record number',
                'The spectrum name',
                'The number of dimensions for the spectrum',
                'The shift list for the spectrum',
                'The full path to the spectrum data file',
               ]
    colHeadings = ['#','Name','Num. Dim','Shift List','Path']
                   
    self.spectrumTable = ScrolledTable(guiFrame, colHeadings, initialRows=5,
                                         grid=(row,0), tipTexts=tipTexts)
    guiFrame.grid_rowconfigure(row, weight=1)
    row += 1
    
    div = LabelDivider(guiFrame, text='PeakLists', grid=(row,0))
    row += 1

    tipTexts = ['The record number',
                'The peak list name',
                'The number of peaks for the peak list',
                'The number assigned in at least one dimension',
                'The percentage assigned in at least one dimension',
                'The number assigned in all dimensions',
                'The percentage assigned in all dimensions',
               ]
    colHeadings = ['#','Name','Peaks','Part Assigned','Part Percent', 'All Assigned', 'All Percent']
                   
    self.peakListTable = ScrolledTable(guiFrame, colHeadings, initialRows=5,
                                       grid=(row,0), tipTexts=tipTexts)
    guiFrame.grid_rowconfigure(row, weight=1)
    row += 1
    
    div = LabelDivider(guiFrame, text='Chains', grid=(row,0))
    row += 1

    tipTexts = ['The chain number',
                'The chain name',
                'The number of residues in the chain',
                'The number of assignable atoms in the chain',
                'The shift list being considered for assignment',
                'The number of assigned atoms in the chain for given shift list',
                'The percentage of atoms assigned in the chain for given shift list',
               ]
    colHeadings = ['#','Name','Residues','Assignable\nAtoms','Shift List', 'Assigned\nAtoms', 'Percent\nAssigned']
                   
    self.chainTable = ScrolledTable(guiFrame, colHeadings, initialRows=3,
                                    grid=(row,0), tipTexts=tipTexts)
    guiFrame.grid_rowconfigure(row, weight=1)
    row += 1
    
    div = LabelDivider(guiFrame, text='Restraint Lists', grid=(row,0))
    row += 1

    tipTexts = ['The restraint list number',
                'The restraint list name',
                'The type of restraints in the restraint list',
                'The number of restraints in the restraint list',
                'The experiments related to the restraint list',
               ]
    colHeadings = ['#','Name','Type','Restraints','Experiments']
                   
    self.restraintListTable = ScrolledTable(guiFrame, colHeadings, initialRows=5,
                                            grid=(row,0), tipTexts=tipTexts)
    guiFrame.grid_rowconfigure(row, weight=1)
    row += 1
    
    div = LabelDivider(guiFrame, text='Structure Ensembles', grid=(row,0))
    row += 1

    tipTexts = ['The structure ensemble number',
                'The structure ensemble name',
                'The number of models in the structure ensemble',
                'The chains in the structure ensemble',
                'The total number of residues in the structure ensemble',
               ]
    colHeadings = ['#','Name','Models','Chains','Residues']
                   
    self.structureTable = ScrolledTable(guiFrame, colHeadings, initialRows=4,
                                        grid=(row,0), tipTexts=tipTexts)
    guiFrame.grid_rowconfigure(row, weight=1)
    row += 1
    
    #tipTexts = ['Print the information to a file']
    #texts    = ['Print to File']
    #commands = [self.printInfo]
    tipTexts = ['Update the information', 'Show the information in a browser (e.g. for printing)']
    texts    = ['Update', 'Show in Browser']
    commands = [self.updateAfter, self.showBrowser]
    self.bottomButtons = ButtonList(guiFrame, texts=texts, commands=commands,
                                    grid=(row,0), tipTexts=tipTexts)

    self.curateNotifiers(self.registerNotify)
    
    self.refresh = False 
    self.updateAfter()
  
  def curateNotifiers(self, notifyFunc):
  
    pass  # have an Update button so probably don't need these
 
  def destroy(self):
  
    self.curateNotifiers(self.unregisterNotify)

    BasePopup.destroy(self) 
  
  def updateAfter(self, obj=None):
    
    if self.refresh:
      return

    self.refresh = True
    self.after_idle(self.update)
    
  def update(self):

    self.updateProject()
    self.updateSpectra()
    self.updatePeakLists()
    self.updateChains()
    self.updateRestraintLists()
    self.updateStructures()

  def updateProject(self):

    project = self.project

    self.projectNameLabel.set('Project name = %s' % project.name)

    textMatrix = []
    objectList = []

    number = 0
    for repository in project.sortedRepositories():
      number += 1
      name = repository.name
      path = repository.url.path

      datum = [number, name, path]
      objectList.append(repository)
      textMatrix.append(datum)

    self.projectTable.update(textMatrix, objectList)


  def updateSpectra(self):

    textMatrix = []
    objectList = []

    nmrProject = self.nmrProject
    number = 0
    for experiment in nmrProject.sortedExperiments():
      for spectrum in experiment.sortedDataSources():
        number += 1
        name = '%s:%s' % (experiment.name, spectrum.name)
        numDim = spectrum.numDim
        shiftList = experiment.shiftList
        if shiftList:
          shiftListText = '%s [%d]' % (shiftList.name or '<No name>', shiftList.serial)
        else:
          shiftListText = ''
        path = getDataSourceFileName(spectrum)

        datum = [number, name, numDim, shiftListText, path]
        objectList.append(spectrum)
        textMatrix.append(datum)

    self.spectrumTable.update(textMatrix, objectList)

  def updatePeakLists(self):

    textMatrix = []
    objectList = []

    nmrProject = self.nmrProject
    number = 0
    for experiment in nmrProject.sortedExperiments():
      for spectrum in experiment.sortedDataSources():
        for peakList in spectrum.sortedPeakLists():
          number += 1
          name = '%s:%s:%s' % (experiment.name, spectrum.name, peakList.serial)
          peaks = peakList.peaks
          numPeaks = len(peaks)
          numPartAssigned = len([peak for peak in peaks if isPeakAssigned(peak, fully=False)])
          percentPartAssigned = int(round((100.0 * numPartAssigned) / max(numPeaks, 1)))
          numAllAssigned = len([peak for peak in peaks if isPeakAssigned(peak, fully=True)])
          percentAllAssigned = int(round((100.0 * numAllAssigned) / max(numPeaks, 1)))

          datum = [number, name, numPeaks, numPartAssigned, percentPartAssigned, numAllAssigned, percentAllAssigned]
          objectList.append(peakList)
          textMatrix.append(datum)

    self.peakListTable.update(textMatrix, objectList)

  def updateChains(self):

    textMatrix = []
    objectList = []

    project = self.project
    nmrProject = self.nmrProject
    number = 0
    for molSystem in project.sortedMolSystems():
      for chain in molSystem.sortedChains():
        number += 1
        name = '%s:%s' % (molSystem.code, chain.code)
        residues = chain.residues
        numResidues = len(residues)
        numAssignableAtoms = 0
        for residue in residues:
          numAssignableAtoms += len([atom for atom in residue.atoms if not atom.chemAtom.waterExchangeable])
        for shiftList in nmrProject.sortedMeasurementLists():
          if shiftList.className != 'ShiftList':
            continue
          shiftListText = '%s [%d]' % (shiftList.name or '<No name>', shiftList.serial)
          numAssignedAtoms = 0
          for residue in residues:
            numAssignedAtoms += len([atom for atom in residue.atoms if atom.atomSet and getAtomSetShifts(atom.atomSet, shiftList)])
          percentAssignedAtoms = int(round((100.0*numAssignedAtoms)/max(numAssignableAtoms, 1)))

          datum = [number, name, numResidues, numAssignableAtoms, shiftListText, numAssignedAtoms, percentAssignedAtoms]
          objectList.append(chain)
          textMatrix.append(datum)

    self.chainTable.update(textMatrix, objectList)

  def updateRestraintLists(self):

    textMatrix = []
    objectList = []

    colHeadings = ['#','Name','Type','Restraints','Experiments']

    project = self.project
    number = 0
    for nmrConstraintStore in project.sortedNmrConstraintStores():
      for constraintList in nmrConstraintStore.sortedConstraintLists():
        number += 1
        name = '%s:%s' % (nmrConstraintStore.serial, constraintList.name or '')
        typ = constraintList.className[:-14]
        numRestraints = len(constraintList.constraints)
        experiments = ', '.join([e.name for e in constraintList.experiments])
        datum = [number, name, typ, numRestraints, experiments]
        objectList.append(constraintList)
        textMatrix.append(datum)

    self.restraintListTable.update(textMatrix, objectList)
        
  def updateStructures(self):

    textMatrix = []
    objectList = []

    project = self.project
    number = 0
    for structureEnsemble in project.sortedStructureEnsembles():
      number += 1
      name = structureEnsemble.ensembleId
      numModels = len(structureEnsemble.models)
      chainCodes = []
      numResidues = 0
      for chain in structureEnsemble.sortedCoordChains():
        chainCodes.append(chain.code)
        numResidues += len(chain.residues)

      chainCodes = ','.join(['%s' % chainCode for chainCode in chainCodes])
      chainCodes = '%s:%s' % (structureEnsemble.molSystem.code, chainCodes)
      datum = [number, name, numModels, chainCodes, numResidues]
      objectList.append(structureEnsemble)
      textMatrix.append(datum)

    self.structureTable.update(textMatrix, objectList)

  def showBrowser(self):

    project = self.project
    projectName = project.name
    url = TABLE_BASE_URL
    webBrowser = WebBrowser(self)

    headerString = '''<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<title>%s</title>
<meta http-equiv="Content-type" content="text/html;charset=UTF-8" />
<meta name="description" content="CCPN project summary" />
<meta name="keywords" content="CCPN, Analysis, NMR" />
<link rel="stylesheet" href="%s/style.css" type="text/css" media="print, projection, screen" />
<style type="text/css">
table.tablesorter tr:nth-child(odd) td {
    background-color: #DDAAAA; color: black;
}
table.tablesorter tr:nth-child(even) td {
    background-color: #AAAADD; color: black;
}
</style>

<script type="text/javascript" src="%s/jquery-latest.js"></script>
<script type="text/javascript" src="%s/jquery.tablesorter.js"></script>
<script type="text/javascript">$(document).ready(function() {        
        $("table").tablesorter();   
}); </script>
</head>
<body>
''' % (projectName, url, url, url)

    htmlStrings = [headerString]
    htmlStrings.append('<h1>Project name = %s</h1>' % projectName)

    for name, table in ( \
        ('Project directories', self.projectTable),
        ('Spectra', self.spectrumTable),
        ('Peak Lists', self.peakListTable),
        ('Chains', self.chainTable),
        ('Restraint Lists', self.restraintListTable),
        ('Structure Ensembles', self.structureTable),
    ):
      htmlStrings.append('<h1>%s</h1>\n' % name)
      htmlStrings.append('<table class=tablesorter>\n')
      htmlStrings.append('<thead>\n')
      htmlStrings.append('<tr>\n')
      for colHeading in table.colHeadings:
        colHeading = colHeading.replace('\n', ' ')
        htmlStrings.append('<th>%s</th>' % colHeading)
      htmlStrings.append('\n</tr>\n')
      htmlStrings.append('</thead>\n')

      htmlStrings.append('<tbody>\n')
      for datum in table.rowTexts:
        htmlStrings.append('<tr>\n')
        for data in datum:
          htmlStrings.append('<td>%s</td>' % data)
        htmlStrings.append('\n</tr>\n')
      htmlStrings.append('</tbody>\n')

      htmlStrings.append('</table>\n')
    htmlStrings.append('</body>\n')
    htmlStrings.append('</html>\n')

    htmlString = ''.join(htmlStrings)

    webBrowser.openHtml(htmlString)

