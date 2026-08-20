"""
============================COPYRIGHT START=============================

DangleFrame.py: Part of the DANGLE package (release v1.1)

DANGLE: Dihedral ANgles from Global Likelihood Estimates 
Copyright (C) 2008 Nicole Cheung, Tim Stevens, Bill Broadhurst (University of Cambridge)

========================================================================


If you make use of the software or any documents in this package, 
please give credit by citing this package, its authors and references
in the literature with the same authors.

We would appreciate hearing of any problems you may encounter, 
but the programs, the documents and any files created by the programs 
are provided WITHOUT ANY WARRANTY and without even the implied warranty of
CORRECTNESS, MERCHANTABILITY or FITNESS FOR A PARTICULAR OR GENERAL USE.

THE RESPONSIBILITY FOR ANY ADVERSE CONSEQUENCES FROM THE USE OF PROGRAMS OR
DOCUMENTS OR ANY FILE OR FILES CREATED BY USE OF THE PROGRAMS OR DOCUMENTS
LIES SOLELY WITH THE USERS OF THE PROGRAMS OR DOCUMENTS OR FILE OR FILES AND
NOT WITH AUTHORS OF THE PROGRAMS OR DOCUMENTS.

You are not permitted to use any pieces of DANGLE in other programs, make
modifications to DANGLE, or make what a lawyer would call a "derived work" 
in any other way without the consent from either author.
 

============================COPYRIGHT END===============================

for further information, please contact the authors:

- Nicole Cheung   : msc51@cam.ac.uk

- Tim Stevens     : tjs23@cam.ac.uk

- Bill Broadhurst : r.w.broadhurst@bioc.cam.ac.uk

========================================================================

If you are using this software for academic purposes, we suggest
quoting the following reference:

===========================REFERENCE START==============================

DANGLE: To be completed.

CCPN:
Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END================================
"""

import Tkinter, os

from memops.gui.Label               import Label
from memops.gui.Frame               import Frame
from memops.gui.LabelFrame          import LabelFrame
from memops.gui.Entry               import Entry
from memops.gui.FloatEntry          import FloatEntry
from memops.gui.Button              import Button
from memops.gui.ButtonList          import ButtonList
from memops.gui.ProgressBar         import ProgressBar
from memops.gui.ScrolledGraph       import ScrolledGraph
from memops.gui.ScrolledMatrix      import ScrolledMatrix
from memops.gui.PulldownList        import PulldownList
from memops.gui.FileSelectPopup     import FileSelectPopup
from memops.gui.MessageReporter     import showError, showYesNo, showInfo, showOkCancel, showWarning
from memops.gui.DataEntry           import askString

from memops.editor.Util             import createDismissHelpButtonList

from ccpnmr.analysis.popups.BasePopup       import BasePopup
from ccpnmr.analysis.core.AssignmentBasic   import getAtomSetShifts, assignAtomsToRes, getShiftLists
from ccpnmr.analysis.core.MoleculeBasic     import getResidueMapping
from ccpnmr.analysis.core.ConstraintBasic   import getFixedResonance
from ccp.lib.StructureLib                   import getResiduePhiPsi

from ccp.lib.MoleculeQuery     import getLinkedResidue


#
from ccpnmr.analysis.core.ChemicalShiftBasic import getChemAtomNmrRef
#

from ccp.gui.ViewRamachandranFrame   import ViewRamachandranFrame

import cambridge.dangle.src.dangle as dangleModule

from cambridge.dangle.src.dangle import Dangle
#from cambridge.dangle.DangleGui import DangleGui

#OUTDIR      = './DanglePred/'
#TEMP_IN     = 'dangle_cs.inp'

whitespace = '\t\n\x0b\x0c\r '

GOOD_COLORS = [None]*10
BAD_COLORS  = ['#FFB0B0']*10
INACTIVE_COLORS = ['#808080']*10
ENSEMBLE_COLOR = '#33FF99'

DEFAULT_MAX_ISLANDS = 2

EDIT_ATTRS = ['phiValue','psiValue',
              'phiUpper','phiLower',
              'psiUpper','psiLower',]

BACKBONE_ATOMS = ('C','CA','CB','HA','HA2','HA3','H','N')

def testMacro(argServer):

  project = argServer.getProject()
  popup = DangleGui(argServer.parent)
  popup.open()
  
  
  
class DangleFrame(Frame):
  
  def __init__(self, parent, dangleGui, project=None, *args, **kw):

    self.guiParent     = parent
    self.dangleGui     = dangleGui
    self.dangleDir     = None
    self.dangleChain   = None
    self.dangleResidue = None
    #self.outDir      = OUTDIR
    self.row           = None
    self.col           = None
    self.project       = project
    self.nmrProject    = None
    self.colorScheme   = 'red'
    
    self.chain         = None
    self.shiftList     = None
    self.dangleStore   = False # Not None
    self.constraintSet = None
    self.ensemble      = None
    
    Frame.__init__(self, parent=parent)

    self.grid_columnconfigure(0, weight=1)
    self.grid_rowconfigure(1, weight=1)
    
    row = 0
    
    # TOP LEFT FRAME
    
    frame = LabelFrame(self, text='Options')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.columnconfigure(5, weight=1)
        
    label = Label(frame, text='Chain')
    label.grid(row=0, column=0, sticky='w')
    self.chainPulldown = PulldownList(frame, callback=self.changeChain,
                                      tipText='Choose the molecular system chain to make predictions for')
    self.chainPulldown.grid(row=0, column=1, sticky='w')
    
    label = Label(frame, text='Shift List')
    label.grid(row=0, column=2, sticky='w')
    self.shiftListPulldown = PulldownList(frame, callback=self.changeShiftList,
                                          tipText='Select the shift list to take input chemical shifts from')
    self.shiftListPulldown.grid(row=0, column=3, sticky='w')
    
    label = Label(frame, text='Max No. of Islands:')
    label.grid(row=0, column=4, sticky='w')
    sizes = range(10)
    texts = [str(s) for s in sizes] + ['Do not reject',]
    self.rejectPulldown = PulldownList(frame, texts=texts, objects=sizes + [None,],
                                       tipText='Select the maximum allowed number of disontinuous prediction islands')
    self.rejectPulldown.set(DEFAULT_MAX_ISLANDS) # Actual value not index
    self.rejectPulldown.grid(row=0, column=5, sticky='w')
    
    label = Label(frame, text='Dangle Run:')
    label.grid(row=1, column=0, sticky='w')
    self.dangleStorePulldown = PulldownList(frame, callback=self.changeDangleStore,
                                            tipText='Select a run number to store DANGLE results within')
    self.dangleStorePulldown.grid(row=1, column=1, sticky='w')
    
    label = Label(frame, text='Ensemble:')
    label.grid(row=1, column=2, sticky='w')
    self.ensemblePulldown = PulldownList(frame, callback=self.changeEnsemble,
                                         tipText='Select the structure ensemble for superimposition of angle values on the GLE plots')
    self.ensemblePulldown.grid(row=1, column=3, sticky='w')
    
    label = Label(frame, text='Restraint Set:')
    label.grid(row=1, column=4, sticky='w')
    self.constrSetPulldown = PulldownList(frame, callback=self.changeConstraintSet,
                                          tipText='Select the CCPN restraint set to store DANGLE dihedral angle restraints in')
    self.constrSetPulldown.grid(row=1, column=5, sticky='w')
    
    # TOP RIGHT FRAME

    outerFrame = Frame(self)
    outerFrame.grid(row=row, column=1, rowspan=2, sticky='nsew')
    outerFrame.rowconfigure(0, weight=1)
    outerFrame.columnconfigure(0, weight=1)

    frame = LabelFrame(outerFrame, text='Global Likelihood Estimates')
    frame.grid(row=0, column=0, sticky='nsew')
    frame.rowconfigure(1, weight=1)
    frame.columnconfigure(2, weight=1)
    
    self.prevPlot = ViewRamachandranFrame(frame, relief='sunken',defaultPlot=False, width=180,height=180,
                                          bgColor=self.cget('bg'), nullColor='#000000', titleText='Previous',
                                          xTicks=False, yTicks=False, xLabel='', yLabel='', showCoords=False)
    self.prevPlot.grid(row=0, column=0, sticky='nsew')
    self.prevPlot.getPlotColor = self.getPlotColor
    
    
    self.nextPlot = ViewRamachandranFrame(frame, relief='sunken',defaultPlot=False, width=180,height=180,
                                          bgColor=self.cget('bg'), nullColor='#000000', titleText='Next',
                                          xTicks=False, yTicks=False, xLabel='', yLabel='', showCoords=False)
    self.nextPlot.grid(row=0, column=1, sticky='nsew')
    self.nextPlot.getPlotColor = self.getPlotColor
    
    
    self.plot = ViewRamachandranFrame(frame, relief='sunken',defaultPlot=False, width=360,height=360,
                                      bgColor=self.cget('bg'), nullColor='#000000')
    self.plot.grid(row=1, column=0, columnspan=2, sticky='nsew')
    self.plot.selectColor  = '#FFB0B0'
    self.plot.getPlotColor = self.getPlotColor
    
    # BOTTOM RIGHT FRAME
    
    frame = Frame(outerFrame)
    frame.grid(row=1, column=0, sticky='nsew')
    frame.rowconfigure(0, weight=1)
    frame.columnconfigure(0, weight=1)

    
    texts = ('Previous', '  Next  ')
    commands = (self.showPrevious, self.showNext)
    tipTexts = ['Show GLE plot of angle predictions for previous residue in chain',
                'Show GLE plot of angle predictions for next residue in chain']
    buttonList = ButtonList(frame, texts, commands, tipTexts=tipTexts)
    buttonList.grid(row=0, column=0, sticky='nsew')
    
    row += 1
    
    # BOTTOM LEFT FRAME
    
    frame = LabelFrame(self, text='Dihedral Angle Predictions')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    
    self.floatEntry = FloatEntry(self,text='', returnCallback=self.setFloatEntry, width=10, formatPlaces=9)

    tipTexts = ['Residue number in chain',
                'Residue type code',
                'Number of high scoring discontinuous angle predictions',
                'Predicted secondary structure code',
                'Predicted phi dihedral angle (CO-N-CA-CO)',
                'Predicted psi dihedral angle (N-CA-CO-N)',
                'Upper bound of phi angle prediction',
                'Lower bound of phi angle prediction',
                'Upper bound of psi angle prediction',
                'Lower bound of phi angle prediction',
                'Chemical shifts used in prediction'] 

    headingList        = ['Res\nNum','Res\nType','No. of\nIslands',
                          'SS','Phi','Psi','Phi\nUpper','Phi\nLower',
                          'Psi\nUpper','Psi\nLower','Chemical Shifts']
                          
    editWidgets        = [None,None,None,None,self.floatEntry,
                          self.floatEntry,self.floatEntry,self.floatEntry,
                          self.floatEntry,self.floatEntry]
                          
    editGetCallbacks   = [None,None,None,None,self.getFloatEntry,
                          self.getFloatEntry,self.getFloatEntry,self.getFloatEntry,
                          self.getFloatEntry,self.getFloatEntry]
                          
    editSetCallbacks   = [None,None,None,None,self.setFloatEntry,
                          self.setFloatEntry,self.setFloatEntry,self.setFloatEntry,
                          self.setFloatEntry,self.setFloatEntry]
                          
    self.predictionMatrix = ScrolledMatrix(frame, 
                                           headingList=headingList, 
                                           multiSelect=True, 
                                           callback=self.selectCell,
                                           tipTexts=tipTexts,
                                           editWidgets=editWidgets,
                                           editGetCallbacks=editGetCallbacks,
                                           editSetCallbacks=editSetCallbacks)
    #                                       doubleCallback=self.loadGLEs)
    self.predictionMatrix.grid(row=0, column=0, sticky='nsew')
    
    row += 1

    tipTexts = ['Remove the predictions for the selected residues',
                'Run the DANGLE method to predict dihedral angles and secondary structure',
                'Delete the DANGLE results stored under the current run number',
                'Store the angle predictions and bounds in a new CCPN dihedral angle restraint list',
                'Store the secondary structure predictions in the CCPN project']
                
    texts = ['Clear\nSelected','Run Prediction!','Delete\nCurrent Run',
             'Commit\nRestraints','Commit\nSecondary Structure']
    commands = [self.clearSelected, self.runDangle, self.deleteRun, 
                self.storeDihedralConstraints, self.storeSecondaryStructure]
    self.buttonList = createDismissHelpButtonList(self, texts=texts, commands=commands, # dismiss_text='Quit',
                                                 dismiss_cmd=self.dangleGui.quit,
                                                 help_url=self.dangleGui.help_url,
                                                 expands=True, tipTexts=tipTexts)
    self.buttonList.grid(row=row, column=0, columnspan=2, sticky='ew')
    self.buttonList.buttons[1].config(bg='#C0FFFF')
    
    self.updateProject(project)
 
    self.notify(dangleGui.registerNotify)


  def destroy(self):
      
    self.notify(self.dangleGui.unregisterNotify)
    Frame.destroy(self)
 
 
  def notify(self, notifyfunc):
     
    for func in ('__init__', 'delete'):
      notifyfunc(self.updateChainPulldown, 'ccp.molecule.MolSystem.Chain', func)
      
    for func in ('__init__', 'delete', 'setName'):
      notifyfunc(self.updateShiftListPulldown, 'ccp.nmr.Nmr.ShiftList', func)
      
    for func in ('__init__', 'delete'):
      notifyfunc(self.updateConstrSetPulldown, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)
      
    for func in ('__init__', 'delete'):
      notifyfunc(self.updateEnsemblePulldown, 'ccp.molecule.MolStructure.StructureEnsemble', func)
    
  
  def updateProject(self, project):
    
    if project:
      self.project    = project
      self.nmrProject = project.currentNmrProject or self.project.newNmrProject(name=project.name)
      self.updateShiftListPulldown()
      self.updateChainPulldown()
      self.updateDangleStorePulldown()
      self.updateConstrSetPulldown()
      self.updateEnsemblePulldown()
      self.updatePredictionMatrixAfter()
   

  
  def makeDangleInput(self, filename, chain, shiftList):
      
    if (not chain) or (not shiftList):
      return
    
    residues = chain.sortedResidues()
    
    seq = ''
    for residue in residues:
      if residue.molResidue.chemComp.code1Letter:
        seq += residue.molResidue.chemComp.code1Letter
      else:
        seq += 'X'
      
    res_0  = residues[0].seqId
    
    fopen = open(filename,'w')
    fopen.write('<entry>\n')
    fopen.write('\t<res_0>%d</res_0>\n' % res_0)
    fopen.write('\t<seq_1>%s</seq_1>\n' % seq)
    fopen.write('\t<chain>%s</chain>\n' % chain.code)
    fopen.write('\t<cs_data>\n')
    fopen.write('\t\t<!-- res_id  res_name  atom_name  chemical_shift -->')
    
    numShift = 0
    
    for residue in residues:
      for atom in residue.atoms:
        atomSet = atom.getAtomSet()
        shifts = getAtomSetShifts(atomSet, shiftList=shiftList)
        if (len(shifts) == 0):
          continue

        # to average ambiguous chemical shifts ???????????????????????????
        value = 0.
        for shift in shifts:
          value += shift.value
        value = value/float(len(shifts))
        
        at_name  = atom.name
        res_name = residue.ccpCode
        res_num  = residue.seqId
        fopen.write('\n\t\t%5s\t%s\t%-4s\t%.3f' % (res_num, res_name, at_name, value))
	
	numShift += 1
    
    fopen.write('\n\t</cs_data>\n')
    fopen.write('</entry>\n')
    fopen.close()
    
    return numShift
    
  def deleteRun(self):
  
    if self.dangleStore:
      msg = 'Really delete DANGLE run "%s"?' % self.dangleStore.name
      
      if showOkCancel('Confirm', msg, parent=self):
        self.dangleStore.delete()
        self.dangleStore = None
        self.dangleChain = None
        self.updatePredictionMatrix()
        self.updateDangleStorePulldown()
   
  def runDangle(self):
    
    chain     = self.chain
    shiftList = self.shiftList
    
    if (not chain) or (not shiftList):
      showError('Cannot Run DANGLE', 'Please specify a chain and a shift list.', parent=self)
      return
      
    # check if there is a DangleChain available
    self.checkDangleStore()
    dangleStore = self.dangleStore
    
    if not dangleStore:
      return
    
    dangleChain = dangleStore.findFirstDangleChain(chain=chain)
    if dangleChain:
      data = (chain.code, dangleChain.shiftList.serial)
      msg = 'Predictions for Chain %s using Shift List %d already exist.\nReplace data?' % data

      if not showYesNo('Replace Data', msg, parent=self):
        return

      else:
        self.dangleChain = dangleChain
        dangleChain.shiftList = shiftList
        
    else:
      self.dangleChain = dangleStore.newDangleChain(chain=chain, shiftList=shiftList)
   
	
    #dangleStore.packageLocator.repositories[0].url.dataLocation = '/home/msc51/ccpn/NexusTestGI'
    #dangleStore.packageName = 'cambridge.dangle'
    
    repository = dangleStore.packageLocator.repositories[0]
    array = dangleStore.packageName.split('.')
    path = os.path.join(repository.url.dataLocation, *array)
    path = os.path.join(path, dangleStore.name, chain.code)  # Dangle_dir/dangleStoreName/chainCode
    if not os.path.exists(path):
      os.makedirs(path)
      
    self.dangleDir = path
    
    inputFile = os.path.join(self.dangleDir, 'dangle_cs.inp')
    if os.path.isfile(inputFile):
      os.unlink(inputFile)
      
    outputFile = os.path.join(self.dangleDir, 'danglePred.txt')
    if os.path.isfile(outputFile):
      os.unlink(outputFile)
      
    numShift = self.makeDangleInput(inputFile, chain, shiftList)
        
    if not os.path.isfile(inputFile):
      msg = 'No DANGLE input has been generated.\nPlease check shift lists.'
      showError('File Does Not Exist', msg, parent=self)
      return
    if numShift == 0:
      msg = 'No shift data in input file.\nPerhaps shifts are not assigned.\nContinue prediction anyway?'
      if not showYesNo('Empty DANGLE input', msg, parent=self):
        return
    
    rejectThresh = self.rejectPulldown.getObject()
    
    # Use the Reference info from the main installation
    # location must be absolute because DANGLE could be run from anywhere
    location = os.path.dirname(dangleModule.__file__)
    
    progressBar = ProgressBar(self)
    self.update_idletasks()
        
    dangle = Dangle(location, inputFile=inputFile, outputDir=self.dangleDir,
                    reject=rejectThresh, angleOnly=False, progressBar=progressBar,
                    writePgm=False)
    
    #self.dangleDir = '/home/msc51/nexus/gItest/DanglePred/'
    #outputFile =  '/home/msc51/nexus/gItest/DanglePred/danglePred.txt'
    
    
    predictions = dangle.predictor.predictions
    gleScores = dangle.predictor.gleScores
    
    self.readPredictions(predictions, gleScores)
    self.updatePredictionMatrix()
   
  def readPredictions(self, predictions, gleScores):
    
    progressBar = ProgressBar(self, text='Reading DANGLE predictions')
    progressBar.total = len(predictions) - 2  # 2 header lines
    
    residues = self.dangleChain.chain.sortedResidues()
    getDangleResidue = self.dangleChain.findFirstDangleResidue
    newDangleResidue = self.dangleChain.newDangleResidue

    for residue in residues:
      seqId = residue.seqId
      prediction = predictions.get(seqId)
      
      if prediction is None:
        continue
      
      gleMatrix = gleScores[seqId]
    
      progressBar.increment()
      
      #resNum, resName = prediction[:2];
      
      numIsland = prediction[2]
      ss = prediction[10]
      angles = [min(179.9999, a) for a in prediction[3:10]]
      
      phi, phiUpper, phiLower, psi ,psiUpper, psiLower, omega = angles

      
      # Normalise to max
      maxVal = max(gleMatrix)
      gleMatrix = [max(0, int(val/maxVal*65535))/65535.0 for val in gleMatrix]

      dangleResidue = getDangleResidue(residue=residue)
      if not dangleResidue:
        dangleResidue = newDangleResidue(phiPsiLikelihoodMatrix=gleMatrix, residue=residue)
      
      else:
        dangleResidue.phiPsiLikelihoodMatrix = gleMatrix

      dangleResidue.numIslands = numIsland  
      dangleResidue.phiValue   = phi
      dangleResidue.phiUpper   = phiUpper
      dangleResidue.phiLower   = phiLower
      dangleResidue.psiValue   = psi
      dangleResidue.psiUpper   = psiUpper
      dangleResidue.psiLower   = psiLower
      dangleResidue.omegaValue = omega
      dangleResidue.secStrucCode = ss
     
    progressBar.destroy()
   
      
  
  def readPredictionFile(self, filename, chain):
    
    try:
      fopen = open(filename,'r')
    except:
      showError('File Reading Error', 'DANGLE prediction file %s cannot be open.' % filename, parent=self)
      return
    
    lines = fopen.readlines()
    progressBar = ProgressBar(self, text='Reading DANGLE predictions')
    progressBar.total = len(lines) - 2  # 2 header lines
    lines = lines[2:]
    
    for line in lines:
      progressBar.increment()
      if (line == '\n'):
        continue
      array     = line.split()   # keep everything as string
      resNum    = int(array[0])
      resName   = array[1]
      numIsland = int(array[2])
      phi       = array[3]
      phiUpper  = array[4]
      phiLower  = array[5]
      psi       = array[6]
      psiUpper  = array[7]
      psiLower  = array[8]
      omega     = array[9]
      ss        = array[10]
      
      if (phi == 'None'):
        phi = None
      else:
        phi = float(phi)
      if (psi == 'None'):
        psi = None
      else:
        psi = float(psi)
      if (omega == 'None'):
        omega = None
      else:
        omega = float(omega)
	if omega == 180:
	  omega = 179.9
      if (phiUpper == 'None'):
        phiUpper = None
      else:
        phiUpper = float(phiUpper)
      if (phiLower == 'None'):
        phiLower = None
      else:
        phiLower = float(phiLower)
      if (psiUpper == 'None'):
        psiUpper = None
      else:
        psiUpper = float(psiUpper)
      if (psiLower == 'None'):
        psiLower = None
      else:
        psiLower = float(psiLower)
      if (ss == 'None'):
        ss = None
  
      
      path = os.path.join(self.dangleDir, 'Res_%d.pgm' % resNum)
      gleMatrix = self.readGLE(path)
      residue = chain.findFirstResidue(seqId=int(resNum))
      
      dangleResidue = self.dangleChain.findFirstDangleResidue(residue=residue)
      if not dangleResidue:
        dangleResidue = self.dangleChain.newDangleResidue(phiPsiLikelihoodMatrix=gleMatrix, residue=residue)
      else:
        dangleResidue.phiPsiLikelihoodMatrix = gleMatrix
	
      dangleResidue.numIslands = numIsland  
      dangleResidue.phiValue   = phi
      dangleResidue.phiUpper   = phiUpper
      dangleResidue.phiLower   = phiLower
      dangleResidue.psiValue   = psi
      dangleResidue.psiUpper   = psiUpper
      dangleResidue.psiLower   = psiLower
      dangleResidue.omegaValue = omega
      dangleResidue.secStrucCode = ss
      
      # Delete temp pgm files to save space once data is in CCPN 
      os.unlink(path)
    
    fopen.close()
    progressBar.destroy()
    
    
  def readGLE(self, gleFile):
    
    if not os.path.isfile(gleFile):
      msg = 'No scorogram Res_%d.pgm\nin directory %s.' % (resNum, self.dangleDir)
      showError('File Reading Error', msg, parent=self)
      return None
      
    fopen = open(gleFile, 'r')
    lines = fopen.readlines()
    dims  = lines[2].split()
    lines = lines[4:]
    fopen.close()

    # only read the top left corner of a 10X10 square bin
    # all readings in the same bin are identical 
    binSize = 10
    matrix = []
    for j in range(36):
      x = j * binSize * binSize * 36
      for i in range(36):
        y = i * binSize
        v = int(lines[x+y].strip())  
	matrix.append(v)
    
    maxVal = float(max(matrix))
    for i in range(len(matrix)):
      matrix[i] = matrix[i] / maxVal
    
    return matrix
    
    
  def getPhiPsiPredictions(self):
    
    #if self.dangleChain:
    # dResidues = self.dangleChain.dangleResidues

    dResidues = self.predictionMatrix.objectList
    
    phiData = []
    psiData = []
    for dResidue in dResidues:
      resNum = dResidue.residue.seqCode
      phi    = dResidue.phiValue
      psi    = dResidue.psiValue
      phiData.append((resNum, phi))
      psiData.append((resNum, psi))
    
    return (phiData, psiData)
    
    
  def clearSelected(self):
    
    for dangleResidue in self.predictionMatrix.currentObjects:
      dangleResidue.numIslands = None
      dangleResidue.phiValue = None
      dangleResidue.psiValue = None
      dangleResidue.omegaValue = None
      dangleResidue.phiUpper = None
      dangleResidue.phiLower = None
      dangleResidue.psiUpper = None
      dangleResidue.psiLower = None
      dangleResidue.secStrucCode = None
    
    self.updatePredictionMatrixAfter()
  
  def storeSecondaryStructure(self):
  
    if not self.dangleChain:
      return
    
    getSpinSystem = self.nmrProject.findFirstResonanceGroup
    newSpinSystem = self.nmrProject.newResonanceGroup
    
    n = 0
    for dangleResidue in self.dangleChain.dangleResidues:
      ssCode = dangleResidue.secStrucCode
      
      if not ssCode:
        continue
      
      residue = dangleResidue.residue
      if not residue:
        continue
      
      spinSystem = getSpinSystem(residue=residue)
      
      if not spinSystem:
        spinSystem = newSpinSystem(residue=residue, ccpCode=residue.ccpCode)
      
      spinSystem.secStrucCode = ssCode
      n += 1
      
    showInfo('Info', 'Stored secondary structure types for %d residues.' % n, parent=self)
   
  def storeDihedralConstraints(self):
    
    if not self.dangleChain:
      return

    # make a new dihedralConstraintList
    head = self.constraintSet
    
    if not head:
      head = self.project.newNmrConstraintStore(nmrProject=self.nmrProject)
      self.constraintSet = head
    
    chain = self.dangleChain.chain
    shiftList = self.dangleChain.shiftList
    name = 'DANGLE Chain %s:%s ShiftList %d' % (chain.molSystem.code, chain.code, shiftList.serial)
    constraintList = head.newDihedralConstraintList(name=name, measureListSerials=[shiftList.serial,])

    # traverse the sequence and make appropriate constraint objects
    residues = chain.sortedResidues()
    for residue in residues:
      # Ensure we have atomSets etc
      getResidueMapping(residue)
    
    residueList = [(dr.residue.seqCode,dr.residue,dr) for dr in  self.dangleChain.dangleResidues]
    residueList.sort()
    
    cnt = 0
    
    for seqCode, residue, dangleResidue in residueList:
      
      phi = dangleResidue.phiValue
      psi = dangleResidue.psiValue
      
      if (phi is None) and (psi is None):
        continue

      # Use below functions because residues may not be sequentially numbered 
      prevRes = getLinkedResidue(residue,'prev')
      nextRes = getLinkedResidue(residue,'next')
      if (prevRes is None) or (nextRes is None):
        continue

      C__1  = prevRes.findFirstAtom(name='C')     # C (i-1)
      N_0   = residue.findFirstAtom(name='N')     # N (i)
      CA_0  = residue.findFirstAtom(name='CA')    # CA(i)
      C_0   = residue.findFirstAtom(name='C')     # N (i)
      N_1   = nextRes.findFirstAtom(name='N')     # C (i+1)

      # get fixedResonances
      fixedResonances = []
      for atom in (C__1, N_0, CA_0, C_0, N_1):
        atomSet = atom.atomSet
        
        if atomSet.resonanceSets:
          resonance = atomSet.findFirstResonanceSet().findFirstResonance()
          
        else:
          # make new resonance
          if not atom.chemAtom:
            print 'no chem atom'
            
          ic = atom.chemAtom.elementSymbol
          if (ic == 'C'):
            ic = '13' + ic
          elif (ic == 'N'):
            ic = '15' + ic
            
          resonance = self.nmrProject.newResonance(isotopeCode=ic)
          assignAtomsToRes([atomSet,], resonance)
        fixedResonances.append(getFixedResonance(head, resonance))
      
      # make dihedralConstraints
      phiResonances = (fixedResonances[0], fixedResonances[1],
                       fixedResonances[2], fixedResonances[3])
      phiConstraint = constraintList.newDihedralConstraint(resonances=phiResonances)

      psiResonances = (fixedResonances[1], fixedResonances[2],
                       fixedResonances[3], fixedResonances[4])
      psiConstraint = constraintList.newDihedralConstraint(resonances=psiResonances)

      # make constraint items
      
      
      if phi is not None:
        phiConstraint.newDihedralConstraintItem(targetValue=phi,
                                                upperLimit =dangleResidue.phiUpper,
                                                lowerLimit =dangleResidue.phiLower)
        cnt += 1
	
      if psi is not None:
        psiConstraint.newDihedralConstraintItem(targetValue=psi,
                                                upperLimit =dangleResidue.psiUpper,
                                                lowerLimit =dangleResidue.psiLower)
        cnt += 1
	
    showInfo('Success', 'DANGLE has generated %d dihedral restraints.' % cnt, parent=self)
      

      
    
  def loadGLEs(self, dRes, row, col):
  
    residue = dRes.residue
    title = '%d %s' % (residue.seqCode, residue.ccpCode)
    
    self.fillGlePlot(self.plot, dRes.phiPsiLikelihoodMatrix, title)
    
    prevDangleRes = self.getDangleResidue(dRes, 'prev')
    if prevDangleRes:
      self.fillGlePlot(self.prevPlot, prevDangleRes.phiPsiLikelihoodMatrix)
    else:
      self.fillGlePlot(self.prevPlot, [0]*1296)    # blank
          
    nextDangleRes = self.getDangleResidue(dRes, 'next')
    if nextDangleRes:
      self.fillGlePlot(self.nextPlot, nextDangleRes.phiPsiLikelihoodMatrix)
    else:
      self.fillGlePlot(self.nextPlot, [0]*1296)    # blank
    
    self.updatePhiPsi(dRes.residue)
        

  def fillGlePlot(self, plot, gleMatrix, title=None):	
    
    scaleCol = plot.scaleColorQuick
    
    if self.colorScheme == 'black':
      plot.nullColor = '#000000'
    else:
      plot.nullColor = '#FFFFFF'
    
    itemconf = plot.canvas.itemconfigure
    matrix   = plot.matrix
    
    for j in range(36):
      for i in range(36):			       
    	v = gleMatrix[j*36+i]			       
    	#if (v < 0.005):				       
    	#  color = plot.nullColor		       
    	#else:					       
    	color = self.getPlotColor(v)		       
    	item = matrix[i][j]
        		       
    	if plot.binWidth < 7:			       
    	  itemconf(item, fill=color, outline=color) 	       
    	elif plot.binWidth < 12:		       
    	  itemconf(item, fill=color, outline=scaleCol(color,0.9))  
    	else:					       
    	  itemconf(item, fill=color, outline=scaleCol(color,0.8))  
    
    if title:
      itemconf(plot.title, text=title)
  
  def getDangleResidue(self, dRes, direction):
    # return a DangleResidue object located offset-residue away from dRes in sequence
    
    # Use below function to guard against non-sequentially numbered residues
    # the below function follows bonds, but uses a cache for speed
    residue = getLinkedResidue(dRes.residue, direction)
    
    if residue and self.dangleChain:
      return self.dangleChain.findFirstDangleResidue(residue=residue)
    
          
  def showPrevious(self):
    
    if not self.dangleResidue:
      return
      
    prevDangleResidue = self.getDangleResidue(self.dangleResidue, 'prev')
    if not prevDangleResidue:
      return
    
    self.predictionMatrix.selectObject(prevDangleResidue)
    #self.dangleResidue = prevDangleResidue
    #self.loadGLEs(self.dangleResidue, None, None)
    #self.predictionMatrix.currentObject = self.dangleResidue
    #self.predictionMatrix.hilightObject(self.predictionMatrix.currentObject)
      
      
  def showNext(self):
  
    if not self.dangleResidue:
      return
      
    nextDangleResidue = self.getDangleResidue(self.dangleResidue, 'next')
    if not nextDangleResidue:
      return
    
    self.predictionMatrix.selectObject(nextDangleResidue)
    #self.dangleResidue = nextDangleResidue
    #self.loadGLEs(self.dangleResidue, None, None)
    #self.predictionMatrix.currentObject = self.dangleResidue
    #self.predictionMatrix.hilightObject(self.predictionMatrix.currentObject)
  
  
  def updatePhiPsi(self, residue):
    
    if self.ensemble:
      phiPsiAccept = []
      plotObjects  = []
      colors       = []
      
      cChain = self.ensemble.findFirstCoordChain(code=residue.chain.code)
      
      if cChain:
        cResidue = cChain.findFirstResidue(residue=residue)
      
        if cResidue:
          for model in self.ensemble.models:
            phiPsiAccept.append(self.getPhiPsi(cResidue, model))
            plotObjects.append((cResidue, model))
            colors.append(ENSEMBLE_COLOR)
      
      if self.colorScheme == 'rainbow':     # default grey circles
        self.plot.updateObjects(phiPsiAccList=phiPsiAccept,
                                objectList=plotObjects)
      else:                                 # bright green circles
        self.plot.updateObjects(phiPsiAccList=phiPsiAccept,
                                objectList=plotObjects,
			        colors=colors)
    
    
  def getPhiPsi(self, residue, model=None):
  
    phi, psi = getResiduePhiPsi(residue, model=model)
    return (phi,psi,1)
    
  
  def getPlotColor(self, i, maxInt=255):
    
    mode = self.colorScheme
    if mode == 'rainbow':
      if (i == 0):
        return '#%02x%02x%02x' % (255,255,255)    # white bg
      elif(i > 0.75):
        red   = 1
        green = (1-i)/0.25
        blue  = 0
      elif (i > 0.5):
        red   = (i-0.5)/0.25
        green = 1
        blue  = 0
      elif (i > 0.25):
        red   = 0
        green = 1
        blue  = (0.5-i)/0.25
      else:
        red   = 0
        green = i/0.25
        blue  = 1
      return '#%02x%02x%02x' % (red*maxInt, green*maxInt, blue*maxInt)
      
    """
    elif mode == 'black':
      if i > 0.5:
        red   = i
        green = 1 - i
        blue  = 1 - i
      else:
        v = 0.1 + (0.9 * i)
        red   = v
        green = v
        blue  = v
        
    elif mode == 'white':
      if i > 0.5:
        red   = i
        green = 1 - i
        blue  = 1 - i
      else:
        v = 1.0 - (0.9 * i)
        red   = v
        green = v
        blue  = v

    return '#%02x%02x%02x' % (red*maxInt, green*maxInt, blue*maxInt)
    """
    
    # default : red to black
    
    if (i == 0):
      return '#%02x%02x%02x' % (255,255,255)    # white bg
  
    return '#%02x%02x%02x' % (((1-i)*255),0,0) 

    
  def updatePredictionMatrixAfter(self, index=None, text=None):
    
    if self.chain and self.shiftList and self.dangleStore:
      self.dangleChain = self.dangleStore.findFirstDangleChain(chain=self.chain, shiftList=self.shiftList)
    else:
      self.dangleChain = None
       
    self.after_idle(self.updatePredictionMatrix)
    
    #if showYesNo('Not Found','No data for Chain %s in Dangle Run %s.\nMake prediction for this chain?' % (self.chain.code, text), parent=self):
    #  self.runDangle()
    
    
  def updatePredictionMatrix(self):
    

    shiftList = self.shiftList
    objectList  = []
    textMatrix  = []
    colorMatrix = []

    if  self.dangleChain:
      residueList = [(dr.residue.seqCode,dr.residue,dr) for dr in  self.dangleChain.dangleResidues]
      residueList.sort()
    else:
      # Chow blank table
      residueList = []
    
    for seqCode, residue, dRes in residueList:
      objectList.append(dRes)
      
      phi = dRes.phiValue
      psi = dRes.psiValue
      ss  = dRes.secStrucCode

      atomNames = []
      for atomName in BACKBONE_ATOMS:
        atom = residue.findFirstAtom(name=atomName)
        if not atom:
          continue
      
        atomSet = atom.atomSet
 
        if atomSet:
          shifts = getAtomSetShifts(atomSet, shiftList=shiftList)
          if shifts:
            atomNames.append(atomName)
      
      atomNames.sort()
      atomNames = ' '.join(atomNames)
      
      textMatrix.append((seqCode, residue.ccpCode, dRes.numIslands,
                         ss, phi, psi,
                         dRes.phiUpper, dRes.phiLower,
                         dRes.psiUpper, dRes.psiLower,
                         atomNames))
      
      if (phi is None) and (psi is None):
        colorMatrix.append(INACTIVE_COLORS)
        
      elif dRes.numIslands >= 5:
        colorMatrix.append(BAD_COLORS)
        
      else:
        colorMatrix.append(GOOD_COLORS)
	
    self.predictionMatrix.update(textMatrix=textMatrix,
                                 objectList=objectList,
                                 colorMatrix=colorMatrix)
                                   
    
  def selectCell(self, dRes, row, col):
 
    self.dangleResidue = dRes
    self.row = row
    self.col = col
    self.loadGLEs(dRes, row, col)
    
  
  def setFloatEntry(self, event):
    
    index = self.col - 4 # index of attribute to set in the EDIT_ATTRS list
    value = self.floatEntry.get()

    if value is not None:
      setattr(self.dangleResidue, EDIT_ATTRS[index], value)
      self.updatePredictionMatrixAfter()
      
  def getFloatEntry(self, dangleResidue):

    if dangleResidue:
      index = self.col - 4 # index of attribute to set in the EDIT_ATTRS list
      self.floatEntry.set(getattr(dangleResidue, EDIT_ATTRS[index]))
    
    
  def changeChain(self, chain):
    
    if chain is not self.chain:
      self.chain = chain
      self.updateEnsemblePulldown() # Ensembles are filtered by chains molSystem
      self.updatePredictionMatrixAfter()
    
  def changeShiftList(self, shiftList):
  
    if shiftList is not self.shiftList:
      self.shiftList = shiftList
      self.updatePredictionMatrixAfter()
      
  def changeEnsemble(self, ensemble):
    
    self.ensemble = ensemble
      
  def changeConstraintSet(self, constraintSet):
  
    if constraintSet is not self.constraintSet:
      self.constraintSet = constraintSet
  
  def changeDangleStore(self, dangleStore):
  
    if self.dangleStore is not dangleStore:
      self.dangleStore = dangleStore
      
      if dangleStore:
        self.dangleChain = dangleStore.findFirstDangleChain()
        self.chain = self.dangleChain.chain
        self.shiftList = self.dangleChain.shiftList
      else:
        self.dangleChain = None
        
      self.updateChainPulldown()
      self.updateShiftListPulldown()
      self.updateEnsemblePulldown() # Ensembles are filtered by chains molSystem
      self.updatePredictionMatrixAfter()
    
  def checkDangleStore(self):

    if not self.dangleStore:
	
      N = len(self.project.dangleStores) + 1
      name = askString('Request','Dangle Run Name:','Run%d' % N,parent=self)
      if not name:
        return None
    
      for character in whitespace:
        if character in name:
          showWarning('Failure','Name cannot contain whitespace',parent=self)
          return None

      if self.project.findFirstDangleStore(name=name):
        showWarning('Failure','Name already used',parent=self)
        return None
      
      self.dangleStore = self.project.newDangleStore(name=name)
      self.updateDangleStorePulldown()
   
  def updateChainPulldown(self, obj=None):
    
    index = 0
    names = []
    chains = []
    chain = self.chain
    
    for molSystem in self.project.molSystems:
      msCode = molSystem.code
    
      for chainA in molSystem.chains:
        residues = chainA.residues
        
        if not residues:
          continue
        
        for residue in residues:
          # Must have at least one protein residue
          if residue.molType == 'protein':
            names.append('%s:%s' % (msCode, chainA.code))
            chains.append(chainA)
            break
    
    if chains:
      if chain not in chains:
        chain = chains[0]
        
      index = chains.index(chain)
        
    else:
      chain = None
    
    if chain is not self.chain:
      self.chain = chain
      self.updatePredictionMatrixAfter()
        
    self.chainPulldown.setup(names, chains, index)
    
    
  def updateShiftListPulldown(self, obj=None):
    
    index = 0
    names = []
    shiftLists = getShiftLists(self.nmrProject)
        
    if shiftLists:
      if self.shiftList not in shiftLists:
        self.shiftList = shiftLists[0]
        
      index = shiftLists.index(self.shiftList)
      names = ['%s:%d' % (sl.name,sl.serial) for sl in shiftLists]  
    
    else:
      self.shiftList = None
      
    self.shiftListPulldown.setup(names, shiftLists, index)
    
    
  def updateDangleStorePulldown(self):
  
    names = ['<New>',]
    dangleStores = [None,]
  
    for dangleStore in self.project.sortedDangleStores():
      names.append(dangleStore.name)
      dangleStores.append(dangleStore)
      
    if self.dangleStore not in dangleStores:
      self.dangleStore = dangleStores[-1]
 
    index = dangleStores.index(self.dangleStore)
 
    self.dangleStorePulldown.setup(names, dangleStores, index)
    
   
  def updateEnsemblePulldown(self, obj=None):

    index = 0
    names = ['<None>',]
    ensembles = [None,]
    
    if self.chain:
      molSystem = self.chain.molSystem
    
      for ensemble in molSystem.sortedStructureEnsembles():
        names.append('%s:%d' % (molSystem.code,ensemble.ensembleId))    
        ensembles.append(ensemble)

      if self.ensemble not in ensembles:
        self.ensemble = ensembles[0]
        
      index = ensembles.index(self.ensemble)
      
    self.ensemblePulldown.setup(names, ensembles, index)
    
    
  def updateConstrSetPulldown(self, obj=None):
       
    names = ['<New>',]
    constraintSets = [None,]
    
    # Use below later, once API speed/loading is improved
    # for constraintSet in self.nmrProject.sortedNmrConstraintStores():
   
    for constraintSet in self.project.sortedNmrConstraintStores():
      names.append('%d' % constraintSet.serial)
      constraintSets.append(constraintSet)
    
    if self.constraintSet not in constraintSets:
      self.constraintSet = constraintSets[0]
    
    index = constraintSets.index(self.constraintSet)
 
    self.constrSetPulldown.setup(names, constraintSets, index)
    
    
  
    
    
  def try1(self):
  
    if not self.project:
      return
  
    ccpCodes = ['Ala','Cys','Asp','Glu','Phe','Gly','His','Ile','Lys','Leu',
                'Met','Asn','Gln','Arg','Ser','Thr','Val','Trp','Tyr','Pro']
		
    atomNames = ['HA','CA','CB','C','N']
    
    molType = 'protein'
  
    for ccpCode in ccpCodes:
      for atomName in atomNames:
      
        chemAtomNmrRef  = getChemAtomNmrRef(self.project, atomName, ccpCode, molType)
	mean = chemAtomNmrRef.meanValue
	sd   = chemAtomNmrRef.stdDev
	
	print '%5s%5s   %.3f   %.3f' % (ccpCode, atomName, mean, sd)



