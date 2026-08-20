
"""
======================COPYRIGHT/LICENSE START==========================

SecStructureGraph.py: Part of the CcpNmr Analysis program

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

from ccpnmr.analysis.core.AssignmentBasic import getAtomSetShifts, getShiftLists
from ccpnmr.analysis.core.ConstraintBasic import getMeanPeakIntensity
from ccpnmr.analysis.popups.BasePopup       import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import getThroughSpacePeakLists
from ccpnmr.analysis.core.MoleculeBasic import getRandomCoilShift, getLinkedResidue
from ccpnmr.analysis.core.CouplingBasic import getCouplingAtomResonances

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FloatEntry import FloatEntry
from memops.gui.FontList import FontList
from memops.gui.Frame import Frame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning
from memops.gui.MultiWidget import MultiWidget 
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledCanvas import ScrolledCanvas
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccp.general.Constants import code1LetterToCcpCodeDict as codeDict

FIXED_HEADINGS = ['Residue','Sec\nStruc',
                    u'd_\u03B1N\ni,i+1',
                    'd_NN\ni,i+1',u'd_\u03B2N\ni,i+1',
                    u'd_\u03B1N\ni,i+3',
                    u'd_\u03B1\u03B2\ni,i+3',
                    u'd_\u03B1N\ni,i+4','d_NN\ni,i+2',
                    u'd_\u03B1N\ni,i+2']

FIXED_TIP_TEXTS = ['The sequence number and residue type',
                   'The stated secondary structure type for the residue',
                   'Indicates residue that show an alpha(i) to amide(i+1) through-space/NOE connection',
                   'Indicates residue that show an amide(i) to amide(i+1) through-space/NOE connection',
                   'Indicates residue that show a beta(i) to amide(i+1) through-space/NOE connection',
                   'Indicates residue that show an alpha(i) to amide(i+3) through-space/NOE connection',
                   'Indicates residue that show an alpha(i) to beta(i+3) through-space/NOE connection',
                   'Indicates residue that show an alpha(i) to amide(i+4) through-space/NOE connectio',
                   'Indicates residue that show an amide(i) to amide(i+2) through-space/NOE connection',
                   'Indicates residue that show an alpha(i) to amide(i+2) through-space/NOE connection']

OPTIONAL_COLUMNS  = [u'\u0394\n\u03B413C\u03B1',u'\u0394\n\u03B413C\u03B2',
                      u"\u0394\n\u03B413C'",u'\u0394\n\u03B41H\u03B1',
                      'CSI',u'3JH-H\u03B1']

OPTIONAL_TIP_TEXTS = {u'\u0394\n\u03B413C\u03B1': 'Carbon alpha secondary chemical shift (difference to random coil value)',
                      u'\u0394\n\u03B413C\u03B2': 'Carbon beta secondary chemical shift (difference to random coil value)',
                      u"\u0394\n\u03B413C'": 'Backbone carbonyl carbon secondary chemical shift (difference to random coil value)',
                      u'\u0394\n\u03B41H\u03B1': 'Hydrogen alpha secondary chemical shift (difference to random coil value)',
                      'CSI': 'Chemical shift index value, using 13C and 1H shifts (Wishart DS, Sykes BD)',
                      u'3JH-H\u03B1': 'The 3-bond amide hydrogen to alpha hydrogen coupling'}

SEC_STRUCT_COLORS = {'H': '#8080F0',
                   'B': '#F0B080',
                   'E': '#F08080',
                   'G': '#80B0F0',
                   'I': '#B080F0',
                   'T': '#80F080',
                   'S': '#80B080',
                   'C': '#F0F080',
                   'U': '#F0F0F0',
                   }

SEC_STRUCT_NAMES = {'H': 'Alpha',
                  'B': 'Bridge',
                  'E': 'Beta',
                  'G': '3-10',
                  'I': 'Pi',
                  'T': 'Turn',
                  'S': 'Bend',
                  'C': 'Coil',
                  'U': 'Undef',
                  }

CSI_REFS = {'Ala':(52.5,177.1,19.0,4.35),
            'Cys':(58.3,174.8,28.6,4.65),
            'Cyss':(58.0,175.1,41.8,4.76),
            'Asp':(54.1,177.2,40.8,4.76),
            'Glu':(56.7,176.1,29.7,4.29),
            'Phe':(57.9,175.8,39.3,4.66),
            'Gly':(45.0,173.6,None,3.97),
            'His':(55.8,175.1,32.0,4.63),
            'Ile':(62.6,176.9,37.5,3.95),
            'Lys':(56.7,176.5,32.3,4.36),
            'Leu':(55.7,177.1,41.9,4.17),
            'Met':(56.6,175.8,32.8,4.52),
            'Asn':(53.6,175.1,39.0,4.75),
            'Pro':(62.9,176.0,31.7,4.44),
            'Gln':(56.2,176.3,30.1,4.37),
            'Arg':(56.3,176.5,30.3,4.38),
            'Ser':(58.3,173.7,62.7,4.50),
            'Thr':(63.1,175.2,68.1,4.35),
            'Val':(63.0,177.1,31.7,3.95),
            'Trp':(57.8,175.8,28.3,4.70),
            'Tyr':(58.6,175.7,38.7,4.60)}
 
CSI_ERRORS = {'Pro':(4.0,4.0,4.0,0.1),
              'Gly':(0.7,0.5,None,0.1),
              None :(0.7,0.5,0.7,0.1)}

DISULPHIDE_CYS_DICT = { 'H'  :8.54,
                      'HA' :4.76,
                      'HB2':3.29,
                      'HB3':3.02,
                      'C'  :175.5,
                      'CA' :55.63,
                      'CB' :41.22,
                      'N'  :118.7,}
 
class SecStructureGraphPopup(BasePopup):
  """
  **Display NMR Evidence for Protein Secondary Structure**
  
  This popup window is designed to collect together several kinds of NMR derived
  information to give a visual representation of evidence to indicate the
  regions of a protein chain with secondary structure, and what that structure
  might be. This system collects data from several sources from through-space
  (e.g. NOE) connectivity; from secondary chemical shifts; from a combined
  chemical shift index; from 3J H-Ha scalar couplings.

  The main chart lists the residue sequence of the selected chain (see
  "Options"), the predicted secondary structure of the residues and various
  types of evidence of the secondary structure. The "Sec Struc" column may be
  set manually based on the contents of the table via the [Set ...] buttons at
  the bottom, or automatically via a program like DANGLE_. In the table the
  eight through-space connectivity columns (d_aN i,i+1 to da_N i,i+2) are always
  displayed, but the other columns relating to chemical shift are optional and
  may be toggled via the top buttons. The connectivity columns are filled with a
  green colour if a given kind of connectivity is observed for a residue,
  additionally for the first three of these there is an extra distinction of
  strong (S) medium (M) and weak (W) categories. The name of a column indicates
  what kind of connectivity is present, for example d_aN i,i+1 represents the
  connection of an alpha position of one amino acid residue to the amide
  position of a residue one position further along (C-terminal) in the chain.

  The Delta-delta columns indicate secondary chemical shifts for Ca, Cb, C' and
  Ha atoms; the difference from recorded chemical shift to (sequence adjusted)
  random coil chemical shift. The random coil chemical shifts are calculated as
  specified in Schwarzinger et  al. as detailed below. These secondary chemical
  shifts are combined to give the CSI values (calculated according to the methods
  in the papers listed below) indicating the secondary structure type; -1
  indicates alpha-helix and +1 beta-strand. The last "3JHNHa" column lists any
  recorded scalar couplings between the amide and alpha hydrogens of a residue.
  These are indicative of the phi backbone angle, according to the Karplus
  relationship. These coupling values may be extracted and analysed in more
  detail using the Data Analysis : `3J H-Ha Coupling`_ system.

  The "Chart" tab shows a graphical representation of the secondary structure
  evidence data. This is the same data that is shown in the first table, it is
  merely presented differently, on a compact graphical form that may be printed
  out by making a PostScript file (right mouse menu). The order of the rows of
  data in the chart is the same order as the columns of the evidence table. The
  grey and white panels of the chart delineate the display into regions of ten
  residues. The text font used in the chart and the number or residues that
  are shown in a row, before the chart wraps back to start at the left again,
  are governed by settings in the last "Options Tab".

  The "Options" tab controls various aspects of how the data in the other tabs
  is displayed. The "Chain" naturally selects the sequence to use and "Shift
  List" is the source of the secondary chemical shift information. The "NOE
  Intensity Classes" section allows the user to set two values which distinguish
  between strong, medium and weak connection strength. Note that these values
  represent a peak intensity relative to the average in the peak list containing
  the connection. The lower table lists the peak lists that are used as sources
  of through-space connectivity information. The use can toggle to consider or
  reject particular peak lists by double clicking in the "Consider?" column. Any
  changes to any of these settings will be reflected one the user selects one
  of the first two tabs.

  **Caveats & Tips**

  The most accurate secondary structure prediction currently available in 
  CCPN comes from the embedded DANGLE_ program. Although DANGLE has a separate
  graphical interface its secondary structure predictions, if committed, will be
  visible in this system via the displayed secondary structure classification. 

  Regions of a polypeptide backbone with stable secondary structure will tend to
  show more near-sequence ("short range") though-space connectivities than
  unstructured, flexible  regions; commonly termini and loops. Naturally
  alpha-helices will tend to show i to i+3 and i to i+4 connections, due to the
  periodicity of the helix. Beta strand regions will often show i to i+2
  connectivity by lack i to i+3 and i to i+4, given that stands are extended
  conformations.

  **References**

  Secondary chemical shifts:

  *Schwarzinger, S., Kroon, G. J. A., Foss, T. R., Chung, J., Wright, P. E., Dyson, H. J.
  "Sequence-Dependent Correlation of Random Coil NMR Chemical Shifts",
  J. Am. Chem. Soc. 123, 2970-2978 (2001)*

  Chemical Shift Index

  *Wishart DS, Sykes BD.
  The 13C chemical-shift index: a simple method for the identification of protein
  secondary structure using 13C chemical-shift data.
  J Biomol NMR. 1994 Mar;4(2):171-80.*

  *Wishart DS, Sykes BD, Richards FM.
  The chemical shift index: a fast and simple method for the assignment of protein
  secondary structure through NMR spectroscopy.
  Biochemistry. 1992 Feb 18;31(6):1647-51.*
  
  .. _DANGLE: DanglePopup.html
  .. _`3J H-Ha Coupling`: CalcHnHaCouplingPopup.html
  """

  def __init__(self, parent, *args, **kw):

    self.waiting   = 0
    self.chain     = None
    self.shiftList = None
    self.peakLists = []
    self.graphData = []
    self.data      = []
    self.randomCoilShifts = {}
     
    BasePopup.__init__(self, parent=parent, title='Chart : Secondary Structure Chart')

  def body(self, guiFrame):

    self.geometry('700x500')
   
    guiFrame.expandGrid(0,0)
    
    tipTexts = ['A table of the evidence used to estimate the secondary structure type of residues',
                'The secondary structure evidence information as a graphical chart, which may then be saved as a PostScript file',
                'The various parameters that control the data display, including which peak lists and molecule chain to use']
    options = ['Residue Data','Chart','Options']
    tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(0,0),
                              callback=self.selectTab, tipTexts=tipTexts)
    self.tabbedFrame = tabbedFrame
    frameA, frameC, frameB = tabbedFrame.frames

    #
    # Residue Table
    #
    
    frameA.expandGrid(1,0)

    frame = Frame(frameA, grid=(0,0), sticky='ew')
    frame.grid_columnconfigure(1, weight=1)    
     
    label = Label(frame, text='Column\nOptions:', grid=(0,0))
    
    tipText = 'Selects which of various secondary chemical shifts, chemical shift index or 3J H-Ha coupling data to display in the table and chart'
    self.optSelector=PartitionedSelector(frame, self.toggleOpt, maxRowObjects=10,
                                         colors=['#C0C0FF']*6, grid=(0,1),
                                         sticky='ew', tipText=tipText)

    self.optSelector.update(objects=OPTIONAL_COLUMNS,
                            labels=OPTIONAL_COLUMNS,
                            selected=OPTIONAL_COLUMNS[:-2])
    
    headingList, tipTexts = self.getHeadingList()
    self.residueMatrix = ScrolledMatrix(frameA, headingList=headingList,
                                        callback=self.selectResidue,
                                        multiSelect=True, grid=(1,0),
                                        tipTexts=tipTexts)
    
    ssCodes = ['H','E','C','U'] # 'B','G','I','T','S',
    texts   = ['Set %s' % SEC_STRUCT_NAMES[code] for code in ssCodes]
    tip = 'Set the secondary structure type for the selected residue to "%s"'
    tipTexts = [tip % SEC_STRUCT_NAMES[code] for code in ssCodes]
    commands = [lambda code=c: self.setSecStructure(code) for c in ssCodes]
    
    tipTexts += ['Read CSI information for the residue from a CSI output file',
                 'Display the secondary structure evidence information as a graphical chart, which may then be saved as a PostScript file']
    texts += ['Read CSI\nFile','Make\nChart']
    commands += [self.readCsi,self.makePlot]
    buttonList = ButtonList(frameA, texts=texts, gridSpan=(1,2),
                            commands=commands, grid=(2,0), tipTexts=tipTexts)

    #
    # Options
    #
    
    frameB.expandGrid(4,0)

    row  = 0
    frame = Frame(frameB, grid=(row,0))
    frame.grid_columnconfigure(8, weight=1)    
     
    label = Label(frame, text='Chain:', grid=(0,0))
    tipText = 'Selects which molecular chain to secondary structure evidence for'
    self.chainPulldown = PulldownList(frame, callback=self.changeChain,
                                       grid=(0,1), tipText=tipText)

    label = Label(frame, text='Shift List:', grid=(0,2))
    tipText = 'Selects which shift list chemical shifts are derived from, for calculating deltas etc.'
    self.shiftListPulldown = PulldownList(frame, callback=self.changeShiftList,
                                          grid=(0,3), tipText=tipText)


    label = Label(frame, text='Chart Font:', grid=(0,4))
    tipText = 'Sets the typeface to use in the main chart display'
    self.fontPulldown = FontList(frame, callback=self.changeFont,
                                 grid=(0,5), tipText=tipText)

    label = Label(frame, text='Residue Wrap Width:', grid=(0,6))
    wraps = [30,40,50,60,70,80,90,100,120,150,200,300,1000]
    tipText = 'Sets how wide, in terms of number residues, the main chart is; will cause sequences to wrap back to the start'
    self.wrapPulldown = PulldownList(frame, callback=self.changeWrap, texts=[str(x) for x in wraps],
                                     index=2, objects=wraps, grid=(0,7), tipText=tipText)


    row += 1
    div = LabelDivider(frameB, text='NOE Intensity classes', grid=(row,0))

    row += 1
    frame = Frame(frameB, grid=(row,0))
    frame.grid_columnconfigure(4, weight=1)
    
 
    label = Label(frame, text='Strong/Medium:', grid=(0,0))
    tipText = 'The intensity value, relative to peak list average, that separates strong and medium NOE strength categories'
    self.strongEntry = FloatEntry(frame, text=1.5, width=8,
                                  grid=(0,1), tipText=tipText)

    label = Label(frame, text='Medium/Weak:', grid=(0,2))
    tipText = 'The intensity value, relative to peak list average, that separates weak and medium NOE strength categories'
    self.weakEntry = FloatEntry(frame, text=0.3, width=8,
                                grid=(0,3), tipText=tipText)
    
    text = ' (Boundaries are relative to mean peak list intensity)'
    label = Label(frame, text=text, grid=(0,4))

    row += 1
    div = LabelDivider(frameB, text='NOE peak lists', grid=(row,0))

    row += 1
    tipTexts = ['The experiment:spectrum to which the through-space peak list belongs',
                'The serial number of the through-space peak list within its spectrum',
                'Whether the through-space peak list should be used to provide connectivity data for the secondary structure evidence table']
    headingList = ['Spectrum','Peak List','Consider?']
    editWidgets      = [None,None,None]
    editGetCallbacks = [None,None,self.togglePeakList]     
    editSetCallbacks = [None,None,None]    
    self.peakListMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                         highlightType=None,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         grid=(row,0), tipTexts=tipTexts)

    #
    # Chart
    #
    
    frameC.expandGrid(0,0)

    self.chartFrame = SecStructureGraph(frameC, grid=(0,0))

    #
    # Main
    #
    

    buttonList = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url,
                                   grid=(0,0), sticky='e')

    self.administerNotifiers(self.registerNotify)
      
    self.updatePeakLists()
    self.updateShiftLists()
    self.updateChains() # will call self.updateAfter()
  
  def changeWrap(self, wrap):
  
    self.chartFrame.wrap = wrap
  
  def changeFont(self, name):
  
    self.chartFrame.font = name
  
  def administerNotifiers(self, notifyFunc):

    for func in ('__init__','delete'):
      notifyFunc(self.updateResidue,'ccp.molecule.MolSystem.Residue',func)
      notifyFunc(self.updateAssignment,'ccp.nmr.Nmr.PeakDimContrib',func)
      notifyFunc(self.updateChains,'ccp.molecule.MolSystem.Chain',func)
      notifyFunc(self.updatePeakLists,'ccp.nmr.Nmr.PeakList',func)
    for func in ('setSecStrucCode','setResidue'):
      notifyFunc(self.updateSpinSystem,'ccp.nmr.Nmr.ResonanceGroup',func)
      
    notifyFunc(self.updateCarbonShift,'ccp.nmr.Nmr.Shift','setValue')
    
    if self.parent.project:
      self.peakLists = getThroughSpacePeakLists(self.parent.project)

  def selectTab(self, i):
  
    if i == 1:
      self.chartFrame.draw(self.data, self.optSelector.state)

  def togglePeakList(self, peakList):
  
    if peakList:
      if peakList in self.peakLists:
        self.peakLists.remove(peakList)
      else:
        self.peakLists.append(peakList)
        
      self.updatePeakLists()
      self.updateAfter()

  def updatePeakLists(self, object=None):
  
    textMatrix = []
    colorMatrix = []
  
    peakLists = []
    if self.parent.project:
      peakLists = getThroughSpacePeakLists(self.parent.project)        
      
      for peakList in peakLists:
        datum = []
        colors = [None,None,'#D08080']
        
        use = 'No'
        if peakList in self.peakLists:
          use = 'Yes'
          colors[2] = '#80D080'
        
        datum.append('%s:%s' % (peakList.dataSource.experiment.name,peakList.dataSource.name))
        datum.append(peakList.serial)
        datum.append(use)
        
        textMatrix.append(datum)
        colorMatrix.append(colors)
    
    self.peakListMatrix.update(objectList=peakLists, textMatrix=textMatrix,colorMatrix=colorMatrix)  

  def selectResidue(self, object, row, col):
  
    pass

  def updateShiftLists(self, *opt):
  
    project = self.parent.project
    shiftLists = getShiftLists(self.nmrProject)
    shiftListNames = self.getShiftListNames(shiftLists)
    
    index = 0
    shiftList = None
    
    if shiftListNames:
      if self.shiftList in shiftLists:
        index = shiftLists.index(self.shiftList)
        shiftList = self.shiftList
      else:
        index = 0
	shiftList = shiftLists[0]
          
    if shiftList is not self.shiftList:
      self.shiftList = shiftList
      self.updateAfter()
    
    self.shiftListPulldown.setup(shiftListNames, shiftLists, index)

  def getShiftListNames(self,shiftLists):
    
    shiftListNames = []

    for shiftList in shiftLists:
      if not hasattr(shiftList, 'name'):
        shiftList.name = "ShiftList "+ str(shiftList.serial)
      elif not shiftList.name:
        shiftList.name = "ShiftList "+ str(shiftList.serial)
      shiftListNames.append(shiftList.name)
 
    return shiftListNames
 
  def changeShiftList(self, shiftList):
    
    if shiftList is not self.shiftList:
      self.shiftList = shiftList
      self.updateAfter()
  
  def setSecStructure(self, code):
    
    for residue, spinSystem in self.residueMatrix.currentObjects:
      spinSystem.secStrucCode = code  
        
    self.updateAfter()

  def readCsi(self):

    fileName  = None
    fileTypes = [  FileType('CSI', ['*.csi']), FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
               title = 'Import CSI file', dismiss_text = 'Cancel',
               selected_file_must_exist = True)

    fileName = fileSelectPopup.getFile() 
    
    if fileName:
      dict = {}
      file = open(fileName,'r')
      file.readline()
      
      line = file.readline()
      while(line):
        fields = line.split()
        if fields and fields[0][0] != '#':
          seqCode = fields[0]
          ccpCode = codeDict['protein'].get(fields[1], '???')
          csi     = int(fields[-2])
          dict['%s%s' % (seqCode,ccpCode)] = csi
 
        line = file.readline()

      missing = 0
      for residue in self.chain.sortedResidues():
        csi = dict.get('%d%s' % (residue.seqCode,residue.ccpCode))
        if csi is not None:
          residue.csi = csi
        else:
          missing += 1

      if missing > 0:
        showWarning('Warning','%d residues are missing data in CSI file' % missing, parent=self)

      self.updateAfter()

  def makePlot(self):
      
    if self.data:
      self.tabbedFrame.select(1)
      self.chartFrame.draw(self.data, self.optSelector.state)

  def updateCarbonShift(self, shift):

    resonance = shift.resonance
    if resonance.isotopeCode != '13C':
      return

    resonanceSet = resonance.resonanceSet
    if resonanceSet:
      atom = resonanceSet.findFirstAtomSet().findFirstAtom()
      if (atom.residue.chain is self.chain) and (atom.name[:2] in ('CA','CB')):
        self.updateAfter()

  def getChains(self):

    chains = []
    for molSystem in self.parent.project.sortedMolSystems():
      for chain in molSystem.sortedChains():
        if chain.molecule.molType in ('protein',None):
          chains.append(chain)

    return chains

  def changeChain(self, chain):
     
     if chain is not self.chain:
       self.chain = chain
       self.updateAfter()
  
  
  def get3jHHa(self, residue):
        
    couplingValues = []
    resonanceA, resonanceB = getCouplingAtomResonances(self.nmrProject, residue, ('H','HA'), makeNew=False)
    if resonanceA and resonanceB:
      for jCoupling in resonanceA.jCouplings:
        if resonanceB in jCoupling.resonances:
          couplingValues.append(jCoupling.value)
    
    if couplingValues:
      return sum(couplingValues)/float(len(couplingValues))
      
  def getCsi(self, residue):
    """
    REFERENCES FOR CSI METHODS AND DATA

    Wishart DS, Sykes BD.
    The 13C chemical-shift index: a simple method for the identification of protein
    secondary structure using 13C chemical-shift data.
    J Biomol NMR. 1994 Mar;4(2):171-80.

    Wishart DS, Sykes BD, Richards FM.
    The chemical shift index: a fast and simple method for the assignment of protein
    secondary structure through NMR spectroscopy.
    Biochemistry. 1992 Feb 18;31(6):1647-51.
    """
    
    ccpCode = residue.ccpCode
    if ccpCode == 'Gly':
      atomNames = ('CA','C','HA2','HA3')
    else:
      atomNames = ('CA','C','CB','HA')
  
    shiftValues = []
    for atomName in atomNames:
      
      atom = residue.findFirstAtom(name=atomName)
      
      resonances = set()
      if atom and atom.atomSet:
        for resonanceSet in atom.atomSet.resonanceSets:
          for resonance in resonanceSet.resonances:
            resonances.add(resonance)

      shiftValue = None
      mean = 0.0
      num  = 0.0
      for resonance in resonances:
        shift = resonance.findFirstShift(parentList=self.shiftList)
        
        if shift:
          mean += shift.value
          num += 1.0
          
      if num:
        shiftValue = mean/num
        
      shiftValues.append(shiftValue)      

    if ccpCode == 'Gly':
      if shiftValues[2]:
        if shiftValues[3]:
          shiftValues[3] = 0.5 * (shiftValues[2]+shiftValues[3])
          shiftValues[2] = None
        else:
          shiftValues[3] = shiftValues[2]
          shiftValues[2] = None
        
      elif shiftValues[3]:
        shiftValues[2] = None
    
    elif ccpCode == 'Cys':
      if shiftValues[3] is not None:
        oxDelta = abs(41.22 - shiftValues[3])
        rdDelta = abs(28.34 - shiftValues[3])
      
        if oxDelta < rdDelta:
          ccpCode = 'Cyss'
  

    csis = [None] * 4
    sign = (-1,-1,1,1)
    refs = CSI_REFS.get(ccpCode)
    errs = CSI_ERRORS.get(ccpCode, CSI_ERRORS.get(None))
        
    csi = None
    if refs and errs:
      for i in range(4):
        shiftValue = shiftValues[i]
 
        if shiftValue is not None:
          if shiftValue > (refs[i]+errs[i]):
            csis[i] = sign[i]
          elif shiftValue < (refs[i]-errs[i]):
            csis[i] = -1  * sign[i]
          else:
            csis[i] = 0
 
      plus  = csis.count(1)
      zero  = csis.count(0)
      minus = csis.count(-1)
      none  = csis.count(None)
 
      if none < 4:
        if (plus > zero) and (plus > minus):
          csi = 1
        elif (minus > zero) and (minus > plus):
          csi = -1
        elif (zero > minus) and (zero > plus):
          csi = 0
        elif (plus == zero) and (zero > minus):
          csi = 1
        elif (minus == zero) and (zero > plus):
          csi = -1
        else:
          csi = 0
  
    return csi
  
  def getRandomCoilShift(self, residue, atomName):

    ccpCode = residue.ccpCode    
    if (ccpCode == 'Gly') and (atomName == 'HA'):
      atomName = 'HA2'
    
    if self.randomCoilShifts.get(residue):
      # could be GLY CB and hence None
      return self.randomCoilShifts[residue].get(atomName)
    
    self.randomCoilShifts[residue] = {}
            
    prev1 = getLinkedResidue(residue,'prev')
    next1 = getLinkedResidue(residue,'next')
    
    if next1:
      next2 = getLinkedResidue(next1,'next')
    else:
      next2 = None
      
    if prev1:
      prev2 = getLinkedResidue(prev1,'prev')
    else:
      prev2 = None
    
    context = [prev2,prev1,residue,next1,next2]
        
    for atom in residue.atoms:
      name = atom.name
      
      if name in ('CA','CB','C','HA','HA2','HA3'):
        chemAtom = atom.chemAtom
        value = getRandomCoilShift(chemAtom, context=context, sourceName='BMRB')
        
        self.randomCoilShifts[residue][name] = value
           
    if self.randomCoilShifts.get(residue):
      return self.randomCoilShifts[residue].get(atomName)
    
    return 


  def updateChains(self, *obj):

    chains = self.getChains()
    
    if len(self.parent.project.molSystems) > 1:
      names  = ['%s:%s' % (c.molSystem.code,c.code) for c in chains]
    else:
      names = [c.code for c in chains]

    index = 0
    chain = None
    if chains:
      if self.chain not in chains:
        chain = chains[0]
        index = 0
      else:
        chain = self.chain
        index = chains.index(chain)

    self.chain = chain
    self.chainPulldown.setup(names, chains, index) 

  def updateSpinSystem(self, spinSystem):
  
    if spinSystem.residue and (spinSystem.residue.chain is self.chain):
      self.updateAfter()

  def updateResidue(self, residue):
    
    if residue.chain is self.chain:
      self.updateAfter()
    else:
      self.updateChains()


  def updateAssignment(self, peakDimContrib):

    if peakDimContrib.resonance.isotopeCode != '1H':
      return

    if peakDimContrib.peakDim.peak.peakList not in self.peakLists:
      return

    resonanceSet = peakDimContrib.resonance.resonanceSet
    if resonanceSet:
      residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
      if residue.chain is self.chain:
        for peakDim in peakDimContrib.peakDim.peak.peakDims:
          if peakDim is peakDimContrib.peakDim:
            continue
          
          for peakDimContrib1 in peakDim.peakDimContribs:
            if peakDimContrib1.resonance.isotopeCode != '1H':
              continue            

            resonanceSet1 = peakDimContrib1.resonance.resonanceSet
            if resonanceSet1:
              residue1 = resonanceSet1.findFirstAtomSet().findFirstAtom().residue
              if residue1.chain is self.chain:
                delta = abs(residue.seqCode - residue1.seqCode)
                if delta < 5:
                  self.updateAfter()
                  return

  def toggleOpt(self, *obj):
  
    self.updateAfter()

  def updateAfter(self, *obj):

    if self.waiting:
      return

    self.waiting = True
    self.after_idle(self.update)

  def queryConnection(self, residue1, atomType1, atomType2, delta):
   
    from ccpnmr.analysis.core.AssignmentBasic import makeResonanceGuiName
    if delta:
      residue2 = residue1.chain.findFirstResidue(seqCode = residue1.seqCode+delta)
    else:
      residue2 = residue1

    if residue2:
      resonances1 = self.getResonancesByType(residue1,atomType1)
      resonances2 = self.getResonancesByType(residue2,atomType2)
 
      if delta != 1:
        for resonance in resonances1:
          for peakDimContrib in resonance.peakDimContribs:
            peak = peakDimContrib.peakDim.peak
            if peak.peakList not in self.peakLists:
              continue
 
            for peakDim in peak.peakDims:
              if peakDim is not peakDimContrib.peakDim:
                for peakDimContrib2 in peakDim.peakDimContribs:
                  if peakDimContrib2.resonance in resonances2:
                    return ' ', '#80C080'

      else:
        text  = ''
        color = None
        maxIntens = (0.0, None)
        for resonance in resonances1:
          for peakDimContrib in resonance.peakDimContribs:
            peak = peakDimContrib.peakDim.peak
            if peak.peakList not in self.peakLists:
              continue
            
            height = 0.0
            intensity = peak.findFirstPeakIntensity(intensityType='height')
            if intensity:
              height = abs(intensity.value)
 
            for peakDim in peak.peakDims:
              if peakDim is not peakDimContrib.peakDim:
                for peakDimContrib2 in peakDim.peakDimContribs:
                  if peakDimContrib2.resonance in resonances2:
                    if height > maxIntens[0]:
                      maxIntens = (height, peak.peakList)
        
        if maxIntens[1]:
          height =  maxIntens[0]/maxIntens[1].meanIntensity
          if height > self.strongEntry.get():
            return 'S', '#80F080'
          elif  height >  self.weakEntry.get():
            return 'M', '#80C080'
          else:
            return 'W', '#A0B0A0'
        
        else:
          return '', None
     
    return '', None
  
  def getResonancesByType(self, residue, atomType):
    
    atoms = []
    if (residue.ccpCode in ('Pro','Hyp')) and atomType in ('H',):
      for atom in residue.atoms:
        if atom.name[:2] == 'HD':
          atoms.append(atom)
      
    elif atomType in ('H',):
      for atom in residue.atoms:
        if atom.name == atomType:
          atoms.append(atom)
 
    else:
      for atom in residue.atoms:
        if atom.name[:2] == atomType:
          atoms.append(atom)
    
    dict = {}
    for atom in atoms:
      if atom.atomSet:
        for resonanceSet in atom.atomSet.resonanceSets:
          for resonance in resonanceSet.resonances:
            dict[resonance] = 1
    
    return dict.keys()
    
  def getAtomShiftDiff(self, residue, atomType):

    delta = None
    randomCoilShift = self.getRandomCoilShift(residue, atomType)

    if self.shiftList and (randomCoilShift is not None):

      atom = residue.findFirstAtom(name=atomType)
      if atom and atom.atomSet:
        shifts = getAtomSetShifts(atom.atomSet, self.shiftList)
        
        if shifts:
          mean = 0.0
          for shift in shifts:
            mean += shift.value
          mean /= float(len(shifts))
 
          delta = mean - randomCoilShift

    return delta

  def getHeadingList(self):
  
    selection = self.optSelector.getSelected()
  
    headingList = FIXED_HEADINGS + selection
    
    tipTexts = FIXED_TIP_TEXTS + [OPTIONAL_TIP_TEXTS[x] for x in selection]
    
    return headingList, tipTexts

  def update(self):

    for peakList in self.peakLists:
      peakList.meanIntensity = getMeanPeakIntensity(peakList.peaks, intensityType='height')

    doCa, doCb, doCO, doHa, doCSI, do3J = self.optSelector.state

    spinSystems = {}
    for spinSystem in self.nmrProject.resonanceGroups:
      residue = spinSystem.residue
      
      if residue:
        spinSystems[residue] = spinSystem

    queryConn = self.queryConnection
    getAtomShiftDiff = self.getAtomShiftDiff
    getCsi = self.getCsi
    get3jHHa = self.get3jHHa

    objectList = []
    textMatrix = []
    colorMatrix = []
    self.data = []
    if self.chain:
      for residue in self.chain.sortedResidues():
        if residue.molResidue.molType != 'protein':
          return
          
        spinSystem = spinSystems.get(residue)
        
        if not spinSystem:
          spinSystem = self.nmrProject.newResonanceGroup(ccpCode=residue.ccpCode)
          spinSystem.residue = residue
        
        secStructure = spinSystem.secStrucCode or 'U'
          
        #if not hasattr(residue, 'csi'):
        residue.csi = getCsi(residue)     

        datum  = []
        colors = [None] * 14
        colors[1] = SEC_STRUCT_COLORS[secStructure]
        
        d_aN0, colors[2] = queryConn(residue,'HA', 'H',1)
        d_NN0, colors[3] = queryConn(residue, 'H', 'H',1)
        d_bN0, colors[4] = queryConn(residue,'HB', 'H',1)
        d_aN3, colors[5] = queryConn(residue,'HA', 'H',3)
        d_ab3, colors[6] = queryConn(residue,'HA','HB',3)
        d_aN4, colors[7] = queryConn(residue,'HA', 'H',4)
        d_NN2, colors[8] = queryConn(residue, 'H', 'H',2)
        d_aN2, colors[9] = queryConn(residue,'HA', 'H',2)

        d_ca = getAtomShiftDiff(residue,'CA')
        d_cb = getAtomShiftDiff(residue,'CB')
        d_co = getAtomShiftDiff(residue,'C')
        d_ha = getAtomShiftDiff(residue,'HA')
        c3j = get3jHHa(residue)

        datum.append('%d %s' % (residue.seqCode, residue.ccpCode)) 
        datum.append(SEC_STRUCT_NAMES[secStructure])
        datum.append(d_aN0)
        datum.append(d_NN0)
        datum.append(d_bN0)
        datum.append(d_aN3)
        datum.append(d_ab3)
        datum.append(d_aN4)
        datum.append(d_NN2)
        datum.append(d_aN2)
        
        if doCa:
          datum.append(d_ca)
        if doCb:
          datum.append(d_cb)
        if doCO:
          datum.append(d_co)
        if doHa:
          datum.append(d_ha)
        if doCSI:
          datum.append(residue.csi)
        if do3J:  
          datum.append(c3j)

        colorMatrix.append(colors)
        objectList.append((residue,spinSystem))
        textMatrix.append(datum)
        self.data.append( [residue,spinSystem,
                           d_aN0,d_NN0,d_bN0,
                           d_aN3,d_ab3,d_aN4,
                           d_NN2,d_aN2,d_ca,
                           d_cb,d_co,d_ha,
                           residue.csi, c3j] )

    headingList, tipTexts = self.getHeadingList()
    
    self.residueMatrix.update(objectList=objectList,
                               headingList=headingList,
                               textMatrix=textMatrix,
                               colorMatrix=colorMatrix,
                               tipTexts=tipTexts)
    self.waiting = False 

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  def open(self):
    
    self.updatePeakLists()
    self.updateShiftLists()
    self.updateChains() # will call self.updateAfter()
    BasePopup.open(self)
  
class SecStructureGraph(Frame):

  def __init__(self, parent, *args, **kw):

    self.font = 'Helvetica 10'
    self.wrap = 50

    Frame.__init__(self, parent=parent, *args, **kw)
      
    self.grid_rowconfigure(0, weight=1)    
    self.grid_columnconfigure(0, weight=1)    
    self.canvasFrame = ScrolledCanvas(self, relief='groove',
                                      borderwidth=2, resizeCallback=None)
    self.canvasFrame.grid(row=0, column=0, sticky='nsew', padx=1, pady=1 )
    self.canvas = self.canvasFrame.canvas
    #self.canvas.bind('<Button-1>', self.toggleResidue)
    
  def draw(self, data, options):
  
    doCa, doCb, doCO, doHa, doCSI, do3J = options
    
    font = self.font
    name, size = font.split()[:2]
    smallFont = '%s %d' % (name, max(int(size)-2, 6))

    # Some Mac Tk implementations use obscure font unless explicitly set fontmap
    self.canvas.tk.call('set', 'fontmap(%s)' % smallFont, smallFont)
    
    if 'italic' not in font:
      font2 = font + ' italic'
    else:
      font2 = font 
      
    spans  = [3,3,4,2,2]
    texts1 = [u'\u03B1N','NN',u'\u03B2N',u'\u03B1N',u'\u03B1\u03B2',u'\u03B1N','NN',u'\u03B1N']
    optLineY = []
    for datum in data:
      if datum[-1] is not None:
        doCsi = True
        break
        
    c = self.canvas
    ctext = c.create_text
    crect = c.create_rectangle
    item = ctext(0,0,text='M',font=font2)
    bbox = c.bbox(item)
    d  = bbox[2]-bbox[0] # X
    th = bbox[3]-bbox[1]
    d2 = th+3 # Y
    ss = th/5
    c.delete('all')
    xmax = 0
    ymax = 0
    deferList = []
    
    resStart = 0
    while resStart < len(data):
      resEnd = min(len(data),resStart+self.wrap)
      dataRange = data[resStart:resEnd]
      resStart += self.wrap
    
      x = 0.0
      y = ymax
 
      i = 0
      y += th
      y += d2
      for text1 in texts1:
        x = 0.0
        y += d2
        item = ctext(x,y,text='d',anchor='sw',font=font2)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]

        x += sdx
        t = text1[0]
        item = ctext(x,y+ss,text=t,font=smallFont,anchor='sw')
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
 
        x += sdx
        t = text1[1]
        item = ctext(x,y+ss,text=t,font=smallFont,anchor='sw')
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
 
        x += sdx + 2
        if i < 3:
          t = '(d'
          item = ctext(x,y,text=t,anchor='sw',font=font2)
          bbox = c.bbox(item)
          sdx = bbox[2]-bbox[0]
          x += sdx
 
          if i == 0:
            item = ctext(x,y+ss,text=u'\u03B1\u03B4',font=smallFont,anchor='sw')
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
            item = ctext(x,y,text=')',anchor='sw',font=font2)
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
 
          elif i == 1:
            item = ctext(x,y+ss,text='N',font=smallFont,anchor='sw')
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
            item = ctext(x,y+ss,text=u'\u03B4',font=smallFont,anchor='sw')
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
            item = ctext(x,y,text=',d',anchor='sw',font=font2)
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
            item = ctext(x,y+ss,text=u'\u03B4',font=smallFont,anchor='sw')
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
            item = ctext(x,y+ss,text='N',font=smallFont,anchor='sw')
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
            item = ctext(x,y,text=')',anchor='sw',font=font2)
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
 
          else:
            item = ctext(x,y+ss,text=u'\u03B2\u03B4',font=smallFont,anchor='sw')
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
            item = ctext(x,y,text=')',anchor='sw',font=font2)
            bbox = c.bbox(item)
            sdx = bbox[2]-bbox[0]
            x += sdx
 
        else:
          t = '(i,i+%d)' % (spans[i-3])
          ctext(x,y,text=t,anchor='sw',font=font2)
 
        i += 1
        xmax = max(x, xmax)

      y = ymax + 9*d2 + d
      yEnd = y
 
      if doCa:
        y += d2
        y += d2
        x = 0.0
        item = ctext(x,y,text=u'\u0394\u03B4(',anchor='sw',font=font2)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y-ss,text='13',anchor='sw',font=smallFont)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text='C',anchor='sw',font=font)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y-ss,text=u'\u03B1',anchor='sw',font=smallFont)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text=')',anchor='sw',font=font2)
        optLineY.append(y)
        yEnd = y

      if doCb:
        y += d2
        y += d2
        x = 0.0
        item = ctext(x,y,text=u'\u0394\u03B4(',anchor='sw',font=font2)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y-ss,text='13',anchor='sw',font=smallFont)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text='C',anchor='sw',font=font)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y-ss,text=u'\u03B2',anchor='sw',font=smallFont)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text=')',anchor='sw',font=font2)
        optLineY.append(y)
        yEnd = y

      if doCO:
        y += d2
        y += d2
        x = 0.0
        item = ctext(x,y,text=u'\u0394\u03B4(',anchor='sw',font=font2)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y-ss,text='13',anchor='sw',font=smallFont)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text="C'",anchor='sw',font=font)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text=')',anchor='sw',font=font2)
        optLineY.append(y)
        yEnd = y

      if doHa:
        y += d2
        y += d2
        x = 0.0
        item = ctext(x,y,text=u'\u0394\u03B4(',anchor='sw',font=font2)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y-4,text='1',anchor='sw',font=smallFont)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text='H',anchor='sw',font=font)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y-4,text=u'\u03B1',anchor='sw',font=smallFont)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text=')',anchor='sw',font=font2)
        optLineY.append(y)
        yEnd = y

      if do3J:
        y += d2
        y += d2
        x = 0.0
        item = ctext(x,y-4,text='3',anchor='sw',font=smallFont)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text='J(',anchor='sw',font=font2)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text=u'H-H',anchor='sw',font=font)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y-4,text=u'\u03B1',anchor='sw',font=smallFont)
        bbox = c.bbox(item)
        sdx = bbox[2]-bbox[0]
        x += sdx
        item = ctext(x,y,text=')',anchor='sw',font=font2)
        optLineY.append(y)
        yEnd = y
 
      if doCSI:
        y += d2
        y += d2
        x = 0.0
        ctext(x,y,text='CSI',anchor='sw',font=font2)
        optLineY.append(y)
        yEnd = y
 
      y = ymax
      x = xmax+d
      n = int( len(dataRange)/10.0 )
      for i in range(n+1):
        if i%2 == 0:
          s = x +(i*10*d)
          e = x +( min((i*10)+10,len(dataRange)) *d)
          crect(s,y,e,yEnd+d2+d2,outline='#EEEEEE', fill='#EEEEEE')
 
      x2 = xmax+d + float(len(dataRange)*d)
      yf = ymax
      y = ymax
      yf3 = 2*d2
      x = xmax
      
      y = ymax + (d2*5)
      for i, w, h in deferList:
        y2 = (d2*i) + y + h
        crect(x+d,y2,x+d+w,y2+max(2, d/6), outline='black', fill='black')
      
      deferList = []
    
    
      levels = [0,0,0,0,0]
      last   = [0,0,0,0,0]
      spinSystems = []
      for residue, spinSystem, d_aN0, d_NN0, d_bN0, d_aN3, d_ab3, d_aN4, \
          d_NN2, d_aN2, d_ca, d_cb, d_co, d_ha, csi, c3j in dataRange:
          
        spinSystems.append(spinSystem)
        y = ymax
        x += d
 
        if residue.seqCode and residue.seqCode % 10 == 0:
          ctext(x,y,text='%d' % residue.seqCode, font=font,anchor='sw')
 
        y += d2
        text = residue.molResidue.chemCompVar.chemComp.code1Letter or 'X'
        
        if not spinSystem:
          fill = '#A0A0A0'
        elif not spinSystem.resonances:
          fill = '#A0A0A0'
        else:
          fill = '#000000'
         
        ctext(x+(d/2),y,text=text,font=font,anchor='s',fill=fill)

 
        for val in (d_aN0,d_NN0,d_bN0):
          y += d2
          if val:
            if val == 'W':
              crect(x,y+(0.8*d),x+d,y+d, outline='black', fill='black')
            elif val == 'M':
              crect(x,y+(0.45*d),x+d,y+d, outline='black', fill='black')
            else: # S
              crect(x,y+(0.1*d),x+d,y+d, outline='black', fill='black')
 
        for i, val in enumerate((d_aN3,d_ab3,d_aN4,d_NN2,d_aN2)):
          y += d2
 
          if val:
            span = spans[i]
            if levels[i] > 4:
              levels[i] = 0
            else:
              if (residue.seqCode - last[i]) > span:
                levels[i] = 0
            y2 = y + d - (max(2, d/6) * (levels[i]))
 
            levels[i] += 1
            last[i] = residue.seqCode
            
            endx = x+(d*span)
            if endx > x2:
              deferList.append((i, endx-x2, d - (max(2, d/6) * (levels[i]))))
              crect(x,y2,x2,y2+max(2, d/6), outline='black', fill='black')
            else:
              crect(x,y2,endx,y2+max(2, d/6), outline='black', fill='black')

        y += d
 
        if doCa:
          y += d2 + d2
          if d_ca is not None:
            crect(x,y,x+d,y-(d2*0.15*d_ca), outline='black', fill='black',width=0.3)

        if doCb:
          y += d2 + d2
          if d_cb is not None:
            crect(x,y,x+d,y-(d2*0.15*d_cb), outline='black', fill='black',width=0.3)

        if doCO:
          y += d2 + d2
          if d_co is not None:
            crect(x,y,x+d,y-(d2*0.15*d_co), outline='black', fill='black',width=0.3)

        if doHa:
          y += d2 + d2
          if d_ha is not None:
            crect(x,y,x+d,y-(d2*d_ha), outline='black', fill='black',width=0.3)
        
        if do3J:
          y += d2 + d2
          if c3j is not None:
            crect(x,y,x+d,y-(0.1*d2*c3j), outline='#000000', fill='#404040',width=0.3)

        if doCSI:
          y += d2 + d2
          if csi is not None:
            crect(x,y,x+d,y-(0.5*d2*csi), outline='#808080', fill='#808080',width=0.3)
 
        y += d2
        y += d2
        yf3 = y

      x  = xmax+d
 
      for yOpt in optLineY:
       c.create_line(x,yOpt,x2,yOpt,fill='#000000',width=0.3)
 
      #c.create_line(x,yf3+40,x2,yf3+40,fill='#FFFFFF',width=0.3)
      crect(-d, ymax-d-d, x2+d, yf3+40,outline='#FFFFFF',width=0.3)
    
      self.plotSecondaryStructure(spinSystems, c, sp=d, xBase=x, yBase=yf3)
      ymax = yf3+80
    
    
 
  def plotSecondaryStructure(self, spinSystems, canvas, sp=6.6, xBase=100, yBase=100):

    X   = float(xBase)
    Y   = float(yBase)
    ps  = ''
    end = []
    i   = 0
    N   = len(spinSystems)
    ssCodes = [ss.secStrucCode for ss in spinSystems]
 
    while i < N:
      if ssCodes[i] in ('H','G','I'):
        X = xBase + (sp*float(i))
        c = 0
        while i < N and ssCodes[i] in ('H','G','I'):
          i += 1
          c += 1
 
        end = self.plotAlphaHelix(canvas,X,Y,c,sp)
 
      elif ssCodes[i] in ('E','B'):
        X = xBase + (sp*float(i))
        if (i>1) and ssCodes[i-1] in ('H','G','I'):
          self.plotCoil(canvas,X-sp,Y,1,sp)
 
        c = 0
        while i < N and ssCodes[i] in ('E','B'):
          i += 1
          c += 1
 
        self.plotBetaSheet(canvas,X,Y,c,sp)
        X = xBase + (sp*i)
        if (i < N-1) and ( ssCodes[i+1] in ('H','G','I') ):
          self.plotCoil(canvas,X,Y,1,sp)
 
        if end:
          canvas.lift(end)
 
      elif ( ssCodes[i] in ('C','U','T','S',None)):
        X = xBase + (sp*float(i))
        c = 0
        if i> 1 and ssCodes[i-1] in ('H','G','I'):
          c += 1
          X -= sp
 
        while i < N and ssCodes[i] in ('C','U','T','S',None):
          i += 1
          c += 1
 
        if (i < N-1) and ssCodes[i+1] in ('H','G','I') :
          c += 1
 
        self.plotCoil(canvas,X,Y,c,sp)

        if end:
          canvas.lift(end)
 

  def plotCoil(self,canvas,X,Y,N,w):

    h   = 0.4 * w
    Y  += h/2.0
    
    item = canvas.create_rectangle(X,Y,X+(N*w),Y+h,outline='#000000',fill='#606060',width=0.3)
    return item


  def plotBetaSheet(self,canvas,X,Y,c,width):

    c = float(c)
    X2     = X+width*c
    height = 0.8 * width
    indent = 0.8 * width
    flare  = 1.2 * height * 0.5
 
    Y  += (height / 2.0)
    Y2  = Y + (height / 2.0)
    Y1  = Y - (height / 2.0)
 
    coords = [(X2,Y),(X2-indent,Y1-flare),(X2-indent,Y1),(X,Y1),(X,Y2),(X2-indent,Y2),(X2-indent,Y2+flare)]
    item = canvas.create_polygon(coords,outline='#000000',fill='#606060',width=0.3)
 
    return item


  def plotAlphaHelix(self,canvas,startX,startY,c,width):

    width  = float(width)
    N      = max(1,int(c/1.8))
    X      = 0.0
    height = 0.8 * width
 
    minor_width = 1.6 * width
    total_width = c * width
 
    width2  = (total_width - minor_width) /float(N)
    startY += height/2
    startX -= width2/2
    s  = 10
    n  = 1
 
    if (N%2 == 1):
      X    = startX + ( N * width2)
      item = self.plotCurve(canvas, X,startY,height,width2,minor_width,s,1,-1,'#B0B0B0')

    while (n < N):
      X    = startX + (n * width2)
      item = self.plotCurve(canvas, X,startY,height,width2,minor_width,s,0,-1,'#B0B0B0')
      n    += 2
 
    n = 2
    while (n < N):
      X    = startX + (n * width2)
      item = self.plotCurve(canvas,X,startY,height,width2,minor_width,s,0,1,'#606060')
      n    += 2
 
    X    = startX + (0 * width2)
    item = self.plotCurve(canvas,X,startY,height,width2,minor_width,s,2,1,'#606060')
 
    end    = None
    if (N%2 == 0):
      X    = startX + ( N * width2)
      end = self.plotCurve(canvas,X,startY,height,width2,minor_width,s,1,1,'#606060')
 
    return end

  def plotCurve(self,canvas,startX,startY,height,width,minor_width,s,half,neg,color):

    from math import cos
    Pi     = 3.14159265359
    startI = 0
    x      = 0.0
    endI   = s
    width  = float(width)
    s      = float(s)
    items  = []
 
    if half == 1:
      endI  = int(s/2)
 
    elif half == 2:
      startI = int(s/2)
      x = int(s/2) * width/s
 
    if neg < 0:
      neg = -1
    else:
      neg = 1

    if neg > 0:
      outline = '#000000'
    else:
      outline = color  
 
    coords = []
    for i in range(startI,endI+1):
      y  = height*cos( x*Pi/width )
      coords.append( [x+startX, startY-(neg*y)] )
      x += width/s
 
    for i in range(startI,endI+1):
      x -= width/s
      y  = height*cos( x*Pi/width )
      coords.append( [x+startX+minor_width, startY-(neg*y)] )

    item = canvas.create_polygon(coords,fill=color,outline=outline,width=0.3)
 
    return item



