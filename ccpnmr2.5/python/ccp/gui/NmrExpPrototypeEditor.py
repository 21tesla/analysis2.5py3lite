""" Editor for NmrExpPrototype. To be usable either as stand-alone
or as part of a program that uses the CCPN DataModel (e.g. CcpNmr Analysis).

======================COPYRIGHT/LICENSE START==========================

NmrExpPrototypeEditor.py: 

Copyright (C) 2005 Rasmus Fogh (University of Cambridge)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in license/LGPL.license
 
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

- email: ccpn@bioc.cam.ac.uk

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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""
import re

from memops.general import Implementation
from memops.general import Util as genUtil

from memops.universal.Io import getTopDirectory

from memops.gui import DataEntry
from memops.gui import Util as guiUtil
from memops.gui import MessageReporter

from memops.gui.Entry import Entry
from memops.gui.IntEntry import IntEntry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.PulldownList import PulldownList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.Spacer import Spacer
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.MultiWidget import MultiWidget
from memops.gui.CheckButton import CheckButton
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.MessageReporter import showOkCancel, showWarning

from memops.editor.BasePopup import BasePopup

from ccp.util import NmrExpPrototype as NmrExpPrototypeUtil

from ccp.api.nmr import NmrExpPrototype


categoryList = list(
 NmrExpPrototype.metaPackage.getElement('ExpCategory').enumeration
)
constantTimeList = list(
 NmrExpPrototype.metaPackage.getElement('ConstantTimeType').enumeration
)
expMeasurementTypeList = list(
 NmrExpPrototype.metaPackage.getElement('ExpMeasurementType').enumeration
)
expTransferTypeList = list(
 NmrExpPrototype.metaPackage.getElement('ExpTransferType').enumeration
)
namingSystemList = list(
 NmrExpPrototype.metaPackage.getElement('NamingSystem').enumeration
)

isotopeList = [
 '1H', '13C', '15N', '19F', '31P', '111Cd', '113Cd', 
 '2H', '14N', '17O', '23Na', '77Se', '79Br'
]

atomSiteTypeList = [
 'H','C','N','P','CA', 'CO', 'CB', 'Cali','Caro','Cmet','H(2)','H(2n)','H(2n+1)'
]

atomSiteTypeData = {
'H':{'name':'H','isotopeCode':'1H'},
'C':{'name':'C','isotopeCode':'13C'},
'N':{'name':'N','isotopeCode':'15N'},
'P':{'name':'P','isotopeCode':'31P'},
'CA':{'name':'CA','isotopeCode':'13C','minShift':40.0,'maxShift':75.0},
'CO':{'name':'CO','isotopeCode':'13C','minShift':155.0,'maxShift':225.0},
'CB':{'name':'CB','isotopeCode':'13C','minShift':0.0,'maxShift':47.0},
'Cali':{'name':'Cali','isotopeCode':'13C','minShift':-10.0,'maxShift':95.0},
'Caro':{'name':'Caro','isotopeCode':'13C','minShift':95.0,'maxShift':155.0},
'Cmet':{'name':'Cmet','isotopeCode':'13C','minShift':-10.0,'maxShift':35.0},
'H(2)':{'name':'H(2)','isotopeCode':'1H','minNumber':2,'maxNumber':2},
'H(2n)':{'name':'H(2n)','isotopeCode':'1H','minNumber':0,'numberStep':2},
'H(2n+1)':{'name':'H(2n+1)','isotopeCode':'1H','minNumber':1,'numberStep':2},
}

expPrototypeLabelFormat = '%3i:  %s'

# TBD DOCUMENTATION MERGE
#
# incorporating
# .../doc/ccpnmr/analysis/ExpGraphPopup.html
# .../doc/ccpnmr/analysis/RefExperimentPopup.html
#
# into
# .../doc/ccpnmr/analysis/NmrExpPrototypePopup.html
#


def nmrExpPrototypeIndexTuple(entries,nmrExpPrototype):
  """ Get index tuple for nmrExpPrototype relative to entries
  """
  
  name = expPrototypeLabelFormat % (nmrExpPrototype.serial,nmrExpPrototype.name)
  for ii,dd1 in enumerate(entries[1:]):
    for jj,dd2 in enumerate(dd1['submenu']):
      if dd2['label'] == name:
        return (ii+1,jj)
  #
  return -1

class NmrExpPrototypePopup(BasePopup):
  """
  **View and Edit NMR Experiment Prototypes**
  
  NmrExpPrototypes describe NMR experiments: The magnetisation transfer flow,
  what is being measured, and how the nuclei on the different axes are connected.
  The NmrExpPrototype groups a series of RefExperiments that share a 
  magnetisation transfer pathway but can have different dimension and axis
  selection. 
  
  **Finding an experiment description**
  
  NmrExpPrototypes and RefExperiments names are complex, but follow a
  systematic  nomenclature; understanding what they mean and what your target
  experiment should be called is a great help. A detailed explanation can be
  found on the documentation Wiki under 'Core Concepts -> NMR Experiment
  Nomenclature', or  (somewhat out of date) at
  http://www.springerlink.com/content/h57q76k245r14558/
  
  The best approach is to zoom in on a limited set of candidates and check the
  names and the detailed descriptions. Sorting on the 'Name', 'Synonym', or 
  'Max Dim' can be a help. The 'Category', 'Synonym' and 'Keywords' columns serve
  to classify experiments, but can be incompletely filled in, or misleading. 
  The Right Mouse -> Filter command is a good way to select sets of experiments;
  you might try filtering on 'NOESY', 'TOCSY', 'CA', 'CO', 'H[N]' or 'H[C] 
  (respectively 15N and 13C one-bond, out-and-back transfer), etc.
  
  The following examples should get you started on the nomenclature:
  
  **HNCAHA**: 4D HNCAHA prototype, starting on HN and ending on HA. The prototype
           includes HNcaHA (3D), HACAnH (3D, note reversed direction), HncaHA (2D) etc.
  
  **H[N[CA[HA]]]**: 4D HNCAHA prototype with out-and-back magnetisation transfer.
  
  **H[N]_H.NOESY**: 3D HSQC-NOESY. Prototype includes NOESY-HSQC (H_H[N])
  
  **H{CA|Cca}NH**: 4D HCBCANNH
  
  **C[h(0)]**: 1D carbon, selecting quaternary carbons only.


  **Viewing experiment prototypes**
  
  You need to look at all the tabs to understand the experiment in detail. The
  Atom & Measurement tab shows which atoms are on the magnetisation transfer 
  pathway, and which shifts etc. are being observed; the Exp Graphs tab show
  the  magnetisation transfers and the flow of magnetisation in the experiment;
  and the RefExperiments tab shows the individual reference experiments and what
  is on the axes.


  **Editing experiment prototypes**

  Editing existing experiment prototypes is discouraged, as these are reference
  information and should be edited only centrally. For this reason existing
  prototypes are normally not editable (see below for instructions). Adding new
  prototypes and new RefExperiments is more likely to be useful. Ultimately it
  is the responsibility of anyone adding reference information to understand
  what he is doing, so the explanation here is focused on the simplest way to 
  proceed.

  You must always start in the Prototypes tab. If there is a closely similar
  experiment prototype, it makes sense to use it  as a starting point. It is
  simple to change atom names, parameters etc., and to remove atom sites. If
  this will do you, select your template experiment and press [Create Copy].
  Then edit the copy, starting with the columns Name, Category, Synonym, and 
  Keywords in the Prototypes tab, and continuing on the other tabs. The 
  numerical columns in the ExpPrototypes table can be left alone, as the values
  will adjust to changes in the other tabs.  Remember to rename and remove or 
  add RefExperiments to fit.

  If you need more complex operations, you are better off creating a new
  prototype. It is still helpful to use an existing one as inspiration.
  
  Start on the Prototypes tab and click [Create New]. Then set the Name, 
  Category, Synonym, and Keywords, ignoring the rest of the columns on this tab.
  Try to make sure you have got the nomenclature right, preferably using another
  prototype as a template, and that the name does not duplicate that of an
  existing prototype.
  
  Now go to the Atoms & Measurements tab, and add Atom Sites and Exp 
  Measurements. For clarity go through the Atom Sites in the order of the
  magnetisation flow of the experiment. Select the atom type in the pulldown
  menu besides the [Create new type] button and press the button to create the
  atom site. If you are measuring the shift (or another parameter) of this atom
  site, press the [Create New] button under the Exp Measurements table while the
  newly created atom site is still selected. This approach will minimise the 
  number of mouse clicks. Finally edit the table contents if needed.
  
  Now go to the Exp Graphs tab. Create a new Exp Graph in the top left table - 
  you only need more than one if you have several distinct magnetisation transfer
  flows in your experiment (essentially if you have curly braces in the systematic
  name). Now add as many magnetisation transfer steps as you need by clicking
  [Create New] under the bottom table, and edit the fields in the table to get the
  contents right. Finally add as many ExpSteps as you need, and edit the contents.
  Note that this table must follow the magnetisation step by step - an HSQC
  experiment would have steps, in order, '1 Shift(H); 2 Shift (N); 1 Shift (H)'.
  
  Finally go to the Ref Experiments tab and add the Ref Experiments. It is 
  easiest if you start by making the RefExperiment with the highest 
  dimensionality. There should always be a RefExperiment with the maximum
  dimensionality of the prototype, and with the same name as the prototype.
  Now add a [New ExpRefDim] in the bottom table for every dimension in your
  RefExperiment, and edit the bottom table. For clarity you should add the
  dimensions in the order of the magnetisation flow; the mapping to actual
  dimensions is done elsewhere. In most cases the only column you may need
  to change is 'expMeasurement'. 
  
  Now start creating the other RefExperiments. It is good practice to have a
  RefExperiment for all combinations of axes that people might sensibly choose
  to measure, even if they are not very likely. That means taking most
  combinations, but avoiding eg. 'H[N[ca]]' (useless), 'hnCAHA' (useless) and
  h[N[CA]] (technically impossible). First select the highest-dimensionality
  RefExperiment in the top table and click [Create Copy]. Change the name so
  that one or more atom sites are in lower rather than upper case, as a sign
  that these are not measured. Now select the corresponding dimensions in the
  lower table and press [Delete RefExpDim]. This will give you a correct
  RefExperiment with lower dimensionality, provided the highest-dimensionality
  RefExperiment was already correct. Most experiments can be reversed; the 
  exceptions are out-and-back experiments like H[n[CA]], or symmetric
  experiments like H[C]_H[C].NOESY. For reversible experiments you should make a
  further copy of your highest-dimensionality RefExperiment. Now change the
  name of the copy so that the atom sites are in the opposite order (e.g. from 
  HNCAHA to HACANH), and set the isReversed column to Yes. Make copies of the
  highest-dimensionality reverse RefExperiment, and use them to create lower 
  dimensionality reversed RefExperiments like you already did for the 
  non-reversed case. Finally add any information you may have for the Synonym
  or Systematic Name columns of individual RefExperiments.
  
  For many high-dimensionality experiments it will be useful to add a 2D
  Projection RefExperiment (and possibly a reversed version as well). For this,
  create a new RefExperiment (not a copy) in the top table, and give it a name
  of the form H[N[CO[CA]]].2D.{N;CO;CA} (for H[N[CO[CA]]]). Now go into the
  bottom table and add two RefExpDim, corresponding to the acquisition dimension
  (here H) and another dimension. At this point select the non-acquisition
  dimension in the bottom table, and click [Add RefExpDimRef] once for each
  dimension in the full prototype that is not already in the table. Finally
  edit the tables as desired.
  
  NOTE: Experiment Prototypes and RefExperiments are editable when created, but
  become frozen when the project has been closed and they are reloaded from 
  disk. If you *really need* to edit these data after the initial creation 
  (remember, this is reference information, and changing data that is used
  in actual projects may introduce errors) this is how to make things editable: 
  Go to the command line, and get hold of the MemopsRoot object (you can get 
  it as top.project). Then select the NmrExpPrototype you want to edit, and
  set myNmrExpPrototype.isEditable = True. Do the same for existing 
  RefExperiments you want to modify. If you wanted to edit the H[n[CA]] 
  RefExperiment, for instance (hint: you do not!), you would type something
  like:
  nmrExpPrototype = top.project.findFirstNmrExpPrototype(serial=20)
  nmrExpPrototype.isEditable = True
  refExperiment = nmrExpPrototype.findFirstRefExperiment(name='H[n[CA]]')
  refExperiment.isEditable = True

  **References**
  
  *A nomenclature and data model to describe NMR experiments.
  Fogh RH, Vranken WF, Boucher W, Stevens TJ, Laue ED.
  J Biomol NMR. 2006 Nov;36(3):147-55* (link_)

  .. _`link`:  http://www.ncbi.nlm.nih.gov/pubmed/17031528

  """

  def __init__(self, parent, project, *args, **kw):

    self.guiParent        = parent
    self.nmrExpPrototype  = None
    self.atomSiteType     = None
    self.help_url =  'file:%s/python/ccpnmr/analysis/doc/build/html/popups/NmrExpPrototypePopup.html' % (
     getTopDirectory()
    )
    self.waiting          = False
    
    if not kw.has_key('popup_name'):
      kw['popup_name'] = "NmrExpPrototypes"
    
    BasePopup.__init__(self, parent, project, **kw)

  def body(self, guiFrame):
  
    self.geometry('600x700')

    guiFrame.grid_columnconfigure(1, weight=1)
    guiFrame.grid_rowconfigure(1, weight=1)

    label = Label(guiFrame, text='Selected Prototype:')
    label.grid(row=0, column=0, sticky='w')
    tipText = 'Select the experiment prototype to operate on'
    
    self.selectedPrototypePulldown = PulldownList(guiFrame,
                                                  tipText=tipText,
                                                  callback=self.selectPulldownPrototype)
    self.selectedPrototypePulldown.grid(row=0, column=1, sticky='w')
   
    options = ['Prototypes','Atom & Measurements',
               'Exp Graphs', 'Ref Experiments',]
    
    tipTexts = [
     'A table of all NMR experiment prototypes',
     'Atoms and NMR measurements used within selected NMR experiment prototypes',
     'Graph of magnetisation transfer types, and magnetisation flow, '
      'for selected experiment prototype',
     'Reference experiments belonging to selected experiment prototype '
      'with description of the experiment axes.',
    ]
    
    tabbedFrame = TabbedFrame(guiFrame, options=options, tipTexts=tipTexts)
    tabbedFrame.grid(row=1, column=0, columnspan=2, sticky='nsew')
    self.tabbedFrame = tabbedFrame
    frameA, frameB, frameC, frameD = tabbedFrame.frames

    frameA.grid_columnconfigure(0, weight=1)
    frameA.grid_rowconfigure(0, weight=1)
    frameB.grid_columnconfigure(0, weight=1)
    frameB.grid_rowconfigure(0, weight=1)
    frameB.grid_rowconfigure(1, weight=1)
    frameC.grid_columnconfigure(0, weight=1)
    frameC.grid_rowconfigure(0, weight=1)
    frameD.grid_columnconfigure(0, weight=1)
    frameD.grid_rowconfigure(0, weight=1)
    
    # edit widgets
    # frame0
    self.nameEntry = Entry(self, 
     text='', returnCallback=self.setName, width=12
    )
    self.categoryPulldown = PulldownList(self, texts=categoryList,
                                         callback=self.setCategory)
    self.synonymEntry = Entry(self, 
     text='', returnCallback=self.setSynonym, width=12
    )
    self.keywordsEntry = MultiWidget(self, Entry, minRows=0,
     callback=self.setKeywords
    )
    self.detailsEntry = Entry(self,
     text='', returnCallback = self.setDetails, width=16
    )

    # frame1
    self.isotopePulldown = PulldownList(self, texts=isotopeList,
                                        callback=self.setIsotope)
 
    self.atomSiteNameEntry = Entry(self, 
     text='', returnCallback=self.setAtomSiteName, width=12
    )
    self.minShiftEntry = FloatEntry(self,
     text='', returnCallback=self.setMinShift, width=10
    )
    self.maxShiftEntry = FloatEntry(self,
     text='', returnCallback=self.setMaxShift, width=10
    )
    self.minNumberEntry = IntEntry(self,
     text='', returnCallback=self.setMinNumber, width=3
    )
    self.maxNumberEntry = IntEntry(self,
     text='', returnCallback=self.setMaxNumber, width=3
    )
    self.numberStepEntry = IntEntry(self,
     text='', returnCallback=self.setNumberStep, width=3
    )
    self.expMeasurementTypePulldown = PulldownList(self, callback=self.setMeasurementType)

    self.atomSitesEntry = MultiWidget(self, CheckButton, minRows=1,
     callback=self.setAtomSites
    )
    self.atomSiteWeightsEntry = MultiWidget(self, FloatEntry, minRows=0,
     callback=self.setAtomSiteWeights
    )
        
    #
    # Prototypes
    #
    
    colHeadings = ['#','Name','Atom\nSites','Graphs','Ref\nExpts',
                   'Max\nDim','Category','Synonym',
                   'Keywords','Details','isEditable?']
                   
    editWidgets = [None, self.nameEntry,
                   None, None, None,
                   None, self.categoryPulldown,
                   self.synonymEntry, self.keywordsEntry,
                   self.detailsEntry, None]
    
    editGetCallbacks = [None, self.getName,
                        None, None, None,
                        None, self.getCategory,
                        self.getSynonym, self.getKeywords,
                        self.getDetails, None]
                        
    editSetCallbacks = [None, self.setName,
                        None, None, None,
                        None, self.setCategory,
                        self.setSynonym, None,
                        self.setDetails, None]
                        
    
    tipTexts = [
     'The serial number of the experiment prototype', 
     'CCPN systematic name for the experiment prototype. '
      'Reflects the magnetisation transfer network', 
     'Number of atom sites in the magnetisation transfer flow(s)', 
     'Number of distinct magnetisation transfer flows', 
     'Number of associated reference experiments', 
     'Maximum reference experiment dimension', 
     'Experiment prototype category. Classification of experiment purpose.', 
     'Human-understandable experiment prototype synonym', 
     'List of keywords, for search and filtering', 
     'A user-editable textual comment for the experiment prototype',  
     'Can experiment prototype be modified? True for newly created experiments', 
    ]
    
    self.prototypeMatrix = ScrolledMatrix(frameA,
                                          multiSelect=False,
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks,
                                          editWidgets=editWidgets,
                                          headingList=colHeadings,
                                          tipTexts=tipTexts,
                                          callback=self.selectNmrExpPrototype)
 
    self.prototypeMatrix.grid(row=0, column=0, sticky='nsew')
    
    nmrExpPrototype = self.project.findFirstNmrExpPrototype()
    self.nmrExpPrototype = nmrExpPrototype
    self.prototypeMatrix.currentObject = nmrExpPrototype
    
    texts    = ['Create New','Create Copy', 'Freeze', 'Delete']
    
    tipTexts = [
     'Make a new blank experiment prototype. Information is added subsequently',
     'Make a new experiment prototype '
      'that is a copy of the currently selected one',
     'Set selected experiment prototype to be no longer modifiable',
     'Delete the selected experiment prototype',
    ]
    
    commands = [self.newNmrExpPrototype, self.newCopyNmrExpPrototype,
                self.freezeNmrExpPrototype, self.deleteNmrExpPrototype]
                
    self.prototypeButtons = ButtonList(frameA, texts=texts, tipTexts=tipTexts,
                                       commands=commands)
    self.prototypeButtons.grid(row=1, column=0, sticky='nsew')

    self.prototypeMatrix.doEditMarkExtraRules = self.prototypeEditMarkExtraRules

    #
    # Atoms and measurements
    #
    
    tipText = 'Table of atom sites'

    frame1 = LabelFrame(frameB, text='Atom Sites', tipText=tipText)
    frame1.grid(row=0, column=0, sticky='nsew')
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=0)
    frame1.grid_columnconfigure(2, weight=1)
    frame1.grid_columnconfigure(3, weight=1)
    
    colHeadings      = [
     '#', 'Isotope', 'Name', 'Min\nShift', 'Max\nShift',
     'Min\nNumber', 'Max\nNumber', 'Number\nStep'
    ]
    editWidgets      = [None, self.isotopePulldown, self.atomSiteNameEntry,
     self.minShiftEntry, self.maxShiftEntry,
     self.minNumberEntry, self.maxNumberEntry, self.numberStepEntry
    ]
    editGetCallbacks = [None, self.getIsotope, self.getAtomSiteName,
     self.getMinShift, self.getMaxShift,
     self.getMinNumber, self.getMaxNumber, self.getNumberStep
    ]
    editSetCallbacks = [None, self.setIsotope, self.setAtomSiteName,
     self.setMinShift, self.setMaxShift,
     self.setMinNumber, self.setMaxNumber, self.setNumberStep
    ]
    
    tipTexts = [
     'The serial number of the atom site',
     'The isotope of the atom site',
     'The name of the atom site',
     'The minimum allowed chemical shift of the atom site',
     'The maximum allowed chemical shift of the atom site',
     'The minimum number of atoms in the site',
     'The maximum number of atoms in the site',
     'The number step - E.g. "a non-zero, even number" would have number_step=2, min_number=2',
    ]
    
    self.atomSiteMatrix = ScrolledMatrix(frame1,
                                         multiSelect=False,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         headingList=colHeadings,
                                         tipTexts=tipTexts,
                                         callback=self.selectAtomSite)
 
    self.atomSiteMatrix.grid(row=0, column=0, columnspan=4, sticky='nsew')

    self.atomSiteMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules
    
    tipText = 'Create new atom site of the selected type'
    
    self.newAtomSiteButton = Button(
     frame1, borderwidth=1, text = 'Create New - Type :', tipText=tipText,
     command=self.newAtomSite
    )
    self.newAtomSiteButton.grid(row=1, column=0,sticky='ew')
    
    
    tipText = 'Select atom site type to create'
    self.atomSiteTypePulldown = PulldownList(frame1, self.changeAtomSiteType,
                                             tipText=tipText,
                                             borderwidth=1, relief='raised')
                                             
    self.atomSiteTypePulldown.grid(row=1, column=1,sticky='w')
 
    tipText = 'Delete atom site'
    self.deleteAtomSiteButton = Button(
     frame1, bd=1, text = 'Delete', tipText=tipText,
     command=self.deleteAtomSite
    )
    self.deleteAtomSiteButton.grid(row=1, column=3,sticky='ew')
    
    tipText = 'Table of Experimental NMR measurements recorded in experiment.'
    
    frame2 = LabelFrame(frameB, text='Exp Measurements', tipText=tipText)
    frame2.grid(row=1, column=0, sticky='nsew')
    frame2.grid_rowconfigure(0, weight=1)
    frame2.grid_columnconfigure(0, weight=1)
    
    colHeadings = ['#', 'Type', 'Atom Sites', 'Atom Site Weights']
    
    tipTexts = [
     'The serial number of the experimental measurement',
     'Type of NMR measurement',
     'Atom site(s) that the measurement relates to',
     'Relative weights for multiple atom sites. '
      'Only for special types, like multiple quantum coherence.',
    ]
    
    editWidgets = [None, self.expMeasurementTypePulldown,
                   self.atomSitesEntry, self.atomSiteWeightsEntry]
                   
    editGetCallbacks = [None, 
     self.getMeasurementType, self.getAtomSites, self.getAtomSiteWeights
    ]
    editSetCallbacks = [None, self.setMeasurementType, None, None
    ]
    self.expMeasurementMatrix = ScrolledMatrix(frame2,
                                               multiSelect=False,
                                               editSetCallbacks=editSetCallbacks,
                                               editGetCallbacks=editGetCallbacks,
                                               editWidgets=editWidgets,
                                               headingList=colHeadings,
                                               tipTexts=tipTexts,
                                               callback=self.selectExpMeasurement)
 
    self.expMeasurementMatrix.grid(row=0, column=0, sticky='nsew')
    
    tipTexts = [
     'Create new experimental measurement, '
      'using the currently selected atom site',
     'Delete selected experimental measurement',
    ]
    
    texts    = ['Create New', 'Delete']
    commands = [self.newExpMeasurement, self.deleteExpMeasurement]
    self.expMeasurementButtons = ButtonList(
     frame2, texts=texts,commands=commands, tipTexts=tipTexts,
    )
    self.expMeasurementButtons.grid(row=1, column=0, sticky='nsew')

    self.expMeasurementMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules
    
    #
    # RefExperiments
    #
    
    self.refExperimentFrame = RefExperimentFrame(frameD, self.project,
                                                 self.nmrExpPrototype,
                                                 popup=self)
    self.refExperimentFrame.grid(row=0, column=0, sticky='nsew')
    
    #
    # Graphs
    #
    
    self.expGraphFrame = ExpGraphFrame(frameC, self.project,
                                       self.nmrExpPrototype,
                                       popup=self)
    self.expGraphFrame.grid(row=0, column=0, sticky='nsew')

    #
    # Main
    #
    
    buttons = createDismissHelpButtonList(tabbedFrame.sideFrame, help_url=self.help_url)
    buttons.grid(row=0, column=2, sticky='e')
    
    if self.atomSiteType not in atomSiteTypeList:
      self.atomSiteType = atomSiteTypeList[0]
      
    ii = atomSiteTypeList.index(self.atomSiteType)
    self.atomSiteTypePulldown.setup(atomSiteTypeList,atomSiteTypeList,ii)
    
    self.updateAfter()
    
    self.administerNotifiers(Implementation.registerNotify)
  
  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete',''):
      for className in ('ccp.nmr.NmrExpPrototype.NmrExpPrototype',
                        'ccp.nmr.NmrExpPrototype.AtomSite',
                        'ccp.nmr.NmrExpPrototype.ExpMeasurement',
                        'ccp.nmr.NmrExpPrototype.RefExperiment',
                        'ccp.nmr.NmrExpPrototype.ExpGraph'):
        notifyFunc(self.updateAfter, className, func)
  
  def destroy(self):
  
    self.administerNotifiers(Implementation.unregisterNotify)
    BasePopup.destroy(self)
  
  def updatePrototypeAll(self):
  
    nmrExpPrototype = self.nmrExpPrototype
    self.updateAfter()
    self.prototypeMatrix.selectObject(nmrExpPrototype)
    self.expGraphFrame.nmrExpPrototype = nmrExpPrototype
    self.refExperimentFrame.nmrExpPrototype = nmrExpPrototype
    self.expGraphFrame.updateAfter()
    self.refExperimentFrame.updateAfter()
   
  
  def selectPulldownPrototype(self, nmrExpPrototype):
  
    self.nmrExpPrototype = nmrExpPrototype
    self.updatePrototypeAll()
  
  def updatePulldownPrototypes(self):
  
    nmrExpPrototype = self.nmrExpPrototype

    index = 0
    names = []
    prototypes = self.project.sortedNmrExpPrototypes()
    cats = []
     
    if prototypes:
      if nmrExpPrototype not in prototypes:
        nmrExpPrototype = prototypes[0]
    
      for prototype in prototypes:
        data = (prototype.serial, prototype.name)
        names.append(expPrototypeLabelFormat % data)
        cats.append(prototype.category)
    
      index = prototypes.index(nmrExpPrototype)
    
    else:
      nmrExpPrototype = None
      
    if self.nmrExpPrototype is not nmrExpPrototype:
      self.nmrExpPrototype = nmrExpPrototype
      self.updatePrototypeAll()

    
    self.selectedPrototypePulldown.setup(names, prototypes, index, None, cats)
    
  def prototypeEditMarkExtraRules(self, nmrExpPrototype, row, col):
  
    if col in (2,3,4):
      return True
  
    if hasattr(nmrExpPrototype,'isEditable'):
      return True
    else:
      return False
      
  def doEditMarkExtraRules(self, obj, row, col):
  
    if self.nmrExpPrototype and hasattr(self.nmrExpPrototype,'isEditable'):
      return True
    else:
      return False
  
  def newNmrExpPrototype(self):
  
    obj = self.project.newNmrExpPrototype(name='dummy', category='other')
    menuString = expPrototypeLabelFormat % (obj.serial,obj.name)
    obj.isEditable = True

    self.update_idletasks()
    self.prototypeMatrix.selectObject(obj)
    self.updatePrototypeAll()
  
  def newCopyNmrExpPrototype(self):
  
    obj = genUtil.copySubTree( self.nmrExpPrototype, self.project,
     topObjectParameters={'name':'copy_of_' + self.nmrExpPrototype.name,
                          'isModifiable':True})
    menuString = expPrototypeLabelFormat % (obj.serial,obj.name)
    #self.nmrExpPrototype = obj
    obj.isEditable = True
    
    for refExperiment in obj.refExperiments:
      refExperiment.isEditable = True

    self.update_idletasks()
    self.prototypeMatrix.selectObject(obj)
    self.updatePrototypeAll()
  
  
  def freezeNmrExpPrototype(self):
  
    if (showOkCancel(
     'Freeze ExpPrototype',
     "Are you sure you want to freeze the ExpPrototype,\nso it cannot be changed?",
    parent=self)):
      obj = self.nmrExpPrototype
      if hasattr(obj,'isEditable'):
        del obj.isEditable
      if self.expGraphFrame:
        self.expGraphFrame.updateAfter()
      if self.refExperimentFrame:
        self.refExperimentFrame.updateAfter()
      self.updateAfter()
          
  def deleteNmrExpPrototype(self):
    
    if self.nmrExpPrototype:
      self.nmrExpPrototype.delete()
      
      if self.nmrExpPrototype is self.prototypeMatrix.currentObject:
        self.prototypeMatrix.currentObject = None
        
      self.nmrExpPrototype = None
 
  def setName(self, event):

    text = self.nameEntry.get()
    if text and text.strip():
      self.nmrExpPrototype.setName( text )

  def getName(self, nmrExpPrototype):
    
    if nmrExpPrototype :
      self.nameEntry.set(nmrExpPrototype.name)
 
  def setSynonym(self, event):

    text = self.synonymEntry.get()
    if text and text.strip():
      self.nmrExpPrototype.setSynonym( text )
    else:
      self.nmrExpPrototype.setSynonym( None )
      

  def getSynonym(self, nmrExpPrototype):
    
    if nmrExpPrototype :
      self.synonymEntry.set(nmrExpPrototype.synonym)

  def setCategory(self, event):
  
    category = self.categoryPulldown.getText()
  
    if self.nmrExpPrototype:
      self.nmrExpPrototype.setCategory(category)

  def getCategory(self, nmrExpPrototype):
  
    names = categoryList
    self.categoryPulldown.setup(names,names,0)
    if nmrExpPrototype:
      name = nmrExpPrototype.category
      self.categoryPulldown.set(name)
  
  def setKeywords(self, values):
    """ values is None if the widget aborts
    """
    if values is not None:
      self.nmrExpPrototype.keywords = values
    self.prototypeMatrix.keyPressEscape()
  
  def getKeywords(self, nmrExpPrototype):
    """ 
    """
    if nmrExpPrototype:
      self.keywordsEntry.setValues(nmrExpPrototype.keywords)  

  def setDetails(self, event):

    text = self.detailsEntry.get()
    if self.nmrExpPrototype:
      if text and text.strip():
        self.nmrExpPrototype.setDetails(text)
      else:
        self.nmrExpPrototype.setDetails(None)
    print self.nmrExpPrototype.details
 
  def getDetails(self, nmrExpPrototype):

    if nmrExpPrototype and nmrExpPrototype.details:
      self.detailsEntry.set(nmrExpPrototype.details)
     
  def selectNmrExpPrototype(self, obj, row, col):
    
    if obj and obj is not self.nmrExpPrototype:
      self.nmrExpPrototype = obj
      self.updateAfter(updateMainTable=False) 
      self.updatePrototypeAll()
 

  def setIsotope(self, event):
  
    isotope = self.isotopePulldown.getText()
    
    obj = self.atomSiteMatrix.currentObject
    if obj:
      obj.setIsotopeCode(isotope)
  
  def getIsotope(self, atomSite):
    if atomSite:
      self.isotopePulldown.set(atomSite.isotopeCode)
  
  def setAtomSiteName(self, event):

    text = self.atomSiteNameEntry.get()
    obj = self.atomSiteMatrix.currentObject
    if obj:
      if text and text.strip():
        obj.setName(text)
      else:
        obj.setName(None)
  
  def getAtomSiteName(self, atomSite):
    
    if atomSite:
      self.atomSiteNameEntry.set(atomSite.name)
  
  def setMinShift(self, event):
  
    x = self.minShiftEntry.get()
    obj = self.atomSiteMatrix.currentObject
    if obj and (obj.maxShift is None or x is None or obj.maxShift > x):
      obj.setMinShift(x)
  
  def getMinShift(self, atomSite):
    
    if atomSite:
      self.minShiftEntry.set(atomSite.minShift)
  
  def setMaxShift(self, event):
  
    x = self.maxShiftEntry.get()
    obj = self.atomSiteMatrix.currentObject
    if obj and (obj.minShift is None or x is None or obj.minShift < x):
      obj.setMaxShift(x)
  
  def getMaxShift(self, atomSite):
    
    if atomSite:
      self.maxShiftEntry.set(atomSite.maxShift)
  
  def setMinNumber(self, event):
  
    x = self.minNumberEntry.get()
    if x is not None:
      obj = self.atomSiteMatrix.currentObject
      if obj and (obj.maxNumber is None or obj.maxNumber >= x):
        obj.setMinNumber(x)
  
  def getMinNumber(self, atomSite):
    
    if atomSite:
      self.minNumberEntry.set(atomSite.minNumber)
  
  def setMaxNumber(self, event):
  
    x = self.maxNumberEntry.get()
    obj = self.atomSiteMatrix.currentObject
    if obj and (x is None or obj.minNumber <= x):
      obj.setMaxNumber(x)
  
  def getMaxNumber(self, atomSite):
    
    if atomSite:
      self.maxNumberEntry.set(atomSite.maxNumber)
  
  def setNumberStep(self, event):
  
    x = self.numberStepEntry.get()
    if x is not None:
      obj = self.atomSiteMatrix.currentObject
      if obj:
        obj.setNumberStep(x)
  
  def getNumberStep(self, atomSite):
    
    if atomSite:
      self.numberStepEntry.set(atomSite.numberStep)
  
  def newAtomSite(self):
  
    dd = atomSiteTypeData.get(self.atomSiteType) or {
     'isotopeCode':'1H','name':'H'
    }
    obj = self.nmrExpPrototype.newAtomSite(**dd)
    self.atomSiteMatrix.currentObject = obj
    self.updateAfter()
      
  def changeAtomSiteType(self,name):
    
    self.atomSiteType = name
          
  def deleteAtomSite(self):
    
    obj = self.atomSiteMatrix.currentObject
    if obj:
      self.atomSiteMatrix.currentObject = None
      obj.delete()
     
  def selectAtomSite(self, obj, row, col):
    if obj:
      self.updateButtons()
    
  def setMeasurementType(self, event):
  
    name = self.expMeasurementTypePulldown.getText()
    
    obj = self.expMeasurementMatrix.currentObject
    if obj:
      obj.setMeasurementType(name)
  
  def getMeasurementType(self, expMeasurement):
    
    if expMeasurement:
      ii = expMeasurementTypeList.index(expMeasurement.measurementType)
      self.expMeasurementTypePulldown.setup(expMeasurementTypeList,
                                            expMeasurementTypeList, ii)
  
  def setAtomSites(self,  values):
    """
    """
    if values is not None:
      # it should not be possible to change the atomSites
      # without the edit widget closing
      atomSites = self.nmrExpPrototype.sortedAtomSites()
      assert len(values) == len(atomSites)
      
      newAtomSites = [x[1] for x in zip(values, atomSites) if x[0]]
 
      obj = self.expMeasurementMatrix.currentObject
      weights = obj.atomSiteWeights
      if weights and len(weights) != len(newAtomSites):
        obj.atomSiteWeights = []
      obj.atomSites = newAtomSites
      
    self.expMeasurementMatrix.keyPressEscape()
  
  def getAtomSites(self, expMeasurement):
    """ 
    """
    if expMeasurement:
      atomSites = self.nmrExpPrototype.sortedAtomSites()
      length = len(atomSites)
      currentSites = expMeasurement.atomSites
 
      options = [ '%s(%s)' % (x.serial,x.name) for x in atomSites]
      values = [False] * length
      for ii in range(length):
        if atomSites[ii] in currentSites:
          values[ii] = True
 
      self.atomSitesEntry.set(options=options, values=values)
      
  def setAtomSiteWeights(self, values):
    """ 
    """
    if values is not None:
      obj = self.expMeasurementMatrix.currentObject
      if not values or len(values) == len(obj.atomSites):
        obj.atomSiteWeights = values
    self.expMeasurementMatrix.keyPressEscape()
  
  def getAtomSiteWeights(self, expMeasurement):
    """
    """
    if expMeasurement:
      nAtomSites = len(expMeasurement.atomSites)
      self.atomSiteWeightsEntry.setValues(expMeasurement.atomSiteWeights)  
  
  def newExpMeasurement(self):
    atomSite = (self.atomSiteMatrix.currentObject
     or self.nmrExpPrototype.findFirstAtomSite()
    )
    obj = self.nmrExpPrototype.newExpMeasurement(
     measurementType='Shift',
     atomSites=(atomSite,)
    )
    self.expMeasurementMatrix.currentObject = obj
    self.updateAfter()
          
  def deleteExpMeasurement(self):
    
    obj = self.expMeasurementMatrix.currentObject
    if obj:
      obj.delete()
      self.expMeasurementMatrix.currentObject = None
     
  def selectExpMeasurement(self, obj, row, col):
    if obj:
      self.updateButtons()  
  
  def editExpGraphs(self):
    """ 
    """
    
    self.tabbedFrame.select(2)
    
    frame = self.expGraphFrame
    if frame.nmrExpPrototype is not self.nmrExpPrototype:
      frame.nmrExpPrototype = self.nmrExpPrototype
      frame.updateAfter()
    

  #def editExpMeasurements(self):
  #  """ 
  #  """
  #  
  #  self.tabbedFrame.select(3)
  #  
  #  frame = self.refExperimentFrame
  #  if frame.nmrExpPrototype is not self.nmrExpPrototype:
  #    frame.nmrExpPrototype = self.nmrExpPrototype
  #    frame.updateAfter()

  def updateAfter(self, obj=None, updateMainTable=True):
    
    if obj is not None and obj.isDeleted:
      ss = obj.qualifiedName
      if ss == 'ccp.nmr.NmrExpPrototype.NmrExpPrototype':
        self.nmrExpPrototype  = None
    
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
 
  def updateButtons(self):
    
    if self.nmrExpPrototype:
      self.prototypeButtons.buttons[1].enable()
    else:
      self.prototypeButtons.buttons[1].disable()

    if self.nmrExpPrototype and hasattr(self.nmrExpPrototype,'isEditable'):
      for ii in (2,3):
        self.prototypeButtons.buttons[ii].enable()
      self.newAtomSiteButton.enable()
      if self.atomSiteMatrix.currentObject:
        self.deleteAtomSiteButton.enable()
      else:
        self.deleteAtomSiteButton.disable()
      if self.nmrExpPrototype.atomSites:
        self.expMeasurementButtons.buttons[0].enable()
      else:
        self.expMeasurementButtons.buttons[0].disable()
      if self.expMeasurementMatrix.currentObject:
        self.expMeasurementButtons.buttons[1].enable()
      else:
        self.expMeasurementButtons.buttons[1].disable()
    else:
      for ii in (2,3):
        self.prototypeButtons.buttons[ii].disable()
      self.newAtomSiteButton.disable()
      self.deleteAtomSiteButton.disable()
      for ii in (0,1):
        self.expMeasurementButtons.buttons[ii].disable()

  def update(self, updateMainTable=True):
    """  Fill tables with values, initial or after change
    """
    
    self.updatePulldownPrototypes()

    
    # Buttons en/disable
    
    self.updateButtons()
    
    # NmrExpPrototype matrix
    
    if updateMainTable:
    
      numCols = 6
 
      objectList = self.project.sortedNmrExpPrototypes()
      textMatrix = len(objectList)*[None]
 
      for ii,expPrototype in enumerate(objectList):
        maxDim = max(0, 0, *(len(x.refExpDims) for x in expPrototype.refExperiments))
        isEditable = hasattr(expPrototype,'isEditable') and 'Yes' or 'No'
        name = expPrototype.name
        synonym = expPrototype.synonym or ''
        
        name = re.sub('(.{20}.+)(\s+[|].+)',r'\1\n\2',name)
        synonym = re.sub('(.{20}.+)(\s+or.+)',r'\1\n\2',synonym)
        
        #atomSites = ','.join([as.name for as in expPrototype.atomSites])
        
        textMatrix[ii] = (expPrototype.serial,
                          name,
                          len(expPrototype.atomSites),
                          len(expPrototype.expGraphs),
                          len(expPrototype.refExperiments),
                          maxDim,
                          expPrototype.category,
                          synonym,
                          expPrototype.keywords,
                          expPrototype.details,
                          isEditable)
 
 
      self.prototypeMatrix.update(objectList=objectList, textMatrix=textMatrix)
    
    #AtomSites Matrix
    
    numCols = 8
    
    expPrototype = self.nmrExpPrototype
    if expPrototype:
      objectList = expPrototype.sortedAtomSites()
      textMatrix = len(objectList)*[None]
      
      for ii,atomSite in enumerate(objectList):
        textMatrix[ii] = (
         atomSite.serial,
         atomSite.isotopeCode,
         atomSite.name,
         atomSite.minShift,
         atomSite.maxShift,
         atomSite.minNumber,
         atomSite.maxNumber,
         atomSite.numberStep,
        )
        
    else:
      objectList = []
      textMatrix = []
      
    self.atomSiteMatrix.update(objectList=objectList, textMatrix=textMatrix)
       
    #ExpMeasurements Matrix
    
    numCols = 4
    
    expPrototype = self.nmrExpPrototype
    if expPrototype:
      objectList = expPrototype.sortedExpMeasurements()
      textMatrix = len(objectList)*[None]
      
      for ii,expMeasurement in enumerate(objectList):
        ll = textMatrix[ii] = numCols*[None]
        ll[0] = expMeasurement.serial
        ll[1] = expMeasurement.measurementType
        ll[2] = ", " .join(
         ['%s(%s)' % (x.serial,x.name) for x in expMeasurement.atomSites]
        )
        ll[3] = ", " .join([str(x) for x in expMeasurement.atomSiteWeights])
    
    else:
      objectList = []
      textMatrix = []
    
    self.expMeasurementMatrix.update(
     objectList=objectList, textMatrix=textMatrix
    )
    
    self.waiting = False


class ExpGraphFrame(Frame):
  """Frame for selecting and editing ExpGraphs. 
  """
  
  peakSignNames = ['<None>','+1','-1']
  peakSignData = {'<None>':None,'+1':1,'-1':-1}

  def __init__(self, parent, project, nmrExpPrototype, popup, *args, **kw):

    self.guiParent        = parent
    self.nmrExpPrototype  = nmrExpPrototype
    self.waiting          = False
    self.project = project
    self.popup = popup
    self.expGraph = None
    
    Frame.__init__(self, parent, **kw)

    self.grid_columnconfigure(0, weight=1)
    self.grid_columnconfigure(1, weight=1)
    self.grid_rowconfigure(1, weight=1)
    self.grid_rowconfigure(2, weight=1)

    row = 0
    
    frame0 = Frame(self)
    frame0.grid(row=row,column=0,columnspan=2,sticky='nsew')
    frame0.grid_columnconfigure(2, weight=1)
    
    self.peakSignPulldown = PulldownList(self.popup, 
                                         texts=self.__class__.peakSignNames,
                                         callback=self.setPeakSign)
     
    self.stepNumberEntry = IntEntry(self.popup,
     text='', returnCallback=self.setStepNumber, width=3
    )
    self.expMeasurementPulldown = PulldownList(self.popup, callback=self.setExpMeasurement)
    self.expTransferTypePulldown = PulldownList(self.popup, callback=self.setTransferType)
    self.atomSite1Pulldown = PulldownList(self.popup, callback=self.setAtomSite1)
    self.atomSite2Pulldown = PulldownList(self.popup, callback=self.setAtomSite2)
    
    row += 1
    self.grid_rowconfigure(row, weight=1)
    tipText = 'Table of magnetisation flow graphs'
    frame1 = LabelFrame(self, text='ExpGraphs', tipText=tipText)
    frame1.grid(row=row, column=0, sticky='nsew')
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)
    colHeadings      = ['#','PeakSign']
    editWidgets      = [None, self.peakSignPulldown]
    editGetCallbacks = [None, self.getPeakSign]
    editSetCallbacks = [None, None]
    tipTexts = [
     'The serial number of the magnetisation flow graph',
     'The relative sign of peaks arising from the graph',
    ]
    self.expGraphMatrix = ScrolledMatrix(frame1,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         initialRows=3, headingList=colHeadings,
                                         tipTexts=tipTexts,
                                         callback=self.selectExpGraph)
    self.expGraphMatrix.grid(row=0, column=0, sticky='nsew')

    self.expGraphMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules
 
    texts    = ['Create New','Delete']
    commands = [self.newExpGraph, self.deleteExpGraph]
    tipTexts = [
     'Create new magnetisation flow graph',
     'Delete selected magnetisation flow graph',
    ]
    self.expGraphButtons = ButtonList(
     frame1, texts=texts,commands=commands, tipTexts=tipTexts,
    )
    self.expGraphButtons.grid(row=1, column=0, sticky='nsew')
    
    tipText = 'Table of magnetisation flow steps for selected graph.'
    frame2 = LabelFrame(self, text='ExpSteps', tipText=tipText)
    frame2.grid(row=row, column=1, sticky='nsew')
    frame2.grid_rowconfigure(0, weight=1)
    frame2.grid_columnconfigure(0, weight=1)
    colHeadings      = ['#','StepNumber', 'Measurement']
    editWidgets      = [
     None,self.stepNumberEntry, self.expMeasurementPulldown
     ]
    editGetCallbacks = [None,self.getStepNumber, self.getExpMeasurement]
    editSetCallbacks = [None,self.setStepNumber, self.setExpMeasurement]
    tipTexts = [
     'Serial number of the magnetisation flow step (corresponds to creation order)',
     'Position of step in magnetisation flow (corresponds to flow order)',
     'Experimental measurement that carries the magnetisation during the step',
    ]
    self.expStepMatrix = ScrolledMatrix(frame2,
     editSetCallbacks=editSetCallbacks, multiSelect=0, tipTexts=tipTexts,
     editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,
     initialRows=3, headingList=colHeadings, callback=self.selectExpStep)
     
    self.expStepMatrix.grid(row=0, column=0, sticky='nsew')

    self.expStepMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules

    texts    = ['Create New', 'Re-number', 'Delete']
    tipTexts = [
     'Create new magnetisation flow step for selected graph. '
      'By default created with next unused experimental measurement',
     'Re-number steps consecutively starting from one',
     'Delete selected magnetisation flow step',
    ]
    commands = [self.newExpStep, self.renumberSteps, self.deleteExpStep]
    self.expStepButtons = ButtonList(
     frame2, texts=texts,commands=commands, tipTexts=tipTexts,
    )
    self.expStepButtons.grid(row=1, column=0, sticky='nsew')

    row += 1
    self.grid_rowconfigure(row, weight=1)
    tipText = ('Table of magnetisation transfer connections between atom sites'
               ' for selected graph')
    frame3 = LabelFrame(self, text='ExpTransfers', tipText=tipText)
    frame3.grid(row=row, column=0, columnspan=2, sticky='nsew')
    frame3.grid_rowconfigure(0, weight=1)
    frame3.grid_columnconfigure(0, weight=1)
    colHeadings      = ['#', 'transfer\nType', 'AtomSite1', 'AtomSite2','transfer\nToSelf']
    editWidgets      = [None, self.expTransferTypePulldown,
     self.atomSite1Pulldown, self.atomSite2Pulldown, None
    ]
    editGetCallbacks = [None, self.getTransferType, 
     self.getAtomSite1, self.getAtomSite2, self.toggleTransferToSelf
    ]
    editSetCallbacks = [None, self.setTransferType, 
     self.setAtomSite1, self.setAtomSite2, None
    ]
    tipTexts = [
     'Serial number of the magnetisation transfer connection',
     'Type of magnetisation transfer',
     'First atom site involved in magnetisation transfer',
     'Second atom site involved in magnetisation transfer',
     'Does the transfer allow magnetisation to remain on the starting site?',
    ]
    
    self.expTransferMatrix = ScrolledMatrix(frame3,
     editSetCallbacks=editSetCallbacks, multiSelect=0,
     editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,
     initialRows=6, headingList=colHeadings, tipTexts=tipTexts,
     callback=self.selectExpTransfer)
     
    self.expTransferMatrix.grid(row=0, column=0, sticky='nsew')

    self.expTransferMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules
    
    texts    = ['Create New', 'Delete']
    commands = [self.newExpTransfer, self.deleteExpTransfer]
    tipTexts = [
     'Create new magnetisation transfer connection for selected graph',
     'Delete selected magnetisation transfer connection',
    ]
    self.expTransferButtons = ButtonList(
     frame3, texts=texts,commands=commands, tipTexts=tipTexts
    )
    self.expTransferButtons.grid(row=1, column=0, sticky='nsew')

    self.updateAfter()
    self.administerNotifiers(Implementation.registerNotify)
  
  def administerNotifiers(self, notifyFunc):
    

    notifyFunc(
     self.updateAfter, 'ccp.nmr.NmrExpPrototype.ExpMeasurement', 'delete'
    )
    for func in ('delete','setName'):
      for className in ('ccp.nmr.NmrExpPrototype.AtomSite',):
        notifyFunc(self.updateAfter, className, func)
    for func in ('delete','setMeasurementType'):
      for className in ('ccp.nmr.NmrExpPrototype.ExpMeasurement',):
        notifyFunc(self.updateAfter, className, func)
    for func in ('__init__', 'delete',''):
      for className in ('ccp.nmr.NmrExpPrototype.ExpGraph',
       'ccp.nmr.NmrExpPrototype.ExpStep','ccp.nmr.NmrExpPrototype.ExpTransfer'
      ):
        notifyFunc(self.updateAfter, className, func)

  def destroy(self):
  
    self.administerNotifiers(Implementation.unregisterNotify)
    Frame.destroy(self)

    
  def doEditMarkExtraRules(self, obj, row, col):
    if self.nmrExpPrototype and hasattr(self.nmrExpPrototype,'isEditable'):
      return True
    else:
      return False

      
  def setPeakSign(self, event):
  
    name = self.peakSignPulldown.getText()
    
    value = self.__class__.peakSignData.get(name)
    
    obj =  self.expGraph
    if obj:
      obj.setPeakSign(value)
  
  def getPeakSign(self, expGraph):
  
    if expGraph:
      self.peakSignPulldown.set(expGraph.peakSign)
  
  def newExpGraph(self):
  
    obj = self.nmrExpPrototype.newExpGraph()
    self.expGraphMatrix.currentObject = self.expGraph = obj
    self.updateAfter()     
    
  def deleteExpGraph(self):
  
    obj =  self.expGraph
    if obj:
      obj.delete()
  
  def selectExpGraph(self, obj, row, col):

    if obj is self.expGraph:
      return
    else:
      self.expGraph = obj
      self.updateExpStep()
      self.updateExpTransfer()
      self.updateButtons()
  
  def setStepNumber(self, event):
    """ Sets stepNumber to desired number. In case of clash, recursively
    adds one to stepNumber of expStep that had the number previously.
    """
  
    num = self.stepNumberEntry.get()
    if num is not None:
      obj1 = self.expStepMatrix.currentObject
      if obj1:
        stepDict = {}
        for expStep in self.expGraph.expSteps:
          stepDict[expStep.stepNumber] = expStep
        del stepDict[obj1.stepNumber]
 
        while True:
          obj2 = stepDict.get(num)
          obj1.stepNumber = num
          if obj2 is None:
            break
          else:
            num = num + 1
            obj1 = obj2
  
  def getStepNumber(self, expStep):
  
    if expStep:
      self.stepNumberEntry.set(expStep.stepNumber)
  
  def renumberSteps(self):
    
    expGraph =  self.expGraph
    if expGraph:
      ll = []
      for expStep in expGraph.expSteps:
        ll.append((expStep.stepNumber,expStep))
      
      ll.sort()
      for ii,xx in enumerate(ll):
        xx[1].stepNumber = ii+1
      
      
  def setExpMeasurement(self, event):
  
    expMeasurement = self.expMeasurementPulldown.getObject()
    
    obj = self.expStepMatrix.currentObject
    if expMeasurement:
      obj.setExpMeasurement(expMeasurement)
  
  def getExpMeasurement(self, expStep):
    expMeasurements = self.nmrExpPrototype.sortedExpMeasurements()
    entries = [
     '%s(%s)' % (x.serial,','.join([y.name for y in x.atomSites])) 
     for x in expMeasurements
    ]
    if entries:
      try:
        ii = expMeasurements.index(expStep.expMeasurement)
      except ValueError:
        ii = 0
      self.expMeasurementPulldown.setup(entries, expMeasurements, ii)
  
  def newExpStep(self):
  
    expMeasurement = NmrExpPrototypeUtil.nextFreeExpMeasurement(
                         self.nmrExpPrototype, 'expSteps')
    if not expMeasurement:
      msg = 'Cannot make step: No measurements in current prototype.'
      showWarning('Failure', msg, parent=self)
      return
    
    
    expGraph = self.expGraph
    if expGraph:
    
      dd = {}
      for expStep in expGraph.expSteps:
        dd[expStep.stepNumber] = None
      num = 1
      while dd.has_key(num):
        num += 1
 
      expGraph.newExpStep(stepNumber=num,expMeasurement=expMeasurement)
          
  def deleteExpStep(self):
    
    obj = self.expStepMatrix.currentObject
    if obj:
      obj.delete()
  
  def selectExpStep(self, obj, row, col):

    if obj:
      self.updateExpStep()
      
  def setTransferType(self, event):
  
    name = self.expTransferTypePulldown.getText()
    
    obj = self.expTransferMatrix.currentObject
    if obj:
      obj.setTransferType(name)
  
  def getTransferType(self, expTransfer):
    
    if expTransfer:
      ii = expTransferTypeList.index(expTransfer.transferType)
      self.expTransferTypePulldown.setup(expTransferTypeList, expTransferTypeList, ii)
      
  def setAtomSite1(self, event):    
  
    obj = self.expTransferMatrix.currentObject
    if obj:
      site1,site2 = obj.atomSites
      if site1.serial > site2.serial:
        site1,site2 = site2,site1
        
      site1 = self.atomSite1Pulldown.getObject()
      if site1.serial > site2.serial:
        site1,site2 = site2,site1
      if site1 is not None and site1 is not site2:
        obj.atomSites = (site1, site2)
  
  def getAtomSite1(self, expTransfer):
  
    atomSites = self.nmrExpPrototype.sortedAtomSites()
    entries = [ '%s(%s)' % (x.serial,x.name) for x in atomSites]
    if expTransfer and atomSites:
      site1,site2 = expTransfer.atomSites
      if site1.serial > site2.serial:
        site1,site2 = site2,site1
      try:
        ii = atomSites.index(site1)
      except ValueError:
        ii = 0
      self.atomSite1Pulldown.setup(entries, atomSites, ii)
  
  def setAtomSite2(self, event):
  
    obj = self.expTransferMatrix.currentObject
    if obj:
      site1,site2 = obj.atomSites
      if site1.serial > site2.serial:
        site1,site2 = site2,site1
        
      site2 =  self.atomSite2Pulldown.getObject()
      if site1.serial > site2.serial:
        site1,site2 = site2,site1
      if site1 is not None and site1 is not site2:
        obj.atomSites = (site1, site2)
 
  def getAtomSite2(self, expTransfer):
  
    atomSites = self.nmrExpPrototype.sortedAtomSites()
    entries = [ '%s(%s)' % (x.serial,x.name) for x in atomSites]
    if expTransfer and atomSites:
      site1,site2 = expTransfer.atomSites
      if site1.serial > site2.serial:
        site1,site2 = site2,site1
      try:
        ii = atomSites.index(site2)
      except ValueError:
        ii = 0
      self.atomSite2Pulldown.setup(entries, atomSites, ii)
    
  def toggleTransferToSelf(self, expTransfer):
    oldValue = expTransfer.transferToSelf
    if oldValue:
      expTransfer.transferToSelf = False
    else:
      # set transferToSelf to True only if atomSites are compatible
      isotopeCodes = set(x.isotopeCode for x in expTransfer.atomSites)
      if len(isotopeCodes) == 1:
        expTransfer.transferToSelf = True
    
    #self.updateExpTransfer()
    
    
  def newExpTransfer(self):
    
    nmrExpPrototype = self.nmrExpPrototype
    atomSites = nmrExpPrototype.sortedAtomSites()
    
    if len(atomSites) < 2:
      msg = 'To make transfers at least two atom sites are required.'
      showWarning('Failure', msg, parent=self)
      return
    
    obj = self.expGraph
    if obj is not None and len(atomSites)>1:
      obj.newExpTransfer( transferType=expTransferTypeList[0],
       atomSites=(atomSites[0],atomSites[-1])
      )
    self.expTransferMatrix.currentObject = obj
    self.updateExpTransfer()
          
  def deleteExpTransfer(self):
    
    obj = self.expTransferMatrix.currentObject
    if obj:
      obj.delete()
  
  def selectExpTransfer(self, obj, row, col):

    if obj:
      self.updateExpTransfer()

  def booleanString(self, boolean):
 
    if boolean:
      return 'Yes'
    else:
      return 'No'
    
  def updateAfter(self, obj=None):
    
    if obj is not None and obj.isDeleted:
      ss = obj.qualifiedName
      if ss == 'ccp.nmr.NmrExpPrototype.NmrExpPrototype':
        self.nmrExpPrototype  = None

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
 
  def update(self):
    """  Fill table with values, initial or after change
    """
    
    
    # select single ExpGraph
    if self.nmrExpPrototype and self.nmrExpPrototype.isDeleted:
      self.nmrExpPrototype = None
      
    if self.nmrExpPrototype:
      expGraphs =  self.nmrExpPrototype.expGraphs
      expGraph = self.expGraph
      if expGraph in expGraphs:
        pass
      elif len(expGraphs) == 1:
        expGraph = iter(expGraphs).next()
      else:
        expGraph = None
    else:
      expGraph = None
      
    self.expGraphMatrix.currentObject = self.expGraph = expGraph 
           
    if self.nmrExpPrototype:
      objectList = self.nmrExpPrototype.sortedExpGraphs()
      textMatrix = len(objectList)*[None]
    else:
      objectList = []
      textMatrix = []
    
    for ii,xx in enumerate(objectList):
      textMatrix[ii] = (xx.serial, xx.peakSign)
      
    self.expGraphMatrix.update(objectList=objectList, textMatrix=textMatrix)
    
    self.updateExpStep()
    self.updateExpTransfer()
    self.updateButtons()
    
    self.waiting = False
  
  def updateButtons(self):
    
    expGraph = self.expGraph
     
    # en/disable buttons
    if self.nmrExpPrototype and hasattr(self.nmrExpPrototype,'isEditable'):
    
      self.expGraphButtons.buttons[0].enable()
      if expGraph:
        self.expGraphButtons.buttons[1].enable()
        self.expStepButtons.buttons[0].enable()
        self.expTransferButtons.buttons[0].enable()
      else:
        self.expGraphButtons.buttons[1].disable()
        self.expStepButtons.buttons[0].disable()
        self.expTransferButtons.buttons[0].disable()
        
    else:
      self.expStepButtons.buttons[2].disable()
      for ii in (0,1):
        for ll in (self.expGraphButtons.buttons,
         self.expStepButtons.buttons, self.expTransferButtons.buttons
        ):
          ll[ii].disable()     
  
  
  
  def updateExpStep(self):
  
    numCols = 3
    
    expGraph = self.expGraph
    
    if expGraph:
      ll = [(x.stepNumber,x) for x in expGraph.expSteps]
      ll.sort()
      objectList = [x[1] for x in ll]
      textMatrix = len(objectList)*[None]
    else:
      objectList = []
      textMatrix = []
    
    for ii,expStep in enumerate(objectList):
      ll = textMatrix[ii] = numCols*[None]
      ll[0] = expStep.serial
      ll[1] = expStep.stepNumber
      x = expStep.expMeasurement
      ll[2] = "%s %s(%s)" % (
       x.serial, x.measurementType, ','.join([y.name for y in x.atomSites])
      )
      
    self.expStepMatrix.update(objectList=objectList, textMatrix=textMatrix)
      
    if expGraph and len(expGraph.expSteps) > 1:
      self.expStepButtons.buttons[1].enable()
    else:
      self.expStepButtons.buttons[1].disable()
      
    if expGraph and self.expStepMatrix.currentObject:
      self.expStepButtons.buttons[2].enable()
    else:
      self.expStepButtons.buttons[2].disable()
  
  
  def updateExpTransfer(self):
  
    numCols = 5
    
    expGraph = self.expGraph
    
    
    if expGraph:
      objectList = expGraph.sortedExpTransfers()
      textMatrix = len(objectList)*[None]
    else:
      objectList = []
      textMatrix = []
    
    for ii,expTransfer in enumerate(objectList):
      ll = textMatrix[ii] = numCols*[None]
      ll[0] = expTransfer.serial
      ll[1] = expTransfer.transferType
      site1,site2 = expTransfer.atomSites
      if site1.serial > site2.serial:
        site1,site2 = site2,site1
      ll[2] = '%s (%s)' % (site1.serial,site1.name)
      ll[3] = '%s (%s)' % (site2.serial,site2.name)
      ll[4] = self.booleanString(expTransfer.transferToSelf)
      
    self.expTransferMatrix.update(objectList=objectList, textMatrix=textMatrix)
    
      
    if self.expTransferMatrix.currentObject:
      self.expTransferButtons.buttons[1].enable()
    else:
      self.expTransferButtons.buttons[1].disable()


class RefExperimentFrame(Frame):
  """Popup for selecting and editing RefExperiments. 
  """

  def __init__(self, parent, project, nmrExpPrototype, popup, *args, **kw):

    self.guiParent        = parent
    self.nmrExpPrototype  = nmrExpPrototype
    self.refExperiment    = None
    self.refExpDimRef     = None
    self.namingSystem     = namingSystemList[0]
    self.waiting          = False
    self.project = project
    self.popup = popup
    
    self.experimentEditable = False
    
    Frame.__init__(self, parent, **kw)

    row = 0
    self.grid_columnconfigure(0, weight=1)
    self.grid_rowconfigure(1, weight=1)
    self.grid_rowconfigure(2, weight=1)
    
    # Naming System Pulldown
    frame0 = Frame(self)
    frame0.grid(row=row,column=0,sticky='nsew')
    frame0.grid_columnconfigure(3, weight=1)
    label = Label(frame0, text='NamingSystem:')
    label.grid(row=0,column=3,sticky='e')
    tipText = 'Naming system for reference experiment systematic name'
    self.namingSystemPulldown = PulldownList(frame0, texts=namingSystemList,
                                             tipText=tipText,
                                             callback=self.setNamingSystem)
    self.namingSystemPulldown.grid(row=0,column=4,sticky='e')
    
    # RefExperiment table
    row += 1
    self.grid_rowconfigure(row, weight=1)
    
    self.nameEntry = Entry(self.popup, 
     text='', returnCallback=self.setName, width=12
    )
    
    self.synonymEntry = Entry(self.popup, 
     text='', returnCallback=self.setSynonym, width=12
    )
    
    self.systematicNameEntry = Entry(self.popup, 
     text='', returnCallback=self.setSystematicName, width=16
    )
    
    self.numDimsEntry = IntEntry(self.popup, 
     text='', returnCallback=self.setNumDims, width=2
    )
    
    tipText = 'Table of reference experiments included in experiment prototype'
    frame1 = LabelFrame(self, text='RefExperiments', tipText=tipText)
    frame1.grid(row=row, column=0, sticky='nsew')
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)
    colHeadings      = ['#', 'Name', 'numDims', 'isReversed', 'Synonym', 
                        'Systematic Name','isEditable']
    editWidgets      = [None, self.nameEntry, self.numDimsEntry, None, 
                        self.synonymEntry, self.systematicNameEntry, None]
    editGetCallbacks = [None, self.getName, self.getNumDims, 
                        self.toggleIsReversed, self.getSynonym, 
                        self.getSystematicName, None]
    editSetCallbacks = [None, self.setName, self.setNumDims, None, 
                        self.setSynonym, self.setSystematicName, None]
    tipTexts = [
     'The serial number of the reference experiment',
     'CCPN systematic name for the reference experiment. '
      'Conforms to the Experiment Prototype name.',
     'Number of dimensions of reference experiment',
     'Is magnetisation transfer flow reversed'
      ' relative to the experiment prototype specification?',
     'Human-understandable reference experiment synonym',
     'Reference experiment systematic name in selected namnig system',
     'Can reference experiment be modified? True for newly created experiments',
    ]
    self.refExperimentMatrix = ScrolledMatrix(frame1,
     editSetCallbacks=editSetCallbacks, multiSelect=0,
     editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,
     initialRows=3, headingList=colHeadings, tipTexts=tipTexts,
     callback=self.selectRefExperiment)
     
    self.refExperimentMatrix.grid(row=0, column=0, sticky='nsew')
    self.refExperimentMatrix.doEditMarkExtraRules = self.refExperimentEditMarkExtraRules
    
    texts    = ['Create New','Create Copy', 'Freeze', 'Delete']
    commands = [
     self.newRefExperiment, self.newCopyRefExperiment,
     self.freezeRefExperiment, self.deleteRefExperiment
    ]
    tipTexts = [
     'Make a new blank reference experiment . Information is added subsequently',
     'Make a new reference experiment '
      'that is a copy of the currently selected one',
     'Set selected reference experiment to be no longer modifiable',
     'Delete the selected reference experiment',
    ]
    self.refExperimentButtons = ButtonList(
     frame1, texts=texts,commands=commands, tipTexts=tipTexts
    )
    self.refExperimentButtons.grid(row=1, column=0, sticky='nsew')

    # RefExpDimRef table
    row += 1
    
    self.constantTimePulldown = PulldownList(self.popup, texts=constantTimeList,
                                             callback=self.setConstantTime)
    
    self.expMeasurementPulldown = PulldownList(self.popup, callback=self.setExpMeasurement)
    self.isotopeEntry = MultiWidget(self.popup, CheckButton, minRows=0,
     callback=self.setCoupledIsotopes
    )
    self.expStepEntry = MultiWidget(self.popup, PulldownMenu, minRows=0,
     callback=self.setExpSteps
    )
    self.scalingFactorsEntry = MultiWidget(self.popup, FloatEntry, minRows=0,
     callback=self.setValidScalingFactors
    )
    self.groupingIdEntry = IntEntry(self.popup,
     text='1', returnCallback=self.setGroupingId, width=3
    )
    
    tipText = 'Table of Reference experiment dimensions and dimension referencings'
    self.grid_rowconfigure(row, weight=1)
    frame2 = LabelFrame(self, text='RefExpDimRefs', tipText=tipText)
    frame2.grid(row=row, column=0, sticky='nsew')
    frame2.grid_rowconfigure(0, weight=1)
    frame2.grid_columnconfigure(0, weight=1)
    colHeadings      = [
     'Dim','#', 'expMeasurement', 'valid scaling\nfactors', 'grouping\nnumber',
     'Coupled\nIsotopes', 'Constant\nTime','expSteps'
    ]
    editWidgets      = [
     None, None, self.expMeasurementPulldown, 
     self.scalingFactorsEntry, self.groupingIdEntry, self.isotopeEntry, 
     self.constantTimePulldown, self.expStepEntry
    ]
    editGetCallbacks = [
     None, None, self.getExpMeasurement, self.getValidScalingFactors,
     self.getGroupingId, self.getCoupledIsotopes, self.getConstantTime, 
     self.getExpSteps
    ]
    editSetCallbacks = [
     None, None, self.setExpMeasurement, None, 
     self.setGroupingId, None, self.setConstantTime, None
    ]
    tipTexts = [
     'Dimension number for reference experiment dimension',
     'Serial number for dimension referencing',
     'Experimental measurement associated with dimension referencing',
     'Scaling factors allowed when converting from observed frequencies to measurement values',
     'Number that classifies dimension referencings from a dimension into consistent groups. '
      'E.g.: Group 1: 13C shift and J(HC); group 2: 15N shift and J(HN)',
     'Isotopes to which couplings are active and observable. '
      'Alternative to specifying couplings as experimental measurements',
     'Acquisition mode: Constant, variable or mixed time.',
     'Experiment steps used for measuring the measurement (in cases where it matters)',
    ]
    self.refExpDimRefMatrix = ScrolledMatrix(frame2,
     editSetCallbacks=editSetCallbacks, multiSelect=0,
     editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,
     initialRows=6, headingList=colHeadings, tipTexts=tipTexts,
     callback=self.selectRefExpDimRef)
     
    self.refExpDimRefMatrix.grid(row=0, column=0, sticky='nsew')

    self.refExpDimRefMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules

    texts    = ['New RefExpDim', 'New RefExpDimRef', 'Delete RefExpDim', 'Delete RefExpDimRef']
    commands = [
     self.newRefExpDim, self.newRefExpDimRef, self.deleteRefExpDim, self.deleteRefExpDimRef
    ]
    tipTexts = [
     'Add new reference experiment dimension',
     'Add new dimension referencing'
      ' to the same dimension as the referencing currently selected',
     'Delete selected reference experiment dimension and all dimension referencings',
     'Delete selected dimension referencing',
    ]
    self.refExpDimRefButtons = ButtonList(
     frame2, texts=texts,commands=commands, tipTexts=tipTexts,
    )
    self.refExpDimRefButtons.grid(row=1, column=0, sticky='nsew')

    row += 1
    texts = ['Purge RefExpDims']
    tipTexts = [
     'Delete reference experiment dimensions that have no dimension '
      'referencings. Used for removing orphan dimensions after deletions elsewhere']
    commands = [self.purgeRefExpDims]
    self.bottomButtons = ButtonList(self,texts=texts, commands=commands, 
                                    tipTexts=tipTexts, expands=0)
    self.bottomButtons.grid(row=row, column=0, columnspan=2, sticky='ew')

    self.updateAfter()

    self.administerNotifiers(Implementation.registerNotify)
  
  def administerNotifiers(self, notifyFunc):
    
    for func in ('delete', 'setName', 'setSynonym'):
      for className in ('ccp.nmr.NmrExpPrototype.AtomSite',):
        notifyFunc(self.updateAfter, className, func)
    for func in ('delete','setMeasurementType'):
      for className in ('ccp.nmr.NmrExpPrototype.ExpMeasurement',):
        notifyFunc(self.updateAfter, className, func)
    for func in ('__init__', 'delete',''):
      for className in ('ccp.nmr.NmrExpPrototype.RefExperiment',
                        'ccp.nmr.NmrExpPrototype.RefExpDim',
                        'ccp.nmr.NmrExpPrototype.RefExpDimRef',
                        'ccp.nmr.NmrExpPrototype.SystematicName'
      ):
        notifyFunc(self.updateAfter, className, func)
    notifyFunc(
     self.updateAfter, 'ccp.nmr.NmrExpPrototype.ExpMeasurement', 'delete'
    )

  def destroy(self):
  
    self.administerNotifiers(Implementation.unregisterNotify)
    Frame.destroy(self)

  
  # Tab level
  
  def setNamingSystem(self, name):
    
    if name != self.namingSystem:
      self.namingSystem = name
      self.updateRefExperiment()
    
  
  def purgeRefExpDims(self):
   
    nmrExpPrototype = self.nmrExpPrototype
    if nmrExpPrototype is not None:
   
      doUpdate = False
      for obj in nmrExpPrototype.sortedRefExperiments():
        for rxd in obj.sortedRefExpDims():
          if not rxd.refExpDimRefs:
            doUpdate = True
            self.deleteRefExpDim(rxd)
         
      if doUpdate:
        self.popup.updateAfter()
      
  
  # RefExperiment Table
  
  def refExperimentEditMarkExtraRules(self, refExperiment, row, col):
    
    if not refExperiment or not hasattr(refExperiment,'isEditable'):
      return False
    
    if (col == 2) and refExperiment:
      if refExperiment.refExpDims:
        return False 
    return True 
    
  def setName(self, event):

    text = self.nameEntry.get().strip()
    
    obj = self.refExperiment
    if obj:
      if text and text.strip():
        obj.setName(text)
      else:
        text = None

  def getName(self, refExperiment):
    
    if refExperiment:
      self.nameEntry.set(refExperiment.name)
 
  def setNumDims(self, *event):
    
    refExperiment = self.refExperiment
    
    if not self.nmrExpPrototype:
      showWarning('Failure','No measurements set for prototype', parent=self)
      return
    
    nn = self.numDimsEntry.get()
    if nn is not None and refExperiment:
      expMeasurements = self.nmrExpPrototype.sortedExpMeasurements()
      nn = min(nn, len(expMeasurements))
      for ii in range(nn):
        if not refExperiment.findFirstRefExpDim(dim=ii+1):
          newDim = refExperiment.newRefExpDim(dim=ii+1)
          obj = newDim.newRefExpDimRef(expMeasurement=expMeasurements[ii])
          self.refExpDimRef = obj
      self.updateRefExpDimRefAfter()

  def getNumDims(self, refExperiment):
    
    if refExperiment :
      self.numDimsEntry.set(len(refExperiment.refExpDims))
  
  def toggleIsReversed(self, refExperiment):
    """ 
    """
    refExperiment.isReversed = not refExperiment.isReversed
    #self.updateAfter()

 
  def setSynonym(self, event):
    
    text = self.synonymEntry.get()
    obj = self.refExperiment
    if obj:
      if text and text.strip():
        obj.setSynonym(text)
      else:
        obj.setSynonym(None)

  def getSynonym(self, refExperiment):
    
    if refExperiment :
      self.synonymEntry.set(refExperiment.synonym)

 
  def setSystematicName(self, event):

    text = self.systematicNameEntry.get()
    obj = self.refExperiment
    if obj:
      curSysNameObj = obj.findFirstSystematicName(
                           namingSystem=self.namingSystem)
      if text and text.strip():
        if curSysNameObj is None:
          obj.newSystematicName(namingSystem=self.namingSystem,
                                name=text)
        else:
          curSysNameObj.name = text
          
      elif curSysNameObj is not None:
        curSysNameObj.delete()
        

  def getSystematicName(self, refExperiment):
    
    if refExperiment :
      sysNameObj = refExperiment.findFirstSystematicName(
                                  namingSystem=self.namingSystem)
      if sysNameObj is None:
        sysName = ''
      else:
        sysName = sysNameObj.name
      self.systematicNameEntry.set(sysName)
  
  def newRefExperiment(self):
  
    obj = self.nmrExpPrototype.newRefExperiment(name='dummy')
    obj.isEditable = True
    self.refExperimentMatrix.currentObject = self.refExperiment = obj
    self.updateAfter()
          
  def newCopyRefExperiment(self):
    
    current = self.refExperiment
    obj = genUtil.copySubTree(
     current, self.nmrExpPrototype,
     topObjectParameters={'name':'copy_of_'+current.name}
    )
    obj.isEditable = True
    self.refExperimentMatrix.currentObject = self.refExperiment = obj
    self.updateAfter()
  
  def freezeRefExperiment(self):
    if (showOkCancel(
     'Freeze RefExperiment',
     "Are you sure you want to freeze the RefExperiment,\nso it cannot be changed?",
    parent=self)):
      obj = self.refExperiment
      if hasattr(obj,'isEditable'):
        del obj.isEditable
      self.updateAfter()
          
  def deleteRefExperiment(self):
    obj = self.refExperiment
    if obj:
      obj.delete()
      self.refExperiment = None
     
  def selectRefExperiment(self, obj, row, col):
    if obj and obj is not self.refExperiment:
      self.refExperiment = obj
      self.refExpDimRef = None
      setEditable = (hasattr(obj, 'isEditable') and obj.isEditable)
      if setEditable != self.experimentEditable:
        self.experimentEditable = setEditable
      self.updateRefExpDimRefAfter()
      self.updateButtons()  
      self.updateRefExpDimRefButtons()
  
  def updateRefExperiment(self):
    
    # select single refExperiment
    if self.nmrExpPrototype:
      refExperiments = self.nmrExpPrototype.refExperiments
      if len(refExperiments) == 1:
        obj = iter(refExperiments).next()
        self.refExperimentMatrix.currentObject = self.refExperiment = obj
        self.isEditable = (hasattr(obj, 'isEditable') and obj.isEditable)
    
    numCols = 6
    if self.nmrExpPrototype and not self.nmrExpPrototype.isDeleted:
      objectList = self.nmrExpPrototype.sortedRefExperiments()
      textMatrix = len(objectList)*[None]
    else:
      objectList = []
      textMatrix = []
      
    for ii,refExperiment in enumerate(objectList):
      sysNameObj = refExperiment.findFirstSystematicName(
                                  namingSystem=self.namingSystem)
      if sysNameObj is None:
        sysName = ''
      else:
        sysName = sysNameObj.name
      textMatrix[ii] = (refExperiment.serial,
                        refExperiment.name,
                        len(refExperiment.refExpDims),
                        refExperiment.isReversed and 'Yes' or 'No',
                        refExperiment.synonym,
                        sysName,
                        hasattr(refExperiment,'isEditable') and 'Yes' or 'No',)
    
    obj = self.refExperiment
    if obj:
      self.experimentEditable = (hasattr(obj, 'isEditable') and obj.isEditable)
      self.refExperimentMatrix.currentObject = obj
    else:
      self.experimentEditable = False
 
    self.refExperimentMatrix.update(objectList=objectList, textMatrix=textMatrix)
    
    # NB this must be called last, 
    # so that the currentObject for the matrices is correctly set
    self.updateButtons()  
    self.updateRefExpDimRefButtons()
  
  
  # RefExpDimRef table
      
  def doEditMarkExtraRules(self, obj, row, col):
    if self.experimentEditable:
      return True
    else:
      return False  
  
      
  def setExpMeasurement(self, event):
  
    expMeasurement = self.expMeasurementPulldown.getObject()
    
    obj = self.refExpDimRefMatrix.currentObject
    if expMeasurement:
      obj.setExpMeasurement(expMeasurement)  
      
  
  def getExpMeasurement(self, refExpDimRef):
    expMeasurements = self.nmrExpPrototype.sortedExpMeasurements()
    entries = [
     '%s(%s)' % (x.serial,','.join([y.name for y in x.atomSites])) 
     for x in expMeasurements
    ]
    if expMeasurements:
      try:
        ii = expMeasurements.index(refExpDimRef.expMeasurement)
      except ValueError:
        ii = 0
      self.expMeasurementPulldown.setup(entries, expMeasurements, ii)
  
  def setValidScalingFactors(self, values):
    """ values is None if the widget aborts
    """
    if values is not None:
      obj = self.refExpDimRefMatrix.currentObject
      obj.validScalingFactors = values
    self.refExpDimRefMatrix.keyPressEscape()
  
  def getValidScalingFactors(self, refExpDimRef):
    """ 
    """
    if refExpDimRef:
      self.scalingFactorsEntry.setValues(refExpDimRef.validScalingFactors) 
  
  def setGroupingId(self, event):
  
    x = self.groupingIdEntry.get()
    if x:
      obj = self.refExpDimRefMatrix.currentObject
      if obj:
        obj.setGroupingId(x)
  
  def getGroupingId(self, refExpDimRef):
    
    if refExpDimRef:
      self.groupingIdEntry.set(refExpDimRef.groupingId) 
    
  def setCoupledIsotopes(self, values):
    
    if values is not None:
      length = len(values)
      assert length == len(isotopeList)
 
      newIsotopeCodes = []
      for ii in range(length):
        if values[ii]:
          newIsotopeCodes.append(isotopeList[ii])
 
      obj = self.refExpDimRefMatrix.currentObject
      obj.setCoupledIsotopeCodes(newIsotopeCodes)
    self.refExpDimRefMatrix.keyPressEscape()
  
  def getCoupledIsotopes(self, refExpDimRef):
    
    if refExpDimRef:
      currentCodes = refExpDimRef.coupledIsotopeCodes
      length = len(isotopeList)
      values = length * [False]
      for ii in range(length):
        if isotopeList[ii] in currentCodes:
          values[ii] = True
      
      self.isotopeEntry.set(options=isotopeList, values=values)

  def setConstantTime(self,event):
  
    name = self.constantTimePulldown.getText()
    
    obj = self.refExpDimRefMatrix.currentObject
    if obj:
      obj.setConstantTime(name)
  
  def getConstantTime(self, refExpDimRef):
    
    if refExpDimRef:
      self.constantTimePulldown.set(refExpDimRef.constantTime)
    
  def setExpSteps(self, values):
    """
    """
    
    if values is not None:
      obj = self.refExpDimRefMatrix.currentObject
      stepdict = {}
      for step in obj.expMeasurement.expSteps:
        stepdict[repr((step.expGraph.serial,step.stepNumber))] = step
      
      newSteps = [stepdict[x] for x in values]
      obj.expSteps = newSteps
    self.refExpDimRefMatrix.keyPressEscape()
  
  def getExpSteps(self, refExpDimRef):
    """
    """
    if refExpDimRef:
      
      steps = refExpDimRef.expMeasurement.expSteps
      currentKeys = set((x.expGraph.serial,x.stepNumber) 
                        for x in refExpDimRef.expSteps)
      options = []
      values = []
      keys = [(x.expGraph.serial,x.stepNumber) for x in steps]
      keys.sort()
      for tt in keys:
        ss = repr(tt)
        options.append(ss)
        if tt in currentKeys:
          values.append(ss)
      
      self.expStepEntry.set(options=options, values=values)
  
  def newRefExpDim(self, expMeasurement=None):
  
    refExperiment = self.refExperiment
    
    if expMeasurement is None:
      expMeasurement = NmrExpPrototypeUtil.nextFreeExpMeasurement(
                           self.nmrExpPrototype, 'refExpDimRefs')
                           
      if expMeasurement is None:
        showWarning('Failure','No measurements set in prototype', parent=self)
        return
    
    dim = len(refExperiment.refExpDims) + 1
    refExpDim = refExperiment.newRefExpDim(dim=dim)
    obj = refExpDim.newRefExpDimRef(expMeasurement=expMeasurement)
    self.refExpDimRefMatrix.currentObject = obj
    self.updateRefExpDimRef()

  def newRefExpDimRef(self, expMeasurement=None):
  
    if expMeasurement is None:
      expMeasurement = NmrExpPrototypeUtil.nextFreeExpMeasurement(
                           self.nmrExpPrototype, 'refExpDimRefs')
                           
      if expMeasurement is None:
        showWarning('Failure','No measurements set in prototype', parent=self)
        return
  
    refExperiment = self.refExperiment
    
    # get refExpDim:
    refExpDimRef = self.refExpDimRefMatrix.currentObject
    
    if refExpDimRef:
      refExpDim = refExpDimRef.refExpDim

      obj = refExpDim.newRefExpDimRef(expMeasurement=expMeasurement)
      self.refExpDimRefMatrix.currentObject = obj
      self.updateRefExpDimRef()
          
  def deleteRefExpDim(self, refExpDim=None):
    """Rasmus Fogh 27/3/08
    NB this is a change in implementation.
    Previously the code deleted the *last* RefExpDim and moved the actual
    RefExpDimRef from i+1 to i to keep the old RefExpDimRefs anyway.
    Now we just delete *this* RefExpDim and move down the values of dim to fit.
    Both are hacks. The new one is simpler and, hopefully, equivalent.
    
    NBNB HACK WARNING
    This code is extremely implementation-dependent.
    Because of one-way links to it, any code that deletes RefExpDims,
    even without hacking, may break existing links. Should only be used for
    pristine RefExperiment copies with no links to them.
    """
    if refExpDim is None:
      if self.refExpDimRef is not None:
        refExpDim = self.refExpDimRef.refExpDim
    if refExpDim and refExpDim.isDeleted:
      return
      
    refExperiment = refExpDim.refExperiment
    parDict = refExperiment.__dict__['refExpDims']
    useDim = refExpDim.dim
    
    refExpDim.delete()
    self.refExpDimRef = None 
    
    for xx in refExperiment.sortedRefExpDims():
      curDim = xx.dim
      if curDim > useDim:
        # this object is after the deleted one in the list
        # give it the dim of one below it.
        del parDict[curDim]
        parDict[useDim] = xx
        xx.__dict__['dim'] = useDim
        useDim = curDim
     
    self.updateRefExpDimRef()
          
  def deleteRefExpDimRef(self):
  
    refExpDimRef = self.refExpDimRefMatrix.currentObject
    if refExpDimRef:
      if len(refExpDimRef.refExpDim.refExpDimRefs) == 1:
        self.deleteRefExpDim()
      else:
        refExpDimRef.delete()
        self.refExpDimRef = None
        
      self.updateRefExpDimRef()
     
  def selectRefExpDimRef(self, obj, row, col):
    if obj and obj is not self.refExpDimRef:
      self.refExpDimRef = obj
      self.updateRefExpDimRef()
      self.updateRefExpDimRefButtons()
      
  def updateRefExpDimRef(self):
  
    numCols = 8
    
    objectList = []
    textMatrix = []
    obj = self.refExperiment
    if obj:
      for refExpDim in obj.sortedRefExpDims():
        objectList.extend(refExpDim.refExpDimRefs)
      textMatrix = len(objectList)*[None]
    
      # inverse order for reversed experiments
      if obj.isReversed:
        objectList.reverse()
    
      for ii,refExpDimRef in enumerate(objectList):
        ll = textMatrix[ii] = numCols*[None]
        ll[0] = refExpDimRef.refExpDim.dim
        ll[1] = refExpDimRef.serial
        x = refExpDimRef.expMeasurement
        ll[2] = "%s %s(%s)" % (
         x.serial, x.measurementType, ','.join([y.name for y in x.atomSites])
        )
        ll[3] = refExpDimRef.validScalingFactors
        ll[4] = refExpDimRef.groupingId
        ll[5] = ', '.join(refExpDimRef.coupledIsotopeCodes)
        ll[6] = refExpDimRef.constantTime
        ll[7] = ', '.join([
         `(x.expGraph.serial,x.stepNumber)`
          for x in refExpDimRef.expSteps
         ])
     
    self.refExpDimRefMatrix.update(
     objectList=objectList, textMatrix=textMatrix
    )
    #
    self.waiting = False
    
  
  def updateRefExpDimRefButtons(self):
      
    obj = self.refExperiment
      
    if self.experimentEditable:
    
      self.refExpDimRefButtons.buttons[0].enable()
      
      curRefExpDimRef = self.refExpDimRef
      if curRefExpDimRef:
        self.refExpDimRefButtons.buttons[1].enable()
        self.refExpDimRefButtons.buttons[2].enable()
        
        if len(curRefExpDimRef.refExpDim.refExpDimRefs) > 1:
          self.refExpDimRefButtons.buttons[3].enable()
        else:
          self.refExpDimRefButtons.buttons[3].disable()
          
      else:
        for ii in (1,2,3):
          self.refExpDimRefButtons.buttons[ii].disable()
    
    else:
      for ii in (0,1,2,3):
        self.refExpDimRefButtons.buttons[ii].disable()
  
    
  def updateRefExpDimRefAfter(self):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.updateRefExpDimRef)
    
    
  def updateAfter(self, obj=None):
    
    if obj is not None and obj.isDeleted:
      ss = obj.qualifiedName
      if ss == 'ccp.nmr.NmrExpPrototype.NmrExpPrototype':
        self.nmrExpPrototype  = None

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
 
  def updateButtons(self):
    
    if self.nmrExpPrototype and self.nmrExpPrototype.isDeleted:
      self.nmrExpPrototype = None
    
    if self.nmrExpPrototype:
    
      self.refExperimentButtons.buttons[0].enable()
      
      obj = self.refExperiment
      
      if obj:
        self.refExperimentButtons.buttons[1].enable()
      else:
        self.refExperimentButtons.buttons[1].disable()
      
      if obj and [x for x in obj.refExpDims if not x.refExpDimRefs]:
        self.bottomButtons.buttons[0].enable()
      else:
        self.bottomButtons.buttons[0].disable()
        
      if self.experimentEditable:
      
        for ii in (2,3):
          self.refExperimentButtons.buttons[ii].enable()
      
      else:
        self.refExperimentButtons.buttons[2].disable()
        self.refExperimentButtons.buttons[3].disable()
    else:
        for ii in (0,1,2):
          self.refExperimentButtons.buttons[ii].disable()
        self.bottomButtons.buttons[0].disable()
  
  def update(self):
    """  Fill table with values, initial or after change
    """
    refExperiment = self.refExperiment
    if (refExperiment is None 
        or refExperiment.nmrExpPrototype is not self.nmrExpPrototype):
      self.refExperiment = self.nmrExpPrototype.findFirstRefExperiment()
  
    self.updateRefExperiment()
    self.updateRefExpDimRef() 
    
    
################################################################

if __name__ == '__main__':

  import Tkinter
  from memops.api.Implementation import MemopsRoot
  
  project = MemopsRoot(name='ccpn', currentUserId='editor')
  
  root = Tkinter.Tk()
  root.withdraw()

  NmrExpPrototypePopup(root, project)

  #root.mainloop()
