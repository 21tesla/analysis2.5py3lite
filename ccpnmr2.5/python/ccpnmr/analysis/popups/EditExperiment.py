
"""
======================COPYRIGHT/LICENSE START==========================

EditExperiment.py: Part of the CcpNmr Analysis program

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
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.DataEntry import askInteger, askString
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelDivider import LabelDivider
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showYesNo, showOkCancel, showWarning
from memops.gui.MultiWidget import MultiWidget
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame


from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import newShiftList, isSpectrum, setExperimentShiftList
from ccpnmr.analysis.core.ExperimentBasic import getPrimaryExpDimRef, getRefExperimentCategories
from ccpnmr.analysis.core.ExperimentBasic import getFilteredRefExperiments, getRefExperiments, setRefExperiment, initExpTransfers
from ccpnmr.analysis.core.ExperimentBasic import cloneExperiment, defaultExperiment
from ccpnmr.analysis.core.AssignmentBasic import getShiftLists
from ccpnmr.analysis.core.MoleculeBasic import STANDARD_ISOTOPES

from ccp.api.nmr import Nmr

VOLUME_UNITS   = ['ul','ml','l']
SAMPLE_STATES  = ['liquid', 'solid', 'ordered', 'powder', 'crystal']

SHIFT_REF_COMPOUNDS = ['DSS','TSP','TMS','HDO','CHCl3','DMSO','DMSO-d5',
                 'dioxane','CH3OH','acetone','acetone-d5','TFE','TFA',
                 'H3PO4','TMP','PCr','liquid NH3']                 

STANDARD_UNITS = ['ppm','ppt','ppb']

NMR_PROBE_TYPES = ['liquid', 'solid', 'nano', 'flow', 'MAS']

class EditExperimentPopup(BasePopup):
  """
  **Curate Experiment Parameters**
  
  This popup window is used to control parameters that relate to the NMR
  experiment entities within a CCPN project. It should be noted that an
  experiment is specifically a record of *what was done*, not a reference to the
  data that was obtained; this is what Analysis refers to as a spectrum. Also,
  an experiment doesn't just refer to what kind of NMR experiment was performed,
  but rather to something that was done on a particular occasion to a particular
  sample. The kinds of information that are stored at the experiment level
  describe how the experiment was performed and how certain data that derives
  from the experiment should be interpreted. Several tabs are used to sub-divide
  the popup window into  several sections to control different aspects of the
  experiments. 

  **Main "Experiments" Tab**

  This table is the main display of all if the NMR experiments described within
  the CCPN project. Each experiment may refer to one or more spectra (and even
  an FID) that resulted from the experimental operation; the names of these are
  listed in the "Data Sources" column. Several experimental parameters may be
  adjusted in this table and, as far as resonance assignment is concerned, the
  most important of these are the "Shift List" and "Mol Systems". It should
  be noted that the "Shift Ref" column can only be configured if chemical
  shift reference details have been entered at the "Shift References" tab.
  
  The shift list of an experiment states how chemical shift measurements should
  be grouped. An experiment's assignments, on the peak of its spectra, only
  contribute to the chemical shift measurenemts in one shift list. In normal
  operation each shift list corresponds to a given set of conditions, where
  resonance positions in spectra a fairly static. Different experiments may use
  different shift lists if their sample conditions are different enough to cause
  peaks to move. Accordingly, a resonance derived from a given atom may have
  several different recorded chemical shift values, each of which resides in a
  different chemical shift list. Because each experiment is associated with a
  single shift list it is thus known which chemical shift average the
  assignments in its spectra contribute to and which chemical shift values to
  look at when suggesting assignments for peaks. The shift list that an
  experiment is linked to may be changed at any time, and when an experiment is
  moved from one shift list to another (which may be new and empty) the
  contributions that its spectrum peaks make to the calculation of average
  chemical shift values will automatically be adjusted; the experiments that are
  linked to a given shift list dictate which peaks to average chemical shift
  values over. 

  The "Mol Systems" for an experiment specify which molecular systems were
  present in the sample. In essence this means which group of molecular chains
  the spectrum peaks, that come from the experiment, could be assigned to (caused
  by). In normal operation the molecular information need not specifically be
  set in this table, because the connection is made automatically when a peak is
  assigned to a resonance that is known to be from a specific atom in a specific
  residue. Once an experiment is associated with a particular molecular system,
  subsequent attempts to assign its peaks to atoms in a different molecular system will
  trigger a warning. The molecular system for an experiment may nonetheless be
  set via this table. Sometimes this is to preemptively associate particular
  experiments, and hence spectra, with specific molecular systems so that there
  is less chance of accidentally assigning a peak to the wrong thing. The
  molecular system, and hence residues, that an experiment is linked to is used
  to facilitate several operations in Analysis. For example, when operating on a
  peak to associate one of its assigned resonances with an atom, the molecular
  system link allows the correct group of chains to be displayed
  automatically. 

  **Experiment Types**
  
  This table is used to specify what kind of NMR experiment was done for each of
  the experiment records within the project. The general idea is that the basic
  magnetisation transfer pathway is given and how this relates to the
  experimental dimensions. Setting such information helps facilitate several
  operations in Analysis. For example when making NMR derived distance
  restraints only peaks from "through-space" experiments like NOESY are listed.
  Also, and perhaps most importantly, the linking of experimental dimensions to
  references, which correspond to particular points of the magnetisation transfer
  pathway, enables the assignment system to have knowledge of experimental
  dimensions that are linked together via covalent "one-bond" connectivity.
  Thus, if an experiment is set to be of 15N HSQC type, then it is known that
  any spectrum peaks represent one-bond correlations between hydrogen and
  nitrogen. This dictates what the isotope type of any assignments must be, and
  if a peak dimension is assigned to a specific atom then the assignment on the
  other peak dimension must be to the covalently bound atom.
  
  Setting the experiment type means setting the 'Full Type' column. The many
  possibilities and difficult names makes this hard to do directly, and we have
  tried to help with the process. The 'Type Synonym' column also sets the 'Full
  Type' column but gives alternative and more human-readable names like 'HNCA'
  instead of 'H[N[CA]]' or '15N HSQC-NOESY' instead of 'H[N]_H.NOESY'. Some
  'Type Synonyms' correspond to more than one Full Type; anything in the 'Alt
  Types' is an alternative Full Type. The possibilities given for these columns
  are filtered, so that experiments that do not fit the nuclei on the axes are
  not shown. The most common experiments are shown at the top of the selection
  menu, the rest can be found in one or more of the categories lower down. The
  'Categories' column divides experiments into a few rough groups and are used
  to filter the possibilities in the other columns. Experiments that could
  belong to more than one group are listed only in the most 'interesting' one.
  The 'use external' category is a special case; it means that the program
  should use the information in the 'External Source' and 'External Name'
  columns to set the Full Type. This option is selected automatically when you
  are loading a spectrum that contains information about the Full Type. At the
  moment the program understands only Bruker pulse program names, and only if
  they follow the standard naming conventions used by the Bruker Applications
  department.

  In normal operation the user specifies the type of NMR experiment that was run
  by selecting options in the "Category", "Type Synonym" and "Full Type" columns.
  If the full CCPN name of the experiment type is known then the user can go
  straight for the "Full Type", but selecting the category and/or synonym first
  allows the number available options to be reduced dramatically; without this
  all possible experiment types that have matching isotopes are shown.
  Setting the category for the NMR experiment gives a rough sub division between
  through-bond, through-space, quantification and other types. Strictly speaking
  an experiment may belong to more than one category, but in this system it is 
  only listed in the least populous. For example a 15N HSQC-NOESY has both
  through-bond and though-space transfers but is categorised as through-space.
  If the category for an experiment is unknown, or not particularly helpful, the
  user may set the synonym in the first instance. The "synonym" of an
  experimental type in Analysis is a common human-readable name, for example
  "HNCA" or "15N HSQC NOESY", but this may still not be sufficient to fully 
  specify the exact NMR experiment that was run. To do this the full CCPN type
  should  be considered.  The External Source and corresponding name columns are
  only used in situations where the loading of a spectrum specifies what kind of
  experiment was run. At present this only occurs for data loaded from Bruker
  files, and then only if the pulse sequence name in the parameters is known to
  the system. Nevertheless, if this data is present the experiment type
  information can be automatically be filled in.  

  The full CCPN type for an experiment uses a special nomenclature that is
  described in the CCPN `experiment nomenclature paper`_ (slightly out of date
  now). In essence the user can distinguish between different magnetisation
  transfer pathways, some of which may have the same common name (synonym). For
  example a 15N HSQC-TOCSY could have either the HSQC step or the TOCSY step
  first. In this instance the system offers a choice between H[N]_H.TOCSY (HSQC
  first) and H_H[N].TOCSY (TOCSY first). The experiment naming system for the
  full CCPN type is fairly complex, and is designed to give a precise
  specification of the magnetisation steps, which atom sites they visit and what
  measurements are made; giving rise to an experimental dimension. It should be
  noted however, that this system does not describe the precise NMR pulse
  sequence that was used. For example no distinction is made between HSQC and
  HMQC. The essential features of the nomenclature are as follows: capital
  letters indicate atom sites that were recorded and result in an experimental
  dimension; lower case letters are atom sites that are part of the pathway but
  not recorded, e.g. carbonyl in H[N[co[CA]]]; square brackets represent
  out-and-back transfers; curly brackets with "|" sub-divisions represent
  alternative pathways; underscores represent transfers that are not one-bond or
  J-coupling (the transfer type is listed after at the end after a dot).

  The lower tables are used to show how the dimensions of the actual experiment
  relate to the reference dimensions that described in the experiment type. Even
  when an experiment type is set it will not always be possible to 
  automatically determine which of the experimental dimensions relates to which
  part of the magnetisation transfer pathway. For example a 3D HCCH TOCSY
  experiment (full type HC_cH.TOCSY) has two hydrogen dimensions; one dimension
  is the hydrogen bound to the measured carbon, and one dimension is the
  hydrogen in the acquisition dimension. Deciding which is which is crucial for
  the correct assignment and interpretation of spectra. The problem only arises
  when there are two or more dimensions with the same nucleus. Sometimes
  Analysis guesses wrong and the user has to check. Changing the dimension
  mapping is a matter of looking in the lower left table and seeing how one
  dimension relates to another. Each dimension that has a direct transfer to
  another recorded dimension is listed. For example, in an HCCH TOCSY dimension
  1 (hydrogen) might be 'onebond' to dimension 3 (carbon), but the user may know
  that it is actually dimension 2 that is really 'onebond' to the carbon. This
  problem may be fixed by double-clicking either the "First Dim" or the "Second
  Dim to to change the transfer pathways so that dimension 3 (carbon) is listed
  as 'onebond' to dimension 2 (hydrogen) (the other rows in the table will
  adjust automatically). The Numbering of dimensions in this table is the same
  as that presented when assigning a peak, i.e. in the `Assignment Panel`_. It
  helps to know that dimension '1' is usually (but not always) the acquisition
  dimension.

  The lower right "Reference Dimension Mapping" is an alternative way of looking
  at the same information and shows how the experimental dimensions have been
  mapped to their reference counterpart in the experiment type database. Here,
  the "Ref Measurement" column can be used to follow the steps in the
  magnetisation transfer pathway by following increasing measurement numbers.
  Changing the "Ref Exp Dim" column in this table is equivalent to making
  changes in the lower left table, but is perhaps more difficult to understand.

  **Experimental Details, Instruments and Shift References**
  
  The "Experimental Details" table is used to list and edit details about the
  recorded experiments in terms of its physical setup. The user may specify
  which instruments  were used and information about the sample and NMR tube. It
  should be noted that in order to specify a spectrometer or probe the
  specification for the instrument must first be entered in the "NMR
  Instruments" tab. Currently, none of the NMR details given in this table have
  any influence on resonance assignment or NMR data analysis, although spinning
  information may be used for solid-state spectra at some point. However, any
  experimental details entered into the CCPN project will be present when
  submissions to the BioMagResBank database are made; initially using the
  `CcpNmr ECI`_.

  The "Shift References" table is use to enter chemical shift reference
  information into the CCPN project. This may them be linked to experiments via
  the first "Experiments" tab, and such information is required for database
  deposition. To add a chemical shift reference specification the user first
  clicks on either "Add Internal Reference" (internal to the sample) or "Add
  External Reference" as appropriate. Then for the rows that appear in the table
  the user double-clicks to edit the columns to specify: which atom in which
  kind of molecule was used, what the reference value and unit is, and whether 
  the reference is direct or indirect. A reference atom allows the direct
  referencing of resonances that have the same kind of isotope, but other
  isotopes may be referenced indirectly by using a shift ratio to relate to a
  direct reference.

  The "NMR Instruments" section contains two table that allows the user to add
  descriptions of the NMR probe and spectrometer that were used during the
  experiments. To achieve this the user adds a new specification for the 
  appropriate kind of instrument then double-clicks to fill in the details for
  each of the rows that appears in the table.

  **Caveats & Tips**

  An experiment may be linked temporarily with a new shift list; selecting
  "<New>" in the Shift List column of the first tab then reseting the shift list
  back to the original one, in order to make a shift list that contains only
  chemical shift value for that experiment at that time. Without any experiment
  links these chemical shift values will not alter as peaks and assignments 
  change.

  **References**
  
  *A nomenclature and data model to describe NMR experiments.
  Fogh RH, Vranken WF, Boucher W, Stevens TJ, Laue ED.
  J Biomol NMR. 2006 Nov;36(3):147-55* (link_)

  .. _`experiment nomenclature paper`:  http://www.ncbi.nlm.nih.gov/pubmed/17031528
  .. _`link`:  http://www.ncbi.nlm.nih.gov/pubmed/17031528
  .. _`Assignment Panel`: EditAssignmentPopup.html
  .. _`CcpNmr ECI`: EntryCompletionPopup.html
  """

  def __init__(self, parent, isModal=False, *args, **kw):

    self.guiParent = parent
    self.detailsExp = None
    self.experiment = None
    self.waiting = False
    self.waitingT = False
    self.waitingD = False
    self.waitingS = False
    self.waitingI = False
    self.isModal = isModal
    self.typeExpt = None
    self.typeExps = []
    self.shiftReference = None
    self.probe = None
    self.spectrometer = None
    self.expDimRef = None
    self.transferExpDimRefs = None
  
    BasePopup.__init__(self, parent, modal=isModal, title="Experiments : Experiments", **kw)

  def body(self, guiFrame):

    self.geometry('750x500')

    guiFrame.expandGrid(0,0)

    tipTexts = ['A table listing the NMR experiments available within the current project including details of which shift lists they use',
                'The experiment type definitions that are associated with the current NMR experiments',
                'A table of ancillary  experimental details, relating to samples and instrumentation',
                'A table of chemical shift reference information that is used by th experiments of the project',
                'Descriptions of NMR spectrometers and probes used in the experiments']
    options = ['Experiments','Experiment Types','Experimental Details',
               'Shift References','NMR Instruments']
    if self.isModal:
      options[1] = 'Specify ' + options[1]
      
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              tipTexts=tipTexts, grid=(0, 0))
    self.tabbedFrame = tabbedFrame
    frameA, frameB, frameC, frameD, frameE = tabbedFrame.frames

    if self.isModal:
      tabbedFrame.select(1)
    
    # Experiment frame
    
    frameA.grid_columnconfigure(0, weight=1)
    frameA.grid_rowconfigure(0, weight=1)

    self.detailsEntry = Entry(self,text='', returnCallback = self.setDetails, width=30)
    self.acqDimPulldown = PulldownList(self, callback=self.setAcqDim)
    self.nameEntry    = Entry(self,text='', returnCallback = self.setName, width=20)
    self.shiftListPulldown = PulldownList(self, callback=self.setShiftList)
    self.shiftRefSelect  = MultiWidget(self, CheckButton, callback=self.setShiftRefs,  minRows=0, useImages=False)
    self.molSystemSelect = MultiWidget(self, CheckButton, callback=self.setMolSystems, minRows=0, useImages=False)

    row = 0
    tipTexts = ['The serial number of the experiment',
                'A short identifying name for the experiment, for graphical display',
                'The number or independent dimensions recorded in the experiment',
                'The name of the reference type for the experiment; to indicate what kind of experiment was run',
                'Which of the recorded experimental dimensions was the final acquisition dimension',
                'Which shift list the experiment uses to record the chemical shifts for it\'s assigned peaks',
                'The names of the spectra (or other data sources) that are derived form the experiment',
                'The molecular systems (groups of chains) that assignments in the experiment are made to',
                'The chemical shift reference records used by the experiment',
                'The number of fully processed spectra resulting from the experiment',
                'Whether there is an unprocessed data source (spectrum) associated with the experiment',
                'A user-editable textual comment for the experiment']
                
    colHeadings      = ['#','Name','Num.\nDim','Exp. Type', 'Acquisition\nDim',
                        'Shift List','Data\nSources','Mol\nSystems',
                        'Shift Ref','Num.\nSpectra','Raw\nData','Details']

    editWidgets      = [None, self.nameEntry, None, None, self.acqDimPulldown,
                        self.shiftListPulldown, None,
                        self.molSystemSelect, self.shiftRefSelect,
                        None, None, self.detailsEntry]
    editGetCallbacks = [None, self.getName, None, None, self.getAcqDim,
                        self.getShiftList, None,
                        self.getMolSystems, self.getShiftRefs,
                        None, None, self.getDetails]
    editSetCallbacks = [None, self.setName, None, None, self.setAcqDim,
                        self.setShiftList, None,
                        self.setMolSystems, self.setShiftRefs,
                        None, None, self.setDetails]
    
    self.scrolledMatrix = ScrolledMatrix(frameA, tipTexts=tipTexts, 
                                         multiSelect=True, 
                                         initialRows=10, grid=(0, 0), 
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks, 
                                         editWidgets=editWidgets,
                                         headingList=colHeadings, 
                                         callback=self.selectCell,
                                         deleteFunc=self.deleteExperiment)
 
    row += 1
    tipTexts = ['Show a table of NMR measurement lists (chemical shifts, T1, J-coupling etc) that relate to the selected experiment',
                'Clone an existing experiment and spectrum (excluding peakLists)',
                'Make a new experiment and spectrum',
                'Delete the selected experiments from the project; any spectrum data on disk is preserved']
    texts    = ['Measurement Lists','Clone Experiment...','New Experiment...','Delete']
    commands = [self.showMeasurementLists, self.cloneExperiment,
                self.newExperiment, self.deleteExperiment]
    self.experimentButtons = ButtonList(frameA, commands=commands, texts=texts,
                                        tipTexts=tipTexts, grid=(row, 0))

    # Experiment Types

    frameB.grid_columnconfigure(0, weight=1)
    frameB.grid_columnconfigure(1, weight=1)
    frameB.grid_rowconfigure(0, weight=1)

    self.categoryPulldown  = PulldownList(self, texts=self.getCategories(), 
                                          callback=self.setCategory )
    self.extSourcePulldown = PulldownList(self, texts=self.getExtSources(), 
                                          callback=self.setExtSource )
    self.extNameEntry = Entry(self, width=16, returnCallback=self.setExtName)
    self.synonymPulldown   = PulldownList(self, callback=self.setSynonym  )
    self.fullTypePulldown  = PulldownList(self, callback=self.setFullType )
    #dimNums = [1,2,3]
    #self.refExpDimPulldown = PulldownList(self, texts=[str(x) for x in dimNums],
    #                                      objects=dimNums, callback=self.setRefExpDim)
    self.refExpDimPulldown = PulldownList(self, callback=self.setRefExpDim)
    self.transferPulldownA = PulldownList(self, callback=self.setTransferExpDimA)
    self.transferPulldownB = PulldownList(self, callback=self.setTransferExpDimB)
    
    row = 0
    tipTexts = ['The serial number of the experiment',
                'The name of the experiment used in graphical displays',
                'Whether the experiment type information comes from an external source, e.g. Bruker pulse programs',
                'The name of the experiment type as specified fom an external source, e.g. from the Bruker pulse program ',
                'A loose, mutually exclusive classification of the experiment to help simplify subsequent setting of experiment type',
                'A common human readable name for the type of experiment',
                'The precise CCPN name for the reference experiment type; indicating what experiment was performed',
                'Alternative reference experiment names which share the same human readable synonym']
    headingList = ['#','Experiment\nName', 'External\nSource', 'External\nName',
                   'Category','Type\nSynonym','Full\nType', 'Alt\nTypes']
    editWidgets      = [None, None, self.extSourcePulldown, self.extNameEntry, 
                        self.categoryPulldown, self.synonymPulldown, 
                        self.fullTypePulldown, None]
    editGetCallbacks = [None, None, self.getExtSource, self.getExtName, 
                        self.getCategory, self.getSynonym, self.getFullType,
                        None]
    editSetCallbacks = [None, None, self.setCategory, self.setExtSource, 
                        self.setExtName, self.setSynonym, self.setFullType,
                        None]
    
    self.expTypeMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                        callback=self.selectTypeExpt,
                                        multiSelect=True, grid=(0,0),
                                        gridSpan=(1,2), tipTexts=tipTexts,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets)
                                               
    row += 1
    tipTexts = ['Spread the experiment type specification, where possible, from the last selected experiment to all highlighted experiments',
                'Open a popup to view and administer the reference experiment type specifications; can be use to describe new experiment types']
    texts    = ['Propagate Experiment Type','Edit Experiment Prototypes']
    commands = [self.propagateExpType,self.guiParent.editExpPrototype]
    self.expTypeButtons = ButtonList(frameB, commands=commands, grid=(row,0),
                                     gridSpan=(1,2),texts=texts, tipTexts=tipTexts)

    row += 1
    frame = LabelFrame(frameB, text='Experiment Dim-Dim Transfers', grid=(row,0))
    frame.expandGrid(0,0)

    tipTexts = ['The first recorded experimental dimension involved in a magnetisation transfer',
                'The type of magnetisation transfer (e.g. NOSY, TOCY, J-coupling, one-bond) between the experimental dimensions',
                'The second recorded experimental dimension involved in a magnetisation transfer']
    headingList = ['First Dim','Transfer Type\nBetween Dims','Second Dim']
    editWidgets      = [self.transferPulldownA, None, self.transferPulldownB]
    editGetCallbacks = [self.getTransferExpDimA, None, self.getTransferExpDimB]
    editSetCallbacks = [self.setTransferExpDimA, None, self.setTransferExpDimB]
    self.transferMatrix = ScrolledMatrix(frame, headingList=headingList,
                                         callback = self.selectExpTransfer,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         grid=(0,0), tipTexts=tipTexts)

    frame = LabelFrame(frameB, text='Reference Dimension Mapping', grid=(row,1))
    frame.expandGrid(0,0)

    tipTexts = ['The number of the experimental dimension',
                'Which dimension in the reference experiment (prototype) specification the real experimental dimension corresponds to',
                'The isotope involved with the reference experimental dimension',
                'What is recorded on the experimental dimension, and which steps of the magnetisation transfer pathway it applies to']
    headingList = ['Exp\nDim','Ref Exp\nDim(.Ref)','Isotope','Ref Measurement']
    editWidgets      = [None, self.refExpDimPulldown, None, None]
    editGetCallbacks = [None, self.getRefExpDim,      None, None]
    editSetCallbacks = [None, self.setRefExpDim,      None, None]
    self.refExpDimMatrix = ScrolledMatrix(frame, headingList=headingList,
                                          callback = self.selectExpDimRef,
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks,
                                          grid=(0,0), editWidgets=editWidgets,
                                          tipTexts=tipTexts)
                                           
    if self.isModal:
      row +=1
      buttons = ButtonList(frameB, commands=[self.close,], grid=(row,0),
                           gridSpan=(1,2), texts=['Close - All Done',],
                           expands=True)
      buttons.buttons[0].config(bg='#B0FFB0')
    
    # Experimental details
    
    frameC.grid_columnconfigure(0, weight=1)
    frameC.grid_rowconfigure(0, weight=1)

    row = 0
    self.scansEntry = IntEntry(self, text=0, width=10,
                               returnCallback=self.setNumScans)

    self.spectrometerPulldown = PulldownList(self, callback=self.setSpectrometer)
    
    self.probePulldown = PulldownList(self, callback=self.setProbe)
    
    self.statePulldown = PulldownList(self, texts=SAMPLE_STATES,
                                      callback=self.setSampleState)

    self.volumeEntry = FloatEntry(self, text=0.0, width=10,
                                  returnCallback=self.setSampleVolume)
                                 
    self.unitPulldown = PulldownList(self, texts=VOLUME_UNITS, index=1,
                                     callback=self.setSampleVolUnit)

    self.tubeEntry = Entry(self, text='', width=10,
                           returnCallback=self.setNmrTube)

    self.spinRateEntry = FloatEntry(self, text=0.0, width=10,
                                    returnCallback=self.setSpinRate)

    self.spinAngleEntry = FloatEntry(self, text=0.0, width=10,
                                    returnCallback=self.setSpinAngle)
    tipTexts = ['The serial number of the experiment',
                'The textual name for the experiment, for graphical displays',
                'The specification of the spectrometer used to record the NMR experiment',
                'The specification of the NMR probe which was used in the experiment',
                'The number of repeated scans made during the experiment',
                'The state that best describes the sample and any molecular ordering; liquid (solution), solid, powder, ordered or crystalline',
                'The total volume of sample used the experiment',
                'The physical units used to describe the sample volume',
                'A description of the type of NMR tube used in the experiment',
                'If the experiment involved a spinning sample, what the angle of spin was',
                'If the experiment involved a spinning sample, what the rate of spin was, in Hz']
                                              
    colHeadings      = ['#','Name','Spectrometer','Probe',
                        'Num.\nScans','Sample\nState',
                        'Sample\nVolume','Volume\nUnit',
                        'NMR Tube\nType','Spinning\nAngle',
                        'Spinning\nRate (Hz)']
                        
    editWidgets      = [None, None,
                        self.spectrometerPulldown,self.probePulldown,
                        self.scansEntry,self.statePulldown,
                        self.volumeEntry,self.unitPulldown,
                        self.tubeEntry,self.spinAngleEntry,
                        self.spinRateEntry]
                        
    editGetCallbacks = [None, None,
                        self.getSpectrometer,self.getProbe,
                        self.getNumScans,self.getSampleState,
                        self.getSampleVolume,self.getSampleVolUnit,
                        self.getNmrTube,self.getSpinAngle,
                        self.getSpinRate]
 
    
    editSetCallbacks =  [None, None,
                         self.setSpectrometer,self.setProbe,
                         self.setNumScans,self.setSampleState,
                         self.setSampleVolume,self.setSampleVolUnit,
                         self.setNmrTube,self.setSpinAngle,
                         self.setSpinRate]
    
    self.detailsMatrix = ScrolledMatrix(frameC, multiSelect=True, 
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        headingList=colHeadings,
                                        callback=self.selectExperiment,
                                        tipTexts=tipTexts, grid=(row, 0))
 

    row += 1
    div = LabelDivider(frameC, text='Sample Conditions')
    #frame4.grid(row=row,column=0,sticky='nsew')
    #frame4.grid_columnconfigure(1, weight=1)
    #frame4.grid_rowconfigure(1, weight=1)

    label = Label(frameC, text='Current set:')
    #label.grid(row=0,column=0,sticky='w')
    self.conditionsPulldown = PulldownList(frameC, callback=self.setConditionSet)
    #self.conditionsPulldown.grid(row=0,column=1,sticky='w')

    tipTexts = ['The type of experimental condition that is described for the sample',
                'The value of the specified condition',
                'The error in the measurement of the condition value',
                'The units of measurement for the condition value and error']
                
    colHeadings      = ['Condition','Value','Error','Unit']
    editWidgets      = [None,None,None,None]
    editGetCallbacks = [None,None,None,None]
    editSetCallbacks = [None,None,None,None]
    self.conditionsMatrix = ScrolledMatrix(frameC, tipTexts=tipTexts,
                                           editSetCallbacks=editSetCallbacks,
                                           editGetCallbacks=editGetCallbacks,
                                           editWidgets=editWidgets,
                                           initialRows=5,
                                           headingList=colHeadings,
                                           callback=None)
                                           
    #self.conditionsMatrix.grid(row=1, column=0, sticky='nsew')
    texts    = ['Edit Conditions',]
    tipTexts = ['Open a table to view and edit the specification of experimental and sample conditions',]
    commands = [self.editConditions,]
    self.conditionsButtons = ButtonList(frameC,texts=texts, commands=commands, tipTexts=tipTexts)
    #self.conditionsButtons.grid(row=2, column=0, sticky='ew')

    row += 1
    tipTexts = ['Display a table to view and administer NMR spectrometer and probe specifications',
                'Spread the experimental details in the table from the last selected experiment to all highlighted experiments']
    texts    = ['Edit NMR Instruments','Propagate Experimental Details']
    commands = [self.editInstruments, self.propagateDetails]
    self.expDetailsButtons = ButtonList(frameC, texts=texts, commands=commands,
                                        tipTexts=tipTexts, grid=(row, 0), gridSpan=(1,2))
    
    # Shift references

    frameD.grid_columnconfigure(0, weight=1)
    frameD.grid_rowconfigure(0, weight=1)
    
    self.isotopePulldown = PulldownList(self, texts=STANDARD_ISOTOPES,
                                        callback=self.setShiftRefIsotope)
                                        
    self.molNamePulldown = PulldownList(self, texts=SHIFT_REF_COMPOUNDS,
                                        callback=self.setShiftRefMolName)
                                        
    self.atomGroupEntry  = Entry(self, text='', width=8, returnCallback=self.setShiftRefAtomGroup)
    self.valueEntry      = FloatEntry(self, text=0.0, width=6, returnCallback=self.setShiftRefValue)
    self.ratioEntry      = FloatEntry(self, text=0.0, width=18, formatPlaces=16,
                                      returnCallback=self.setShiftRefRatio)
    self.unitPulldown2    = PulldownList(self,  texts=STANDARD_UNITS,
                                        callback=self.setShiftRefUnit)
                                        
    self.geometryEntry   = Entry(self, text='', returnCallback=self.setShiftRefGeometry)
    self.locationEntry   = Entry(self, text='', returnCallback=self.setShiftRefLocation)
    self.axisEntry       = Entry(self, text='', returnCallback=self.setShiftRefAxis)
    
    tipTexts = ['The serial number of the chemical shift reference specification',
                'Whether the chemical shift reference is internal or external to the sample',
                'The kind of nuclear isotope to which the reference applies',
                'The number of experiments in the project which use the shift reference specification',
                'The name of the molecule used to give a reference value to chemical shifts',
                'Which atom of the reference molecule provides the reference chemical shift value',
                'The reference value of the chemical shift for the specified atom',
                'Which measurement unit the chemical shift reference value is in; ppm, ppb or ppt',
                'Whether the chemical shift referencing is direct or indirect (and thus uses a shift ratio)',
                'The precise numeric ratio used to indirectly get the reference shift value of an isotope, given the direct measurement of a different isotope',
                'For external references, a description of the geometry of the container used to hold the reference compound, e.g. cylindrical or spherical',
                'For external references, a description of the location of the reference',
                'For external references, orientation of the reference container with respect to external magnetic field, e.g. parallel or perpendicular']
    colHeadings      = ['#','Class','Isotope','Experiments',
                        'Mol. Name','Atom','Value','Unit',
                        'Ref Type','Indirect\nShift Ratio',
                        'Sample\nGeometry','Location','Axis']
                        
    editWidgets      = [None, None, self.isotopePulldown, None, self.molNamePulldown,
                        self.atomGroupEntry,self.valueEntry, self.unitPulldown2,None,
                        self.ratioEntry,self.geometryEntry,self.locationEntry,self.axisEntry]
                        
    editGetCallbacks = [None, None, self.getShiftRefIsotope, self.editExperiments, self.getShiftRefMolName,
                        self.getShiftRefAtomGroup,self.getShiftRefValue,self.getShiftRefUnit,self.toggleShiftRefType,
                        self.getShiftRefRatio ,self.getShiftRefGeometry,self.getShiftRefLocation,self.getShiftRefAxis]
                        
    editSetCallbacks = [None, None, self.setShiftRefIsotope, None, self.setShiftRefMolName,
                        self.setShiftRefAtomGroup,self.setShiftRefValue,self.setShiftRefUnit,None,
                        self.setShiftRefRatio ,self.setShiftRefGeometry,self.setShiftRefLocation,self.setShiftRefAxis]
                        
    self.shiftRefMatrix = ScrolledMatrix(frameD, multiSelect=True,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         headingList=colHeadings,
                                         callback=self.selectShiftRef,
                                         tipTexts=tipTexts, grid=(0, 0))
    self.shiftRefMatrix.doEditMarkExtraRules = self.doShiftRefEditMarkExtraRules

    tipTexts = ['Add a new record for a chemical shift reference that is internal to the sample',
                'Add a new record for a chemical shift reference that is external to the sample',
                'Delete the selected chemical shift reference records']
    texts    = ['Add Internal\nReference','Add External\nReference','Delete\nSelected']
    commands = [self.addInternalShiftRef,self.addExternalShiftRef,self.removeShiftRefs]
    self.shiftRefButtons = ButtonList(frameD, texts=texts, tipTexts=tipTexts,
                                      commands=commands, grid=(1, 0))
    
    # NMR Instruments
    
    frameE.grid_columnconfigure(0, weight=1)
    frameE.grid_rowconfigure(1, weight=1)

    row = 0
    div = LabelDivider(frameE, text='NMR Probe')
    div.grid(row=row,column=0,sticky='ew')

    row += 1
    self.probeNameEntry = Entry(self, text='', width=10,
                                returnCallback=self.setProbeName)

    self.probeTypePulldown = PulldownList(self, texts=NMR_PROBE_TYPES,
                                          callback=self.setProbeType)
    
    self.probeModelEntry = Entry(self, text='', width=10,
                                 returnCallback=self.setProbeModel)
                                 
    self.probeSerialEntry = Entry(self, text='', width=10,
                                  returnCallback=self.setProbeSerial)
                                       
    self.probeDiameterEntry = FloatEntry(self, text=0.0, width=10,
                                         returnCallback=self.setProbeDiameter)
                                   
    self.probeDetailsEntry = Entry(self, text='', width=10,
                                   returnCallback=self.setProbeDetails)

    tipTexts = ['Serial number of NMR probe specification',
                'The name of the probe for graphical representation',
                'A classification for the kind of probe used, e.g. liquid, solid, nano, flow or MAS',
                'The manufacturer\'s definition of the probe model',
                'The manufacturer\'s serial number for the specific NMR probe',
                'The probe diameter in cm',
                'A user-specified textual comment about the probe']
    colHeadings      = ['#','Name','Type','Model','Serial #','Diameter (cm)','Details']
    editWidgets      = [None,self.probeNameEntry, self.probeTypePulldown,
                        self.probeModelEntry, self.probeSerialEntry,
                        self.probeDiameterEntry, self.probeDetailsEntry]
                        
    editGetCallbacks = [None,self.getProbeName,self.getProbeType,
                        self.getProbeModel,self.getProbeSerial,
                        self.getProbeDiameter,self.getProbeDetails]
    
    editSetCallbacks = [None,self.setProbeName,self.setProbeType,
                        self.setProbeModel,self.setProbeSerial,
                        self.setProbeDiameter,self.setProbeDetails]
                        
    self.probeMatrix = ScrolledMatrix(frameE, tipTexts=tipTexts,
                                      editSetCallbacks=editSetCallbacks,
                                      editGetCallbacks=editGetCallbacks,
                                      editWidgets=editWidgets,
                                      headingList=colHeadings,
                                      callback=self.selectProbe,
                                      grid=(row, 0))

    row += 1
    tipTexts = ['Add a new NMR probe specification to the CCPN project',
                'Delete the selected NMR probe specification']
    texts    = ['New Probe Specification','Delete Probe Specification']
    commands = [self.newProbe,self.deleteProbe]
    self.probeButtons = ButtonList(frameE,texts=texts, tipTexts=tipTexts,
                                   commands=commands, grid=(row, 0))

    row += 1
    div = LabelDivider(frameE, text='Spectrometer', grid=(row, 0))

    row += 1
    self.spectrometerNameEntry = Entry(self, text='', width=10,
                                       returnCallback=self.setSpectrometerName)

    self.spectrometerFreqEntry = FloatEntry(self, text='', width=12,
                                       returnCallback=self.setSpectrometerFreq)
    
    self.spectrometerModelEntry = Entry(self, text='', width=10,
                                       returnCallback=self.setSpectrometerModel)
                                 
    self.spectrometerSerialEntry = Entry(self, text='', width=10,
                                         returnCallback=self.setSpectrometerSerial)
                                                                          
    self.spectrometerDetailsEntry = Entry(self, text='', width=10,
                                         returnCallback=self.setSpectrometerDetails)


    tipTexts = ['Serial number of the NMR spectrometer specification',
                'A name for the spectrometer, for graphical displays',
                'The rounded spectrometer frequency, from the 1H resonance frequency in MHz, used in textual description, e.g. "500", "900"',
                'The actual numeric magnetic field strength expressed as a 1H resonance frequency in MHz, e.g. "500.013"',
                'The manufacturer\'s definition of the spectrometer model',
                'The manufacturer\'s serial number for the specific NMR spectrometer',
                'A user-specified textual comment about the NMR spectrometer']
    colHeadings      = ['#','Name','Nominal Freq.','Proton Freq. (MHz)','Model','Serial #','Details']
    editWidgets      = [None,self.spectrometerNameEntry,
                        None,self.spectrometerFreqEntry,
                        self.spectrometerModelEntry, self.spectrometerSerialEntry,
                        self.spectrometerDetailsEntry]
                        
    editGetCallbacks = [None,self.getSpectrometerName,
                        None,self.getSpectrometerFreq,
                        self.getSpectrometerModel,self.getSpectrometerSerial,
                        self.getSpectrometerDetails]
    
    editSetCallbacks = [None,self.setSpectrometerName,
                        None,self.setSpectrometerFreq,
                        self.setSpectrometerModel,self.setSpectrometerSerial,
                        self.setSpectrometerDetails]
                        
    self.spectrometerMatrix = ScrolledMatrix(frameE, tipTexts=tipTexts,
                                             editSetCallbacks=editSetCallbacks,
                                             editGetCallbacks=editGetCallbacks,
                                             editWidgets=editWidgets,
                                             headingList=colHeadings,
                                             callback=self.selectSpectrometer,
                                             grid=(row, 0))

    row += 1
    tipTexts = ['Add a new NMR spectrometer specification to the CCPN project',
                'Delete the selected NMR spectrometer specification']
    texts    = ['New Spectrometer Specification','Delete Spectrometer Specification']
    commands = [self.newSpectrometer,self.deleteSpectrometer]
    self.spectrometerButtons = ButtonList(frameE,texts=texts, tipTexts=tipTexts,
                                          commands=commands, grid=(row, 0))

   # Main window

    dismissText = None
    if self.isModal:
      dismissText = 'Done'
    bottomButtons = UtilityButtonList(self.tabbedFrame.sideFrame, helpUrl=self.help_url,
                                      closeText=dismissText)
    bottomButtons.grid(row=0, column=0, sticky='e')

    self.updateShiftRefsAfter()
    self.updateExpDetailsAfter()
    self.updateExpTypesAfter()
    self.updateInstrumentsAfter()
    self.update()
  
    self.administerNotifiers(self.registerNotify)
  
  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete','setName'):
      for clazz in ('ccp.nmr.Nmr.Experiment','ccp.nmr.Nmr.DataSource'):
        notifyFunc(self.updateAfter,clazz, func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateAfter,'ccp.nmr.Nmr.ShiftList', func)
	
    for func in ('setDetails', 'setName', 'setExperimentType','setShiftList',
                 'setShiftReferences','addShiftReference','removeShiftReference',
                 'setMolSystems','addMolSystem','removeMolSystem', 'setVolumeUnit'):
      notifyFunc(self.updateAfter,'ccp.nmr.Nmr.Experiment', func)

    for clazz in ('ccp.nmr.Nmr.DataSource','ccp.nmr.Nmr.ShiftList'):
      notifyFunc(self.updateAfter, clazz, 'setName')
 
    for func in ('setIsotopeCode','setMolName'):
      for clazz in ('ccp.nmr.Nmr.ExternalShiftReference','ccp.nmr.Nmr.InternalShiftReference'):
        notifyFunc(self.updateAfter,clazz, func)
    
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.ExpDim', 'setIsAcquisition')

    # Experiment Types

    for func in ('setRefExperiment','__init__','delete','setName'):
      for clazz in ('ccp.nmr.Nmr.Experiment',):
        notifyFunc(self.updateExpTypesAfter,clazz, func)

    for func in ('setRefExpDimRef',):
      for clazz in ('ccp.nmr.Nmr.ExpDimRef',):
        notifyFunc(self.updateExpTypesAfter,clazz, func)
 
    # Experiment Details
     
    for func in ('__init__', 'delete','setName',
                  'setProbe','setSampleConditionSet',
                  'setSpectrometer',
                  'setName','setNmrTubeType',
                  'setNumScans','setSampleState',
                  'setSampleVolume','setSpinningAngle',
                  'setSpinningRate','setVolumeUnit'):
      notifyFunc(self.updateExpDetailsAfter,'ccp.nmr.Nmr.Experiment', func)

    for func in ('__init__','delete','setName'):
      notifyFunc(self.updateConditionSets,'ccp.nmr.Nmr.SampleConditionSet', func)
    
    # Shift References
      
    for func in ('__init__','delete','setAtomGroup','setIndirectShiftRatio',
                 'setValue','setIsotopeCode','setMolName','setReferenceType',
                 'setUnit','setExperiments','addExperiment','removeExperiment'):
      for clazz in ('ccp.nmr.Nmr.ExternalShiftReference','ccp.nmr.Nmr.InternalShiftReference'):
        notifyFunc(self.updateShiftRefsAfter,clazz, func)

    for func in ('setShiftReferences',
                 'addShiftReference','removeShiftReference'):
      notifyFunc(self.updateShiftRefsAfter,'ccp.nmr.Nmr.Experiment', func) 
    
    # NMR Instruments
    
    for func in ('__init__','delete','setName',
                 'setSerialNumber','setDetails',
                 'setModel','setNominalFreq'
                 'setProtonFreq','setExperiments',
                 'addExperiment','removeExperiment'):
      notifyFunc(self.updateInstrumentsAfter,'ccp.general.Instrument.NmrSpectrometer',func)

    for func in ('__init__','delete','setName',
                 'setSerialNumber','setDetails',
                 'setModel','setProbeType'
                 'setDiameter','setExperiments',
                 'addExperiment','removeExperiment'):
      notifyFunc(self.updateInstrumentsAfter, 'ccp.general.Instrument.NmrProbe', func)
    
      
  def open(self):
    
    self.updateAfter()
    self.updateExpTypesAfter()
    self.updateExpDetailsAfter()
    self.updateShiftRefsAfter()
    self.updateInstrumentsAfter()
    BasePopup.open(self)

  def close(self):
  
    if self.isModal:
      names = []
      for experiment in self.typeExps:
        if not experiment.refExperiment:
          names.append(experiment.name)
  
      if names:
        if len(names) == 1:
          msg = 'Experiment %s does not have a reference experiment type set. Try again?' \
                 % names[0]
        else:
          msg = 'Experiments %s and %s do not have reference experiment types set. Try again?' \
                 % (','.join(names[:-1]), names[-1])
      
        if showYesNo('Warning', msg, parent=self):
          return
  
    BasePopup.close(self)
  
  def editInstruments(self):
  
    self.tabbedFrame.select(4)
  
  """
  def editExperimentTypes(self):
  
    self.guiParent.editExpType()
  """
  """
  def readMdd(self):
    from gothenburg import Usf3Io
    from memops.gui.FileSelectPopup import FileSelectPopup

    popup = FileSelectPopup(self)
    file = popup.getFile()
    popup.destroy()
    
    if file:
      Usf3Io.readDataSource(self.nmrProject, file)
  """
  
  def getName(self, experiment):

    if experiment :
      width = max(20, len(experiment.name))
      self.nameEntry.config(width=width)
      self.nameEntry.set(experiment.name)

  def setName(self, event):

    text = self.nameEntry.get()
    if text and text != ' ':
      
      if text != self.experiment.name:
        if self.experiment.nmrProject.findFirstExperiment(name=text):
          showWarning('Failure','Name %s already in use' % text, parent=self)
          return
    
        self.experiment.setName( text )
  
  def getDetails(self, experiment):

    if experiment and experiment.details:
      self.detailsEntry.set(experiment.details)
  
  def setDetails(self, event):

    text = self.detailsEntry.get()
    if text and text != ' ':
      self.experiment.setDetails( text )

  def getMolSystems(self, experiment):
  
    molSystems = self.project.sortedMolSystems()
    names  = []
    values = []
    for molSystem in molSystems:
      names.append(molSystem.code)
      
      if molSystem in experiment.molSystems:
        values.append(True)
      else:
        values.append(False)  
  
    self.molSystemSelect.set(values=values,options=names)
  
  def setMolSystems(self, obj):
    
    if self.experiment:
      if obj is None:
        self.scrolledMatrix.keyPressEscape()
      else:
        molSystems = self.project.sortedMolSystems()
        values = self.molSystemSelect.get()
        selectedMolSystems = [molSystems[i] for i in range(len(values)) if values[i]]
        self.experiment.setMolSystems(selectedMolSystems)
        self.scrolledMatrix.keyPressEscape()

  def getShiftRefs(self, experiment):
  
    shiftRefs = self.nmrProject.sortedShiftReferences()
    names  = []
    values = []
    for shiftReference in shiftRefs:
      data = (shiftReference.serial,
              shiftReference.isotopeCode,
              shiftReference.molName)
      names.append('%d:%s:%s' % data)
      
      if shiftReference in experiment.shiftReferences:
        values.append(True)
      else:
        values.append(False)  
        
    self.shiftRefSelect.set(values=values,options=names)

  def setShiftRefs(self, obj):
  
    if self.experiment:
      if obj is None:
        self.scrolledMatrix.keyPressEscape()
      else:
        shiftRefs = self.nmrProject.sortedShiftReferences()
        values = self.shiftRefSelect.get()
        selectedRefs = [shiftRefs[i] for i in range(len(values)) if values[i]]
        self.experiment.setShiftReferences(selectedRefs)
        self.scrolledMatrix.keyPressEscape()

  def getShiftListNames(self, shiftLists):
  
    names = []
    for shiftList in shiftLists:
      if shiftList is None:
        names.append('<New>')
        continue
        
      if not shiftList.name:
        shiftList.name = 'ShiftList %d' % shiftList.serial
      names.append('%s [%d]' % (shiftList.name, shiftList.serial))
    
    return names
  
  def getAcqDim(self, experiment):

    expDims = experiment.sortedExpDims()
    objects = [None] + expDims
    names = ['None']
    for expDim in expDims:
      expDimRef = getPrimaryExpDimRef(expDim)
      if expDimRef and expDimRef.isotopeCodes:
        name = '%d (%s)' % (expDim.dim, ','.join(expDimRef.isotopeCodes))
      else:
        name = '%d' % expDim.dim
      names.append(name)

    expDim = experiment.findFirstExpDim(isAcquisition=True)
    index = objects.index(expDim)

    self.acqDimPulldown.setup(names, objects, index)

  def setAcqDim(self, *extra):

    if self.experiment:
      acqExpDim = self.acqDimPulldown.getObject() # could be None
      for expDim in self.experiment.sortedExpDims():
        expDim.isAcquisition = (expDim is acqExpDim)

  def getShiftList(self, experiment):

    index = 0
    shiftLists = getShiftLists(self.nmrProject) + [None,]
    names = self.getShiftListNames(shiftLists)
    shiftList = experiment.shiftList
    
    if shiftList and (shiftList in shiftLists):
      index = shiftLists.index(shiftList)

    self.shiftListPulldown.setup(names, shiftLists, index)

  def setShiftList(self, null):
  
    shiftList = self.shiftListPulldown.getObject()
    
    if self.experiment:
      project = self.experiment.root
      if shiftList is None:
        shiftList = newShiftList(project, unit='ppm')
      
      if shiftList and (shiftList is not self.experiment.shiftList):
        setExperimentShiftList(self.experiment, shiftList)
  """
  def showExperimentalDetails(self):

    self.guiParent.editExpDetails(experiment=self.experiment)
  """
  def showMeasurementLists(self):

    if self.experiment:
      self.guiParent.editMeasurementLists(experiment=self.experiment)
    else:
      self.guiParent.editMeasurementLists(experiment='<Any>')
      
  def cloneExperiment(self):

    if self.experiment:
      popup = NewExperimentPopup(self, self.experiment)
      popup.destroy()
    else:
      showWarning('Warning','Need to select experiment to clone', parent=self)

  def newExperiment(self):

    popup = NewExperimentPopup(self)
    popup.destroy()
    """
    numDim = int(askInteger('Experiment Dimensions','Number of dimensions?',2,parent=self) or 0)
    
    if numDim:
      n = len(self.nmrProject.experiments)+1
      name = 'New Exp.' + str(n)
      while self.nmrProject.findFirstExperiment(name=name):
        n   += 1
        name = 'New Exp.' + str(n)
    
      if numDim < 6:
        experiment = Nmr.Experiment(self.nmrProject,name=name,numDim=numDim)
      else:
        showWarning('Error','Experiment dimensionality greater\nthan 5 not implemented', parent=self)
"""

  def deleteExperiment(self, *event):

    experiments = self.scrolledMatrix.currentObjects[:]
    if len(experiments) == 1:
      if showOkCancel('Delete Experiment','Really delete experiment and\nany spectra, lists, peaks etc... ?', parent=self):
        self.experiment.delete()
        self.experiment = None
    elif len(experiments) > 1:
      if showOkCancel('Delete Experiment','Really delete %d experiments and\nany spectra, lists, peaks etc... ?' % (len(experiments)), parent=self):
        for experiment in experiments:
          experiment.delete()
        self.experiment = None

  def selectCell(self, object, row, col):

    self.experiment = object
    if self.experiment:
      for n in (1, 3):
        self.experimentButtons.buttons[n].enable()
 
  def updateAfter(self, *opt):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
 
  def update(self):

    if self.experiment:
      for n in (1, 3):
        self.experimentButtons.buttons[n].enable()
    else:
      for n in (1, 3):
        self.experimentButtons.buttons[n].disable()
 
    allMolSys = len(self.project.molSystems)
    objectList  = []
    textMatrix  = []
    for experiment in self.nmrProject.sortedExperiments():
      objectList.append(experiment)
      shiftList = experiment.shiftList
      if shiftList:
        shiftListText = '%s [%d]' % (shiftList.name or '<No name>', shiftList.serial)
      else:
        shiftListText = ''

      experimentType = None
      if experiment.refExperiment:
        experimentType = experiment.refExperiment.name
        
      dataSourcesText = ''
      numSpectra = 0
      
      dataSources = experiment.sortedDataSources()
      
      for dataSource in dataSources:
        dataSourcesText += dataSource.name
        if isSpectrum(dataSource):
          numSpectra += 1
        if dataSource is not dataSources[-1]:
          dataSourcesText += ' '

      rawData = 'No'
      if experiment.rawData:
        rawdata = 'Yes'
      
      shiftRefs = ['%d:%s:%s' % (sr.serial,sr.isotopeCode,sr.molName) for sr in experiment.shiftReferences]
        
      msCodesSorted = [(ms.code, ms) for ms in experiment.molSystems]
      msCodesSorted.sort()
      numMolSys = len(msCodesSorted)
      msCodes = []
      
      if (numMolSys == allMolSys) and (numMolSys > 5):
        msCodes = ['** ALL **',]
      
      elif numMolSys > 5:
        msCodes = ['** %d **' % numMolSys,]
      
      else:
        for i, (code, ms) in enumerate(msCodesSorted):
          if i > 0:
            if i%5 == 0:
              msCodes.append('\n%s' % code)
            else:
              msCodes.append(' %s' % code)
          else:
            msCodes.append(code)
        
      expDims = experiment.findAllExpDims(isAcquisition=True)
      names = []
      for expDim in expDims:
        expDimRef = getPrimaryExpDimRef(expDim)
        if expDimRef and expDimRef.isotopeCodes:
          name = '%d (%s)' % (expDim.dim, ','.join(expDimRef.isotopeCodes))
        else:
          name = '%d' % expDim.dim
        names.append(name)
      acqDim = ', '.join(names)

      datum   = [experiment.serial,
                 experiment.name,
                 experiment.numDim,
                 experimentType,
                 acqDim,
                 shiftListText,
                 dataSourcesText,
                 ''.join(msCodes),
                 ','.join(shiftRefs),
                 numSpectra,
                 rawData,
                 experiment.details or '']
       
      textMatrix.append(datum)

    if not objectList:
      textMatrix.append([])
 
    self.scrolledMatrix.update(objectList=objectList, textMatrix=textMatrix)
 
    self.waiting = False

  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
        
    BasePopup.destroy(self)
    
  # Experiment Type functionss


  def propagateExpType(self):
  
    experiments = self.expTypeMatrix.currentObjects
    if self.typeExpt and (len(experiments) > 1):
      refExpDict = {}
      refExperiment = None
      
      for experiment in experiments:
        refExperiment0 = experiment.refExperiment
        
        if refExperiment0:
          refExpDict[refExperiment0] = True
          refExperiment = refExperiment0
    
      if len(refExpDict.keys()) > 1:
        showWarning('Propagate type failed',
                    'Ambiguous experiment types in current selection', parent=self)
        return            
    
      if not refExperiment:
        showWarning('Warning',
                    'No experiment type to propagate', parent=self)
        return            
    
      name = refExperiment.name
      text = 'Really propagate type %s to selected experiments?' % name
      
      if showWarning('Confirm',text, parent=self):
        for experiment in experiments:
          if experiment.refExperiment:
            continue
        
          refExperiments = getFilteredRefExperiments(experiment)
          if refExperiment not in refExperiments:
             data = (name, experiment.name)
             showWarning('Propagate type error',
                         'Experiment type %s not possible for experiment %s' % data, parent=self)
             continue
          
          setRefExperiment(experiment, refExperiment)
          initExpTransfers(experiment)

    if self.isModal:
      self.updateExpTypesAfter() # Modal version doesn't let the notifier calls get through...
  
  def selectExpTransfer(self, obj, row, col):
  
    self.transferExpDimRefs = obj
  
  def selectExpDimRef(self, obj, row, col):
  
    self.expDimRef = object
  
  def getCategories(self):

    names = ['<None>', 'through-bond', 'through-space',
             'quantification', 'use external', 'other']
    
    return names
  
  def getExtSources(self):

    names = ['<None>','ccpn', 'bruker']
    
    return names
   
  def getSynonyms(self):
   
    names = []
    refExperimentSets = []
    
    if self.typeExpt:
    
      refExperiments = getRefExperiments(self.typeExpt)
      namesDict = {}
      for refExperiment in refExperiments:
        nmrExpPrototype = refExperiment.nmrExpPrototype
        if nmrExpPrototype.synonym:
          name = nmrExpPrototype.synonym
        else:
          name = nmrExpPrototype.name
        
        if name not in namesDict:
          namesDict[name] = set()
          
        namesDict[name].add(refExperiment)
    
      names = namesDict.keys()
      names.sort()
      refExperimentSets = [namesDict[name] for name in names]
    
    return names, refExperimentSets
   
  def getFullTypes(self, experiment=None):
    
    if experiment is None:
      experiment = self.typeExpt
    
    names = []
    refExperiments = []    
    
    if experiment:
      if experiment.refExperiment:
        if not (hasattr(experiment, 'pulProgName') and 
                hasattr(experiment, 'pulProgType') and
                hasattr(experiment, 'category') and 
                experiment.category == 'use external'):
          # no external source, try to use nmrExpPrototype synonym
          synonym = experiment.refExperiment.nmrExpPrototype.synonym
          if synonym:
            names = [(rx.name, rx) for rx in getRefExperiments(experiment)
                     if rx.nmrExpPrototype.synonym == synonym]
                 
      if not names:
        names = [(rx.name,rx) for rx in getRefExperiments(experiment)]
    
      names.sort()
      
      refExperiments = [x[1] for x in names]
      names = [x[0] for x in names]
      
    #names = ['<None>',] + names
    return names, refExperiments
  
  def getCategory(self, experiment):

    names = self.getCategories()
    
    if names:
      index = 0
      if experiment:
        if experiment.refExperiment:
          name = experiment.refExperiment.nmrExpPrototype.category or '<None>'
          if name in names:
            index = names.index(name)
        
        elif hasattr(experiment, 'category'):
          name = experiment.category or '<None>'
          if name in names:
            index = names.index(name)
        
        elif (hasattr(experiment, 'pulProgName') and 
              hasattr(experiment, 'pulProgType')):
          if experiment.pulProgName and experiment.pulProgType:
            index = names.index('use external')
        
      self.categoryPulldown.setup(names,names,index)
    else:
      self.categoryPulldown.setup([],[], 0)
  
  def getExtSource(self, experiment):

    names = self.getExtSources()
    if names:
      index = 0
      if experiment and hasattr(experiment, 'pulProgType'):
        name = experiment.pulProgType or '<None>'
        if name in names:
          index = names.index(name)
        
      self.extSourcePulldown.setup(names,names,index)
    else:
      self.extSourcePulldown.setup([],[], 0)
   
  def getSynonym(self, experiment):
   
    synonyms, refExperimentSets = self.getSynonyms()
    if synonyms:
      cats = [None,]
      names = ['<None>',] + synonyms
      
      for refExperiments in refExperimentSets:
        cat = set()
        
        nn = 0
        for refExperiment in refExperiments:
          xx = getRefExperimentCategories(refExperiment)
          if 'Projected' in xx:
            nn += 1
          cat.update(xx)
        
        if 'Projected' in cat and nn == len(refExperiments):
          # This means that if a synonym matches only because of
          # projection experiments, then the synonym appears only 
          # in the 'Projected' category
          cat = set(['Projected'])
        
        cats.append(cat or None) 
           
      index = 0
      if experiment and experiment.refExperiment:
        name = experiment.refExperiment.nmrExpPrototype.synonym
        if name in names:
          index = names.index(name)
    
      self.synonymPulldown.setup(names,names, index, None, cats)
    else:
      self.synonymPulldown.setup([],[], 0)

  def getExtName(self, experiment):

    if experiment and hasattr(experiment, 'pulProgName'):
      self.extNameEntry.set(experiment.pulProgName)

  def getFullType(self, experiment):

    names, refExperiments = self.getFullTypes()
    
    if names:
      
      if len(refExperiments) > 5:
        cats = [None,]
        
        for refExperiment in refExperiments:
          cat = getRefExperimentCategories(refExperiment)
          
          if 'Projected' in cat:
            cat = set(['Projected'])
          
          cats.append(cat or None)  
        
      else:
        cats = None
            
      names.insert(0,'<None>')
      index = 0
      if experiment and experiment.refExperiment:
        
        name = experiment.refExperiment.name
        if name in names:
          index = names.index(name)
    
    
      self.fullTypePulldown.setup(names, names, index, None, cats)
    else:
      self.fullTypePulldown.setup([],[], 0)

  def getExpDimRefName(self, expDimRef):
  
    expDim = expDimRef.expDim
    name =  '%d' % expDim.dim
    if len(expDim.expDimRefs) > 1:
      name += '.%d' % (expDimRef.serial)
    
    return name

  def getRefExpDimRefName(self, refExpDimRef):
  
    name =  '%d' % refExpDimRef.refExpDim.dim
    if len(refExpDimRef.refExpDim.refExpDimRefs) > 1:
      name += '.%d' % (refExpDimRef.serial)
    
    return name
      
  def getRefExpDimRefs(self, expDimRef):
  
    isotopeCodes = list(expDimRef.isotopeCodes)
    isotopeCodes.sort()
    
    refExperiment = expDimRef.expDim.experiment.refExperiment
    refExpDimRefs = []
    if refExperiment:
      for refExpDim in refExperiment.sortedRefExpDims():
        for refExpDimRef in refExpDim.sortedRefExpDimRefs():
          refIsotopeCodes = [atomSite.isotopeCode for atomSite in refExpDimRef.expMeasurement.atomSites]
          refIsotopeCodes.sort()
          
          if refIsotopeCodes == isotopeCodes:
            refExpDimRefs.append(refExpDimRef)
        
    return refExpDimRefs
  
  def getRefExpDim(self, expDimRef):

    refExpDimRefs = self.getRefExpDimRefs(expDimRef)
    if refExpDimRefs:
      names = [self.getRefExpDimRefName(x) for x in refExpDimRefs]
      index = 0
      
      if expDimRef.refExpDimRef:
        name = self.getRefExpDimRefName(expDimRef.refExpDimRef)
        if name in names:
          index = names.index(name)
      
      self.refExpDimPulldown.setup(names, refExpDimRefs, index)
    else:
      self.refExpDimPulldown.setup([],[],0)

  def getTransferExpDimA(self, transferExpDimRefs):
    
    expDimRefA, expDimRefB = transferExpDimRefs
    
    expDimRefs = []
    for expDim in self.typeExpt.expDims:
      for expDimRef in expDim.expDimRefs:
        if expDimRef is expDimRefB:
          continue
      
        if expDimRef.isotopeCodes == expDimRefA.isotopeCodes:
          expDimRefs.append(expDimRef)
    
    names = [self.getExpDimRefName(x) for x in expDimRefs]
    index = expDimRefs.index(expDimRefA)
    
    self.transferPulldownA.setup(names, expDimRefs, index)

  def getTransferExpDimB(self, transferExpDimRefs):
    
    expDimRefA, expDimRefB = transferExpDimRefs
    
    expDimRefs = []
    for expDim in self.typeExpt.expDims:
      for expDimRef in expDim.expDimRefs:
        if expDimRef is expDimRefA:
          continue
      
        if expDimRef.isotopeCodes == expDimRefB.isotopeCodes:
          expDimRefs.append(expDimRef)
    
    names = [self.getExpDimRefName(x) for x in expDimRefs]
    index = expDimRefs.index(expDimRefB)
    
    self.transferPulldownB.setup(names, expDimRefs, index)
  
  def swapRefExpDimRefs(self, expDimRef1, expDimRef2):
  
    # Old ref mappings
    expDim1 = expDimRef1.expDim
    expDim2 = expDimRef2.expDim
    refExpDim1 = expDim1.refExpDim
    refExpDim2 = expDim2.refExpDim
    refExpDimRef1 = expDimRef1.refExpDimRef
    refExpDimRef2 = expDimRef2.refExpDimRef
 
    # Swap expDimRef refMappings
 
    expDim1.setRefExpDim(None)
    expDimRef1.setRefExpDimRef(None)
    expDim2.setRefExpDim(None)
    expDimRef2.setRefExpDimRef(None)

    expDim1.setRefExpDim(refExpDim2)
    expDimRef1.setRefExpDimRef(refExpDimRef2)
    expDim2.setRefExpDim(refExpDim1)
    expDimRef2.setRefExpDimRef(refExpDimRef1)
    
  def setTransferExpDimA(self, null):
      
    if self.transferExpDimRefs:
      expDimRefA, expDimRefB = self.transferExpDimRefs
      expDimRefC = self.transferPulldownA.getObject()
 
      if expDimRefC is not expDimRefA:
        self.swapRefExpDimRefs(expDimRefA, expDimRefC)
        initExpTransfers(self.typeExpt)   
    
  def setTransferExpDimB(self, null):
      
    if self.transferExpDimRefs:
      expDimRefA, expDimRefB = self.transferExpDimRefs
      expDimRefC = self.transferPulldownB.getObject()
 
      if expDimRefC is not expDimRefB:
        self.swapRefExpDimRefs(expDimRefB, expDimRefC)
        initExpTransfers(self.typeExpt)   
      
  def setCategory(self, null):
    
    name = self.categoryPulldown.getText()
    
    if self.typeExpt:
      if name == '<None>':
        self.typeExpt.category = None
        if self.typeExpt.refExperiment:
          setRefExperiment(self.typeExpt,None)
         
      else:
        self.typeExpt.category = name
        if (self.typeExpt.refExperiment and 
            self.typeExpt.refExperiment.nmrExpPrototype.category != name):
          setRefExperiment(self.typeExpt,None)
        
      self.updateExpTypesAfter()
   
  def setExtSource(self, null):
    
    name = self.extSourcePulldown.getText()
   
    if self.typeExpt:
      if name == '<None>':
        self.typeExpt.pulProgType = None
         
      else:
        self.typeExpt.pulProgType = name
        
      self.updateExpTypesAfter()
   
  def setSynonym(self, null):

    name = self.synonymPulldown.getText()
   
    if self.typeExpt:
      if not self.typeExpt.refExperiment or (self.typeExpt.refExperiment.nmrExpPrototype.synonym != name ):
   
        if name == '<None>':
          if self.typeExpt.refExperiment:
            setRefExperiment(self.typeExpt,None)
            initExpTransfers(self.typeExpt)
        else:
          refExperiments = getRefExperiments(self.typeExpt)
          for refExperiment in refExperiments:
            if refExperiment.nmrExpPrototype.synonym == name:
              setRefExperiment(self.typeExpt, refExperiment)
              initExpTransfers(self.typeExpt)
              break
          else: 
            for refExperiment in refExperiments:
              if refExperiment.nmrExpPrototype.name == name:
                setRefExperiment(self.typeExpt, refExperiment)
                initExpTransfers(self.typeExpt)
                break
 
     
      if self.isModal:
        self.updateExpTypesAfter()
      

  def setExtName(self, null):

    if self.typeExpt:
      text = self.extNameEntry.get()
      if text and text != ' ':
        self.typeExpt.pulProgName = text
      else:
        self.typeExpt.pulProgName = None
      if self.isModal:
        self.updateExpTypesAfter()
   
  def setFullType(self, null):

    name = self.fullTypePulldown.getText()
   
    if self.typeExpt:
      if not self.typeExpt.refExperiment or (self.typeExpt.refExperiment.name != name ):
        
        if name == '<None>':
          if self.typeExpt.refExperiment:
            setRefExperiment(self.typeExpt,None)
            initExpTransfers(self.typeExpt)
        else:
          refExperiments = getRefExperiments(self.typeExpt)
          for refExperiment in refExperiments:
            if refExperiment.name == name:
              setRefExperiment(self.typeExpt, refExperiment)
              initExpTransfers(self.typeExpt)
              break

      if self.isModal:
        self.updateExpTypesAfter()

  def setRefExpDim(self, null):
    expDimRef = self.refExpDimMatrix.currentObject
    if expDimRef:
      refExpDimRef = self.refExpDimPulldown.getObject()
 
      if not expDimRef.refExpDimRef or \
        (expDimRef.refExpDimRef is not refExpDimRef):
 
        expDimRef2 = None
        refExpDimRef2 = expDimRef.refExpDimRef
        for expDim in self.typeExpt.expDims:
          for expDimRef3 in expDim.expDimRefs:
            if expDimRef3.refExpDimRef and \
              (expDimRef3.refExpDimRef is refExpDimRef):
              
              expDimRef2 = expDimRef3
              break
          else:
            continue
          break
 
        expDimRef.setRefExpDimRef(None)
        expDimRef.expDim.setRefExpDim(None) 
        
        if refExpDimRef:
 
          if refExpDimRef2 and expDimRef2:
            expDimRef2.expDim.setRefExpDim(refExpDimRef2.refExpDim)
            expDimRef2.setRefExpDimRef(refExpDimRef2)
 
          expDimRef.expDim.setRefExpDim(refExpDimRef.refExpDim)
          expDimRef.setRefExpDimRef(refExpDimRef)

        initExpTransfers(expDimRef.expDim.experiment)

  def selectTypeExpt(self, object, row, col):
  
    if object is not self.typeExpt:
      self.typeExpt = object
      self.updateExpTransfers()
      self.updateRefDims()

      
  def updateExpTypes(self, experiments=None):

    textMatrix = []
    objectList = []
    
    if experiments:
      self.typeExps = experiments
    else:
      self.typeExps = self.nmrProject.sortedExperiments()
    
    
    for experiment in self.typeExps:
      
      # get extSource and extName,
      extSource = None
      if hasattr(experiment,'pulProgName'):
        extName = experiment.pulProgName
        
        if hasattr(experiment,'pulProgType'):
          # expSource is irrelevant without an extName
          extSource = experiment.pulProgType
          
      else:
        hasExtName = False
        extName = None
      
      # initialise: use external source, if available, to set refExperiment
      refExperiment = experiment.refExperiment
      if extSource is not None and not hasattr(experiment, 'category'):
        # this is first time we get here, and we have external name and source 
        # use external source to set fullType
        experiment.category = 'use external'
        refExperiments = getRefExperiments(experiment)
        if refExperiments:
          if len(refExperiments) < 8:
            # set to first refExp found. NB len() check should not be necessary.
            refExperiment = refExperiments[0]
            setRefExperiment(experiment, refExperiment)
            initExpTransfers(experiment)
        else:
          # no match found. use category None instead.
          del experiment.category
 
      # get list of allowed type strings (NB must be after initialisation)
      altTypes, altRefExps = self.getFullTypes(experiment)
      altTypeStr = None
      
      if refExperiment:
        # refExperiment was set - set parameters
        category = refExperiment.nmrExpPrototype.category
        synonym  = refExperiment.nmrExpPrototype.synonym
        fullType = refExperiment.name
        if fullType in altTypes:
          # type fits possibilities - remove from altTypes
          altTypes.remove(fullType)
        else:
          # fullType not in possibilities - add warning
          altTypeStr = 'NB! Inconsistent full type'
      else:
        # no refExperiment - set up without it
        if hasattr(experiment, 'category'):
          category = experiment.category
        else:
          category = None
        synonym  = None
        fullType = None
      
      if altTypeStr is None:
        # set string for altTypes, if not pre-set
        if len(altTypes) > 4:
          altTypeStr = ' ... '
        elif altTypes:
          altTypeStr = '; '.join(altTypes)
        else:
          altTypeStr = '<None>'
      
      datum = [experiment.serial,
               experiment.name,
               extSource,
               extName,
               category,
               synonym,
               fullType,
               altTypeStr]
       
      objectList.append(experiment)
      textMatrix.append(datum)

    self.expTypeMatrix.update(textMatrix=textMatrix, objectList=objectList)
    
    if self.typeExpt not in objectList:
      self.typeExpt = None
    
    self.updateExpTransfers()
    self.updateRefDims()
    self.waitingT = False
    
  def updateRefDims(self):

    textMatrix = []
    objectList = []
    if self.typeExpt:
      for expDim in self.typeExpt.sortedExpDims():
        for expDimRef in expDim.sortedExpDimRefs():
        
          measurementText = None
          redr = expDimRef.refExpDimRef
          if redr:
            if redr.expSteps:
              steps   = ','.join([str(step.stepNumber) for step in redr.expSteps])
            else:
              steps   = ','.join([str(step.stepNumber) for step in redr.expMeasurement.expSteps])
            measure = redr.expMeasurement.measurementType
            atoms   = '(%s)' % ( ','.join([site.name for site in redr.expMeasurement.atomSites]) )
            time    = ' %s timing' % redr.constantTime
            couple  = ''
            
            if redr.coupledIsotopeCodes:
              couple = ' Coupling:%s' % ( ','.join(redr.coupledIsotopeCodes))
          
            measurementText = '%s %s%s%s%s' % (steps,measure,atoms,couple,time)
        
          expDimText = '%d' % expDim.dim
          if len(expDim.expDimRefs) > 1:
            expDimText += '.%d' % expDimRef.serial
          ss = expDimRef.displayName
          if ss:
            expDimText += ' (%s)' % ss
        
          refExpDimText = ''
          if redr:
            refExpDimText = self.getRefExpDimRefName(redr)
              
          datum = [ expDimText,
                    refExpDimText,
                    ','.join(expDimRef.isotopeCodes),
                    measurementText]
            
          objectList.append(expDimRef)
          textMatrix.append(datum)
    
    self.refExpDimMatrix.update(textMatrix=textMatrix, objectList=objectList)

  def updateExpTransfers(self):

    textMatrix = []
    objectList = []
    if self.typeExpt:
      expTransfers = set()
      
      for expDim in self.typeExpt.expDims:
        for expDimRef in expDim.expDimRefs:
          for expTransfer in expDimRef.expTransfers:
            expTransfers.add(expTransfer)
    
      for expTransfer in expTransfers:
        expDimRefA, expDimRefB = list(expTransfer.expDimRefs)
        isotopesA = ','.join(expDimRefA.isotopeCodes)
        isotopesB = ','.join(expDimRefB.isotopeCodes)
        
        if expDimRefA.serial > 1:
          dimA = '%d.%d (%s)' % (expDimRefA.expDim.dim,expDimRefA.serial,isotopesA)
        else:
          dimA = '%d (%s)' % (expDimRefA.expDim.dim,isotopesA)
        
        if expDimRefB.serial > 1:
          dimB = '%d.%d (%s)' % (expDimRefB.expDim.dim,expDimRefB.serial,isotopesB)
        else:
          dimB = '%d (%s)' % (expDimRefB.expDim.dim,isotopesB)
      
        (dimA, expDimRefA), (dimB, expDimRefB) = sorted([(dimA, expDimRefA), (dimB, expDimRefB)])
    
        datum = [dimA, expTransfer.transferType, dimB]
 
        objectList.append( (expDimRefA, expDimRefB) )
        textMatrix.append(datum)
    
    self.transferMatrix.update(textMatrix=textMatrix,
                               objectList=objectList)

  def updateExpTypesAfter(self, *obj):

    if self.waitingT:
      return
    else:
      self.waitingT = True
      self.after_idle(self.updateExpTypes)
  

  # Experiment Details Functions
  
  def propagateDetails(self):
  
    experiments = self.detailsMatrix.currentObjects
    if len(experiments) < 2:
      return
    
    expt  = self.detailsExp
    experiments.remove(expt)
    names = ','.join([e.name for e in experiments])
    msg   = 'Propagate details from experiment %s to %s?' % (expt.name,names)
    if showOkCancel('Confirm',msg, parent=self):
      spectrometer  = expt.spectrometer
      probe         = expt.probe
      numScans      = expt.numScans
      sampleState   = expt.sampleState
      sampleVolume  = expt.sampleVolume
      volumeUnit    = expt.volumeUnit
      nmrTubeType   = expt.nmrTubeType
      spinningAngle = expt.spinningAngle
      spinningRate  = expt.spinningRate
 
      for expt2 in experiments:
        expt2.spectrometer  = spectrometer
        expt2.probe         = probe
        expt2.numScans      = numScans
        expt2.sampleState   = sampleState
        expt2.sampleVolume  = sampleVolume
        expt2.volumeUnit    = volumeUnit
        expt2.nmrTubeType   = nmrTubeType
        expt2.spinningAngle = spinningAngle
        expt2.spinningRate  = spinningRate
  
  def selectExperiment(self, obj, row, col):
  
    self.detailsExp = obj
      
  def setNumScans(self, event):

    if self.detailsExp:
      value = self.scansEntry.get()
      
      if (value is not None) and (value < 1):
        value = None 
      
      self.detailsExp.numScans = value

  def getSpectrometers(self):

    store          = self.project.currentInstrumentStore
    spectrometers  = [None,]
    spectrometers += store.findAllInstruments(className='NmrSpectrometer')
      
    return spectrometers
  
  def getProbes(self):
  
    store   = self.project.currentInstrumentStore
    probes  = [None, ]
    probes += store.findAllInstruments(className='NmrProbe')
      
    return probes

  def setSpectrometer(self, null):

    spectrometer = self.spectrometerPulldown.getObject()
    
    if self.detailsExp:
      self.detailsExp.spectrometer = spectrometer

  def setProbe(self, null):

    probe = self.probePulldown.getObject()
    
    if self.detailsExp:
      self.detailsExp.probe = probe

  def setSampleState(self, null):

    sampleState = self.statePulldown.getText()
    
    if self.detailsExp:
      self.detailsExp.sampleState = sampleState

  def setSampleVolume(self, event):

    if self.detailsExp:
      value = self.volumeEntry.get()
      
      if (value is not None) and (value <= 0.0):
        value = None
        
      self.detailsExp.sampleVolume = value

  def setSampleVolUnit(self, index, name=None):

    volumeUnit = self.unitPulldown.getText()

    if self.detailsExp:
      self.detailsExp.volumeUnit = volumeUnit

  def setNmrTube(self, event):

    if self.detailsExp:
      self.detailsExp.nmrTubeType = self.tubeEntry.get() or None

  def setSpinRate(self, event):

    if self.detailsExp:
      value = self.spinRateEntry.get()
      
      if (value is not None) and (value < 0.0):
        value = None
      
      self.detailsExp.spinningRate = value

  def setSpinAngle(self, event):

    if self.detailsExp:
      self.detailsExp.spinningAngle = self.spinAngleEntry.get()

  def getNumScans(self, experiment):
  
    value = 0
    if experiment:
      value = experiment.numScans
      
    self.scansEntry.set(value)
    
  def getSampleState(self, experiment):

    if experiment and experiment.sampleState:
      self.statePulldown.set(experiment.sampleState)

  def getSpectrometer(self, experiment):

    index = 0
    names = ['<None>',]
    
    if experiment:
      spectrometers = self.getSpectrometers()
      names += ['%d:%s' % (s.serial,s.name) for s in spectrometers[1:]]
      index  = spectrometers.index(experiment.spectrometer)

    self.spectrometerPulldown.setup(names, spectrometers, index)

  def getProbe(self, experiment):

    index = 0
    names = ['<None>',]
    
    if experiment:
      probes = self.getProbes()
      names += ['%d:%s' % (p.serial,p.name) for p in probes[1:]]
      index  = probes.index(experiment.probe)

    self.probePulldown.setup(names, probes, index)

  def getSampleVolume(self, experiment):

    value = 0.0
    if experiment:
      value = experiment.sampleVolume
  
    self.volumeEntry.set(value)


  def getSampleVolUnit(self, experiment):

    index = 1
    if experiment.volumeUnit in VOLUME_UNITS:
      index = VOLUME_UNITS.index(experiment.volumeUnit)

    self.unitPulldown.setup(VOLUME_UNITS, VOLUME_UNITS, index)


  def getNmrTube(self, experiment):
  
    text = ''
    if experiment:
      text = experiment.nmrTubeType or ''
  
    self.tubeEntry.set(text)
  

  def getSpinRate(self, experiment):
  
    value = 0.0
    if experiment:
      value = experiment.spinningRate

    self.spinRateEntry.set(value)


  def getSpinAngle(self, experiment):
    
    value = 0.0
    if experiment:
      value = experiment.spinningAngle 
      
      if value is None:
        value = 54.74
 
    self.spinAngleEntry.set(value)
    

  def setConditionSet(self, null):

    sampleConditionSet = self.conditionsPulldown.getObject()

    if self.detailsExp:
      self.detailsExp.setSampleConditionSet(sampleConditionSet)
      # Could be none


  def updateConditionSets(self, sampleConditionSet=None):
  
    index = 0
    sampleConditionSets = self.nmrProject.sortedSampleConditionSets()
    names = ['%d' % scs.serial for scs in sampleConditionSets]
    
    sampleConditionSets.append(None)
    names.append('<None>')
    
    if self.detailsExp and self.detailsExp.sampleConditionSet:
      index = sampleConditionSets.index(self.detailsExp.sampleConditionSet)
    
    self.conditionsPulldown.setup(names, sampleConditionSets, index)

  def editConditions(self):
  
    self.parent.editSampleConditionSets(sampleConditionSet=self.detailsExp.sampleConditionSet)

  def setExperiment(self, index, name):
  
    experiment = self.nmrProject.findFirstExperiment(name=name)
    if experiment is not self.detailsExp:
      self.detailsExp = experiment
      self.updateExpDetailsAfter()

  def updateExpDetailsAfter(self, *opt):

    if self.waitingD:
      return
    else:
      self.waitingD = True
      self.after_idle(self.updateExpDetails)
 
  def updateExpDetails(self):

    
    # Exp details
    objectList = []
    textMatrix = []
    for expt in self.nmrProject.sortedExperiments():
    
      probe = expt.probe
      if probe:
        probe = '%d:%s' % (probe.serial,probe.name)
    
      spectrometer = expt.spectrometer
      if spectrometer:
        spectrometer = '%d:%s' % (spectrometer.serial,spectrometer.name)
    
      datum = [expt.serial,
               expt.name,
               spectrometer,
               probe,
               expt.numScans,
               expt.sampleState,
               expt.sampleVolume,
               expt.volumeUnit,                                         
               expt.nmrTubeType,                                        
               expt.spinningAngle,                                      
               expt.spinningRate]                                      
                                                                                   
      objectList.append(expt)
      textMatrix.append(datum)                                                     
                                                                                   
    self.detailsMatrix.update(objectList=objectList, textMatrix=textMatrix)

    # Conditions
    objectList  = []
    textMatrix  = []
    if self.detailsExp and self.detailsExp.sampleConditionSet:
      for sampleCondition in self.detailsExp.sampleConditionSet.sampleConditions:
        datum = [sampleCondition.condition,
                 sampleCondition.value,
                 sampleCondition.error,
                 sampleCondition.unit]
 
        textMatrix.append(datum)
        objectList.append(sampleCondition)
    
    self.conditionsMatrix.update(objectList=objectList, textMatrix=textMatrix)
    
    self.waitingD = False
  
    # Shift Reference Functions

  def doShiftRefEditMarkExtraRules(self, obj, row, col):

    if (col > 9) and (obj.className != 'ExternalShiftReference'):
      return False
      
    return True  
  
  def editExperiments(self, obj):
  
    self.parent.editExperiment()
    
  def toggleShiftRefType(self, shiftReference):
  
    if shiftReference:
      if shiftReference.referenceType == 'direct':
        shiftReference.setReferenceType('indirect')
      else:
        shiftReference.setReferenceType('direct')
    
  def getShiftRefValue(self, shiftReference):
  
    value = 0.0
    if shiftReference:
      value = shiftReference.value
     
    self.valueEntry.set(value)
    
  def getShiftRefRatio(self, shiftReference):
  
    value = 0.0
    if shiftReference:
      value = shiftReference.indirectShiftRatio
     
    self.ratioEntry.set(value)
    
  def getShiftRefGeometry(self, shiftReference):
  
    text = ''
    if shiftReference and (shiftReference.className == 'ExternalShiftReference'):
      text = shiftReference.sampleGeometry
  
    self.geometryEntry.set(text)
  
  def getShiftRefLocation(self, shiftReference):
  
    text = ''
    if shiftReference and (shiftReference.className == 'ExternalShiftReference'):
      text = shiftReference.location
  
    self.locationEntry.set(text)
  
  def getShiftRefAxis(self, shiftReference):
  
    text = ''
    if shiftReference and (shiftReference.className == 'ExternalShiftReference'):
      text = shiftReference.axis
  
    self.axisEntry.set(text)
  
  def getShiftRefAtomGroup(self, shiftReference):
    
    text = ''
    if shiftReference:
      text = shiftReference.atomGroup
  
    self.atomGroupEntry.set(text)
  
  def getShiftRefIsotope(self, shiftReference):
  
    self.isotopePulldown.set(shiftReference.isotopeCode)
    
    
  def getShiftRefMolName(self, shiftReference):
  
    self.molNamePulldown.set(shiftReference.molName)


  def getShiftRefUnit(self, shiftReference):
  
    self.unitPulldown2.set(shiftReference.unit)
    
    
  def setShiftRefValue(self, event):
    
    if self.shiftReference:
      value = self.valueEntry.get() or 0.0
    
      self.shiftReference.value = value
    
  def setShiftRefRatio(self, event):
  
    if self.shiftReference:
      value = self.ratioEntry.get() or None
    
      self.shiftReference.indirectShiftRatio = value
    
  def setShiftRefGeometry(self, event):
  
    if self.shiftReference:
      text = self.geometryEntry.get() or None
  
      self.shiftReference.sampleGeometry = text
    
  def setShiftRefLocation(self, event):
  
    if self.shiftReference:
      text = self.locationEntry.get() or None
  
      self.shiftReference.location = text
    
  def setShiftRefAxis(self, event):
  
    if self.shiftReference:
      text = self.axisEntry.get() or None
  
      self.shiftReference.axis = text
    
  def setShiftRefAtomGroup(self, event):
    
    if self.shiftReference:
      text = self.atomGroupEntry.get() or None
  
      self.shiftReference.atomGroup = text
    
  def setShiftRefIsotope(self, null):
  
    isotopeCode = self.isotopePulldown.getText()
     
    self.shiftReference.isotopeCode = isotopeCode
    
  def setShiftRefMolName(self, null):
  
    molName = self.molNamePulldown.getText() 
     
    self.shiftReference.molName = molName
    
  def setShiftRefUnit(self, null):
  
    unit = self.unitPulldown2.getText()
     
    self.shiftReference.unit = unit
  
  def addInternalShiftRef(self):
  
    if self.nmrProject:
      newRef = self.nmrProject.newInternalShiftReference
      self.shiftReference = newRef(isotopeCode='1H', molName='TSP',
                                   atomGroup='H', value=0.000,
                                   referenceType='direct')

  def addExternalShiftRef(self):
  
    if self.nmrProject:
      newRef = self.nmrProject.newExternalShiftReference
      self.shiftReference = newRef(isotopeCode='1H', molName='TSP',
                                   atomGroup='H', value=0.000,
                                   referenceType='direct')

  def removeShiftRefs(self):
  
    haveExpts = False
    for shiftReference in self.shiftRefMatrix.currentObjects:
      if shiftReference.experiments:
        haveExpts = True
        break
    
    if haveExpts and not showOkCancel('Confirm','Really delete shift references with links to experiments?'):
      return
 
    for shiftReference in self.shiftRefMatrix.currentObjects:
      shiftReference.delete()

  def selectShiftRef(self, object, row, col):

    if object:
      self.shiftReference = object

  def updateShiftRefsAfter(self, *opt):

    if self.waitingS:
      return
    else:
      self.waitingS = True
      self.after_idle(self.updateShiftRefs)

  def updateShiftRefs(self):

    objectList  = []
    textMatrix  = []
    if self.nmrProject:
      for shiftReference in self.nmrProject.sortedShiftReferences():
        
        refClass = shiftReference.className[:8]
        
        if refClass == 'External':
          geometry = shiftReference.sampleGeometry
          location = shiftReference.location
          axis     = shiftReference.axis
        else:
          geometry = location = axis = None
        
        #' '.join([e.name for e in shiftReference.experiments]),
        datum = [shiftReference.serial,
                 refClass,
                 shiftReference.isotopeCode,
                 len(shiftReference.experiments),
                 shiftReference.molName,
                 shiftReference.atomGroup,
                 shiftReference.value,
                 shiftReference.unit,
                 shiftReference.referenceType,
                 shiftReference.indirectShiftRatio,
                 geometry,location,axis]
                 
        textMatrix.append(datum)
        objectList.append(shiftReference)

    self.shiftRefMatrix.update(objectList=objectList, textMatrix=textMatrix)

    self.waitingS = False
  
  # NMR Instrument functions
  

  def selectProbe(self, obj, row, col):
  
    self.probe = obj

  def selectSpectrometer(self, obj, row, col):
  
    self.spectrometer = obj
  
  def newProbe(self):
  
    instrumentStore = None
  
    if instrumentStore is None:
      instrumentStore = self.project.currentInstrumentStore  
        
    if instrumentStore is None:
      instrumentStore = self.project.findFirstInstrumentStore()
      
    if instrumentStore is None:
      instrumentStore = self.project.newInstrumentStore(name='Default')
  
    name = ''
    while not name:
      name = askString('Text Entry','Enter NMR Probe Name')
  
    probe = instrumentStore.newNmrProbe(name=name)
  
  def newSpectrometer(self):
  
    instrumentStore = None
  
    if instrumentStore is None:
      instrumentStore = self.project.currentInstrumentStore  
        
    if instrumentStore is None:
      instrumentStore = self.project.findFirstInstrumentStore()
      
    if instrumentStore is None:
      instrumentStore = self.project.newInstrumentStore(name='Default')
  
    name = ''
    while not name:
      name = askString('Text Entry','Enter NMR Spectrometer Name')
  
    spectrometer = instrumentStore.newNmrSpectrometer(name=name)
 
  def deleteProbe(self):
    
    if self.probe:
      msg = 'Really remove specification of NMR probe %s?' % self.probe.name
    
      if showOkCancel('query',msg):
        self.probe.delete()
        self.probe = None
  
  def deleteSpectrometer(self):
    
    if self.spectrometer:
 
      msg = 'Really remove specification of NMR Spectrometer %s?' %  self.spectrometer.name
    
      if showOkCancel('query',msg):
        self.spectrometer.delete()
        self.spectrometer = None
    
    
  def setProbeName(self, event):
      
    if self.probe:
      default = 'Probe%d' % self.probe.serial
      self.probe.name = self.probeNameEntry.get() or default


  def setProbeType(self, index, name=None):
  
    if name is None:
      index = self.probeTypePulldown.getSelectedIndex()
    
    if self.probe:
      self.probe.setProbeType( NMR_PROBE_TYPES[index] )
      self.updateInstrumentsAfter()


  def setProbeModel(self, event):
      
    if self.probe:
      self.probe.setModel( self.probeModelEntry.get() or None )
   
        
  def setProbeSerial(self, event):

    if self.probe:
      self.probe.serialNumber = self.probeSerialEntry.get() or None
 
  
  def setProbeDiameter(self, event):

    if self.probe:
      self.probe.setDiameter( self.probeDiameterEntry.get() or None )
      self.updateInstrumentsAfter()

  
  def setProbeDetails(self, event):

    if self.probe:
      self.probe.setDetails( self.probeDetailsEntry.get() or None )
    
    
  def getProbeName(self,probe):
  
    text = ''
    if probe:
      text = probe.name

    self.probeNameEntry.set(text)
  
  def getProbeType(self,probe):
  
    index = -1
    if probe and probe.probeType:
      index = NMR_PROBE_TYPES.index(probe.probeType)

    self.probeTypePulldown.setIndex(index)
    
  def getProbeModel(self,probe):
  
    text = ''
    if probe:
      text = probe.model

    self.probeModelEntry.set(text)
  
  def getProbeSerial(self,probe):

    text = ''
    if probe:
      text = probe.serialNumber

    self.probeSerialEntry.set(text)
  
  def getProbeDiameter(self,probe):
  
    value = 0.0
    if probe:
      value = probe.diameter

    self.probeDiameterEntry.set(value)
  
  def getProbeDetails(self,probe):
  
    text = ''
    if probe:
      text = probe.details

    self.probeDetailsEntry.set(text)


  def setSpectrometerName(self, event):
      
    if self.spectrometer:
      default = 'self.spectrometer%d' % self.spectrometer.serial
      self.spectrometer.name = self.spectrometerNameEntry.get() or default


  def setSpectrometerFreq(self, event):
   
    if self.spectrometer:
      value = self.spectrometerFreqEntry.get() or None
      self.spectrometer.setProtonFreq( value )
      
      if value is not None:
        value = '%d' % round(value)
      
      self.spectrometer.setNominalFreq( value )
      self.updateInstrumentsAfter()
     
        
  def setSpectrometerModel(self, event):

    if self.spectrometer:
      self.spectrometer.setModel( self.spectrometerModelEntry.get() or None )
       
        
  def setSpectrometerSerial(self, event):

    if self.spectrometer:
      self.spectrometer.serialNumber = self.spectrometerSerialEntry.get() or None
    
    
  def setSpectrometerDetails(self, event):

    if self.spectrometer:
      self.spectrometer.setDetails( self.spectrometerDetailsEntry.get() or None )
    
    
  def getSpectrometerName(self, spectrometer):
  
    text = ''
    if spectrometer:
      text = spectrometer.name

    self.spectrometerNameEntry.set(text)
  
  def getSpectrometerFreq(self, spectrometer):
  
    value = 0.0
    if spectrometer:
       value = spectrometer.protonFreq

    self.spectrometerFreqEntry.set(value)  
      
  def getSpectrometerModel(self, spectrometer):
  
    text = ''
    if spectrometer:
      text = spectrometer.model

    self.spectrometerModelEntry.set(text)
  
  def getSpectrometerSerial(self, spectrometer):

    text = ''
    if spectrometer:
      text = spectrometer.serialNumber

    self.spectrometerSerialEntry.set(text)
  
  def getSpectrometerDetails(self, spectrometer):
  
    text = ''
    if spectrometer:
      text = spectrometer.details

    self.spectrometerDetailsEntry.set(text)
  
  def updateInstrumentsAfter(self, *opt):

    if self.waitingI:
      return
    else:
      self.waitingI = True
      self.after_idle(self.updateInstruments)
 
  def updateInstruments(self):
    
    store = self.project.currentInstrumentStore
    getInstruments = store.findAllInstruments
    
    # Probe
    objectList = []
    textMatrix = []
    
    probes = [(p.serial, p) for p in getInstruments(className='NmrProbe')]
    probes.sort()
    for serial, probe in probes:
      datum = [serial,
               probe.name,
               probe.probeType,
               probe.model,
               probe.serialNumber,
               probe.diameter,
               probe.details]
 
      objectList.append(probe)
      textMatrix.append(datum)

    self.probeMatrix.update(objectList=objectList, textMatrix=textMatrix)

    # Spectrometer
    objectList = []
    textMatrix = []
    spectrometers = [(s.serial, s) for s in getInstruments(className='NmrSpectrometer')]
    spectrometers.sort()
    for serial, spectrometer in spectrometers:
      datum = [serial,
               spectrometer.name,
               spectrometer.nominalFreq,
               spectrometer.protonFreq,
               spectrometer.model,
               spectrometer.serialNumber,
               spectrometer.details]
 
      objectList.append(spectrometer)
      textMatrix.append(datum)

    self.spectrometerMatrix.update(objectList=objectList, textMatrix=textMatrix)

    self.waitingI = False

class NewExperimentPopup(BasePopup):

  def __init__(self, parent, experiment=None):

    self.guiParent = parent
    self.experiment = experiment
  
    BasePopup.__init__(self, parent, modal=True, title="Experiment detail")

  def body(self, guiFrame):

    experiment = self.experiment

    self.geometry('400x100')

    guiFrame.expandGrid(0,1)

    row = 0

    if not experiment:
      Label(guiFrame, text='Number of dimensions: ', grid=(row, 0))
      objects = range(1,7)
      texts = ['%s' % object for object in objects]
      self.numDimList = PulldownList(guiFrame, texts=texts, objects=objects, index=1,
                                     tipText='The number of dimensions of the experiment',
                                     callback=self.updateIsotopeLists, grid=(row,1), sticky='w')
      row += 1

      Label(guiFrame, text='Isotopes: ', grid=(row, 0))
      self.isotopeFrame = Frame(guiFrame, grid=(row, 1), sticky='w')
      self.isotopeLists = []
      self.updateIsotopeLists()
      row += 1

    Label(guiFrame, text='Experiment name: ', grid=(row, 0))
    name = ''
    self.exptNameEntry = Entry(guiFrame, text=name, tipText='The name of the experiment', width=30, grid=(row,1), sticky='ew')
    row += 1

    self.fileCheckButton = CheckButton(guiFrame, text='Include link for data file?',
                                       tipText='Whether the new spectrum should have a link to a data file',
                                       selected=True, grid=(row,0), gridSpan=(1,2), sticky='w')

    if experiment:
      text = 'Clone %s' % experiment.name
      ss = ''
    else:
      text = 'Create new'
      ss = ' and number of dimensions'
    tipTexts = ['%s experiment and spectrum using this name%s' % (text, ss)]
    texts    = [text]
    commands = [self.ok]
    self.buttons = UtilityButtonList(guiFrame, commands=commands, texts=texts,
                                     tipTexts=tipTexts, doClone=False, grid=(row, 1), sticky='e')

  def updateIsotopeLists(self, *extra):

    numDim = self.numDimList.getObject()
    n = 0
    for n, isotopeList in enumerate(self.isotopeLists):
      if n < numDim:
        isotopeList.grid(row=0, column=n)
      else:
        isotopeList.grid_forget()

    for m in range(n, numDim):
      texts = ('1H', '2H', '13C', '15N', '19F', '31P', '79Br')
      isotopeList = PulldownList(self.isotopeFrame, texts=texts,
                                 tipText='The isotope for dimension %s' % (m+1), grid=(0,m))
      self.isotopeLists.append(isotopeList)
 
  def apply(self):

    experimentOld = self.experiment
    if experimentOld:
      numDim = experimentOld.numDim
    else:
      numDim = self.numDimList.getObject()
      isotopeCodes = [isotopeList.getText() for isotopeList in self.isotopeLists[:numDim]]

    exptNameNew = self.exptNameEntry.get()
    spectrumNameNew = None

    if not exptNameNew:
      showWarning('Experiment name', 'Experiment name not specified', parent=self)
      return False

    names = [expt.name for expt in self.nmrProject.experiments]
    if exptNameNew in names:
      showWarning('Name used', 'Name %s already used' % exptNameNew, parent=self)
      return False

    includeDataFile = self.fileCheckButton.get()
    if experimentOld:
      experimentNew = cloneExperiment(experimentOld, exptNameNew, cloneDataFile=includeDataFile)
      dataSourceNew = experimentNew.findFirstDataSource()
    else:
      if includeDataFile:
        import os
        dataPath = os.path.join(os.getcwd(), 'dataFile.spc')
      else:
        dataPath = None
      experimentNew = defaultExperiment(self.nmrProject, exptNameNew, isotopeCodes, dataPath=dataPath)
      dataSourceNew = experimentNew.findFirstDataSource()
      from ccpnmr.analysis.popups.EditSpectrum import EditSpectrumPopup
      popup = EditSpectrumPopup(self.guiParent, transient=False, modal=True,
                                useReducedDim=False, spectrum=dataSourceNew)
      cancelled = popup.cancelled
      popup.destroy()
      if cancelled:
        dataSourceNew.delete()

    self.parent.parent.finishInitSpectrum(dataSourceNew)

    return True
