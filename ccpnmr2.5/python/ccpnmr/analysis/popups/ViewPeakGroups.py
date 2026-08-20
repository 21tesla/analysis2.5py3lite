
"""
======================COPYRIGHT/LICENSE START==========================

ViewPeakGroups.py: Part of the CcpNmr Analysis program

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
from memops.general import Implementation

from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.ButtonList import UtilityButtonList

from ccpnmr.analysis.popups.BasePopup import BasePopup

class ViewPeakGroupsPopup(BasePopup):
  
  """
  **Present Data About Different Peak Categories**

  This popup is used to allow the user to see information about different
  categories of  peak. At present this system is only used by the [Test Shift
  Match] function  of the `Make Distance Restraints`_ popup.

  Textual data about the different groups of peaks is presented in the upper
  panel. The lower buttons allow the user to open tables which lists the peaks
  in each of the categories. Specifically these open the `Selected Peaks`_
  table, and from here the user can perform many more functions, including peak
  assignment, deletion and locating peaks in spectrum windows.

  .. _`Make Distance Restraints`: CalcDistConstraintsPopup.html
  .. _`Selected Peaks`: SelectedPeaksPopup.html

  """ 

  def __init__(self, parent, title, message, groups, groupNames, tipTexts, *args, **kw):

    self.guiParent = parent
    self.message = message
    self.groups = groups
    self.groupNames = groupNames
    self.text = title
    self.tipTexts = tipTexts
      
    if len(self.groupNames) < len(self.groups):
      for i in range(len(self.groupNames),len(self.groups)):
        self.groupNames.append('Group %d' % i)
    
    BasePopup.__init__(self, parent=parent, title='View Peak Groups', **kw)
  
  def body(self, guiFrame):
       
    guiFrame.grid_columnconfigure(0, weight=1)
    
    subFrame0 = LabelFrame(guiFrame, text=self.text, grid=(0,0))
    label = Label(subFrame0, text=self.message, anchor='w', grid=(0,0))
    
    subFrame1 = LabelDivider(guiFrame, text='View Peaks', grid=(1,0))
      
    texts = []
    commands = []
    for i in range(len(self.groups)):
      texts.append(self.groupNames[i])
      commands.append(lambda n=i: self.viewGroup(n))
     
    buttonList = UtilityButtonList(guiFrame, commands=commands,
                                   texts=texts, grid=(2,0),
                                   tipTexts=self.tipTexts,
                                   doClone=False)
    
    
  def viewGroup(self,i):
  
    peaks = self.groups[i]
    self.guiParent.viewPeaks(peaks)
