from memops.editor.BasePopup import BasePopup

from nijmegen.cing.CingFrame import CingFrame

def testCingPopup(argServer):

  project = argServer.getProject()
  popup = CingPopup(argServer.parent)
  popup.open()


class CingPopup(BasePopup):
  """
  **Setup and Submit NMR Structure Validation Analyses to CING**
  
  This popup is used to connect the data in a CCPN project to the structure
  validation analyses that are available at the iCING server. When submitting to
  this service the CCPN project (without spectrum data files) will be sent over
  the Internet for analysis. Many different NMR and structure based analyses
  will be performed under the CING platform. These include: PROCHECK, ShiftX,
  DSSP, WHAT_CHECK, Wattos etc.

  For more information detailed information see the `CING home page
  <https://nmr.le.ac.uk/icing/>`_.
  
  **Use**

  To setup and run a CING job from Analysis the user first has to create a
  "run", which contains all the settings. A run can be made with 
  [New Run] or [Copy Run].

  The data that relates to a job/run is divided into three sub-tabs to indicate
  what the input data from the CCPN project is, what settings are used during
  the calculation, and what the output CCPN data is. The input data is further
  sub-divided into various categories of data, and the user can select  from
  pulldown menus to dictate which CCPN entities will be selected for the
  calculation.

  **Input Data**

  Typically to run a CING job the user will select the appropriate Molecular
  System from the "Input Data" panel, then select an structure ensemble in the
  "Structures sub-tab; these are minimum requirements for CING. The user may
  also a add shift lists, peak lists and restraint lists  to the analysis by
  going into the relevant sub-tab and clicking   the [Add ...] and [Remove ...]
  buttons to add the CCPN data entity that is selected at the bottom right of
  the panel.

  **Run Settings**

  The "Run Settings" tab allows the user to adjust settings that control the
  calculation job and actually submit all the data to the iCING server for
  analysis. Usually the user simply clicks [Submit Project!], assuming a run has
  been setup and the computer has an Internet connection. After submission the
  user can [Check Run Status] to get a short report, but otherwise the user will
  be informed when the job is complete in the CCPN interface. The CING analysis
  data will be available to download from the specified "Results URL" to the
  "Results File" by using [Download Results] Also, the data will be available
  for a time  on the CING web site ([View Results HTML]). To specifically force
  the iCING  server to remove all of the user's data from its system [Purge
  Server Result] can be clicked. Otherwise the results will be available on the
  CING website until the server does its next data clean-up.

  **Output Data**

  The "Output Data" table lists all of the CCPN entities that were generated or
  modified by the selected calculation job (if it were run). However, this data
  is not currently filled in by the CING server, although it will be at some
  stage in the future.

  """

  def __init__(self, parent, **kw):

    self.parent = parent

    BasePopup.__init__(self, parent=parent, title='CING Setup', **kw)
    
  def body(self, guiFrame):

    self.geometry('740x600')
  
    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    self.frame = CingFrame(guiFrame, self.parent, closeButton=True)
    self.frame.grid(row=0, column=0, sticky='nsew')
    self.frame.updateAll()
    # wb104 9 Feb 2015
    # not sure why, but need below as well as grid above or otherwise
    # just get blank dialog until it is resized
    self.after_idle(lambda: self.frame.grid(row=0, column=0, sticky='nsew'))
  
  def open(self):
  
    BasePopup.open(self)
   
   
  def close(self):
  
    BasePopup.close(self)
    
    
  def destroy(self):
  
    BasePopup.destroy(self)     
     
if __name__ == '__main__':

  print "Run testCingPopup() as a CcpNmr Analysis macro"
