from memops.gui.Frame           import Frame
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.Util            import createDismissHelpButtonList
from memops.editor.BasePopup    import BasePopup

from IsdFrame import IsdFrame

#
# Funtion to load in as CCPN Macro
#   
def testIsdPopup(argServer):

  popup = IsdPopup(argServer.parent, argServer.getProject())
  popup.open()

   
class IsdPopup(BasePopup):

  def __init__(self, parent, ccpnProject):

    self.parent      = parent
    self.ccpnProject = ccpnProject
    
    BasePopup.__init__(self, parent=parent, title='ISD Parameter Setup')
                       
  def body(self, guiFrame):

    # Ensure that the first row and column in popup expand
    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1, minsize=300)

    self.isdFrame = IsdFrame(guiFrame, self.ccpnProject)
    self.isdFrame.grid(row=0, column=0, sticky='nsew')
    self.parent.isdFrame = self.isdFrame
      
  def destroy(self):
  
    BasePopup.destroy(self)
