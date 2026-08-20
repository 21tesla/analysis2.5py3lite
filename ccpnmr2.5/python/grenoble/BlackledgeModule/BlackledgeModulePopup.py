
from ccpnmr.analysis.popups.BasePopup import BasePopup
from grenoble.BlackledgeModule.BlackledgeModuleFrame import BlackledgeModuleFrame

class BlackledgeModulePopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.calcStore  = None
    self.waiting    = False

    BasePopup.__init__(self, parent, title="Data Analysis : MODULE ", **kw)

  def body(self, guiFrame):

    self.geometry('700x600')

    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)
    
    self.frame = BlackledgeModuleFrame(self, self.project, grid=(0,0))

def launchModulePopup():

  import Tkinter

  global top

  root = Tkinter.Tk()
  root.withdraw()
  top  = BlackledgeModulePopup(root)

  top.update_idletasks()


if __name__ == '__main__':

  launchModulePopup()

