

from ccpnmr.analysis.popups.BasePopup import BasePopup
from gottingen.PalesFrame import PalesFrame

class PalesPopup(BasePopup):


  def __init__(self, parent, *args, **kw):

    self.calcStore = None
    self.waiting = False

    BasePopup.__init__(self, parent, title="Data Analysis : PALES ", **kw)

  def body(self, guiFrame):

    self.geometry('700x600')

    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)
    
    self.frame = PalesFrame(self, self.project, grid=(0,0))
