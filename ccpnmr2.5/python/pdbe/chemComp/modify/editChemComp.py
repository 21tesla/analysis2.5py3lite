
from ccp.general.Util import setCurrentStore
from ccp.gui.ChemCompEditor import ChemCompEditPopup
from ccp.gui.ChemCompFrame import ChemCompFrame
from memops.api import Implementation
from memops.gui.Button import Button

if (__name__ == '__main__'):


  project = Implementation.MemopsRoot(name='edit')
  setCurrentStore(project,'ChemElementStore')

  root = Tkinter.Tk()
  root.top = root

  root.grid_rowconfigure(0, weight=1)
  root.grid_columnconfigure(0, weight=1)

  frame = ChemCompFrame(root, project, chemCompEntries = ('ChemComp',)) #path = editChemCompDataDir,
  frame.grid(sticky=Tkinter.NSEW)

  def getSelected():

    chemComp = frame.getSelectedChemComp()

    ChemCompEditPopup(root, chemComp)

  button = Button(root, text='Edit selected', command=getSelected)
  button.grid()

  root.mainloop()

