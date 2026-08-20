from sys import platform

from os import listdir, path, R_OK, access

from ccpnmr.analysis.popups.BasePopup import BasePopup
from Tkinter import PhotoImage

import mimetypes

from memops.gui.Tree import Tree
from memops.gui.Frame import Frame
from memops.gui.Button import Button


from memops.universal.Io import getTopDirectory


class FileBrowserPopup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
    BasePopup.__init__(self, parent=parent, title='File Selection')

                       

  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = FileBrowserFrame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('300x450')
    
    self.update_idletasks()
 
class FileBrowserFrame(Frame):

  def __init__(self, guiParent, basePopup, callback=None):

    # Base popup required to handle notification of data model changes
    # e.g. new peak lists, so that the GUI can update to the latest
    # state
    self.basePopup = basePopup
    self.guiParent = guiParent
 
    # we need to be able to define a callback before exit
    self.callback = callback

    Frame.__init__(self, guiParent)

    # Basic principle is that we want to supply a root directory
    # and we will get a tree that will allow us to browse anything
    # below the directory. We will also have a button for going
    # up a directory.

    self.rootDir = '/home/jionides/work/CCPN/test_WMS_archive'
  
    self.icon_dict = {'none':'text-x-generic-template',
                      'text':'text-x-generic',
                      'text/plain':'text-x-generic',
                      'text/html':'text-html',
                      'text/script':'text-x-script',
                      'text/x-python':'text-x-script',
                      'text/x-perl':'text-x-script',
                      'image':'image-x-generic',
                      'application':'applications-other',
                      'video':'video-x-generic',
                      'message':'internet-mail',
                      'audio':'audio-x-generic',
                      'chemical':'applications-internet',
                      'application/x-tar':'package-x-generic',
                      'application/msword':'x-office-document',
                      'application/postscript':'x-office-drawing',
                      'application/pdf':'x-office-drawing',
                      'application/xml':'text-x-script',
                      'application/x-python-code':'application-x-executable',}

    self.grid_rowconfigure(1, weight=1)
    self.grid_columnconfigure(1, weight=1)
    
    self.fileTree = Tree(self)
    self.fileTree.grid(row=1, column=1, padx=10, pady=10, sticky='nsew')

    self.ok_button = Button(self, width=10, text='OK', command=self.isOk)
    #self.ok_button = Button(self, text='OK', command=self.isOk)
    self.ok_button.grid(row=2, column=1, padx=10, pady=10, sticky='e')
    
    self.cancel_button = Button(self, width=10, text='Cancel', command=self.isCancel)
    #self.cancel_button = Button(self,text='Cancel', command=self.isCancel)
    self.cancel_button.grid(row=2, column=1, padx=10, pady=10, sticky='w')

    files = listdir(self.rootDir)   
    parents = [None] * len(files)
    icons = []
    objects = []
    texts = []
    callbacks = []
    for file in files:
      file = path.join(self.rootDir, file)
    
      text, icon, isDir = self.getData(file)
      icons.append(icon)
      texts.append(text)    
      objects.append(file)
    
      if isDir:
        callbacks.append(self.openDir)
      else:
        callbacks.append(None)

    print 'updating', texts
    
    self.fileTree.update(parents, objects, texts, icons, callbacks)
               

  def getData(self, file):
    
    dirr, text = path.split(file)
    
    readable = access(file, R_OK)
    
    isDir = False
    if path.isdir(file):
      
      if not readable:
        icon = 'emblem-unreadable'
      elif path.islink(file):
        isDir = True
        icon = 'folder-remote'
      else:
        isDir = True
        icon = 'folder'
      text += '/'  
        
    elif path.isfile(file):
      
      if not readable:
        icon = 'emblem-unreadable'
      #elif path.islink(file):
      #  icon = 'text-html'
      else:
        mimeType = mimetypes.guess_type(file)[0]
        
        if mimeType:
          if not self.icon_dict.get(mimeType):
            print mimeType
        
          icon = self.icon_dict.get(mimeType) \
                  or self.icon_dict.get(mimeType.split('/')[0]) \
                  or self.icon_dict['none']
        else:
          icon = self.icon_dict['none']
    
    elif not readable:
      icon = 'emblem-unreadable'
      
    else:
      icon = 'emblem-system'
    
    return text, icon, isDir


  def openDir(self, node):
    
    tree = node.tree
    filePath = node.object
    files = listdir(filePath)
        
    for file in files:
       
       file = path.join(filePath, file)
       text, icon, isDir = self.getData(file)
       parent, null = path.split(file)
  
       if isDir and access(file, R_OK):
         callback = self.openDir
       else:
         callback = None
           
       tree.add(node, file, text, icon, callback)

  def isCancel(self):

    print 'should be exiting'
    self.quit()

  def isOk(self):

    print 'setting selected ', self.fileTree.getSelected()

  def destroy(self):

    print 'in main destroy method'
    Frame.destroy(self)
    
  def quit(self):
  
    print 'in main quit method'
    self.guiParent.parent.destroy()

if __name__ == '__main__':

  import sys
  import Tkinter

  root = Tkinter.Tk()
  root.withdraw()
   
  popup = FileBrowserPopup(root)
 
  root.mainloop()
  

