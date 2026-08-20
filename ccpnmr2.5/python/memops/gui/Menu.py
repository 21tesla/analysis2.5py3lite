
"""
======================COPYRIGHT/LICENSE START==========================

Menu.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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
import Tkinter

from memops.gui.Base import Base
from memops.gui.ToolTip import ToolTip

class SimpleMenuItem:

  def __init__(self, parent, kind, index=None, **kw):
    
    # Added insert option (with inde) RHF Feb 2010

    self.parent = parent
    self.kind = kind
    self.label = kw.get('label')
    
    if index is None:
      Tkinter.Menu.add(parent, kind, **kw)
    else:
      Tkinter.Menu.insert(parent, kind, index, **kw)
      

class Menu(Tkinter.Menu, Base):

  def __init__(self, parent, menu_items=None, include_event=False, *args, **kw):

    if (menu_items is None):
      menu_items = []
    
    self.menu_items = []
    self.submenus = []
    self.submenuDict = {}
    self.tipTexts = {}
    self.activeEntry = None

    self.font = kw.get('font')

    # cannot set any colours here otherwise overrides the MacOS defaults

    # if not kw.get('bg'):
    #   kw['bg'] = 'grey85'
    #
    # if not kw.get('borderwidth'):
    #   kw['borderwidth'] = 1
    #
    # if not kw.get('activebackground'):
    #   kw['activebackground'] = '#D0B0A0'
    #
    # if not kw.get('activeforeground'):
    #   kw['activeforeground'] = '#400000'
    #
    # if not kw.get('activeborderwidth'):
    #   kw['activeborderwidth'] = 1
      
    self.label = kw.get('label')
    self.kind = 'menu'

    Tkinter.Menu.__init__(self, parent, *args, **kw)

    self.kw = kw
    self.parent = parent

    self.include_event = include_event
    self.menu_event = None

    self.setMenuItems(menu_items)

    self.toolTip = ToolTip(self, text='')
    
    self.bind('<Motion>', self.hover)
    self.bind('<Leave>', self.clearToolTip)
    self.bind('<ButtonPress>', self.clearToolTip)
    self.bind('<Unmap>', self.clearToolTip)
    self.bind('<Destroy>', self.clearToolTip)
    #self.shortcuts = []

  # setMenuItems() argument menu_items is list of dicts
  # each dict has to have key 'kind' with value one of menu item types:
  #   'command', 'cascade', 'separator', 'checkbutton' or 'radiobutton'
  # if kind == cascade:
  #   dict has to have key 'label' with value the label for the menu
  #   dict has to have key 'submenu' with value a list of submenu_items
  #   dict optionally has key 'tearoff' (default value is 0)
  #   dict optionally has key 'shortcut' (default value is '')
  # else:
  #   the rest of the dict should be the normal options for that menu type
  #   (e.g. 'label', 'command', ...)
  #   and also optionally key 'shortcut' (default value is '')

  
  def clearToolTip(self, event):
  
    self.toolTip.text = ''
    self.toolTip.deactivate()
    self.activeEntry = None

  def hover(self, event):
    
    # Tk bug, this does not work: self.index("active")
    
    index = self.index("@%d" % (event.y))
    if index is None:
      self.toolTip.deactivate()
      self.activeEntry = None
    
    elif self.type(index) in ('command','cascade'):
      label = self.entrycget(index, 'label')
      
      if label != self.activeEntry:
        self.activeEntry = label
        tipText = self.tipTexts.get(label)
        self.toolTip.deactivate()
  
        if tipText:
          self.toolTip.text = tipText
          self.toolTip.activate()
    
    else:
      self.activeEntry = None
      self.toolTip.deactivate()
    
  def setMenuItems(self, menu_items = None):

    if (menu_items is None):
      menu_items = []

    self.deleteMenuItems()

    for item in menu_items:
      self.addMenuItem(item)

  def deleteMenuItems(self):

    self.delete(0, Tkinter.END)
    self.activeEntry = None
    self.toolTips = {}

  def delete(self, index1, *args):

    # nasty Tkinter bug which got fixed in late 2008
    # it causes a serious memory leak if you update menus a lot
    # see: http://bugs.python.org/issue1342811
    # add equivalent patch here just in case people
    # are using broken Tkinter rather than patched code
    cmds = set()

    def deleteIndex(ind):
      if ind >= 0 and ind < len(self.menu_items):
        if self.entryconfig(ind).has_key('command'):
          cmd = str(self.entrycget(ind, 'command'))
          if cmd and cmd in self._tclCommands:
            cmds.add(cmd)
        item = self.menu_items[ind]
        if item in self.submenus:
          self.submenus.remove(item)
          label = item.label
          
          del self.submenuDict[label]
          if label in self.tipTexts:
            del self.tipTexts[label]
            if label == self.activeEntry:
              self.activeEntry = None
          
          ##item.delete(0, Tkinter.END)
          item.destroy() # otherwise get memory leak
        del self.menu_items[ind]

    n = len(args)
    if n == 0:
      deleteIndex(index1)
    elif n == 1:
      index2 = args[0]
      if index2 == Tkinter.END:
        index2 = len(self.menu_items) - 1
      for ind in range(index2, index1-1, -1):
        deleteIndex(ind)

    Tkinter.Menu.delete(self, index1, *args)

    # final bit of patch for Tkinter bug
    for cmd in cmds:
      if cmd in self._tclCommands:
        self.deletecommand(cmd)

  def addMenuItem(self, item):

    kind = item['kind']
    rest = item.copy()
    del rest['kind']

    if (kind == 'cascade'):
      tearoff = rest.get('tearoff', 0)
      shortcut = rest.get('shortcut', '')
      
      kw = self.kw.copy()
      kw['tearoff'] = tearoff
      
      label = rest['label']
      items = rest['submenu']
      tipText = rest.get('tipText')
      menu = Menu(self, include_event=self.include_event, **kw)
      self.add(kind, label=label, menu=menu, shortcut=shortcut, tipText=tipText)
      menu.setMenuItems(items)

    else:

      self.add(kind, **rest)

  def setMenuEvent(self, event):

    self.menu_event = event
    for submenu in self.submenus:
      submenu.setMenuEvent(event)

  def popupMenu(self, event):

    self.setMenuEvent(event)
    self.post(event.x_root, event.y_root)

  def popdownMenu(self, *event):

    self.unpost()

  def add_command(self, shortcut='', index=None, **options):

    self.add('command', shortcut, index, **options)

  def add_cascade(self, shortcut='', index=None, **options):

    self.add('cascade', shortcut, index, **options)

  def add_radiobutton(self, shortcut='', index=None, **options):

    self.add('radiobutton', shortcut, index, **options)

  def add(self, kind, shortcut = '', index=None, **options):
    """ Added insert option, with 'index' RHF Feb 2010.
    """

    label = options.get('label')

    if (shortcut):

      assert len(shortcut) == 1, 'len(shortcut) = %d' % len(shortcut)

      n = label.find(shortcut)
      assert n >= 0, 'illegal shortcut: label = %s, shortcut = %s' % (label, shortcut)

      s = shortcut.lower()
      #assert s not in self.shortcuts, 'duplicate shortcut: label = %s, shortcut = %s' % (label, shortcut)

      #self.shortcuts.append(s)
      options['underline'] = n
 
    if (self.include_event and options.has_key('command')):

      options_copy = options.copy()
      func = options['command']
      options_copy['command'] = lambda: func(self.menu_event)
      options = options_copy

    if 'tipText' in options:
      self.tipTexts[label] = options['tipText']
      del options['tipText']


    if kind == 'cascade':
      menu = options.get('menu')
      menu.label = label
      menu.kind = 'menu'
      menu.parent = self
      self.menu_items.append(menu)
      self.submenus.append(menu)
      self.submenuDict[label] = menu
      if index is None:
        Tkinter.Menu.add(self, kind, **options)
      else:
        Tkinter.Menu.insert(self, index, kind, **options)
        
    else:
      menu_item = SimpleMenuItem(self, kind, index, **options)
      self.menu_items.append(menu_item)

if __name__ == '__main__':

  def new():
    print 'new'

  def pick():
    print 'pick'

  root = Tkinter.Tk()
 
  fontNames = ('Courier', 'Helvetica', 'Lucida', 'System', 'Times')
  BOLD = 'bold'
  ITALIC = 'italic'
  UNDERLINE = 'underline'
  fontSpec = '%s %d %s' % (fontNames[2], 24, ITALIC)
  fontSpecUnderline = '%s %d %s %s' % (fontNames[2], 24, ITALIC, UNDERLINE)
  mediumFontSpec = '%s %d %s' % (fontNames[2], 16, ITALIC)
  tinyFontSpec = '%s %d %s' % (fontNames[2], 1, ITALIC)

  menubar = Menu(root)
  menu = Menu(menubar, tearoff=0, font=mediumFontSpec)

  # cannot set any colours here otherwise overrides the MacOS defaults

  from PIL import ImageDraw, Image, ImageTk, ImageFont
  import tkFont

  tkFontNames = tkFont.families()
  menuImages = []
  x, y = 0, 0
  for ii, fontName in enumerate(sorted(tkFontNames[:20])):
    _ul = ii % 6
    _numUl = 1
    # try to load fontName
    tkF = tkFont.Font(name=fontName)
    actual = tkF.actual()
    try:
      fnt = ImageFont.truetype(fontName, 24)
    except Exception as es:
      try:
        fnt = ImageFont.truetype(actual['family'].lower(), 24)
      except Exception as es:
        print(str(es))
        continue

    # split font into three for underlining the middle group
    nameGroup = (fontName[0:_ul], fontName[_ul:_ul+_numUl], fontName[_ul+_numUl:])
    print(nameGroup)

    # make two different images for pos/neg - can't catch mouse hover on the menu though so may be redundant
    _imgW, _imgH = fnt.getsize(fontName)
    img = Image.new(mode="RGBA", size=(_imgW, _imgH+2), color=(0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    imgNeg = Image.new(mode="RGBA", size=(_imgW, _imgH+2), color=(255, 255, 255, 0))
    dNeg = ImageDraw.Draw(imgNeg)

    # get metrics to place the underline
    ascend, descend = fnt.getmetrics()

    x = 0
    for gr, txt in enumerate(nameGroup):
      if not txt:
        # skip the first group if empty
        continue

      _w, _h = fnt.getsize(txt)
      ascend, descend = fnt.getmetrics()

      d.text((x, 0), txt, fill=(0, 0, 0, 255), font=fnt)
      dNeg.text((x, 0), txt, fill=(255, 255, 255, 255), font=fnt)
      if gr == 1:
        d.line((x, min(ascend + 1, _h + 1), x+_w-1, min(ascend + 1, _h + 1)), fill=(0, 0, 0, 255))
        dNeg.line((x, min(ascend + 1, _h + 1), x + _w - 1, min(ascend + 1, _h + 1)), fill=(255, 255, 255, 255))
      x += _w

    imgTk = ImageTk.PhotoImage(img)         # was (img.crop(img.getbbox())) but have created the correct box size
    imgTkNeg = ImageTk.PhotoImage(imgNeg)
    menuImages.append((fontName, imgTk, imgTkNeg, _imgW, _imgH, ascend, descend))

  menu.add_command(label='New', shortcut='N', command=new, compound="left")
  menuItem = menu.add_command(label='Pick', shortcut='P', command=pick, font=fontSpec)
  menu.entryconfig(1, underline=0, accelerator='P')

  for ii, (item, img, imgNeg, _, _, _, _) in enumerate(menuImages[:20]):
    # alternate between black/white - testing contrast
    # spacing seem very large on Windows
    if ii % 2:
      menuItem = menu.add_command(image=img, compound='top')
    else:
      menuItem = menu.add_command(image=imgNeg, compound='top')

  menu.add_command(label='Pick2', shortcut='2', command=pick, tipText='Tip B')
  menu.add_command(label='Pick3', shortcut='3', command=pick, tipText='Tip C')
  menubar.add_cascade(label='Project', shortcut='P', menu=menu, font=fontSpec)
  menu.entryconfig(3, state = Tkinter.DISABLED)
  menu.entryconfig(4, state = Tkinter.DISABLED)

  menu = Menu(menubar, tearoff=0, font=fontSpec)
  menubar.add_cascade(label='View', menu=menu)
  root.config(menu=menubar)

  root.mainloop()
