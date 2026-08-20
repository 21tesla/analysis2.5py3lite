"""
======================COPYRIGHT/LICENSE START==========================

makeAnalysisDocRst.py: Part of the CcpNmr Analysis program

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
import cgi  # just to escape characters in output
import os
import shutil
import sys

RST_PATH = 'ccpnmr/analysis/doc/source'

# TTD
# Disable UtilityButtonList clone/help/close
# Some entries get value not nearest label
# Pulldown sub-types like fontlist
# Make popup help buttons work
# Combine with Sphix doc gen

# characters defining what is above and below text
headingDecoration = {
  1: ('=', '='),
  2: (None, '='),
  3: (None, '-'),
  4: (None, '~'),
}

# looks like Tkinter.Canvas is used sometimes instead of memops.gui.Canvas
from Tkinter import Canvas, Radiobutton

from memops.universal.Io import getPythonDirectory

from memops.gui.BasePopup import BasePopup
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.CheckButtons import CheckButtons
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileSelect
from memops.gui.Frame import Frame
#from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.LabeledEntry import LabeledEntry
from memops.gui.Menu import Menu
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.PulldownList import PulldownList
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.RadioButtons import RadioButtons
from memops.gui.Scale import Scale
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.Scrollbar import Scrollbar
from memops.gui.TabbedFrame import TabbedFrame

from memops.editor.SaveProjectFrame import SaveProjectFrame

from ccp.gui.DataLocationFrame import DataLocationFrame

from ccpnmr.analysis.frames.PeakTableFrame import PeakTableFrame
from ccpnmr.analysis.frames.ResonanceFrame import ResonanceFrame

from ccpnmr.analysis.popups.WindowPopup import WindowPopup

unknownName = '???'
UNICODE = type(u'')
FLOAT = type(0.0)
INT = type(1)
MISSING = '*Documentation missing*'
EDITABLE = ' *(Editable)*'

# TBD: need to keep doc for non-English languages

def makeAnalysisDocMacro(argServer):

  top = argServer.parent
  makeAnalysisDoc(top)

def getModuleName(popup):

  moduleName = popup.__class__.__module__
  n = moduleName.rfind('.')
  if n >= 0:
    moduleName = moduleName[n+1:]
  return moduleName

def getUnitName(popup):

  return popup.__class__.__name__

def getTitle(popup):

  title = popup.title()
  n = title.find('(')
  if n >= 0:
    title = title[:n].strip()

  return title

def getWidgetParent(widget):

  if hasattr(widget, 'parent'):
    return widget.parent
  elif hasattr(widget, 'master'):
    return widget.master
  else:
    return None

def getWidgetName(widget):

  parent = getWidgetParent(widget)
  while parent:
    for key in parent.__dict__.keys():
      if getattr(parent, key) is widget:
        return key
    parent = getWidgetParent(parent)

  return unknownName

def findWidgetOfClass(widget, clazz):

  parent = getWidgetParent(widget)
  while parent and not isinstance(parent, clazz):
    parent = getWidgetParent(parent)

  return parent

def findWidgetOfParentClass(widget, clazz):

  parent = getWidgetParent(widget)
  while parent and not isinstance(parent, clazz):
    widget = parent
    parent = getWidgetParent(parent)

  if not parent:
    widget = None

  return widget

def guessScrolledMatrixText(widget):
  
  
  n = 0
  words = []
  for i, heading in enumerate(widget.headingList):
    heading = heading.replace('\n',' ')
    words.append(heading)
    n += len(heading)

  text = ' | '.join(words)
  
  try:
    text = text.encode("utf-8")
  except UnicodeDecodeError:
    text = text.encode("utf-16")
  
  return text

def getScrolledMatrixTable(text, widget, index):

  headingList = widget.headingList
  editGetCallbacks = widget.editGetCallbacks
  ncols = len(headingList)
  headings = [' '.join(heading.split()) for heading in headingList]

  if hasattr(widget, 'tipTexts'):
    tipTexts = widget.tipTexts
  else:
    tipTexts = [MISSING] * ncols
    
  title = '**Table %d**\n' % (index) #, text)
  w3 = len(EDITABLE)

  w1 = 11
  w2 = len(MISSING)
  for i, heading in enumerate(headings):
    heading = '*' + heading + '*'
    n = len(heading)
    if n > w1:
      w1 = n
    
    tipText = widget.tipTexts[i]
    if not tipText:
      continue
    
    m = len(tipText) + w3
    if m > w2:
      w2 = m
  
  w2 = max(w2, len(title)-w1)
  
  div = '=' * w1 + '  ' + '=' * w2 + '\n' 
    
  format = '%%%d.%ds' % (w1,w1)
  
  line = div
  line += title
  line += '-' * (max(len(title)-1, len(div))-1)
  line += '\n'
  
  for n in range(ncols):
    documentation = tipTexts[n] or MISSING
    heading = headings[n]
    
    if (not heading) and documentation == MISSING:
      continue
    
    
    heading = '*' + heading + '*'
    heading = format % heading
    
    try:
      heading = heading.encode("utf-8")
    except UnicodeDecodeError:
      heading = heading.encode("utf-16")
    
    try:
      documentation = documentation.encode("utf-8")
    except UnicodeDecodeError:
      documentation = documentation.encode("utf-16")
    
    isEditable = ((len(editGetCallbacks) > n) and editGetCallbacks[n] is not None) and EDITABLE or ''
    #writeString(fp, ":%s: %s %s" % (name or '?', heading, isEditable), indent)
     
    line += '%-s  %s %s\n' % (heading or '?', documentation, isEditable)
  
  line += div
  line += '\n'
 
  return line


def getWidgetText(widget):
  
  # WHY DOESN'T THIS WORK NORMALLY ??!!
  if isinstance(widget, LabelDivider):
    return widget.label.cget('text')
  

  parent = getWidgetParent(widget)
  if isinstance(parent, TabbedFrame) and widget in parent.frames:
    n = parent.frames.index(widget)
    text = parent.options[n].strip()
  else:
    try:
      text = widget.cget('text')
      if type(text) == type(()):
        text = ' '.join(text)
      text = text.strip()
    except Exception, e:
      text = ''

    if hasattr(widget, 'text'):
      text = widget.text

    if not text and hasattr(widget, 'docKey'):
      text = widget.docKey
  
    if not text:
      for clazz in (Radiobutton, RadioButtons, CheckButton, CheckButtons,
                    PulldownList, PulldownMenu, Entry, Scale, PartitionedSelector):
        if isinstance(widget, clazz):
          ww = findWidgetToLeft(widget)
          if not ww:
            ww = findWidgetToRight(widget)
          if ww:
            if isinstance(ww, Label):
              text = getWidgetText(ww)
            elif isinstance(ww, Button):
              text = getWidgetText(ww) + ' entry'
          break
      else:
        if isinstance(widget, LabeledEntry):
          text = getWidgetText(widget.label)
        elif isinstance(widget, LabelFrame) and hasattr(widget, 'label'):
          text = getWidgetText(widget.label)
        elif isinstance(widget, LabelDivider) and hasattr(widget, 'label'):
          text = getWidgetText(widget.label)
          print "QQQ"
        elif isinstance(widget, ScrolledMatrix):
          text = guessScrolledMatrixText(widget)
  
  
  if type(text) is UNICODE:
    text = text.encode( "utf-8" )
  elif type(text) is FLOAT:
    text = str(text)
  elif type(text) is INT:
    text = str(text)
  elif text is None:
    text = '*None*'
   
  if text.endswith(':'):
    text = text[:-1]

  text = text.replace('\n', ' ')
  text = text.strip()

  widget.widgetText = text

  return text

def getWidgetClass(widget):

  parent = getWidgetParent(widget)
  if isinstance(parent, TabbedFrame) and widget in parent.frames:
    name = 'Tab'
  elif isinstance(widget, ScrolledMatrix):
    name = 'Table'
  else:
    name = widget.__class__.__name__

  if name == 'Radiobutton':
    # this is the one class we use the Tkinter widget rather than our own
    # but want to keep our convention for the spelling
    name = 'RadioButton'

  return name

def writeString(fp, s, indent = 0):

  if type(s) is UNICODE:
    s = s.encode( "utf-8" )

  fp.write('%s%s\n' % (indent*' ', s))

def writeHeading(fp, text, level=1):

  (top, bottom) = headingDecoration[level]

  if top:
    ss = len(text) * top
    writeString(fp, ss)

  writeString(fp, text)

  if bottom:
    ss = len(text) * bottom
    writeString(fp, ss)

  fp.write('\n')

def cmpChildren(child1, child2):

  d1 = child1.grid_info()
  d2 = child2.grid_info()

  if d1 and d2:
    for key in ('row', 'column'):
      v1 = d1.get(key)
      v2 = d2.get(key)
      if v1 is not None and v2 is not None:
        c = cmp(int(v1), int(v2))
        if c:
          return c

  for key in ('winfo_y', 'winfo_x'):
    f1 = getattr(child1, key)
    f2 = getattr(child2, key)
    c = cmp(f1(), f2())
    if c:
      return c

  return cmp(child1, child2)

def getWidgetChildren(widget):

  if isinstance(widget, LabeledEntry):
    return []

  children = widget.children.values()
  
  # below condition excludes popups
  children = [child for child in children if hasattr(child, 'grid_info')]
  children.sort(cmpChildren)

  return children

def getGriddedWidgetChildren(widget):

  if isinstance(widget, LabeledEntry):
    return []

  children = widget.children.values()
  # below condition excludes popups
  children = [child for child in children if hasattr(child, 'grid_info') and child.grid_info().get('row') and child.grid_info().get('column')]
  children.sort(cmpChildren)

  return children

def findWidgetToLeft(widget):

  parent = getWidgetParent(widget)
  children = getGriddedWidgetChildren(parent)
  n = children.index(widget)
  toLeft = None
  if n > 0:
    ww = children[n-1]
    # now check that to left and not something else
    d1 = widget.grid_info()
    d2 = ww.grid_info()

    if d1 and d2:
      v1 = d1.get('row')
      v2 = d2.get('row')
      if v1 is not None and v2 is not None and v1 == v2:
        toLeft = ww

  return toLeft

def findWidgetToRight(widget):

  parent = getWidgetParent(widget)
  children = getGriddedWidgetChildren(parent)
  n = children.index(widget)
  toRight = None
  if n >= 0 and n < len(children)-1:
    ww = children[n+1]
    # now check that to right and not something else
    d1 = widget.grid_info()
    d2 = ww.grid_info()

    if d1 and d2:
      v1 = d1.get('row')
      v2 = d2.get('row')
      if v1 is not None and v2 is not None and v1 == v2:
        toRight = ww

  return toRight

def findWidgetAbove(widget):

  above = None
  d1 = widget.grid_info()
  if d1:
    r1 = d1.get('row')
    if r1 is not None:
      r1 = int(r1)
      if r1 > 0:
        c1 = int(d1.get('column'))
        parent = getWidgetParent(widget)
        children = getGriddedWidgetChildren(parent)
        for child in children:
          if child is not widget:
            d2 = child.grid_info()

            if d2:
              r2 = d2.get('row')
              if r2 is not None:
                r2 = int(r2)
                if r2 == r1-1:
                  c2 = int(d2.get('column'))
                  if c1 == c2: # is this sensible?
                    above = child
                    break

  return above

def widgetHasEntry(widget):

  if not widget.winfo_viewable():
    return False

  for clazz in (ButtonList, Canvas, FileSelect, Label, Menu, Scrollbar,
                SaveProjectFrame, TabbedFrame, DataLocationFrame,
                PeakTableFrame, ResonanceFrame):
    if isinstance(widget, clazz):
      return False

  # cannot do isinstance check because that removes too much (e.g. ScrolledMatrix)
  # but below might be too weak
  if widget.__class__.__name__ in ('Frame',):
    parent = getWidgetParent(widget)
    if not isinstance(parent, TabbedFrame) or not widget in parent.frames:
      return False

  # the below is aimed at buttons in UtilityButtonList which just have images, no text
  name = getWidgetName(widget)
  text = getWidgetText(widget)
  if name == unknownName and not text:
    return False

  return True

def processChildren(fp, widget, widgetsDone, indent, tableCount):
  
  children = getWidgetChildren(widget)

  for child in children:

    tableCount = makeWidgetDoc(fp, child, widgetsDone, indent, tableCount, widget)

  return tableCount

def findDocText(obj):

  doc = obj.findDoc()
  if doc:
    doc.used = True
    text = doc.text
  else:
    text = ''

  return text


def makeWidgetDoc(fp, widget, widgetsDone, indent=6, tableCount=1, parent=None):
  

  if widget in widgetsDone:
    return tableCount

  hasEntry = widgetHasEntry(widget)

  hasIndent = False
  if hasEntry:
    widgetClass = getWidgetClass(widget)
    widgetsDone.add(widget)
    name = getWidgetName(widget)
    text = getWidgetText(widget)
    text = text.replace('\n', ' ')
    
    if (name == unknownName) and text:
      name = text
    if hasattr(widget, 'tipText'):
      if widgetClass in ('LabelFrame','LabelDivider'):
        documentation = widget.label.tipText or MISSING
      else:
        documentation = widget.tipText or MISSING
    else:
      documentation = '' # Not capable of being documented
      
    if widgetClass == 'Tab':
      writeHeading(fp, name, level=2)
      writeString(fp, '%s' % (documentation))
    
    elif widgetClass == 'Table':
      lines = getScrolledMatrixTable(text, widget, tableCount)
      writeString(fp, lines)
      tableCount += 1

    elif widgetClass  == 'RadioButtons':
      if documentation != MISSING:
        writeString(fp, '|radio| |radio| **%s**: %s' % (text or widgetClass, documentation))
    elif widgetClass == 'CheckButtons':
      if documentation != MISSING:
        writeString(fp, '|check| |check| **%s**: %s' % (text or widgetClass, documentation))
    elif widgetClass == 'PartitionedSelector':
      if documentation != MISSING:
        writeString(fp, '|selector| **%s**: %s' % (text or widgetClass, documentation))

    elif widgetClass == 'ValueRamp':
      writeString(fp, '|ramp| %s' % (documentation))
    elif widgetClass == 'PulldownList':
      writeString(fp, '|pulldown| **%s**: %s' % (text or widgetClass, documentation))
    elif widgetClass == 'PulldownMenu':
      writeString(fp, '|pulldown| **%s**: %s' % (text or widgetClass, documentation))
    elif widgetClass == 'RadioButton':
      if parent and (documentation == MISSING) and (getWidgetClass(parent) == 'RadioButtons'):
        pass

      else:
        writeString(fp, '|radio| **%s**: %s' % (text or widgetClass, documentation))
    
    elif widgetClass == 'CheckButton':
      if parent and (documentation == MISSING) and (getWidgetClass(parent) == 'CheckButtons'):
        pass

      else:
        writeString(fp, '|check| **%s**: %s' % (text or widgetClass, documentation))
    
    elif widgetClass == 'Spacer':
      pass
    elif widgetClass == 'FloatEntry':
      writeString(fp, '|float| **%s**: %s' % (text or widgetClass, documentation))
    elif widgetClass == 'IntEntry':
      writeString(fp, '|int| **%s**: %s' % (text or widgetClass, documentation))
    elif widgetClass == 'Entry':
      writeString(fp, '|entry| **%s**: %s' % (text or widgetClass, documentation))
    elif widgetClass == 'Button':
      writeString(fp, '|button| **%s**: %s' % (text or widgetClass, documentation))
    elif widgetClass == 'Text':
      writeString(fp, '|entry| %s' % (documentation))
    elif widgetClass == 'LabelDivider':
      writeHeading(fp, text, level=4)
      if documentation != MISSING:
        writeString(fp, '%s\n' % ( documentation))
    elif widgetClass == 'LabelFrame':
      writeHeading(fp, text, level=4)
      if documentation != MISSING:
        writeString(fp, '%s\n' % ( documentation))
    elif widgetClass == 'MultiWidget':
      if widget.widgetType == 'CheckButton':
        writeString(fp, '|check| **%s**: %s' % (text or widgetClass, documentation))
      elif widget.widgetType == 'PulldownList':
        writeString(fp, '|pulldown| **%s**: %s' % (text or widgetClass, documentation))
      elif widget.widgetType == 'FloatEntry':
        writeString(fp, '|float| **%s**: %s' % (text or widgetClass, documentation))
      elif widget.widgetType == 'IntEntry':
        writeString(fp, '|int| **%s**: %s' % (text or widgetClass, documentation))
      elif widget.widgetType == 'Entry':
        writeString(fp, '|entry| **%s**: %s' % (text or widgetClass, documentation))
    elif isinstance(widget, Frame):
      if documentation != MISSING:
        writeString(fp, '**%s**: %s' % (text or widgetClass, documentation))
 
    else:  
      writeString(fp, '%s: %s: %s' % (widgetClass, text, documentation))
    
    fp.write('\n')
  
  if isinstance(widget, TabbedFrame):
    tableCount = processChildren(fp, widget.sideFrame, widgetsDone, indent, tableCount)
    for n in range(widget.numTabs):
      widget.select(n)
      widget.update_idletasks()
      tableCount = processChildren(fp, widget, widgetsDone, indent, tableCount)
  elif not isinstance(widget, PartitionedSelector):
    tableCount = processChildren(fp, widget, widgetsDone, indent, tableCount)

  return tableCount

def startFile(title, fileName, subdirectory):

  docDir = getDocDir()

  path = os.path.join(docDir, subdirectory, fileName + '.rst')
  fp = open(path, 'w')
  writeHeading(fp, title, level=1)
  writeString(fp, '\n.. |pulldown| image:: ../images/pulldown.png\n   :align: bottom\n') 
  writeString(fp, '\n.. |check| image:: ../images/check.png\n   :align: bottom\n') 
  writeString(fp, '\n.. |radio| image:: ../images/radio.png\n   :align: bottom\n') 
  writeString(fp, '\n.. |float| image:: ../images/float.png\n   :align: bottom\n') 
  writeString(fp, '\n.. |int| image:: ../images/int.png\n   :align: bottom\n') 
  writeString(fp, '\n.. |entry| image:: ../images/entry.png\n   :align: bottom\n') 
  writeString(fp, '\n.. |button| image:: ../images/button.png\n   :align: bottom\n') 
  writeString(fp, '\n.. |ramp| image:: ../images/ramp.png\n   :align: bottom\n') 
  writeString(fp, '\n.. |selector| image:: ../images/selector.png\n   :align: bottom\n') 

  return fp

def endFile(fp):

  fp.close()

def getDocDir():

  docDir = os.path.join(getPythonDirectory(), RST_PATH)
  if not os.path.exists(docDir):
    os.makedirs(docDir)

  return docDir

from re import match

def makePopupDoc(popup):

  tableCounter = 1

  fileName = getUnitName(popup)
  title = getTitle(popup)
  fp = startFile(title, fileName, 'popups')
  
  docString = popup.__doc__ or '*Main Documentation Missing*'
  
  indent = 8
  lines = docString.split('\n')
  for line in lines:
    found = match(r'(\s+)\S+', line)
    if found and (len(found.group(1)) < indent):
      indent = len(found.group(1))
   
  for i, line in enumerate(lines):
    lines[i] = line[indent:].rstrip()
  
  docString = '\n'.join(lines)
  
  writeString(fp, docString+'\n')
  writeHeading(fp, 'Main Panel', level=2)
  
  widgetsDone = set()
  children = getWidgetChildren(popup)
  for child in children:
    tableCounter = makeWidgetDoc(fp, child, widgetsDone, 6, tableCounter)

  endFile(fp)

def makeAnalysisPopupDoc(top):

  # TBD: if there is more than one popup in module then have problem because
  # documentation would end up in same file so would trample each other
  popupDict = top.popups
  moduleDict = {}
  doneWindowPopup = False
  for popup in popupDict.values():
    isWindowPopup = isinstance(popup, WindowPopup)
    
    if (isWindowPopup and not doneWindowPopup) or (not isWindowPopup):
      moduleName = getModuleName(popup)
      if not moduleDict.has_key(moduleName):
        moduleDict[moduleName] = []
      moduleDict[moduleName].append(popup)

  for moduleName in sorted(moduleDict.keys()):
    print 'working on module %s' % moduleName
    print moduleDict[moduleName]
    for popup in moduleDict[moduleName]:
      ###if not isinstance(popup, WindowPopup):
      makePopupDoc(popup)

  print "makeAnalysisPopupDoc ALL DONE"

def menuItemHasEntry(menu_item):

  if not menu_item.label:
    return False

  if menu_item.kind == 'separator':
    return False

  parent = menu_item.parent
  label = parent.label
  menu_items = parent.menu_items
  n = len(menu_items)
  ind = menu_items.index(menu_item)

  # need to remove windows from Windows menu
  # below a hack

  NUM_WINDOW_ITEMS = 5
  WINDOW_LABEL = 'Windows'
  if label == WINDOW_LABEL:
    if ind >= NUM_WINDOW_ITEMS:
      return False

  # need to remove windows-specific items from Navigate menu
  # below a hack

  NUM_NAVIGATE_ITEMS = 2
  NAVIGATE_LABEL = 'Navigate'
  if label == NAVIGATE_LABEL:
    if (n-ind) > NUM_NAVIGATE_ITEMS:
      return False

  return True

def makeMenuItemDoc(fp, menu_item, level):

  hasEntry = menuItemHasEntry(menu_item)

  if hasEntry:
    if isinstance(menu_item, Menu):
      makeMenuDoc(fp, menu_item, level+1)
    else:
      tipText = menu_item.parent.tipTexts.get(menu_item.label) or '*Menu Documentation Missing*'
      writeHeading(fp, 'Menu %s' % menu_item.label, level=level)
      writeString(fp, tipText+'\n')

def makeMenuDoc(fp, menu, level=2):

  for menu_item in menu.menu_items:
    makeMenuItemDoc(fp, menu_item, level)

def makeAnalysisWindowMenuDoc(top):

  popups = top.getWindowPopups()
  if popups:
    popup = popups[0]
    menus = popup.windowFrames[0].scrolled_window.menu.menu_items
    title = fileName = 'WindowMenus'
    fp = startFile(title, fileName, 'menu')
    for menu in menus:
      if isinstance(menu, Menu): # might be a Separator so skip unless a Menu
        makeMenuDoc(fp, menu)
    endFile(fp)

    print "makeAnalysisWindowMenuDoc ALL DONE"
  else:
    print 'cannot do window menus because no window'

def makeAnalysisDoc(top):

  makeAnalysisPopupDoc(top)
  makeAnalysisWindowMenuDoc(top)

