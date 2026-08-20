
"""
======================COPYRIGHT/LICENSE START==========================

ScrolledHtml.py: <write function here>

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
import base64
import formatter
import htmllib
import os
import string
try:
  import urllib
except:
  print 'Warning: ScrolledHtml will only be able to open local files'
import urlparse

import Tkinter

from memops.universal.Io import splitPath

from memops.gui.MessageReporter import showError
from memops.gui.ScrolledText import ScrolledText

font_family = 'Helvetica'
font_sizes = { 'h1': 24, 'h2' : 20, 'h3': 16 }
default_size = 12

class ScrolledHtml(formatter.NullWriter, htmllib.HTMLParser, ScrolledText):

  def __init__(self, parent, startUrlCallback = None, endUrlCallback = None,
               enterLinkCallback = None, leaveLinkCallback = None, *args, **kw):

    self.startUrlCallback = startUrlCallback
    self.endUrlCallback = endUrlCallback
    self.enterLinkCallback = enterLinkCallback
    self.leaveLinkCallback = leaveLinkCallback

    self.initTags()

    self.protocol = ''
    self.location = ''
    self.path = ''
    self.dir = ''

    formatter.NullWriter.__init__(self)
    self.formatter = formatter.AbstractFormatter(self)
    htmllib.HTMLParser.__init__(self, self.formatter)
    ScrolledText.__init__(self, parent, *args, **kw)
    self.text_area.config(font=(font_family, default_size))

    self.bind('<Enter>', lambda event: self.focus())
    self.bind('<KeyPress-Prior>', lambda event: self.pageChange(-1))
    self.bind('<KeyPress-Next>', lambda event: self.pageChange(1))

  def openUrl(self, url, forceLoad = False, allowModifyPath = True):

    #print 'openUrl "%s" "%s" "%s" "%s" "%s"' % (url, self.protocol, self.location, self.path, self.dir)

    if (self.startUrlCallback):
      self.startUrlCallback(self.yview())

    self.setState(Tkinter.NORMAL)

    (protocol, location, path, query, fragment) = urlparse.urlsplit(url)
    #print 'openUrl1 "%s" "%s" "%s" "%s" "%s"' % (protocol, location, path, query, fragment)

    if (not protocol):
      protocol = self.protocol
    if (protocol == self.protocol):
      if (not location):
        location = self.location
      if (location == self.location):
        if (allowModifyPath and not path):
          allowModifyPath = False
          path = self.path

    #print 'openUrl2 "%s" "%s" "%s"' % (protocol, location, path)
    if (allowModifyPath and self.dir and path and (path[0] != '/')):
      path = self.dir + '/' + path

    #print 'openUrl3 "%s" "%s" "%s"' % (protocol, location, path)

    if (forceLoad or (protocol != self.protocol) or \
        (location != self.location) or (path != self.path)):

      (dir, rest) = splitPath(path)
      if (rest and (rest[0] == '~')):
        dir = dir + rest
      elif (dir == '/'):
        dir = ''

      self.protocol = protocol
      self.location = location
      self.path = path
      self.dir = dir

      self.fragmentIndices = {}
      self.clearTags()
      self.clear()
      self.images = {} # need to keep reference to images otherwise they get gc'ed

      try:
        if (self.path[-4:].lower() == '.gif'):
          self.loadImage(self.path)
        else:
          if (protocol == 'file'):
            u = open(path)
          else:
            nonFragmentUrl = urlparse.urlunsplit((protocol, location, path, query, ''))
            u = urllib.urlopen(nonFragmentUrl)
          t = u.read()
          u.close()
          self.feed(t)
          self.close()
      except:
        showError('Opening url', 'Could not open URL "%s"' % url)

    self.moveToFragment(fragment)
    self.setDefaultCursor()
    self.setState(Tkinter.DISABLED)

    if (self.endUrlCallback):
      fullUrl = urlparse.urlunsplit((protocol, location, path, query, fragment))
      self.endUrlCallback(fullUrl)

  def pageChange(self, n):
 
    self.text_area.yview('scroll', str(n), 'pages')

  def initTags(self):

    self.tags = 0
    self.font_tag = None
    self.link_tag = None
    self.margin_tag = None

  def setDefaultCursor(self):

    self.setCursor('left_ptr')

  def clearTags(self):

    for tag in range(self.tags):
      self.tag_delete(str(tag))

    self.initTags()

  def stripSpaceTagAdd(self, tag, start, end):

    while (1):
      #print 'stripSpaceTagAdd1', start
      c = self.text_area.get(start)
      if (c not in string.whitespace):
        break
      if (self.text_area.compare(start, '>=', end)):
        return
      start = self.text_area.index('%s+1c' % start)

    while (1):
      #print 'stripSpaceTagAdd2', end
      c = self.text_area.get('%s-1c' % end)
      if (c not in string.whitespace):
        break
      if (self.text_area.compare(start, '>=', end)):
        return
      end = self.text_area.index('%s-1c' % end)

    self.tag_add(tag, start, end)

  def enterLink(self, url):

    self.setCursor('hand2')
    if (self.enterLinkCallback):
      self.enterLinkCallback(url)

  def leaveLink(self):

    self.setDefaultCursor()
    if (self.leaveLinkCallback):
      self.leaveLinkCallback()

  def loadImage(self, src, width = 0, height = 0):

    url = self.protocol + '://' + self.location + src
    if (self.images.has_key(url)):
      image = self.images[url]
    elif (self.protocol == 'file'):
      try:
        image = Tkinter.PhotoImage(file=src)
      except:
        return
    else:
      try:
        u = urllib.urlopen(url)
        data = u.read()
        u.close()
        data = base64.encodestring(data)
        image = Tkinter.PhotoImage(data=data)
      except:
        return

    if (not self.images.has_key(url)):
      if (width and height):
        w = image.width()
        h = image.height()
        if ((width < w) and (height < h)):
          xscale = int(0.5 + float(w)/width)
          yscale = int(0.5 + float(h)/height)
          image = image.subsample(xscale, yscale)
      self.images[url] = image

    self.image_create(Tkinter.CURRENT, image=image, align=Tkinter.CENTER)

  def moveToFragment(self, fragment):

    if (not fragment):
      return

    ind = self.fragmentIndices.get(fragment)
    if (ind is not None):
      self.yview(ind)

  ##### implementation for parser #####

  def anchor_bgn(self, href, name, type):

    #print 'anchor_bgn1 "%s" "%s" "%s"' % (href, name, type)

    htmllib.HTMLParser.anchor_bgn(self, href, name, type)

    if (name):
      self.fragmentIndices[name] = self.index(Tkinter.CURRENT)
    elif (href):
      tag = self.tags
      self.tags = tag + 1
      self.link_tag = tag
      self.tag_config(tag, underline=1, foreground='blue')
      command = lambda event: self.openUrl(href)
      self.tag_bind(tag, '<Button-1>', command)
      command = lambda event: self.enterLink(href)
      self.tag_bind(tag, '<Enter>', command)
      command = lambda event: self.leaveLink()
      self.tag_bind(tag, '<Leave>', command)
 
  def anchor_end(self):

    #print 'anchor_end'
    #htmllib.HTMLParser.anchor_end(self)
    self.link_tag = None

  def handle_image(self, src, alt, *args):

    #print 'handle_image1 "%s" "%s"' % (src, alt), args

    if (args):
      (width, height) = args[2:4]
    else:
      width = height = 0

    if (self.link_tag is not None):
      start = self.index(Tkinter.CURRENT)

    try:
      if (src[-4:].lower() != '.gif'):
        raise
      if (src[0] != '/'):
        src = self.dir + '/' + src
      self.loadImage(src, width, height)
    except:
      htmllib.HTMLParser.handle_image(self, src, alt, *args)

    if (self.link_tag is not None):
      end = self.index(Tkinter.CURRENT)
      #print 'handle_image2', start, end
      #self.stripSpaceTagAdd(self.link_tag, start, end)
      self.tag_add(self.link_tag, start, end)

    # TBD: below is here because seem to have bug which means
    # that image does not break new paragraph settings
    # which means that <P> after image was being ignored
    # if there was a <P> just before the image
    self.formatter.parskip = 0
    self.formatter.hard_break = 0
    self.formatter.para_end = 0
    self.formatter.nospace = 0

  ##### implementation for writer #####

  def new_alignment(self, align):

    #print 'new_alignment', align
    pass

  def new_font(self, font):

    #print 'new_font', font
    if (font):
      (size, italic, bold, teletype) = font
      size = font_sizes.get(size, default_size)
      s = []
      if (bold):
        s.append('bold')
      if (italic):
        s.append('italic')
      style = ' '.join(s)
      #print 'new_font1', (font_family, size, style)
      tag = self.tags
      self.tags = tag + 1
      self.font_tag = tag
      self.tag_config(tag, font=(font_family, size, style))
    else:
      self.font_tag = None

  def new_margin(self, margin, level):

    #print 'new_margin', margin, level

    if (margin is None):
      self.margin_tag = None
    else:
      # TBD: assume margin == 'ul' for now
      self.margin_level = level
      tag = self.tags
      self.tags = tag + 1
      self.margin_tag = tag
      
  def new_spacing(self, spacing):

    #print 'new_spacing', spacing
    pass

  def send_line_break(self):

    #print 'send_line_break'
    self.append('\n')

  def send_paragraph(self, blankline):

    #print 'send_paragraph', blankline
    self.append(blankline*'\n')

  def send_hor_rule(self, *args, **kw):

    #print 'send_hor_rule', args, kw
    pass

  def send_flowing_data(self, data):

    #print 'send_flowing_data1: "%s"' % data

    start = self.index(Tkinter.CURRENT)
    self.append(data)
    end = self.index(Tkinter.CURRENT)

    #print 'send_flowing_data2', start, end
    # link tag has to be done differently because there seems to be
    # an unfortunate bug in the parser and we don't want white
    # space at the beginning or end of the highlighted section
    if (self.link_tag is not None):
      self.stripSpaceTagAdd(self.link_tag, start, end)
    for tag in (self.font_tag, self.margin_tag):
      if (tag is not None):
        self.tag_add(tag, start, end)

  def send_literal_data(self, data):

    #print 'send_literal_data', data
    self.append(data)

  def send_label_data(self, data):

    #print 'send_label_data', data
    if (self.margin_tag is not None):
      s = (2*self.margin_level)*' ' + data + ' '
      start = self.index(Tkinter.CURRENT)
      self.append(s)
      #print 'send_label_data1', start
      #for n in range(len(s)):
      #  print 'send_label_data2', n, s[n], self.bbox('%s+%d*c' % (start, n))
      #self.tag_config(tag, lmargin2=30)

if (__name__ == '__main__'):

  r = Tkinter.Tk()
  s = ScrolledHtml(r, xscroll=False, width=80, height=30)
  r.grid_rowconfigure(0, weight=1)
  r.grid_columnconfigure(0, weight=1)
  s.grid(sticky=Tkinter.NSEW)

  url = 'http://www.bio.cam.ac.uk'
  #url = 'file:/home/wb104/think_in_java/TIJ3_t.htm'
  s.openUrl(url)

  r.mainloop()
