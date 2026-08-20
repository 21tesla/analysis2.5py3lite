"""
======================COPYRIGHT/LICENSE START==========================

Ticks.py: Utility code for CCPN code generation framework

Copyright (C) 2005  (CCPN Project)

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

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""
import math
ceil = math.ceil
floor = math.floor
log = math.log

EPS = 1.0e-6

class Ticks:

  # width and height should be in terms of screen pixels even for PS output
  def __init__(self, width, height, isHorizontal=True, region=None,
               font=None, fill='black', maxMajor=5, maxMinor=25,
               tickScale=None, doRuleLines=True, haveYInversion=True,
               haveTicksOnZeroSide=True, minDecimalPlaces=0, isTextDrawn=True,
               deltaMajor=None, deltaMinor=None, numberDecimals=None,
               tickMajor=None, tickMinor=None):

    if not font:
      font = ('Helvetica', 10)

    self.width = width
    self.height = height
    self.isHorizontal = isHorizontal
    self.region = region
    self.font = font
    self.fill = fill
    self.maxMajor = maxMajor
    self.maxMinor = maxMinor
    self.tickScale = tickScale
    self.doRuleLines = doRuleLines
    self.haveYInversion = haveYInversion
    self.haveTicksOnZeroSide = haveTicksOnZeroSide
    self.minDecimalPlaces = minDecimalPlaces
    self.isTextDrawn = isTextDrawn
    self.deltaMajor = deltaMajor
    self.deltaMinor = deltaMinor
    self.numberDecimals = numberDecimals
    self.tickMajor = tickMajor
    self.tickMinor = tickMinor

    self.updateTicks()

  # drawLine should be overridden in subclass
  def drawLine(self, x0, y0, x1, y1):

    pass

  # drawText should be overridden in subclass
  def drawText(self, text, x, y, a, b):

    pass

  # convenience function:
  # should use same length for canvas and PS to be consistent
  # so always use canvas length, which is what is stored in API
  def setMaxMajorMinor(self, length):

    self.maxMajor = length / 100
    self.maxMinor = length / 20

  def setRegion(self, region):

    self.region = region
    self.updateTicks()

  def updateTicks(self):

    if self.region:
      (r0, r1) = self.region
      if abs(r1-r0) > EPS:
        self.findTicks()
        self.createTicks()
        self.drawTicks()

  def convertCoord(self, r, r0, r1, d, haveInversion = False):
 
    c = d * (r-r0) / (r1-r0)
 
    if (haveInversion):
      c = d - c
 
    return c

  def createTicks(self):
 
    w = self.width
    h = self.height
    (r0, r1) = self.region
    isTextDrawn = self.isTextDrawn
 
    #print 'createTicks', w, h, self.minorTicks, self.majorTicks, self.majorFormat
 
    # 30 Mar 2010: complete hack to get PDF output to look ok
    isPdf = False
    if hasattr(self, 'outputHandler'):
      if self.outputHandler.__class__.__name__ == 'Pdf':
        isPdf = True

    if (self.isHorizontal):
 
      if (self.doRuleLines): # meant for canvases
        self.ruleLine1 = (0, 0, w, 0)
        self.ruleLine2 = (0, h-1, w, h-1)
 
      s = self.tickMinor or ((self.tickScale or h) // 5)
      if (self.haveTicksOnZeroSide):
        y0 = 0
        y1 = s
      else:
        y0 = h
        y1 = h - s
      self.minorLines = []
      for tick in self.minorTicks:
        x = self.convertCoord(tick, r0, r1, w)
        self.minorLines.append((x, y0, x, y1))
 
      s = self.tickMajor or (2 * s)
      if (self.haveTicksOnZeroSide):
        y0 = 0
        y1 = s
        yt = y1 + 2
        if isPdf:
          at = 0.5
        else:
          at = 0.4
        if (self.haveYInversion):
          bt = 1.0
        else:
          bt = 0.0
      else:
        y0 = h
        y1 = h - s
        if isPdf:
          yt = y1
          at = 0.5
        else:
          yt = y1 + 2
          at = 0.4
        if (self.haveYInversion):
          bt = 0.0
        else:
          bt = 1.0
      self.majorLines = []
      self.majorTexts = []
      for tick in self.majorTicks:
        x = self.convertCoord(tick, r0, r1, w)
        self.majorLines.append((x, y0, x, y1))
        if isTextDrawn:
          text = self.doMajorFormat(tick)
          self.majorTexts.append((text, x, yt, at, bt))
 
    else: # not self.isHorizontal
 
      if (self.doRuleLines): # meant for canvases
        self.ruleLine1 = (0, 0, 0, h)
        self.ruleLine2 = (w-1, 0, w-1, h)
 
      s = self.tickMinor or ((self.tickScale or w) // 5)
      if (self.haveTicksOnZeroSide):
        x0 = 0
        x1 = s
      else:
        x0 = w
        x1 = w - s
      self.minorLines = []
      for tick in self.minorTicks:
        y = self.convertCoord(tick, r0, r1, h, self.haveYInversion)
        self.minorLines.append((x0, y, x1, y))
 
      s = self.tickMajor or (2 * s)
      if (self.haveTicksOnZeroSide):
        x0 = 0
        x1 = s
        xt = x1 + 2
        at = 0.0
        bt = 0.3
      else:
        x0 = w
        x1 = w - s
        if isPdf:
          xt = x1 - 2
          at = 1.2
        else:
          xt = x1 + 2
          at = 1.0
        bt = 0.3
      self.majorLines = []
      self.majorTexts = []
      for tick in self.majorTicks:
        y = self.convertCoord(tick, r0, r1, h, self.haveYInversion)
        self.majorLines.append((x0, y, x1, y))
        if isTextDrawn:
          text = self.doMajorFormat(tick)
          self.majorTexts.append((text, xt, y, at, bt))
 
  def doMajorFormat(self, tick):

    text = self.majorFormat % tick

    n = text.find('e')
    if (n >= 0):
      # strip out 0's in exponent
      n = n + 1
      m = text.rfind('0', n)
      if (m >= 0):
        if (m == len(text)-1):
          n = n - 1 # exponent is 0 so remove e as well
        text = text[:n] + text[m+1:]

    return text

  def drawTicks(self):
 
    if (self.doRuleLines):
      self.drawLine(*self.ruleLine1)
      self.drawLine(*self.ruleLine2)
 
    for line in self.minorLines:
      self.drawLine(*line)
 
    for line in self.majorLines:
      self.drawLine(*line)
 
    for data in self.majorTexts:
      self.drawText(*data)

  def findTicks(self):
 
    (r0, r1) = self.region
    s0 = float(min(r0, r1))
    s1 = float(max(r0, r1))
    d = s1 - s0
    dd = float(max(abs(r0), abs(r1)))
 
    A = 0.999
    B = 0.001
    w = B * d
 
    deltaMajor = self.deltaMajor
    if deltaMajor is None:
      n = int(floor(log(A*d)/log(10)))
      deltaMajor = pow(10.0, n)
      t = 0
      if 5*deltaMajor < d:
        deltaMajor = 5 * deltaMajor
        t = 2
      elif (2*deltaMajor < d):
        deltaMajor = 2 * deltaMajor
        t = 1
 
      nticks = self.numberTicks(s0, s1, w, deltaMajor)
      d = deltaMajor
      #print 'findTicks0', nticks, deltaMajor, t
      while 1:
        if t == 2:
          d = (2 * d) / 5
        else:
          d = d / 2
        nticks = self.numberTicks(s0, s1, w, d)
        #print 'findTicks1', nticks, deltaMajor, d, t
        if nticks > self.maxMajor:
          break
        if t == 0:
          n = n - 1
        deltaMajor = d
        t = (t - 1) % 3
    else:
      if deltaMajor > 0:
        n = log(deltaMajor) / log(10)
      else:
        n = 0

    self.majorTicks = self.calcTicks(s0, s1, w, deltaMajor)

    numberDecimals = self.numberDecimals
    if dd >= 1000:
      if numberDecimals is None:
        nn = int(floor(log(A*dd)/log(10)))
        numberDecimals = max(0, nn-n)
      self.majorFormat = '%%.%de' % numberDecimals
    else:
      if numberDecimals is None:
        numberDecimals = max(self.minDecimalPlaces, -n)
      self.majorFormat = '%%.%df' % numberDecimals
 
    #print 'findTicks2', self.majorTicks, nticks, deltaMajor, t

    deltaMinor = self.deltaMinor
    if deltaMinor is None and self.deltaMajor is None:
      d = deltaMinor = deltaMajor
      if t == 2:
        d = d / 5
      else:
        d = d / 2
      while 1:
        nticks = self.numberTicks(s0, s1, w, d)
        if nticks > self.maxMinor:
          break
        deltaMinor = d
        t = (t - 1) % 3
        if (t == 2):
          d = (2 * d) / 5
        else:
          d = d / 2

    self.minorTicks = []
    if self.deltaMinor is not None or self.deltaMajor is None:
      self.minorTicks = self.calcTicks(s0, s1, w, deltaMinor)
#      minorTicks = self.calcTicks(s0, s1, w, deltaMinor)
#      # remove minor ticks that are also major (slow code but shouldn't matter)
#      for minorTick in minorTicks:
#        for majorTick in self.majorTicks:
#          if abs(majorTick-minorTick) < 0.1*deltaMinor:
#            break
#        else:
#          self.minorTicks.append(minorTick)

  def numberTicks(self, s0, s1, w, delta):
 
    #n0 = int(math.ceil((s0+w)/delta))
    #n1 = int(math.floor((s1-w)/delta))
 
    #return (n1 - n0 + 1)
 
    # below might under-estimate by 1 but that is ok
    return int(floor((s1-s0-2*w) / delta))
 
  def calcTicks(self, s0, s1, w, delta):
 
    if delta <= 0:
      return []

    n0 = int(ceil((s0+w)/delta))
    n1 = int(floor((s1-w)/delta))
 
#    ticks = [ n * delta for n in xrange(n0, n1+1) ]
    ticks = map(lambda x: x*delta,  xrange(n0, n1+1))
 
    return ticks
 
if (__name__ == '__main__'):

  class TestTicks(Ticks):

    def drawLine(self, x0, y0, x1, y1):

      print 'drawLine', x0, y0, x1, y1

    def drawText(self, text, x, y, a, b):

      print 'drawText', text, x, y, a, b

  t = TestTicks(width = 200.0, height = 20.0, isHorizontal = True,
                region = (4.6, 12.3), doRuleLines = False, haveYInversion = False)

  print 'minorTicks', t.minorTicks
  print 'majorTicks', t.majorTicks
  print 'majorFormat', t.majorFormat
  print 'minorLines', t.minorLines
  print 'majorLines', t.majorLines
  print 'majorTexts', t.majorTexts
