
"""
======================COPYRIGHT/LICENSE START==========================

Color.py: <write function here>

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
from math import sqrt

NAMES = []
colors = {}


class Color:

  def __init__(self, name, hexString):

    assert name not in NAMES, 'name = %s' % name

    self.name = name
    self.hex  = hexString.upper()

    NAMES.append(name)
    colors[name] = self

  def rgb(self):

    return hexToRgb(self.hex)

  def hsb(self):
  
    return hexToHsb(self.hex)
    
white =      Color('white',     "#FFFFFF")
midgrey =    Color('midgrey',   "#C0C0C0")
grey =       Color('grey',      "#808080")
darkgrey =   Color('darkgrey',  "#404040")
black =      Color('black',     "#000000")
maroon =     Color('maroon',    "#800000")
red =        Color('red',       "#FF0000")
midpink =    Color('midpink',   "#FF8080")
lightpink =  Color('lightpink', "#FFC0C0")
orange =     Color('orange',    "#FF8000")
yellow =     Color('yellow',    "#FFFF00")
chartreuse = Color('chartreuse',"#80FF00")
lime =       Color('lime',      "#00FF00")
midgreen =   Color('midgreen',  "#00C000")
green =      Color('green',     "#008000")
lightgreen = Color('lightgreen',"#80FF80")
spring =     Color('spring',    "#00FF80")
cyan =       Color('cyan',      "#00FFFF")
teal =       Color('teal',      "#008080")
skyblue =    Color('skyblue',   "#0080FF")
lightblue =  Color('lightblue', "#C0C0FF")
blue =       Color('blue',      "#0000FF")
midblue =    Color('midblue',   "#8080FF")
navy =       Color('navy',      "#000080")
mauve =      Color('mauve',     "#8000FF")
magenta =    Color('magenta',   "#FF00FF")
purple =     Color('purple',    "#800080")
pink =       Color('pink',      "#FF0080")

standardColors = (white,midgrey,grey,darkgrey,black,
                  maroon,red,midpink,lightpink,
                  orange,yellow,chartreuse,lime,
                  midgreen,green,lightgreen,spring,
                  cyan,teal,skyblue,lightblue,
                  blue,midblue,navy,mauve,
                  magenta,purple,pink)
    
# r, g, b expected to be between 0 and 1
def hexRepr(r, g, b):

  r = max(0, min(int(255 * r), 255))
  g = max(0, min(int(255 * g), 255))
  b = max(0, min(int(255 * b), 255))

  hexRep = '#%02x%02x%02x' % (r, g, b)

  return hexRep.upper()

def hexToHsb(hexColor):
  
  rh = hexColor[1:3]
  gh = hexColor[3:5]
  bh = hexColor[5:7]

  r = int(rh, 16)/255.0
  g = int(gh, 16)/255.0
  b = int(bh, 16)/255.0

  return rgbToHsb(r, g, b)

def hexToRgb(hexColor):
  
  rh = hexColor[1:3]
  gh = hexColor[3:5]
  bh = hexColor[5:7]

  r = int(rh, 16)/255.0
  g = int(gh, 16)/255.0
  b = int(bh, 16)/255.0

  return r, g, b

# r, g, b expected to be between 0 and 1
def rgbToHsb(r, g, b):

    Min = min( r, g, b )
    Max = max( r, g, b )
    v = Max
    delta = float(Max - Min)
    if Max != 0:
      s = delta / Max 
    else:
      s = 0.0
      h = 0.0 # arbitrary
      return (h, s, v)
    
    if s == 0.:
      h = 0 # arbitrary
      return (h, s, v)
 
    if   r == Max:
      h = ( g - b ) / delta       
    elif g == Max:
      h = 2 + ( b - r ) / delta   
    else:
      h = 4 + ( r - g ) / delta   
    
    h /= 6.0                       
    if h < 0:
      h += 1.0
      
    return (h, s, v)

# H, S, V expected to be between 0 and 1
def hsbToRgb(H, S, V):

    H = H * 6
    i = int(H)
    f = H - i
    p = V*(1-S)
    q = V*(1-(S*f))
    t = V*(1 - (S * (1-f)))

    if i == 0:
      R = V
      G = t
      B = p
    elif i == 1:
      R = q
      G = V
      B = p
    elif i == 2:
      R = p
      G = V
      B = t
    elif i == 3:
      R = p
      G = q
      B = V
    elif i == 4:
      R = t
      G = p
      B = V
    else:
      #elif i = 5:
      R = V
      G = p
      B = q
   
    return (R,G,B)

def hexXor(hexColor1, hexColor2):
  
  r1 = int('0x'+hexColor1[1:3], 16)
  g1 = int('0x'+hexColor1[3:5], 16)
  b1 = int('0x'+hexColor1[5:7], 16)
  r2 = int('0x'+hexColor2[1:3], 16)
  g2 = int('0x'+hexColor2[3:5], 16)
  b2 = int('0x'+hexColor2[5:7], 16)

  return '#%02x%02x%02x' % (r1^r2, g1^g2, b1^b2)

def hexInvert(hexColor):
  
  r = int('0x'+hexColor[1:3], 16)
  g = int('0x'+hexColor[3:5], 16)
  b = int('0x'+hexColor[5:7], 16)

  return '#%02x%02x%02x' % (255-r, 255-g, 255-b)

# returns value between 0 and 255
def getIntRgb(widget, color):

  rgb = widget.winfo_rgb(color)
  (r, g, b) = (rgb[0]/256, rgb[1]/256, rgb[2]/256)
  r = min(255, r)
  g = min(255, g)
  b = min(255, b)

  return (r, g, b)

def getFloatRgb(widget, color):

  rgb = widget.winfo_rgb(color)
  n = float(256 * 256)
  (r, g, b) = (rgb[0]/n, rgb[1]/n, rgb[2]/n)

  return (r, g, b)

def invertColor(widget, color):

  [r0,g0,b0] = getFloatRgb(widget, color)
  m = (0.35*r0)+(0.50*g0)+(0.15*b0)
  if   m<0.25:
    m = 1-m
  elif m<0.50:
    m = m+0.5
  elif m<0.75:
    m = m-0.5
  else:
    m = 1-m
  
  r = m
  g = m
  b = m

  return '#%02x%02x%02x' % (r*255, g*255, b*255)

def inverseGrey(color):

  r,g,b = hexToRgb(color)
  
  m = (0.35*r)+(0.50*g)+(0.15*b)
  if   m<0.25:
    m = 1-m
  elif m<0.50:
    m = m+0.5
  elif m<0.75:
    m = m-0.5
  else:
    m = 1-m
  
  m *=255
  
  return '#%02x%02x%02x' % (m,m,m)
  

def inverseRgb(color):

  r,g,b = hexToRgb(color)
  
  r = 255*(1-r)
  g = 255*(1-g)
  b = 255*(1-b)
  
  return '#%02x%02x%02x' % (r,g,b)
  
  
def invertColorRgb(widget, color):

  rgb = widget.winfo_rgb(color)
  (r, g, b) = (255-(rgb[0]/256.),255-(rgb[1]/256.),255-(rgb[2]/256.))

  return '#%02x%02x%02x' % (r, g, b)

def scaleColor(widget, color, scale):

  assert scale >= 0, 'scale = %s' % scale
  rgb = widget.winfo_rgb(color)
  (r, g, b) = (rgb[0]/256, rgb[1]/256, rgb[2]/256)

  r = min(255, int(round(scale * r)))
  g = min(255, int(round(scale * g)))
  b = min(255, int(round(scale * b)))

  return '#%02x%02x%02x' % (r, g, b)

if __name__ == '__main__':

  c = Color('red2', '#FF0000')

  print c.hex, c.rgb()

  c = Color('my_color', '#5419B8')

  print c.hex, c.rgb()

  for c in standardColors:
    print c.name, c.hex, c.rgb(), c.hsb()
