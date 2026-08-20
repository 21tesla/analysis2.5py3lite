"""
======================COPYRIGHT/LICENSE START==========================

formats.py: license header generation code for CCPN framework

Copyright (C) 2004  (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/GPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
 
You should have received a copy of the GNU General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

======================COPYRIGHT/LICENSE START==========================

licenses.py: license header generation code for CCPN framework

Copyright (C) 2004  (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/GPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
 
You should have received a copy of the GNU General Public
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

===========================REFERENCE END===============================


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

===========================REFERENCE END===============================

"""
import re

# information for differnt file type formats

formats = {
 'python' : {
  'commentStart': '\"\"\"',
  'commentEnd': '\"\"\"',
  'includeFileMatches': [r'.*\.py$'],
  'ignoreStartLines': ['^#!/']
 },
          
 'xml' : {
  'commentStart': '<!--',
  'commentEnd': '-->',
  'includeFileMatches': [r'.*\.xml$'],
  'ignoreStartLines': ['^\<\?xml']
 },      
 'html' : {
  'commentStart': '<!--',
  'commentEnd': '-->',
  'includeFileMatches': [r'.*\.html$', r'.*\.htm$'],
  'ignoreStartLines': ['^\<\!DOCTYPE', '^\<html']
 },
 
 'css' : {
  'commentStart': '/*',
  'commentEnd': '*/',
  'includeFileMatches': [r'.*\.css'],
  'ignoreStartLines': []
 },

 'c' : {
  'commentStart': '/*',
  'commentEnd': '*/',
  'includeFileMatches': [r'.*\.c$', r'.*\.h$'],
  'ignoreStartLines': []
 },
 
 'java' : {
  'commentStart': '/*',
  'commentEnd': '*/',
  'includeFileMatches': [r'.*\.java$'],
  'ignoreStartLines': []
 },
 
 'text' : {
  'commentStart': '<!--',
  'commentEnd': '-->',
  'includeFileMatches': [],
  'ignoreStartLines': []
 }
}


# preprocess match expressions

emptyLineExpr = re.compile("^\s*$")

for dd in formats.values():

  ll = dd['includeFileMatches']
  for ii in range(len(ll)):
    ll[ii] = re.compile(ll[ii])

  ll = dd['ignoreStartLines']
  for ii in range(len(ll)):
    ll[ii] = re.compile(ll[ii])
    
  dd['ignoreStartLines'].append(emptyLineExpr)
