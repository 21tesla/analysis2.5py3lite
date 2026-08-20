"""Basic data type Boolean true/false, implemented as (0,1) Int enum

Correct procedure is that values of type in compatibleTypes must be 
accepted as input, and coerced to PythonType. isValid should be applied only
after coercion. Note that the coercion may well cause an error
(e.g. out-of-range)


======================COPYRIGHT/LICENSE START==========================

Boolean.py: Basic data type implementation for CCPN Python data model API

Copyright (C) 2005 Rasmus Fogh (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
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
Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================
"""

# corresponding python type
from types import IntType as PythonType

from memops.general.Constants import True, False
from memops.general.Constants import trueString, falseString

# Python types acceptable as input instead of main type
compatibleTypes = ()

# Special isValid function
# None

try:
  # conversion to String
  toString  = repr

  # conversion from string
  fromString = bool

  # casting/creation function.
  create = bool

except:
  # for old python versions, like Jython
  def toString(value):
    return value and trueString or falseString
  
  def fromString(value):
    if value == trueString:
      return True
    elif value == falseString:
      return False
    else:
      raise ValueError("invalid literal for fromString(): %s" % value)
  
  def create(value):
    return not not value
