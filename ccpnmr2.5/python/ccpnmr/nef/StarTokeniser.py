# encoding: utf-8

"""

STAR file tokenizer

# Copyright © 2011, 2013 Global Phasing Ltd. All rights reserved.
#
# Author: Peter Keller
#
# This file forms part of the GPhL StarTools library.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#  Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the
#  distribution.
#
#  If the regular expression used to match STAR/CIF data in the
#  redistribution is not identical to that in the original version,
#  this fact must be stated wherever the copyright notice is
#  reproduced.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.

'''
Created on 25 Nov 2013

@author: pkeller
'''

#
# Modified by Rasmus Fogh, CCPN project, 5/2/2016
#

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2019"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: CCPN $"
__dateModified__ = "$dateModified: 2017-07-07 16:33:02 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import re
import collections


# STAR parsing REGEX, following International Tables for Crystallography volume G section 2.1
_REGEX = r"""(?xmi) # $Revision$  # No 'u' flag for perl 5.8.8/RHEL5 compatibility
^;([\S\s]*?)(?:\r\n|\s)^;(?:(?=\s)|$)  # Multi-line string
|(?:^|(?<=\s))(\#.*?)\r?$              # Comment
|(?:^|(?<=\s))(?:
  (global_)                            # STAR global block
  |(save_\S*)                          # STAR save frame header or terminator
  |(\$\S+)                             # STAR save frame reference
  |(stop_)                             # STAR nested loop terminator
  |(data_\S+)                          # Data block header
  |(loop_)                             # Loop header
  |((?:global_\S+)|(?:stop_\S+)|(?:data_)|(?:loop_\S+))  # Invalid privileged construct
  |(_\S+)                              # Data name
  |'(.*?)'                             # Single-quoted string
  |"(.*?)"                             # Double-quoted string
  |(\.)                                # CIF null
  |(\?)                                # CIF unknown/missing
  |([\[\]]\S*)                         # Square bracketed constructs (reserved)
  |((?:[^'";_$\s]|(?<!^);)\S*)         # Non-quoted string
  |(\S+)                               # Catch-all bad token
)
(?:(?=\s)|$)"""

# Compiled form of _REGEX
_star_pattern = re.compile(_REGEX, re.UNICODE)

# Token types. NB numbers must be synced to regex - these are used directly!!!
TOKEN_MULTILINE         = 1
TOKEN_COMMENT           = 2
TOKEN_GLOBAL            = 3
TOKEN_SAVE_FRAME        = 4
TOKEN_SAVE_FRAME_REF    = 5
TOKEN_LOOP_STOP         = 6
TOKEN_DATA_BLOCK        = 7
TOKEN_LOOP              = 8
TOKEN_BAD_CONSTRUCT     = 9
TOKEN_DATA_NAME         = 10
TOKEN_SQUOTE_STRING     = 11
TOKEN_DQUOTE_STRING     = 12
TOKEN_NULL              = 13
TOKEN_UNKNOWN           = 14
TOKEN_SQUARE_BRACKET    = 15
TOKEN_STRING            = 16
TOKEN_BAD_TOKEN         = 17

# Rasmus Fogh, CCPN project 5/2/2016
# # Modified Tokeniser to
# - use namedtuples instead of custom objects
# - to use a string input instead of a memory map (which gave string/byte conflict errors)
# - to wrap the regex iterator without a wrapping class.

#
StarToken = collections.namedtuple('StarToken', ('type', 'value'))


# Removed for Python 2 compatibility:
# StarToken.__doc__ = "StarToken named tuple (with fields 'type', 'value')"
# "returned by the STAR token iterator"

def getTokenIterator(text):
    """Iterator that returns an iterator over all STAR tokens in a generic STAR file"""
    return (StarToken(x.lastindex, x.group(x.lastindex))
            for x in _star_pattern.finditer(text))
