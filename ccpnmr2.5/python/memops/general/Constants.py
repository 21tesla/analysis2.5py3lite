"""Constants for initialisation of global variables,
used both at meta level and data level
NB Implementation Constants (used only at metalevel) are set in ImpConstants


======================COPYRIGHT/LICENSE START==========================

Constants.py: code for CCPN data model and code generation framework

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

# NBNB TBD clear out unneeded

from memops.universal.Constants import *

from memops.general import Version


# Current version of data model.
# Used by generation scripts to mark generated code.
# Main way of tracking IO code and IO mappings for compatibiilty.
# Incremented by hand when model (or I/O generators) changes

currentModelVersion = Version.Version(2, 1, '', 2,'', 'DataModel')

# valid characters for file names 
# NB string.ascii_letters and string.digits are not compatible 
# with Python 2.1 (used in ObjectDomain)
defaultFileNameChar = '_'
separatorFileNameChar = '+'
validFileNamePartChars = ('abcdefghijklmnopqrstuvwxyz'
                          'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789' 
                          + defaultFileNameChar)
validCcpnFileNameChars  = validFileNamePartChars + '-.' + separatorFileNameChar


# Repository tag prefix
repositoryTagPrefix = 'branch_'

# NB the ''.join() construction is necessary to avoid the string being overwritten on check-in
emptyRepositoryId = ''.join(('$','I','d',':  $'))
emptyRepositoryTag = ''.join(('$','N','a','m','e',':  $'))

# Data Type constants    

# Special constants
infinity = -1

# names of basic data types
int_code = 'Int'
long_code = 'Long'
float_code = 'Float'
double_code = 'Double'
string_code = 'String'
boolean_code = 'Boolean'
datetime_code = 'DateTime'
dict_code = 'Dict'
typeCode_enumeration = [int_code, long_code, float_code, double_code,
 string_code, boolean_code, datetime_code, dict_code,
]
typeCode_enumeration.sort()

baseDataTypeModule = __import__(
 'memops.general.baseDataTypes', globals(), locals(), typeCode_enumeration
)

# NBNB TBD the codes for 'Dict and List are short-term ahcks

# names of jdbc data types.
jdbcTypeCode_enumeration = [
 'CLOB', 'VARCHAR', 'BIT', 'INTEGER', 'BIGINT', 'FLOAT', 'DOUBLE', 'TIMESTAMP',
 'Dict','List'
]
jdbcTypeCode_enumeration.sort()

# names of java data types. 
javaTypeCode_enumeration = ['java.lang.Integer','java.lang.Long',
'java.lang.Float','java.lang.Double', 'java.lang.String','java.lang.Boolean',
'java.util.Map', 'java.util.Collection'
]
javaTypeCode_enumeration.sort()

# names of java simple data types. 
# NB the simpleType for strings is the empty string
javaSimpleTypeCode_enumeration = [
 'int','long','float','double','boolean', None
]
javaSimpleTypeCode_enumeration.sort()

# names of xml data types.
xmlTypeCode_enumeration = [
 'string', 'boolean', 'int', 'long', 'float', 'double',
 'Dict','List'
]
xmlTypeCode_enumeration.sort()

# dictionary of all typeCodes and their enumeration, to simplify processing
typeCodes = {
 'typeCode':typeCode_enumeration,
 'jdbcTypeCode':jdbcTypeCode_enumeration,
 'javaTypeCode':javaTypeCode_enumeration,
 'javaSimpleTypeCode':javaSimpleTypeCode_enumeration,
 'xmlTypeCode':xmlTypeCode_enumeration
}

# file, string constants
keyStringSeparator = ','

# directory name for python directory
pythonDir = 'python'

# directory name for model directory (parallel with Python directory)
modelDir = 'model'

# directory name for java directory
javaDir = 'java'

# directory name for java sql api
javaSqlDir = 'sqlApi'

# directory name for java xml api
javaXmlDir = 'xmlApi'

# directory name for java model
javaModelDir = 'model'

# directory name for sql table creation scripts
sqlDir = 'sql'

# Directory names for Python code model subdirectories
modelSubDirs = ['model']

# Directory names for Python code api subdirectories
apiCodeDir = 'api'
apiSubDirs = [apiCodeDir]

# Directory names for Python code xml subdirectories
xmlCodeDir = 'xml'
xmlSubDirs = [xmlCodeDir]

# Directory names for Python code xml Schema subdirectories (in model)
xmlSchemaDir = 'xmlSchema'
xmlSchemaSubDirs = [xmlSchemaDir]
