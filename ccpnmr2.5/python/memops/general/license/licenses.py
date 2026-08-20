gnuLicense = """
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU %sGeneral Public
License as published by the Free Software Foundation; either
version 2%s of the License, or (at your option) any later version.
 
A copy of this license can be found in %%(licenseLocation)s/%s.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%sGeneral Public License for more details.
 
You should have received a copy of the GNU %sGeneral Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

""" 

licenses = {

'GPL'  : gnuLicense % ('','','GPL','',''),

'LGPL' : gnuLicense % ('Lesser ','.1','LGPL','Lesser ','Lesser '),

 'restricted' : """
This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of %(organization)s
and in accordance with the guidelines of %(organization)s.
""",

 'ccpn' :"""
This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in %(licenseLocation)s/CCPN.license.
""",

 'brief' : 
 """The license for this file is in %(licenseLocation)s/%(licenseFileName)s.license"""
}
