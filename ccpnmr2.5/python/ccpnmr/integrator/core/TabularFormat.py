
"""
======================COPYRIGHT/LICENSE START==========================

Io.py: code for CCPN data model and code generation framework

Copyright (C) 2011  (CCPN Project)

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

class TabularFormat(object):
  """generic table and keyword file format
  """
  
  # format settings.
  # NB should be called as self.xyz, so they can be overridden
  format = 'Missing'
  
  lineSep = '\n'
  fieldSep = ' '
  dummyValue = '.'
  dummiesAllowed=True
  
  startCommentLine = '# '
  startDataLine = ''
  startComment = ' # '   # End-of-line comment. If None, skip such comments.
  
  def __init__(self, **kw):
  
    self.lines = []
    self.numTableColumns = None
    self.columnFormats = None
    
    # Set keywords - used to set and override format settings
    for key,val in kw.items():
      setattr(self, key, val)
    
    self.addComment("Written by CcpNmr Integrator. Format is %s" % self.format)
    self.addLine()
  
  def getText(self):
    """ get text as block
    """
    return lineSep.join(self.lines)
  
  def export(self, target):
    """ Write text to file
    Input: target can be a string or a stream
    """
    if isinstance(target, basestring):
      stream = open(target,'w')
    else:
      stream = target
    for line in self.lines:
      stream.write(line)
      stream.write('\n')
  
  def clear(self):
    """ Clear data and make ready for new use
    """
    del self.lines[:]
    self.numTableColumns = None
    self.columnFormats = None
  
  def addLine(self):
    """Add empty line
    """
    
    self.lines.append('')
      
  def addData(self, keyword, value, extraKeywords=None, comment=None, 
              multipleValue=False):
    """ Add data line
    """
    if self.startDataLine:
      keyword = self.startDataLine + keyword
    ll = [keyword]
    
    if extraKeywords:
      ll.extend(extraKeywords)
      
    if multipleValue:
      if None in value:
        if self.dummiesAllowed:
          for ii,val in enumerate(value):
            if val is None:
              val[ii] = self.dummyValue
      else:
        raise ValueError("value 'None' not allowed")
      ll.extend(value)
      
    else:
      if value is None:
        if self.dummiesAllowed:
          value = self.dummyValue
      else:
        raise ValueError("value 'None' not allowed")
      ll.append(value)
    
    text = self.fieldSep.join(ll)
    
    if comment:
      if self.startComment is not None:
        text = self.startComment.join((text, comment))
    #
    self.lines.append(text)
      

  def addComment(self, comment):
    """ Add single line comment
    """
    self.lines.append(self.startCommentLine + comment)
  
  
  def addMultiComment(self,comment):
    """Add multi-line comment
    NB should be overwritten in formats with proper multiline comments
    """
    for ss in comment.splitlines():
      self.addComment(ss)
  
  def addSequence(self, sequence):
    """ Add sequence record 
    Input: a list of MolResidues or Residues in correct order
    """
    raise NotImplementedError("addSequence not implemented in format %s" 
                              % self.format)
  
  def _initTable(self, colNames, formats=None):
    """ Initialize table start record 
    Input:
    colNames: Column names
    formats (optional) Format expressions, in C/Python '%' syntax 
    
    Must set self.numTableColumns and self.columnFormats
    """
    self.numTableColumns = len(colNames)
    self.columnFormats = formats
  
  def startTable(self, colNames, formats=None):
    """ Add table start record 
    Input:
    colNames: Column names
    formats (optional) Format expressions, in C/Python '%' syntax 
    
    Must set self.numTableColumns and self.columnFormats
    """
    self._initTable(colNames, formats)
    self.addLine()
    self.lines.append(self.fieldSep.join(colNames))
    self.addLine()
    
  
  def addTableLine(self, values):
    """ Add line to existing table
    """
    
    formats = self.columnFormats
    
    if self.numTableColumns is None:
      raise exception("No table currenlty open")
    
    if None in values and not self.dummiesAllowed:
      raise ValueError("value 'None' not allowed")
    
    if len(values) != self.numTableColumns:
      raise ValueError("%s values passed to %s-column table" 
                       % (len(values),self.numTableColumns))
      
    ll = []
    for ii,val in enumerate(values):
      if val is None:
        ss = self.dummyValue
      elif formats:
        ss = formats[ii] % val
      else:
        ss = repr(val)
      ll.append(ss)
    
    self.lines.append(self.fieldSep.join(ll))
    
  
  def _closeTable(self):
    """ Close table
    """
    self.numTableColumns = None
    self.columnFormats = None
  
  def endTable(self):
    """ Add table end record 
    """
    self._closeTable()
    self.addLine()
   
