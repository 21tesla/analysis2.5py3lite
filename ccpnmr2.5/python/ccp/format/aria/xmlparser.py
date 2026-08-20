"""
Authors: Wolfgang Rieping
         Michael Habeck
         Institut Pasteur, Paris
        
         Copyright (C) 2003 Michael Habeck and Wolfgang Rieping
         No warranty implied or expressed.
         All rights reserved.

$Author: wimvranken $
$Revision: 1.3.12.1 $
$Date: 2010-03-05 11:34:30 $

======================COPYRIGHT/LICENSE START==========================

xmlparser.py: 

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

- contact Wim Vranken (wim@ebi.ac.uk)
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

import re, string

def quote(pattern):

    if type(pattern).__name__ == 'SRE_Pattern':
        return '(\'%s\'|"%s")' % (pattern.pattern, pattern.pattern)
    else:
        return '(\'%s\'|"%s")' % (pattern, pattern)

CHAR            = re.compile('(.|\n)')
CHAR_DATA       = re.compile('[^<&]*')
NOT_CHAR_DATA   = re.compile('[^<&]*\]\]>[^<&]*')
LETTER          = re.compile('[a-zA-Z]')
DIGIT           = re.compile('[0-9]')
SPACE           = re.compile('\s+')
OPTIONAL_SPACE  = re.compile('\s*')
ATTRIBUTE_VALUE = re.compile('("[^"]*"|\'[^\']*\')')
COMMENT         = re.compile('<!--([^-]|-[^-])*-->')
ENCODING_NAME   = re.compile('[A-Za-z]([A-Za-z0-9\._]|\-)*')
SYSTEM_LITERAL  = re.compile('("[^"]*")|(\'[^\']*\')')
VERSION_NUMBER  = re.compile('([a-zA-Z0-9_.:]|\-)+')
ILLEGALCHARS    = re.compile('[^\t\r\n -\176\240-\377]') 
EQUALS          = re.compile(OPTIONAL_SPACE.pattern + '=' + \
                             OPTIONAL_SPACE.pattern)

VERSION_INFO = re.compile(SPACE.pattern + 'version' + \
                          EQUALS.pattern + quote(VERSION_NUMBER))

ENCODING_DECLARATION = re.compile(SPACE.pattern + 'encoding' + \
                                  EQUALS.pattern + quote(ENCODING_NAME))

STANDALONE_DECLARATION = re.compile(SPACE.pattern + 'standalone' + \
                                    EQUALS.pattern + quote('(yes|no)'))

NAMECHAR = re.compile(LETTER.pattern + '|' + DIGIT.pattern + '|\.|\-|_|:')

NAME = re.compile('(' + LETTER.pattern   + '|_|:)' + \
                  '(' + NAMECHAR.pattern + ')*')

ATTRIBUTE = re.compile(NAME.pattern + EQUALS.pattern + ATTRIBUTE_VALUE.pattern)

START_TAG = re.compile('<' + NAME.pattern + \
                       '(' + SPACE.pattern + ATTRIBUTE.pattern + ')*' + \
                       OPTIONAL_SPACE.pattern + '>') 

END_TAG = re.compile('</' + NAME.pattern + OPTIONAL_SPACE.pattern + '>')

EMPTY_ELEMENT = re.compile(START_TAG.pattern[:-1] + '/>')

XML_DECLARATION = re.compile('<\?xml' + VERSION_INFO.pattern + \
                             '(' + ENCODING_DECLARATION.pattern + ')?'+ \
                             '(' + STANDALONE_DECLARATION.pattern + ')?' + \
                             OPTIONAL_SPACE.pattern + '\?>')

EXTERNAL_ID = re.compile('(SYSTEM' + SPACE.pattern + \
                         SYSTEM_LITERAL.pattern + '|' + \
                         'PUBLIC' + SPACE.pattern + \
                         SPACE.pattern + \
                         SYSTEM_LITERAL.pattern + ')')

MARKUP_DECLARATION = re.compile(CHAR.pattern)
DECLARATION_SEPARATOR = re.compile(SPACE.pattern)
                         
DOCTYPE_DECLARATION = re.compile('<!DOCTYPE' + SPACE.pattern + NAME.pattern +\
                                 '(' + SPACE.pattern + EXTERNAL_ID.pattern + \
                                 ')?' + OPTIONAL_SPACE.pattern + \
                                 '(\[' + '(' + MARKUP_DECLARATION.pattern + \
                                 '|' + DECLARATION_SEPARATOR.pattern + ')*'+\
                                 '\]' + OPTIONAL_SPACE.pattern + ')?>')

MISCELLANEOUS = re.compile('(' + COMMENT.pattern + \
                           '|' + SPACE.pattern + ')')

PROLOG = re.compile('(' + XML_DECLARATION.pattern + ')?' + \
                    '(' + MISCELLANEOUS.pattern + ')*' + \
                    '(' + DOCTYPE_DECLARATION.pattern + \
                    MISCELLANEOUS.pattern + '*)?')

class SelfmadeXMLParser:

    is_self_made = 1

    def __init__(self):

        self.setup()

    def setup(self):

        opt_space = OPTIONAL_SPACE.pattern
        
        self.junk = re.compile('%s(%s)%s' % (opt_space, COMMENT.pattern,
                                             opt_space))
        
        name = '(?P<name>%s)' % NAME.pattern
        attr = '(?P<attr>(%s%s)*)' % (SPACE.pattern, ATTRIBUTE.pattern)
        
        self.start = re.compile(opt_space + '<' + name + attr + opt_space + \
                                '(?P<end>/?)>' + opt_space) 
        self.end   = re.compile('%s</%s%s>%s' % (opt_space, name, opt_space,
                                                 opt_space))

        self.attr = re.compile(SPACE.pattern + \
                               '(?P<attr>%s)' % NAME.pattern + \
                               EQUALS.pattern + \
                               '(?P<value>%s)' % ATTRIBUTE_VALUE.pattern)

    def parse_prolog(self, xmlfile):

        match = PROLOG.match(xmlfile)
            
        if match is None:
            return 0
        else:
            return match.end()

    def handle_prolog(self, prolog):

        match = DOCTYPE_DECLARATION.search(prolog)

    def parse_junk(self, xmlfile):
        return self.junk.match(xmlfile)

    def go(self, xmlfile):

        ## parse for junk
        
        match = self.parse_junk(xmlfile)
        if match is not None:
            return xmlfile[match.end():]
            
        ## parse for start

        match = self.start.match(xmlfile)
        if match is not None:
            tag, attrs = self.handle_element_match(match)
            self.StartElementHandler(tag, attrs)

            ## empty element
            if match.group('end'):
                self.EndElementHandler(tag)                

            return xmlfile[match.end():]

        ## parse for end

        match = self.end.match(xmlfile)
        if match is not None:
            tag = self.handle_end_match(match)
            self.EndElementHandler(tag)
            return xmlfile[match.end():]

        raise StandardError, 'mis-formatted xml file'

    def handle_element_match(self, match):

        name = match.group('name')
        attributes = match.group('attr')

        ## extract attributes and values

        attr_dict = {}

        match = self.attr.match(attributes)

        while match:

            ## get rid of the quotation marks therefore 'eval'
            attr_dict[match.group('attr')] = eval(match.group('value'))
            attributes = attributes[match.end():]
            match = self.attr.match(attributes)

        return name, attr_dict

    def handle_end_match(self, match):
        
        return match.group('name')

    def StartElementHandler(self, tag, attrs):
        raise NotImplementedError

    def EndElementHandler(self, tag):
        raise NotImplementedError
            
    def ParseFile(self, xmlfile):

        file = xmlfile.read()

        ## remove prolog

        end = self.parse_prolog(file)
        self.handle_prolog(file[:end])
        file = file[end:]

        ## parse loop

        while file:

            file = self.go(file)

if __name__ == '__main__':

    p = SelfmadeXMLParser()

    def start_handler(x, y):
        return
        print x, y

    def end_handler(x):
        return
        print x


    p.StartElementHandler = start_handler
    p.EndElementHandler = end_handler

