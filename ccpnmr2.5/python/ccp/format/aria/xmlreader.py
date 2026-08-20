"""
Authors: Wolfgang Rieping
         Michael Habeck
         Institut Pasteur, Paris
        
         Copyright (C) 2003 Michael Habeck and Wolfgang Rieping
         No warranty implied or expressed.
         All rights reserved.

$Author: wimvranken $
$Revision: 1.4.12.2 $
$Date: 2011-03-15 14:45:14 $

======================COPYRIGHT/LICENSE START==========================

xmlreader.py: 

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

class XMLReaderError(Exception):
    pass

class XMLTagError(XMLReaderError):

    def __str__(self):
        return 'XML Element/Attribute "%s" not known or missing.' % \
               Exception.__str__(self)

class XMLDocument:

    def __init__(self, doc_type, dtd):
        self.__doc_type = doc_type
        self.__dtd = dtd

    def get_doc_type(self):
        return self.__doc_type

    def get_dtd(self):
        return self.__dtd

class XMLElement:

    def __init__(self, name = None, attr = None, tag_order = None):

        self.__cdata = None

        if name is not None:
            self.__name = name

        if attr is not None:
            for name, value in attr.items():
                setattr(self, name, value)

    def set_cdata(self, data):
        if self.__cdata is None:
            self.__cdata = data
        else:
            self.__cdata += data

    def get_cdata(self):
        return self.__cdata

    def get_name(self):
        return self.__name

    def __getattr__(self, name):
        raise XMLTagError(name)

class ContentHandler:

    """
    Abstract base class for content handler. Event-driven parsers
    must subclass from that one.
    """
    
    def startElementHandler(self, name, attr):
        raise NotImplementedError

    def endElementHandler(self, name):
        raise NotImplementedError

    def charDataHandler(self, data):
        raise NotImplementedError

class BaseReader:

    """
    Abstract reader class
    """

    def __init__(self, handler):
        self.__handler = handler

    def createParser(self):
        raise NotImplementedError

    def getContentHandler(self):
        return self.__handler

class XMLContentHandler(ContentHandler):

    def __init__(self):
        self.elements = []
        self.current = None

    def setDocument(self, d):
        """
        sets the main document of the xml file, that is the
        root of the hierarchy of XMLElement instances.
        """
        self.elements = [d]
        self.current = d

    def startElementHandler(self, name, attr):

        try:
            name = str(name).split(':')[1]
        except:
            pass

        for key, value in attr.items():
          attr[key] = str(value)

        new_element = XMLElement(name, attr)
        self.elements.append(self.current)
        self.current = new_element

    def endElementHandler(self, name):

        try:
            name = str(name).split(':')[1]
        except:
            pass

        if not len(self.elements):
            return

        e = self.load_from_element(name, self.current)
	
        self.current = self.elements.pop()
	
        try:
            value = getattr(self.current, name)
            try:
                value.append(e)
            except:
                setattr(self.current, name, [value, e])

        except:

            ## if single elements shall not be python
            ## lists, change the following line into:
            ## setattr(self.current, name, e)
            ## is setattr(self.current, name, [e]) for lists...
            
            setattr(self.current, name, e)

    def charDataHandler(self, data):

        data = data.strip()
        
        if not len(data):
            return

        self.current.set_cdata(str(data))

    def load_from_element(self, name, e):
        return e

class XMLReader(BaseReader):
    
    """
    is able to read xml files and converts it into
    a hierarchy of XMLElements.
    """

    def __init__(self, content_handler):
        BaseReader.__init__(self, content_handler)
        
    def createParser(self):
        """
        creates the parser
        """

        content_handler = self.getContentHandler()
        
        if content_handler is None:
            raise 'No content handler set.'

        try:
            import xml.parsers.expat as expat
            parser = expat.ParserCreate()
            
        except:
            from xmlparser import SelfmadeXMLParser as Parser
            parser = Parser()

        parser.StartElementHandler = content_handler.startElementHandler
        parser.EndElementHandler = content_handler.endElementHandler
        parser.CharacterDataHandler = content_handler.charDataHandler

        return parser

    def parse_doc_type(self, file):

        seek = file.tell()

        doctype_found = 0

        seek_previous = seek
        
        while not doctype_found:

            line = file.readline()
            doctype_found = line.find('<!DOCTYPE') > -1

            if seek_previous == file.tell(): break

            seek_previous = file.tell()
            
        file.seek(seek)

        if not doctype_found:
            tag, dtd = None, None

        else:
            tag = line.split()[1]
            dtd = line.split()[3][1:-2]

        try:
            tag = str(tag).split(':')[1]
        except:
            pass

        return tag, dtd

    def load(self, filename, gzip = 0):
        import os
        
        filename = os.path.expanduser(filename)

        if gzip:
            from tools import gzip_open as open_func
        else:
            open_func = open

        f = open_func(filename)

        document = XMLDocument(*self.parse_doc_type(f))

        handler = self.getContentHandler()
        handler.setDocument(document)

        parser = self.createParser()
        parser.ParseFile(f)
        f.close()

        failed = 0

        if len(handler.elements) <> 1:
            failed = 1
        elif handler.elements[0] <> document:
            failed = 1

        if failed:
            raise XMLReaderError, 'XML document misformatted.'

        return document
