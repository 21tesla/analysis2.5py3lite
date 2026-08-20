#!/usr/bin/python

"""This module provides entry, saveframe, and loop objects. Use python's built in help function for documentation.

There are two variables you can set to control our behavior. Setting bmrb.verbose to True will print some of what is going on to the terminal. Setting raise_parse_warnings to true will raise an exception if the parser encounters something problematic. Normally warnings are suppressed.

Some errors will be detected and exceptions raised, but this does not implement a full validator (at least at present).

Call directly (rather than importing) to run a self-test.
"""

#############################################
#                 Imports                   #
#############################################

import os
import re
import sys
import csv
import bmrb
import copy
import gzip
import shutil
import urllib2
import itertools
from cStringIO import StringIO

# Import our clone of an ordered dict if we are using a low version of python
if sys.version_info < (2,7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from sans import STARLexer
from sans import SansParser
from sans import ErrorHandler, ContentHandler
from sans import quote

#############################################
#            Global Variables               #
#############################################

verbose = False
raise_parse_warnings = False

# Rasmus Fogh October 2014:
skip_empty_loops = False

#############################################
#             Module methods                #
#############################################

def diff(entry1,entry2):
    """Prints the differences between two entries. Non-equal entries will always be detected, but specific differences detected depends on order of entries."""
    diffs = entry1.compare(entry2)
    if len(diffs) == 0:
        print "Identical entries."
    for difference in diffs:
        print difference

def validate(entry,schema=None):
    """Prints a validation report of an entry."""
    validation = entry.validate(schema)
    if len(validation) == 0:
        print "No problems found during validation."
    for err in validation:
        print err

def __cleanValue__(value):
    """Automatically quotes the value in the appropriate way. Don't quote values you send to this method or they will show up in another set of quotes as part of the actual data. E.g.:

    __cleanValue__('"e. coli"') returns '\'"e. coli"\''

    while

    __cleanValue__("e. coli") returns "'e. coli'"

    In fact, you probably will never have to call this method directly. (All print calls automatically use it.)
    """

    # NBNB None/True/False handling added by Rasmus Fogh 20141007
    if value is None:
      value = '.'
    elif value is True:
      value = 'true'
    elif value is False:
      value = 'false'
    else:
      # stings pass through this unchanged
        value = str(value)

    # This was the original code:
    # if type(value) != str:
    #   value = str(value)

    # If it's going on it's own line, don't touch it
    if "\n" in value:
        if value[-1] != "\n":
            return value + "\n"
        return value

    if value == "":
        raise ValueError("Empty strings are not allowed as values. Use a '.' or a '?' if needed.")

    # Put tags with newlines on their own ; delimited line
    if " " in value or "\t" in value or value[0] == "_" or "#" in value \
     or (len(value) > 4 and (value[:5] == "data_" or value[:5] == "save_" or value[:5] == "loop_" or value[:5] == "stop_")):
        if '"' in value:
            return "'" + value + "'"
        elif "'" in value:
            return '"' + value + '"'
        elif '"' in value and "'" in value:
            return value + '\n'
        else:
            return "'" + value + "'"

    # It's good to go
    return value

def __formatCategory__(value):
    """Adds a '_' to the front of a tag (if not present) and strips out anything after a '.'"""
    if value:
        if value[:1] != "_":
            value = "_" + value
        if "." in value:
            value = value[:value.index(".")]
    return value

def __formatTag__(value):
    """Strips anything before the '.'"""
    if '.' in value:
        value = value[value.index('.')+1:]
    return value

def __getSchema__(passed_schema=None):
    """If passed a schema (not None) it returns it. If passed none, it checks if the default schema has been initialized. If not initialzed, it initializes it. Then it returns the default schema."""

    global standard_schema
    if passed_schema is None:
        passed_schema = standard_schema
    if passed_schema is None:
        standard_schema = schema()
        passed_schema = standard_schema

    return passed_schema

def __interpretFile__(the_file):
    """Helper method returns some sort of object with a read() method. the_file could be a URL, a file location, a file object, or a gzipped version of any of the above."""

    if hasattr(the_file, 'read') and hasattr(the_file,'readline'):
        star_buffer = the_file
    elif type(the_file) is str:
        if the_file[0:7] == "http://" or the_file[0:8] == "https://" or the_file[0:6] == "ftp://":
            star_buffer = urllib2.urlopen(the_file)
            return star_buffer
        else:
            star_buffer = open(the_file, 'r')
    else:
        raise ValueError("Cannot figure out how to interpret the file you passed.")

    # Decompress the buffer if we are looking at a gzipped file
    try:
        gzip_buffer = gzip.GzipFile(fileobj=star_buffer)
        gzip_buffer.readline()
        gzip_buffer.seek(0)
        return gzip_buffer
    # Apparently we are not looking at a gzipped file
    except (IOError,AttributeError) as e:
        star_buffer.seek(0)
        return star_buffer

#############################################
#                Classes                    #
#############################################

class __entryParser__(ContentHandler, ErrorHandler):
    """Parses an entry directly using a SANS parser. You should not ever use this class directly."""
    ent = None
    curframe = None
    curloop = None
    mode = None

    def __init__(self, entry=None):
        if entry is None:
            raise ValueError("You must provide an entry to parse into. Also, why are you using this class?")
        self.ent = entry
        self.curframe = None
        self.curloop = None

    def comment(self, line, text):
        if verbose: print "Comment '%s' on line: %d" % (text, line)
        pass

    def startData(self, line, name):
        if verbose: print "Data '%s' started on line: %d" % (name, line)
        self.ent.bmrb_id = name

    def endData(self, line, name):
        if verbose: print "Data '%s' ended on line: %d" % (name, line)
        pass

    def startSaveFrame(self, line, name):
        if verbose: print "Saveframe '%s' started on line: %d" % (name, line)
        self.curframe = saveframe.fromScratch(saveframe_name=name)
        self.ent.addSaveframe(self.curframe)

    def endSaveFrame(self, line, name):
        if verbose: print "Saveframe '%s' ended on line: %d" % (name, line)
        self.curframe = None

    def startLoop(self, line):
        if verbose: print "Loop started on line: %d" % (line)
        self.curloop = loop.fromScratch()
        self.curframe.addLoop(self.curloop)

    def endLoop(self, line):
        if verbose: print "Loop ended on line: %d" % (line)
        self.curloop = None

    def data(self, tag, tagline, val, valline, delim, inloop):

        if verbose: print "Tag / value: %s : %s ( %d : %d ) d %s" % (tag, val, tagline, valline, delim)

        if delim == 13:
                val = "$"+str(val)

        if inloop :
            # Update the columns and then add the data
            self.curloop.addColumn(tag,ignore_duplicates=True)
            self.curloop.addDataByColumn(tag,val)
        else :
            self.curframe.addTag(tag,val)

    def error(self, line, msg):
        raise ValueError("Parse error: " + str(msg),line)

    def fatalError(self, line, msg):
        raise ValueError("Fatal parse error: " + str(msg),line)

    def warning(self, line, msg):
        if raise_parse_warnings:
            raise Warning("Parse warning: " + str(msg),line)
        if verbose:
            print "Parse warning: " + str(msg) + " " + str(line)


class schema:
    """A BMRB schema."""

    headers = []
    schema = OrderedDict()
    types = OrderedDict()

    def __init__(self, schema_file=None):
        """Initialize a BMRB schema. With no arguments the most up-to-date schema will be fetched from the BMRB FTP site. Otherwise pass a URL or a file to load a schema from using the schema_url or schema_file optional arguments."""

        self.headers = []
        self.schema = OrderedDict()
        self.types = OrderedDict()

        if schema_file is None:
            schema_file = 'https://svn.bmrb.wisc.edu/svn/nmr-star-dictionary/bmrb_star_v3_files/adit_input/xlschem_ann.csv'
        self.schema_file = schema_file

        schem_stream = __interpretFile__(schema_file)
        fix_newlines = StringIO('\n'.join(schem_stream.read().splitlines()))

        csv_reader = csv.reader(fix_newlines)
        self.headers = csv_reader.next()

        # Skip the header descriptions and header index values and anything else before the real data starts
        while csv_reader.next()[0] != "TBL_BEGIN":
            pass

        for line in csv_reader:

            # Stop at the end
            if line[0] == "TBL_END":
                break

            if line[8].count(".") == 1:
                self.schema[line[8]] = (line[27],line[28],line[1])
                self.types[line[8][:line[8].index(".")]] = (line[1],line[42])
            else:
                if verbose:
                    print "Detected invalid tag in schema: %s" % str(line)

    def __repr__(self):
        """Return how we can be initialized."""
        if self.schema_file:
            return "bmrb.schema(schema_file='"+str(self.schema_file)+"')"
        else:
            return "Invalid BMRB schema."

    def __str__(self):
        """Print the schema that we are adhering to."""
        if self.schema_file:
            return "BMRB schema loaded from: '"+str(self.schema_file)+"'"
        else:
            return "BMRB schema from ??? (If you see this you did something wrong.)"

    def valType(self, tag, value, category=None,linenum=None):

        if not tag in self.schema:
            return ["Tag '" + str(tag) + "' not found in schema. Line " + str(linenum) + "."]

        valtype,null_allowed,allowed_category = self.schema[tag]

        if category != None:
            if category != allowed_category:
                return ["The tag '" + str(tag) + "' in category '" + str(category) + "' should be in category '" + str(allowed_category) + "'."]

        if value == ".":
            if null_allowed == "NOT NULL":
                return ["Value cannot be NULL but is: '" + tag + "':'" + value + "' on line " + str(linenum) + "."]
            return []

        if "VARCHAR" in valtype:
            length = valtype[valtype.index("(")+1:valtype.index(")")]
            if len(str(value)) > length:
                return ["Length of value (" + str(len(value)) + ") is too long for VARCHAR(" + length + "): '" + tag + "':'" + value + "' on line " + str(linenum) + "."]
        elif "CHAR" in valtype:
            length = valtype[valtype.index("(")+1:valtype.index(")")]
            if len(str(value)) > length:
                return ["Length of value (" + str(len(value)) + ") is too long for CHAR(" + length + "): '" + tag + "':'" + value + "' on line " + str(linenum) + "."]
        elif "FLOAT" in valtype:
            try:
                a = float(value)
            except Exception as e:
                return ["Value is not of type FLOAT.:'" + tag + "':'" + value + "' on line " + str(linenum) + "."]
        elif "INTEGER" in valtype:
            try:
                a = int(value)
            except Exception as e:
                return ["Value is not of type INTEGER: '" + tag + "':'" + value + "' on line " + str(linenum) + "."]
        return []

class entry:
    """An OO representation of a BMRB entry. You can initialize this object several ways; (e.g. from a file, from the official database, from scratch) see the classmethods."""

    # Put these here for reference
    bmrb_id = 0
    frame_list = []

    def __cmp__(self, other):
        """Returns 1 if this entry is not equal to another entry, 0 if it is equal."""
        return len(self.compare(other))

    def __delitem__(self, item):
        """Remove the indicated saveframe."""

        if isinstance(item,saveframe):
            del self.frame_list[self.frame_list.index(item)]
            return
        else:
            self.__delitem__(self.__getitem__(item))

    def __getitem__(self, item):
        """Get the indicated saveframe."""
        try:
            return self.frame_list[item]
        except TypeError:
            return self.getSaveframeByName(item)

    def __init__(self, **kargs):
        """Don't use this directly, use fromFile, fromScratch, fromString, or fromDatabase to construct."""

        # They initialized us wrong
        if len(kargs) == 0:
            raise ValueError("You must provide either a BMRB ID, a file name, an entry number, or a string to initialize. Use the class methods.")
        elif len(kargs) > 1:
            raise ValueError("You cannot provide multiple optional arguments. Use the class methods instead of initializing directly.")

        # Initialize our local variables
        self.frame_list = []

        if 'the_string' in kargs:
            # Parse from a string by wrapping it in StringIO
            star_buffer = StringIO(kargs['the_string'])
        elif 'file_name' in kargs:
            star_buffer = __interpretFile__(kargs['file_name'])
        elif 'entry_num' in kargs:
            # Parse from the official BMRB library
            try:
                star_buffer = urllib2.urlopen('http://rest.bmrb.wisc.edu/bmrb/NMR-STAR3/' + str(kargs['entry_num']))
            except urllib2.HTTPError:
                raise IOError("Entry " + str(kargs['entry_num']) + " does not exist in the public database.")
        else:
            # Initialize a blank entry
            self.bmrb_id = kargs['bmrb_id']
            return

        # Load the BMRB entry from the file
        t = __entryParser__(entry=self)
        SansParser.parser( STARLexer(star_buffer), t, t ).parse()

    def __repr__(self):
        """Returns a description of the entry."""
        return "<bmrb.entry '" + str(self.bmrb_id) + "'>"

    def __setitem__(self, key, item):
        """Set the indicated saveframe."""

        # It is a saveframe
        if isinstance(item, saveframe):
            # Add by ordinal
            try:
                self.frame_list[key] = item
            except TypeError:
                # Add by key
                if key in self.frameDict():
                    dict((frame.name,frame) for frame in self.frame_list)
                    for pos,frame in enumerate(self.frame_list):
                        if frame.name == key:
                            self.frame_list[pos] = item
                else:
                    raise KeyError("Saveframe with name '" + str(key) + "' does not exist and therefore cannot be written to. Use the addSaveframe method to add new saveframes.")
        else:
            raise ValueError("You can only assign an entry to a saveframe splice.")

    def __str__(self):
        """Returns the entire entry in STAR format as a string."""
        ret_string = "data_" + str(self.bmrb_id) + "\n\n"
        for frame in self.frame_list:
            ret_string += str(frame) + "\n"
        return ret_string

    @classmethod
    def fromFile(cls, the_file):
        """Create an entry by loading in a file. If the_file starts with http://, https://, or ftp:// then we will use those protocols to attempt to open the file."""
        return cls(file_name=the_file)

    @classmethod
    def fromDatabase(cls, entry_num):
        """Create an entry corresponding to the most up to date entry on the public BMRB server. (Requires ability to initiate outbound HTTP connections.)"""
        return cls(entry_num=entry_num)

    @classmethod
    def fromString(cls, the_string):
        """Create an entry by parsing a string."""
        return cls(the_string=the_string)

    @classmethod
    def fromScratch(cls, bmrb_id):
        """Create an empty entry that you can programatically add to. You must pass a number corresponding to the BMRB ID. If this is not a "real" BMRB entry, use 0 as the BMRB ID."""
        return cls(bmrb_id=bmrb_id)

    def addSaveframe(self, saveframe):
        """Add a saveframe to the entry."""
        self.frame_list.append(saveframe)

    def compare(self, other):
        """Returns the differences between two entries as a list. Otherwise returns 1 if different and 0 if equal. Non-equal entries will always be detected, but specific differences detected depends on order of entries."""
        diffs = []
        if self is other:
            return diffs
        if type(other) is str:
            if str(self) == other:
                return diffs
            else:
                return ['String was not exactly equal to entry.']
        try:
            if str(self.bmrb_id) != str(other.bmrb_id):
                diffs.append("BMRB ID does not match between entries: '"+str(self.bmrb_id)+"' vs '"+str(other.bmrb_id) + "'")
            if len(self.frame_list) != len(other.frame_list):
                diffs.append("The number of saveframes in the entries are not equal: "+str(len(self.frame_list))+" vs "+str(len(other.frame_list)))
            for frame in self.frameDict():
                if other.frameDict().get(frame,None) is None:
                    diffs.append("No saveframe with name '"+str(self.frameDict()[frame].name) + "' in other entry.")
                else:
                    comp = self.frameDict()[frame].compare(other.frameDict()[frame])
                    if len(comp) > 0:
                        diffs.append("Saveframes do not match: '"+str(self.frameDict()[frame].name) + "'")
                        diffs.extend(comp)
        except Exception as e:
            diffs.append("An exception occured while comparing: '" + str(e) + "'")

        return diffs

    def frameDict(self):
        """Returns a dictionary of saveframe name -> saveframe object"""
        return dict((frame.name,frame) for frame in self.frame_list)

    def getLoopsByCategory(self, value):
        """Allows fetching loops by category."""
        results = []
        for frame in self.frame_list:
            for loop in frame.loops:
                if loop.category == __formatCategory__(value):
                    results.append(loop)
        return results

    def getSaveframeByName(self, frame):
        """Allows fetching a saveframe by name."""
        frames = self.frameDict()
        if frame in frames:
            return frames[frame]
        else:
            raise KeyError("No saveframe with name '" + str(frame) + "'")

    def getSaveframesByCategory(self, value):
        """Allows fetching saveframes by category."""
        return self.getSaveframesByTagAndValue("Sf_category", value)

    def getSaveframesByTagAndValue(self, tag_name, value):
        """Allows fetching saveframe(s) by tag and tag value."""

        ret_frames = []

        for frame in self.frame_list:
            try:
                if frame.getTag(tag_name) == value:
                    ret_frames.append(frame)
            except KeyError:
                pass

        return ret_frames

    def getTagValue(self, tag):
        """ Given a tag (E.g. _Assigned_chem_shift_list.Data_file_name) return a list of all values for that tag. """
        category,tag = __formatCategory__(tag), __formatTag__(tag)

        results = []
        for frame in self.frame_list:
            if frame.tag_prefix == category:
                results.append(frame.getTag(tag))
        return results

    def printTree(self):
        """Prints a summary, tree style, of the frames and loops in the entry."""
        print repr(self)
        for pos,frame in enumerate(self):
            print "\t[" + str(pos) + "] " + repr(frame)
            for pos2,loop in enumerate(frame):
                print "\t\t[" + str(pos2) + "] " + repr(loop)

    def validate(self,validation_schema=None):
        """Validate an entry against a STAR schema. You can pass your own custom schema if desired, otherwise the schema will be fetched from the BMRB servers. Returns a list of errors found. 0-length list indicates no errors found."""

        errors = []
        for frame in self:
            errors.extend(frame.validate(validation_schema=validation_schema))
        return errors

    def exportString(self):
      """Return string version with empty loops suppressed"""
      global skip_empty_loops
      skip_empty_loops = True
      try:
        return str(self)
      finally:
        skip_empty_loops = False


class saveframe:
    """A saveframe. Use the classmethod fromScratch to create one."""

    tags = OrderedDict()
    loops = []
    name = ""
    tag_prefix = None

    def __cmp__(self, other):
        """Returns 1 if this saveframe is not equal to another saveframe, 0 if it is equal."""
        return len(self.compare(other))

    def __delitem__(self, item):
        """Remove the indicated tag or loop."""

        if isinstance(item,loop):
            del self.loops[self.loops.index(item)]
            return

        else:
            self.__delitem__(self.__getitem__(item))

    def __getitem__(self, item):
        """Get the indicated loop or tag."""
        try:
            return self.loops[item]
        except TypeError:
            try:
                return self.getTag(item)
            except KeyError as e:
                try:
                    return self.loopDict()[item]
                except KeyError:
                    raise KeyError("No tag or loop matching '" + str(item) + "'")

    def __len__(self):
        """Return the number of loops in this saveframe."""
        return len(self.loops)

    def __init__(self, **kargs):
        """Don't use this directly. Use the class methods to construct."""

        # They initialized us wrong
        if len(kargs) == 0:
            raise ValueError("Use the class methods to initialize.")

        # Initialize our local variables
        self.tags = OrderedDict()
        self.loops = []

        if 'the_string' in kargs:
            # Parse from a string by wrapping it in StringIO
            star_buffer = StringIO(kargs['the_string'])
        elif 'file_name' in kargs:
            star_buffer = __interpretFile__(kargs['file_name'])
        elif 'saveframe_name' in kargs:
            # If they are creating from scratch, just get the saveframe name
            self.name = kargs['saveframe_name']
            if 'tag_prefix' in kargs:
                self.tag_prefix = __formatCategory__(kargs['tag_prefix'])
            return

        # If we are reading from a CSV file, go ahead and parse it
        if 'csv' in kargs and kargs['csv'] is True:
            csvreader = csv.reader(star_buffer)
            tags = csvreader.next()
            values = csvreader.next()
            if len(tags) != len(values):
                raise ValueError("Your CSV data is invalid. The header length does not match the data length.")
            for x in range(0,len(tags)):
                self.addTag(tags[x], values[x])
            return

        star_buffer = StringIO("data_1\n"+star_buffer.read())
        tmp_entry = entry.fromScratch(0)

        # Load the BMRB entry from the file
        t = __entryParser__(entry=tmp_entry)
        SansParser.parser( STARLexer(star_buffer), t, t ).parse()

        # Copy the first parsed saveframe into ourself
        self.tags = tmp_entry[0].tags
        self.loops = tmp_entry[0].loops
        self.name = tmp_entry[0].name
        self.tag_prefix = tmp_entry[0].tag_prefix

    @classmethod
    def fromScratch(cls, saveframe_name, tag_prefix=None):
        """Create an empty saveframe that you can programatically add to. You may also pass the tag prefix as the second argument. If you do not pass the tag prefix it will be set the first time you add a tag."""
        return cls(saveframe_name=saveframe_name, tag_prefix=tag_prefix)

    @classmethod
    def fromFile(cls, the_file, csv=False):
        """Create a saveframe by loading in a file. Specify csv=True is the file is a CSV file. If the_file starts with http://, https://, or ftp:// then we will use those protocols to attempt to open the file."""
        return cls(file_name=the_file, csv=csv)

    @classmethod
    def fromString(cls, the_string, csv=False):
        """Create a saveframe by parsing a string. Specify csv=True is the string is in CSV format and not NMR-STAR format."""
        return cls(the_string=the_string, csv=csv)

    def __repr__(self):
        """Returns a description of the saveframe."""
        return "<bmrb.saveframe '" + str(self.name) + "'>"

    def __setitem__(self, key, item):
        """Set the indicated loop or tag."""

        # It's a loop
        if isinstance(item,loop):
            try:
                integer = int(str(key))
                self.loops[integer] = item
            except ValueError:
                if key in self.loopDict():
                    for pos,tmp_loop in enumerate(self.loops):
                        if tmp_loop.category == key:
                            self.loops[pos] = item
                else:
                    raise KeyError("Loop with category '" + str(key) + "' does not exist and therefore cannot be written to. Use addLoop instead.")
        else:
            # If the tag already exists, set its value
            self.addTag(key,item,update=True)

    def __str__(self):
        """Returns the saveframe in STAR format as a string."""

        if self.tag_prefix == None:
            raise ValueError("The tag prefix was never set!")

        # Make sure this isn't a dummy saveframe before proceeding
        try:
            width = max(map(lambda x:len(self.tag_prefix+"."+x), self.tags))
        except ValueError:
            return "\nsave_" + str(self.name) + "\n\nsave_\n"

        # Print the saveframe
        ret_string = "save_" + str(self.name) + "\n"
        pstring = "   %-" + str(width) + "s  %s\n"
        mstring = "   %-" + str(width) + "s\n;\n%s;\n"

        # Print the tags
        for tag in self.tags:
            cleanTag = __cleanValue__(self.tags[tag])
            if "\n" in cleanTag:
                ret_string +=  mstring % (self.tag_prefix+"."+tag, cleanTag)
            else:
                ret_string +=  pstring % (self.tag_prefix+"."+tag, cleanTag)

        # Print any loops
        for loop in self.loops:
            ret_string +=  str(loop)

        # Close 'er up
        ret_string += "save_\n"
        return ret_string

    def addLoop(self, loop):
        """Add a loop to the saveframe loops."""

        if loop.category in self.loopDict():
            if loop.category is None:
                raise ValueError("You cannot have two loops with the same category in one saveframe. You are getting this error because you haven't yet set your loop categories.")
            else:
                raise ValueError("You cannot have two loops with the same category in one saveframe. Category: '" + str(loop.category) + "'.")

        self.loops.append(loop)

    def addTag(self, name, value, update=False):
        """Add a tag to the tag list. Does a bit of validation and parsing. Set update to true to update a tag if it exists rather than raise an exception."""

        if "." in name:
            if name[0] != ".":
                prefix = __formatCategory__(name)
                if self.tag_prefix == None:
                    self.tag_prefix = prefix
                elif self.tag_prefix != prefix:
                    raise ValueError("One saveframe cannot have tags with different prefixes (or tags that don't match the set prefix)! " + str(self.tag_prefix) + " vs " + str(prefix))
                name = name[name.index(".")+1:]
            else:
                name = name[1:]

        # No duplicate tags
        try:
            self.getTag(name)
            if not update:
                raise ValueError("There is already a tag with the name '" + name + "'.")
            else:
                self.tags[name] = value
                return
        except KeyError:
            pass

        if "." in name:
            raise ValueError("There cannot be more than one '.' in a tag name.")
        if " " in name:
            raise ValueError("Column names can not contain spaces.")

        if verbose:
            print "Adding tag: ("+name+") with value ("+value+")"

        self.tags[name] = value

    def addTags(self, tag_list, update=False):
        """Adds multiple tags to the list. Input should be a list of tuples that are either [key,value] or [key]. In the latter case the value will be set to ".".  Set update to true to update a tag if it exists rather than raise an exception."""
        for tag_pair in tag_list:
            if len(tag_pair) == 2:
                self.addTag(tag_pair[0],tag_pair[1],update=update)
            elif len(tag_pair) == 1:
                self.addTag(tag_pair[0],".",update=update)
            else:
                raise ValueError("You provided an invalid tag/value or tag to add: %s" % str(tag_pair))

    def compare(self, other):
        """Returns the differences between two saveframes as a list. Non-equal saveframes will always be detected, but specific differences detected depends on order of saveframes."""
        diffs = []

        # Check if this is literally the same object
        if self is other:
            return diffs
        # Check if the other object is our string representation
        if type(other) is str:
            if str(self) == other:
                return diffs
            else:
                return ['String was not exactly equal to saveframe.']
        # Do STAR comparison
        try:
            if str(self.name) != str(other.name):
                diffs.append("\tSaveframe names do not match: "+str(self.name)+"' vs '"+str(other.name) + "'")
            if str(self.tag_prefix) != str(other.tag_prefix):
                diffs.append("\tTag prefix does not match: "+str(self.tag_prefix)+"' vs '"+str(other.tag_prefix) + "'")
            if len(self.tags) != len(other.tags):
                diffs.append("\tNumber of tags does not match: "+str(len(self.tags))+"' vs '"+str(len(other.tags)) + "'")
            for tag in self.tags:
                if self.tags[tag] != other.tags.get(tag,None):
                    diffs.append("\tMismatched tag values for tag '" + str(tag) + "': '" + str(self.tags[tag]).replace("\n","\\n")+"' vs '"+str(other.tags.get(tag,"[not present]")).replace("\n","\\n") + "'")
            if len(self.loops) != len(other.loops):
                diffs.append("\tNumber of children loops does not match: "+str(len(self.loops))+" vs "+str(len(other.loops)))
            for loop in self.loopDict():
                if other.loopDict().get(loop,None) is None:
                    diffs.append("\tNo loop with category '"+str(self.loopDict()[loop].category) + "' in other entry.")
                else:
                    comp = self.loopDict()[loop].compare(other.loopDict()[loop])
                    if len(comp) > 0:
                        diffs.append("\tLoops do not match: '"+str(self.loopDict()[loop].category) + "'")
                        diffs.extend(comp)
        except Exception as e:
            diffs.append("\tAn exception occured while comparing: '" + str(e) + "'")

        return diffs

    def deleteTag(self, tag):
        """Deletes a tag from the saveframe based on tag name."""
        tag = __formatTag__(tag)
        if tag in self.tags:
            del self.tags[tag]
        else:
            raise KeyError("There is no tag with name " + str(tag) + " to remove!")

    def getDataAsCSV(self, header=True, show_category=True):
        """Return the data contained in the loops, properly CSVd, as a string. Set header to False omit the header. Set show_category to False to omit the loop category from the headers."""
        csv_buffer = StringIO()
        cwriter = csv.writer(csv_buffer)
        if header:
            if show_category:
                cwriter.writerow([str(self.tag_prefix)+"."+str(x) for x in self.tags])
            else:
                cwriter.writerow([str(x) for x in self.tags])

        data = []
        for piece in [self.tags[x] for x in self.tags]:
            #if type(piece) is str and piece[0] == "$":
            #    piece = piece[1:]
            data.append(piece)

        cwriter.writerow(data)

        csv_buffer.seek(0)
        return csv_buffer.read().replace('\r\n', '\n')

    def getLoopByCategory(self, name):
        """Return a loop based on the loop name (category)."""
        name = __formatCategory__(name)
        for loop in self.loops:
            if loop.category == name:
                return loop
        raise KeyError("No loop with category " + name)

    def getTag(self, query):
        """Allows fetching a tag by name."""

        query = __formatTag__(query)

        if query in self.tags:
            return self.tags[query]
        if self.tag_prefix and (self.tag_prefix + "." + query) in self.tags:
            return self.tags[self.tag_prefix + "." + query]

        raise KeyError("No tag with name '" + str(query) + "'")

    def loopDict(self):
        """Returns a hash of loop category -> loop."""
        return dict((loop.category,loop) for loop in self.loops)

    def loopIterator(self):
        """Returns an iterator for saveframe loops."""
        for loop in self.loops:
            yield loop

    def setTagPrefix(self, tag_prefix):
        """Set the tag prefix for this saveframe."""
        self.tag_prefix = __formatCategory__(tag_prefix)

    def sortTags(self, validation_schema=None):
        """ Sort the tags and loops so they are in the same order as a BMRB schema. Will automatically use the standard schema if none is provided."""

        my_schema = __getSchema__(validation_schema)

        # Sort the tag list
        new_tag_list = OrderedDict()
        for item in my_schema.schema:
            if item[:item.index(".")] == self.tag_prefix:
                tag = item[item.index(".")+1:]
                if tag in self.tags:
                    new_tag_list[tag] = self.tags[tag]
        self.tags = new_tag_list

        # Now sort the loops
        tag_order = my_schema.types.keys()
        self.loops.sort(key=lambda x: tag_order.index(x.category))

        # Sort the loops
        full_ordering = my_schema.schema.keys()

        for tmp_loop in self.loops:
            new_data = []
            new_order = sorted(tmp_loop.columns, key=lambda x:full_ordering.index(tmp_loop.category+"."+x))

            for row in tmp_loop.data:
                temp_data = [None]*len(row)
                for pos in range(0,len(row)):
                    temp_data[new_order.index(tmp_loop.columns[pos])] = row[pos]
                    new_data.append([row[tmp_loop.columns.index(x)] for row in tmp_loop.data])

            tmp_loop.data = new_data
            tmp_loop.columns = new_order

            #for column in tmp_loop.columns:
            #    print column,full_ordering.index(tmp_loop.category+"."+column)
            # Here is loop sorting code - TODO: move to loops

    def tagIterator(self):
        """Returns an iterator for saveframe tags."""
        for tag in self.tags:
            yield [tag,self.tags[tag]]

    def printTree(self):
        """Prints a summary, tree style, of the loops in the saveframe."""
        print repr(self)
        for pos,loop in enumerate(self):
            print "\t[" + str(pos) + "] " + repr(loop)

    def validate(self,validation_schema=None):
        """Validate a saveframe against a STAR schema. You can pass your own custom schema if desired, otherwise the schema will be fetched from the BMRB servers. Returns a list of errors found. 0-length list indicates no errors found."""

        # Get the default schema if we are not passed a schema
        my_schema = __getSchema__(validation_schema)

        errors = []
        for pos,tag in enumerate(self.tags):
            errors.extend(my_schema.valType(self.tag_prefix + "." + tag, self.tags[tag], category=self.getTag("Sf_category"), linenum=pos))

        for loop in self.loops:
            errors.extend(loop.validate(validation_schema=validation_schema, category=self.getTag("Sf_category")))

        return errors


class loop:
    """A loop."""

    category = None
    columns = []
    data = []

    def __cmp__(self, other):
        """Returns 1 if this loop is not equal to another loop, 0 if it is equal."""
        return len(self.compare(other))

    def __getitem__(self, item):
        """Get the indicated row from the data array."""
        try:
            return self.data[item]
        except TypeError:
            if type(item) is tuple:
                item = list(item)
            return self.getDataByTag(tags=item)

    def __len__(self):
        """Return the number of rows of data."""
        return len(self.data)

    def __init__(self, **kargs):
        """Use the classmethods to initialize."""

        # Initialize our local variables
        self.columns = []
        self.data = []
        self.category = None

        # Update our category if provided
        if 'category' in kargs:
            self.category = __formatCategory__(kargs['category'])
            return

        # They initialized us wrong
        if len(kargs) == 0:
            raise ValueError("Use the class methods to initialize.")

        if 'the_string' in kargs:
            # Parse from a string by wrapping it in StringIO
            star_buffer = StringIO(kargs['the_string'])
        elif 'file_name' in kargs:
            star_buffer = __interpretFile__(kargs['file_name'])
        elif 'saveframe_name' in kargs:
            # If they are creating from scratch, just get the saveframe name
            self.name = kargs['saveframe_name']
            if 'tag_prefix' in kargs:
                self.tag_prefix = __formatCategory__(kargs['tag_prefix'])
            return

        # If we are reading from a CSV file, go ahead and parse it
        if 'csv' in kargs and kargs['csv'] is True:
            csvreader = csv.reader(star_buffer)
            self.addColumn(csvreader.next())
            for row in csvreader:
                self.addData(row)
            return

        star_buffer = StringIO("data_0\nsave_dummy_frame\n" + star_buffer.read() + "save_")
        tmp_entry = entry.fromScratch(0)

        # Load the BMRB entry from the file
        t = __entryParser__(entry=tmp_entry)
        SansParser.parser( STARLexer(star_buffer), t, t ).parse()

        # Copy the first parsed saveframe into ourself
        self.columns = tmp_entry[0][0].columns
        self.data = tmp_entry[0][0].data
        self.category = tmp_entry[0][0].category

    @classmethod
    def fromFile(cls, the_file, csv=False):
        """Create a saveframe by loading in a file. Specify csv=True is the file is a CSV file. If the_file starts with http://, https://, or ftp:// then we will use those protocols to attempt to open the file."""
        return cls(file_name=the_file, csv=csv)

    @classmethod
    def fromScratch(cls, category=None):
        """Create an empty saveframe that you can programatically add to. You may also pass the tag prefix as the second argument. If you do not pass the tag prefix it will be set the first time you add a tag."""
        return cls(category=category)

    @classmethod
    def fromString(cls, the_string, csv=False):
        """Create a saveframe by parsing a string. Specify csv=True is the string is in CSV format and not NMR-STAR format."""
        return cls(the_string=the_string, csv=csv)

    def __repr__(self):
        """Returns a description of the loop."""
        return "<bmrb.loop '" + str(self.category) + "'>"

    def __str__(self):
        """Returns the loop in STAR format as a string."""

        if skip_empty_loops and not self.data:
          # Rasmus Fogh October 2014. Suppress empty loops
          return ''

        # If we have no columns than return without printing
        if len(self.columns) == 0:
            if len(self.data) == 0:
                return "\n   loop_\n\n   stop_\n"
            else:
                raise ValueError("Impossible to print data if there are no associated tags. Loop: " + str(self.category))

        # Make sure that if there is data, it is the same width as the column tags
        if len(self.data) > 0:
            for row in self.data:
                if len(self.columns) != len(row):
                    raise ValueError("The number of column tags must match width of the data. Loop: " + str(self.category))

        # Start the loop
        ret_string = "\n   loop_\n"
        # Print the columns
        pstring = "      %-s\n"
        # Check to make sure our category is set
        if self.category == None:
            raise ValueError("The category was never set for this loop. Either add a column with the category intact, specify it when generating the loop, or set it using setCategory.")
        # Print the categories
        for column in self.columns:
            ret_string += pstring % (self.category + "." + column)
        ret_string += "\n"

        if len(self.data) != 0:
            # The nightmare below creates a list of the maximum length of elements in each column in the self.data matrix
            #  Don't try to understand it. It's an imcomprehensible list comprehension.
            title_widths = [max(map(lambda x:len(str(x))+3, col)) for col in [[row[x] for row in self.data] for x in range(0,len(self.data[0]))]]
            # Generate the format string
            pstring = "     " + "%-*s"*len(self.columns) + " \n"

            # Print the data, with the columns sized appropriately
            for datum in copy.deepcopy(self.data):

                # Put quotes as needed on the data
                datum = map(lambda x:__cleanValue__(x), datum)
                for pos,item in enumerate(datum):
                    if "\n" in item:
                        datum[pos] = "\n;\n" + item + ";\n"

                # Print the data (combine the columns widths with their data)
                ret_string += pstring % tuple(itertools.chain.from_iterable([d for d in zip(title_widths,datum)]))

        # Close the loop
        ret_string += "   stop_\n"
        return ret_string

    def addColumn(self, name, ignore_duplicates=False):
        """Add a column to the column list. Does a bit of validation and parsing. Set ignore_duplicates to true to ignore attempts to add the same tag more than once rather than raise an exception.

        You can also pass a list of column names to add more than one column at a time."""

        # If they have passed multiple columns to add, call ourself on each of them in succession
        if isinstance(name, (list, tuple)):
            for x in name:
                self.addColumn(x, ignore_duplicates=ignore_duplicates)
            return

        name = name.strip()

        if "." in name:
            if name[0] != ".":
                category = name[0:name.index(".")]
                if category[:1] != "_":
                    category = "_" + category

                if self.category == None:
                    self.category = category
                elif self.category != category:
                    raise ValueError("One loop cannot have columns with different categories (or columns that don't match the set prefix)!")
                name = name[name.index(".")+1:]
            else:
                name = name[1:]

        # Ignore duplicate tags
        if name in self.columns:
            if ignore_duplicates:
                return
            else:
                raise ValueError("There is already a column with the name '" + name + "'.")
        if "." in name:
            raise ValueError("There cannot be more than one '.' in a tag name.")
        if " " in name:
            raise ValueError("Column names can not contain spaces.")
        self.columns.append(name)

    def addData(self, the_list):
        """Add a list to the data field. Items in list can be any type, they will be converted to string and formatted correctly. The list must have the same cardinality as the column names."""
        if len(the_list) != len(self.columns):
            raise ValueError("The list must have the same number of elements as the number of columns! Insert column names first.")
        # Add the user data
        self.data.append(the_list)

    def addDataByColumn(self, column_id, value):
        """Add data to the loop one element at a time, based on column. Useful when adding data from SANS parsers."""
        column_id = __formatTag__(column_id)
        if not column_id in self.columns:
            raise ValueError("The column tag (" + column_id + ") to which you are attempting to add data does not yet exist. Create the columns before adding data.")
        pos = self.columns.index(column_id)
        if len(self.data) == 0:
            self.data.append([])
        if len(self.data[-1]) == len(self.columns):
            self.data.append([])
        if len(self.data[-1]) != pos:
            raise ValueError("You cannot add data out of column order.")
        self.data[-1].append(value)

    def clearData(self):
        """Erases all data in this loop. (But not the data columns.)"""
        self.data = []

    def compare(self, other):
        """Returns the differences between two loops as a list. Otherwise returns 1 if different and 0 if equal. Order of loops being compared does not make a difference on the specific errors detected."""
        diffs = []

        # Check if this is literally the same object
        if self is other:
            return diffs
        # Check if the other object is our string representation
        if type(other) is str:
            if str(self) == other:
                return diffs
            else:
                return ['String was not exactly equal to loop.']
        # Do STAR comparison
        try:
            if str(self.category) != str(other.category):
                diffs.append("\t\tCategory of loops does not match: '"+str(self.category)+"' vs '"+str(other.category) + "'")
            if self.data != other.data or self.columns != other.columns:
                diffs.append("\t\tLoop data does not match for loop with category '"+str(self.category) + "'")
        except Exception as e:
            diffs.append("\t\tAn exception occured while comparing: '" + str(e) + "'")

        return diffs

    def deleteDataByTagValue(self, tag, value, index_tag=None):
        """Deletes all rows which contain the provided value in the provided column. If index_tag is provided, that column is renumbered starting with 1."""

        # Parse the category and tag
        full = str(self.category) + "." + __formatTag__(str(tag))
        tag,category = __formatTag__(full),__formatCategory__(full)

        if category != self.category:
            raise ValueError("Category provided in your column '" + category + "' does not match this loop's category '" + str(self.category) + "'.")

        # Can't delete a row if the tag they want to delete isn't there!
        if tag not in self.columns:
            raise ValueError("The tag you provided '" + tag + "' isn't in this loop!")

        search_column = self.columns.index(tag)

        # Delete all rows in which the user-provided tag matched
        cur_row = 0
        while cur_row < len(self.data):
            if self.data[cur_row][search_column] == value:
                self.data.pop(cur_row)
                continue
            cur_row += 1

        # Re-number if they so desire
        if index_tag is not None:
            self.renumberRows(index_tag)

    def getDataByTag(self, tags=None):
        """Provided a list of tag names, or ordinals corresponding to columns, return the selected tags by row as a list of lists."""

        # All tags
        if tags is None:
            return self.data
        # Turn single elements into lists
        if type(tags) != list:
            tags = [tags]

        # Strip the category if they provide it
        for pos,item in enumerate(map(lambda x:str(x),tags)):
            if "." in item and __formatCategory__(item) != self.category:
                raise ValueError("Cannot fetch data with column (" + str(item) + ") because the category does not match the category of this loop (" + self.category + ").")
            tags[pos] = __formatTag__(item)

        # Map column name to column position in list
        column_mapping = dict(itertools.izip(reversed(self.columns), reversed(xrange(len(self.columns)))))

        # Make sure their fields are actually present in the entry
        column_ids = []
        for query in tags:
            if str(query) in column_mapping:
                column_ids.append(column_mapping[query])
            elif type(query) == int:
                column_ids.append(query)
            else:
                raise ValueError("Your column '" + str(query) + "' was not found in this loop.")

        # Use a list comprehension to pull the correct tags out of the rows
        return [[row[col_id] for col_id in column_ids] for row in self.data]

    def getDataAsCSV(self, header=True, show_category=True):
        """Return the data contained in the loops, properly CSVd, as a string. Set header to False to omit the header. Set show_category to false to omit the loop category from the headers."""
        csv_buffer = StringIO()
        cwriter = csv.writer(csv_buffer)

        if header:
            if show_category:
                cwriter.writerow([str(self.category)+"."+str(x) for x in self.columns])
            else:
                cwriter.writerow([str(x) for x in self.columns])

        for row in self.data:

            data = []
            for piece in row:
                data.append(piece)

            cwriter.writerow(data)

        csv_buffer.seek(0)
        return csv_buffer.read().replace('\r\n', '\n')

    def printTree(self):
        """Prints a summary, tree style, of the loop."""
        print repr(self)

    def renumberRows(self, index_tag, start_value=1, maintain_ordering=False):
        """Renumber a given column incrementally. Set start_value to initial value if 1 is not acceptable. Set maintain_ordering to preserve sequence with offset. E.g. 2,3,3,5 would become 1,2,2,4."""

        full = str(self.category) + "." + __formatTag__(str(index_tag))
        tag,category = __formatTag__(full),__formatCategory__(full)

        # The column to replace in is the column they specify
        try:
            renumber_column = self.columns.index(tag)
        except ValueError:
            # Or, perhaps they specified an integer to represent the column?
            try:
                renumber_column = int(tag)
            except ValueError:
                raise ValueError("The renumbering column you provided '%s' isn't in this loop!" % tag)

        # Verify the renumbering column ID
        if renumber_column >= len(self.columns) or renumber_column < 0:
            raise ValueError("The renumbering column ID you provided '%s' is too large or too small! Value column ids are 0-%d." % (tag, len(self.columns)-1))

        # Do nothing if we have no data
        if len(self.data) == 0:
            return

        if maintain_ordering:
            # If they have a string buried somewhere in the row, we'll have to restore the original values
            data_copy = copy.deepcopy(self.data)

            for x in range(0,len(self.data)):
                try:
                    if x == 0:
                        offset = start_value - int(self.data[0][renumber_column])
                    self.data[x][renumber_column] = int(self.data[x][renumber_column]) + offset
                except ValueError:
                    self.data = data_copy
                    raise ValueError("You can't renumber a row containing anything that can't be coerced into an integer using maintain_ordering. I.e. what am I suppose to renumber '%s' to?" % self.data[x][renumber_column])

        # Simple renumbering algorithm if we don't need to maintain the ordering
        else:
            for x in range(0,len(self.data)):
                self.data[x][renumber_column] = x + start_value

    def setCategory(self, category):
        """ Set the category of the loop. Usefull if you didn't know the category at loop creation time."""
        self.category = __formatCategory__(category)

    def sortRows(self, tags, cmp=None):
        """ Sort the data in the rows by their values for a given column or columns. Specify the columns using their names or ordinals. Accepts a list or an int/float. By default we will sort numerically. If that fails we do a string sort. Supply a function as cmp and we will sort based on that function. See the help for sorted() for more details. If you provide multiple columns to sort by, they are interpreted as increasing order of sort priority."""

        # Do nothing if we have no data
        if len(self.data) == 0:
            return

        # This will determine how we sort
        sort_ordinals = []

        processing_list = []
        if type(tags) is list:
            processing_list = tags
        else:
            processing_list = [tags]

        # Process their input to determine which columns to operate on
        for cur_tag in processing_list:
            full = str(self.category) + "." + __formatTag__(str(cur_tag))
            tag,category = __formatTag__(full),__formatCategory__(full)

            # The column to replace in is the column they specify
            try:
                renumber_column = self.columns.index(tag)
            except ValueError:
                # Or, perhaps they specified an integer to represent the column?
                try:
                    renumber_column = int(tag)
                except ValueError:
                    raise ValueError("The sorting column you provided '%s' isn't in this loop!" % tag)

            # Verify the renumbering column ID
            if renumber_column >= len(self.columns) or renumber_column < 0:
                raise ValueError("The sorting column ID you provided '%s' is too large or too small! Value column ids are 0-%d." % (tag, len(self.columns)-1))

            sort_ordinals.append(renumber_column)

        # Do the sort(s)
        for column in sort_ordinals:
            # Going through each column, first attempt to sort as integer. Then fallback to string sort.
            try:
                tmp_data = sorted(self.data,key=lambda x:float(x[column]),cmp=cmp)
            except ValueError:
                tmp_data = sorted(self.data,key=lambda x:x[column],cmp=cmp)
            self.data = tmp_data

    def validate(self,validation_schema=None, category=None):
        """Validate a loop against a STAR schema. You can pass your own custom schema if desired, otherwise the schema will be fetched from the BMRB servers. Returns a list of errors found. 0-length list indicates no errors found."""

        # Get the default schema if we are not passed a schema
        my_schema = __getSchema__(validation_schema)

        errors = []

        # Check the data
        for rownum,row in enumerate(self.data):
            # Make sure the width matches
            if len(row) != len(self.columns):
                errors.append("Loop '" + str(self.category) + "' data width does not match it's column tag width on row " + str(rownum) + ".")
            for pos,datum in enumerate(row):
                errors.extend(my_schema.valType(self.category + "." + self.columns[pos], datum, category=category, linenum=str(rownum)+" column "+str(pos)))

        return errors

try:
    from json import JSONEncoder

    class JSON_Encoder(JSONEncoder):
        """ Used with the python json library to convert NMR-STAR to JSON. """

        def default(self, o):

            if isinstance(o,entry):
                return { "data_" + o.bmrb_id : o.frame_list }
            elif isinstance(o,saveframe):
                saveframe_hash = OrderedDict()
                saveframe_hash[o.tag_prefix] = o.tags
                for each_loop in o:
                    saveframe_hash[each_loop.category] = each_loop

                return saveframe_hash
            elif isinstance(o,loop):
                return [dict(zip(o.columns,row)) for row in o.data]
            else:
                return JSONEncoder.default(self, o)



except:
    class JSON_Encoder():
        def __init__(self, *args, **kwargs):
            raise RuntimeError("No json module! Cannot encode to JSON. Consider updating python.")


standard_schema = None

# Do unit tests if we are ran directly
if __name__ == '__main__':

    import optparse
    # Specify some basic information about our command
    parser = optparse.OptionParser(usage="usage: %prog",version="%prog",description="STAR handling python module. Usually you'll want to import this. When ran directly a self test is performed.")
    parser.add_option("--diff", action="store_true", dest="diff", default=False, help="Compare two entries.")
    # Options, parse 'em
    (options, cmd_input) = parser.parse_args()

    if options.diff:
        if len(cmd_input) < 2:
            raise ValueError("You must supply two file names as arguments.")
        diff(entry.fromFile(cmd_input[0]), entry.fromFile(cmd_input[1]))
        sys.exit(0)

    print "Running unit tests..."

    errors = 0

    def printError(e,x,ent_str):
        if len(e.args) == 1:
            print str(x)+": ",str(e)
        else:
            print str(x)+": ",str(e.args[0]),"on line",(e.args[1]-1)
            splitted = ent_str.split("\n")
            with open("/tmp/"+str(x),"w") as tmp:
                tmp.write(ent_str)
            for x in range(e.args[1]-5,e.args[1]+2):
                try:
                    print "\t%-5d: %s" % (x,splitted[x])
                except IndexError:
                    print "\t%-5d: %s" % (x,"EOF")
                    return

    myrange = (15000,15200)
    if len(sys.argv) == 3: myrange = (int(sys.argv[1]),int(sys.argv[2]))
    if len(sys.argv) == 2: myrange = (int(sys.argv[1]), (int(sys.argv[1])+1))

    use_stardiff = False
    if os.path.exists("/bmrb/linux/bin/stardiff"):
        import subprocess
        use_stardiff = True
        print "External stardiff detected. Will use to verify results."

    for x in xrange(*myrange):
        try:
            orig_str = urllib2.urlopen('http://rest.bmrb.wisc.edu/bmrb/NMR-STAR3/' + str(x)).read()
        except urllib2.HTTPError:
            continue
        try:
            ent = entry.fromString(orig_str)
            ent_str = str(ent)
        except IOError:
            continue
        except Exception as e:
            printError(e,x,orig_str)
            errors += 1
            continue
        try:
            reent = entry.fromString(ent_str)
            if ent_str != str(reent):
                print str(x)+": Inconsisent output when re-parsed."
                errors += 1
        except Exception as e:
            printError(e,x,ent_str)
            errors += 1
            continue

        if use_stardiff:
            # Write our data to a pipe for stardiff to read from
            with open("/tmp/comparator1","wb") as tmp1:
                tmp1.write(str(ent))
            with open("/tmp/comparator2","wb") as tmp2:
                tmp2.write(str(orig_str))

            compare = subprocess.Popen(["/bmrb/linux/bin/stardiff","-ignore-tag","_Spectral_peak_list.Text_data","/tmp/comparator1","/tmp/comparator2"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

            # Wait for stardiff to complete
            compare.wait()
            results = compare.stdout.read()
            if not "NO DIFFERENCES REPORTED" in results:
                print str(x)+": Output inconsistent with original: " + results.strip()
                open("/tmp/" + str(x),"wb").write(str(ent_str))
                errors += 1

        comp = ent.compare(reent)
        if len(comp) > 0:
            print str(x)+": Internal entry comparator detects difference(s):"
            diff(ent,reent)

    if errors == 0:
        print "If you didn't see any errors, than everything is working!"
    else:
        print "At least %d errors were found." % (errors)

    sys.exit(0)
