"""
Classes for dealing with STAR syntax
"""

__author__    = "$Author: jurgenfd $"
___revision__ = "$Revision: 1.8.4.2 $"
___date__     = "$Date: 2011-04-06 09:00:38 $"

"""
$Log: not supported by cvs2svn $
Revision 1.8.4.1  2010/02/11 12:44:33  wimvranken
Fix to deal with loop_ only chemical shift files - these are now modified on-the-fly in order to read them.

Revision 1.8  2007/07/09 12:22:27  wimvranken
Merge with branch 4.

Revision 1.7.4.1.4.2  2007/06/13 08:35:52  wimvranken
Removed copyright - is not CCPN

Revision 1.7.4.1.4.1  2006/10/25 10:57:09  wimvranken
Fixed problem that shows up in windows only (\ to continue end of line)

Revision 1.7.4.1  2006/07/21 17:29:01  wimvranken
Now only contain full path imports

Revision 1.7  2005/11/03 14:06:39  wimvranken
Replaced var text_lenght by text_length, added extra error message printout

Revision 1.6  2005/06/27 16:43:41  wb104
Updated licenses.

Revision 1.5  2004/12/15 17:57:20  tjs23
TJS: Updated licenses.

Revision 1.4  2003/08/08 12:45:54  wfv20
Changed preferred quote, added comment writing support

Revision 1.3  2003/07/14 09:34:32  wfv20
Modified so list references are not carried through class initialization

Revision 1.2  2003/07/10 16:13:31  wfv20
Changed to new universal setup

Revision 1.1  2003/07/01 12:56:18  wfv20
Jurgen Doreleijers modified Python nmrStar reader. Bug fixes and additions (notably for reading nmrView star files) have been made.

Revision 1.1.1.1  2001/11/02 20:16:40  jurgen
Initial package capable of read/write access to STAR files without nested loops

"""
## Standard modules
import os

## BMRB modules
#import ccp.format.nmrStar.bmrb as STAR # TODO: check
from ccp.format.nmrStar.bmrb.SaveFrame      import * #@UnusedWildImport
from ccp.format.nmrStar.bmrb.TagTable       import * #@UnusedWildImport
from ccp.format.nmrStar.bmrb.Text           import * #@UnusedWildImport
from ccp.format.nmrStar.bmrb.Utils          import * #@UnusedWildImport

## CCP stuff
from memops.universal.Io import getTopDirectory

if os.name != 'posix': 
    tempdir     = None # The module tempfile will pick one
else:
    tempdir     = '/tmp'

# No changes required after this line
###############################################################################    

"""
STAR file
Only methods for reading and writing are currently implemented.
datanodes is a list of possibly mixed saveframes and tagtables
"""
class File (Lister):
    def __init__( self,
                    title                   = 'general_star_file_title',
                    filename                = '',
                    datanodes               = None,
                    flavor                  = None, # Call set_flavor when changing
                    preferred_quote         = '"',  # Put somewhere else?
                    verbosity   = 2
                  ):
        self.title      = title
        self.filename   = filename
        
        if datanodes:
          self.datanodes  = datanodes
        else:
          self.datanodes = []
          
        self.flavor     = flavor
        self.verbosity  = verbosity
        
    "Simple checks on integrity"
    def check_integrity( self,  recursive = 1  ):
        if recursive:
            for datanode in self.datanodes:
                if datanode.check_integrity( recursive = 1):
                    print "ERROR: integrity check failed for Saveframe"
                    return 1
        if self.verbosity >= 9:
            print 'Checked integrity of File    (%2s datanodes,  recurs.=%s)  : OK [%s]' % (
                len(self.datanodes), recursive, self.title )

    "Returns the STAR text representation"
    def star_text(self, flavor = None):
        if flavor == None:
            flavor = self.flavor
        str = 'data_%s\n' % self.title
        # Data node objects can be of type SaveFrame OR TagTable only
        # Data node object can now also contain comment information
        #      these comments are printed before the saveframe (Wim 2003/08/05)
        for datanode in self.datanodes:
            str = str + datanode.comment
            str = str + datanode.star_text( flavor = flavor)
        return str


    """
    Reads a NMR-STAR formatted file using
    the filename attribute.
    
    Added option to pass text in straight away - this is for handling crappy loop_ only files so
    I can fix them in memory (Wim 11/02/10).
    """
    
    def read (self, strip_comments=1, nmrView_type = 0, text = ""):

        if not text:
          if not self.filename:
              print 'ERROR: no filename in STARFile with title:', self.title
              return 1
          text = open(self.filename,'r').read()
          
        if self.parse( text=text, strip_comments=strip_comments, nmrView_type = nmrView_type):
            print "ERROR: couldn't parse file"
            return 1
         
        return 0

    """
    - Parses text into save frames and tagtables.
    - Input text should start at position given with non-white space character
    - Appends a list of datanodes(save frames or tagtables)
    """
    def parse (self, text='', strip_comments=1, nmrView_type = 0):

        if self.verbosity > 1:        
            print 'Parsing STAR file:', self.filename

        """
        '"Begin at the beginning," the King said, gravely,
        "and go on till you come to the end; then stop."' (LC)
        """

        ## Collapse the semicolon block for ease of parsing
        ## Very expensive to do
        ## Timed at: xx seconds for a xx Mb file with xx semicolon blocks
        text = semicolon_block_collapse( text )

        ## Now it's easy to strip comments
        if strip_comments:
            text = comments_strip( text )
        
        
        ## For nmrView 'nmrStar' also compress {  } into {}
        ## Wim 05/03/2003
        
        if nmrView_type:
            text = nmrView_compress( text )
        
        ## TITLE
        match_data_tag = re.search(r'\s*data_(\S+)\s+', text, 0 )
        if not match_data_tag:
            print "Warning: Found no 'data_title' string in file's text."
            print "Warning: Your file is not valid NMR-STAR - to attempt reading" 
            print "Warning: this file a data_title tag was added automatically."
            text = "data_autoTitle\n" + text
            match_data_tag = re.search(r'\s*data_(\S+)\s+', text, 0 )

        self.title = match_data_tag.group(1)
        pos = match_data_tag.end()

        ## Four quick searches for possible continuations
        next_sf_begin   = None      # SAVE FRAME BEGIN
        next_sf_end     = None      # SAVE FRAME END
        next_free_tt    = None      # FREE TAGTABLE
        next_loop_tt    = None      # LOOP TAGTABLE
        sf_open         = None      # When a saveframe is open
        text_length     = len(text)

        ## Only break when parsed to the eof
        while pos < text_length:
            if self.verbosity >= 9:
                print 'Parse text from position:%s : [%s]' % (
                    pos, text[pos:pos+10] )
            
            match_save_begin_nws = pattern_save_begin_nws.search(text,pos,pos+len('save_1'))
            if match_save_begin_nws:
                if match_save_begin_nws.start() == pos:
                    next_sf_begin = 1
            if not (next_sf_begin):
                match_save_end_nws = pattern_save_end_nws.search(text,pos,pos+len('save_ '))
                if match_save_end_nws:
                    if match_save_end_nws.start() == pos:
                        next_sf_end = 1
            if not (next_sf_begin or next_sf_end):
                match_tag_name_nws = pattern_tag_name_nws.search(text,pos,pos+len(' _X'))
                if match_tag_name_nws:
                    if match_tag_name_nws.start() == pos:
                        next_free_tt = 1
            if not (next_sf_begin or next_sf_end or next_free_tt):
                match_tagtable_loop_nws = pattern_tagtable_loop_nws.search(text,pos,pos+len('loop_ '))
                if match_tagtable_loop_nws:
                    if match_tagtable_loop_nws.start() == pos:
                        next_loop_tt = 1

            ## Just checking
            if not ( next_sf_begin or next_sf_end or next_free_tt or next_loop_tt ):
                print 'ERROR: No new item found in data_nodes_parse.'
                print 'Items looked for are a begin or end of a saveframe, or'
                print 'a begin of a tagtable(free or looped).'
                print 
                print "At text:"
                print text[pos:pos+70]
                print "Preceded by:"
                print text[pos-200:pos]
                return None
            
            ## SAVE FRAME BEGIN
            if next_sf_begin:
                if sf_open:
                    print "ERROR: Found the beginning of a saveframe but"
                    print "ERROR: saveframe before is still open(not closed;-)"
                    return None
                match_save_begin = pattern_save_begin.search( text, pos )
                if not match_save_begin:
                    print "ERROR: Code error (no second match on sf begin)";
                    return None
                if match_save_begin.start() != pos:
                    print "ERROR: Code error (wrong second match on sf begin)";
                    return None
                self.datanodes.append( SaveFrame(  tagtables    = [], # Need resetting
                                            verbosity    = self.verbosity ) )
                self.datanodes[-1].title = match_save_begin.group(1)
                sf_open         = 1
                next_sf_begin   = None
                pos             = match_save_begin.end()
                continue

            ## SAVE FRAME END
            if next_sf_end:
                if not sf_open:
                    print "ERROR: Found the end of a saveframe but"
                    print "ERROR: saveframe was not open"
                    return None
                match_save_end = pattern_save_end.search( text, pos )
                if not match_save_end:
                    print "ERROR: Code error (no second match on sf end)";
                    return None
                if match_save_end.start() != pos:
                    print "ERROR: Code error (wrong second match on sf end)";
                    return None
                sf_open     = None
                next_sf_end = None
                pos         = match_save_end.end()
                continue

            ## FREE or LOOP TAGTABLE
            if next_free_tt:
                free            = 1
                next_free_tt    = None
            else: # next_loop_tt must be true as this was checked before
                if not next_loop_tt:
                    print 'ERROR: code bug in File.parse()'
                    return None
                free            = None
                next_loop_tt    = None

                match_tagtable_loop = pattern_tagtable_loop.search( text, pos )
                if not match_tagtable_loop:
                    print 'ERROR: Code error, no second match on tagtable_loop'
                    return None
                if match_tagtable_loop.start() != pos:
                    print "ERROR: Code error (wrong second match on tagtable_loop)"
                    return None
                pos = match_tagtable_loop.end()

            if sf_open:
                dn = self.datanodes[-1].tagtables # Insert in last saveframes' tagtables
            else:
                dn = self.datanodes
                
            dn.append(      
                    TagTable(   free      = free,
                                tagnames  = [],
                                tagvalues = [],
                                verbosity = verbosity ) )
            tt = dn[-1] # Just to be verbose for the beloved reader
            pos = tt.parse( text=text, pos=pos )
            
            if pos ==  None:
                print "ERROR: In parsing tagtable"
                return None
            if self.verbosity >=9:                
                print 'Parsed tagtable up to pos: [%s]' % pos
            
        if self.verbosity >= 9:
            print 'Parsed: [%s] datanodes (top level count only)' % len( self.datanodes )
            
        if self.check_integrity( recursive = 0):
            print "ERROR: integrity not ok"
            return 1

        # Save some memory
        text = ''
        return 0



    """
    Writes the object to a STAR formatted file using
    the filename attribute.
    """
    def write (self):
        if not self.filename:
            print 'ERROR: no filename in STARFile with title:', self.title
            return 1
        open(self.filename,'w').write( self.star_text() )
        if self.verbosity > 1:
            print 'Written STAR file:', self.filename


    """
    Tries to reformat a file on disk with the filename given in the
    attribute of this object.
    Running Steve Madings (BMRB) formatNMRSTAR program if available    
    NOTE: this does NOT do anything with the datanodes of this object!
    """
    def formatNMRSTAR( self,
                    comment_file_str_dir    = '/bmrb/lib',
                    ):

        if self.verbosity >= 9:
            print "Attempting to reformat STAR file using external program if available"
        
        if os.name != 'posix':
            print "WARNING: No external program available on non-posix systems for reformatting STAR files"
            return 1

        ##  Try command and check for non-zero exit status
        ##  Note that these commands are only valid on Unix 
        ##  Standard error is thrown on the bit bucket.
        cmd = "%s < %s 2>/dev/null" % ('formatNMRSTAR', self.filename)
        pipe = os.popen( cmd )
        output = pipe.read()
        
        ##  The program exit status is available by the following construct
        ##  The status will be the exit number (in one of the bytes)
        ##  unless the program executed successfully in which case it will
        ##  be None.
        status = pipe.close()
        if self.verbosity >= 9:
            print "Got status:", status

        ## Success
        if ( status == None ):
            try:
                open(self.filename, 'w').write(output)
            except IOError:
                print 'ERROR: Could not open the file for writing', self.filename
                return 1            
            if self.verbosity >= 9:
                print "Reformatted STAR file:", self.filename
            return 0
        else:
            if self.verbosity :
                print "WARNING: Not pretty printing STAR file", self.filename
            return 1


################################################################################
if __name__ == '__main__':
    if 1:
        # 0 is only errors, 1 is warnings as well, 2 is normal and 9 is debug
        STAR.verbosity              = 2
        strf                        = File( verbosity=STAR.verbosity)
        
        if os.name == 'posix':
            pathname = '.'
        else:
            pathname = r'C:\Temp'

##        filename = 'local/test_nopound.txt'
##        filename = 'local/block_21921.txt'
        
        topDir = getTopDirectory()
        #filename = os.path.join(topDir, '../rawdata/8/1/1/info.general/bmr5106.str')
        filename = os.path.join(topDir, '../reference/nmrStar/1b4c.restraints')
        #filename = os.path.join(topDir, '../rawdata/2/2/1/info.general/bmr5307.str')

        strf.filename = filename

##        def myfunc():
        if strf.read():
            print "ERROR: In read. Exiting program"
            
            
        # strf.filename
        # strf.formatNMRSTAR: REFORMATS FILE! Don't use!!
        # strf.attrnames:     bulk output of everything in file
        # strf.flavor:        mmCIF, ... (in principle only nmrStar works?)
        # strf.title:         pdb code?
        #
        # strf.datanodes:     list of saveframes
        #               
        #               .title:       saveframe name
        #               .tagtables:   list o/ebi/msd/nmrqual/rawdata/8/1/1/info.general/bmr5106.stf tags in saveframe
        #                         .tagnames:  list of tagnames
        #                         .tagvalues: list of list of values for tagnames
        #                         .free:      if None, is a LOOP tagtable! If 1, is a list of simple tags
        # 
        
        # time.time()
        # strf.read()
        # time.time()

##        profile.run( 'myfunc()' )

##        print "Exiting program (test done)"
##        sys.exit(0)
                    

##            strf.flavor = 'mmCIF'
        #strf.filename = strf.filename + '_new1.str'
        #if strf.write():
        #    print "ERROR: In write. Exiting program"


