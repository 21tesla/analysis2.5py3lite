"""
Classes for dealing with STAR syntax
"""

__author__    = "$Author: jurgenfd $"
___revision__ = "$Revision: 1.6.4.1 $"
___date__     = "$Date: 2011-04-06 09:00:38 $"

"""
$Log: not supported by cvs2svn $
Revision 1.6  2007/07/09 12:22:27  wimvranken
Merge with branch 4.

Revision 1.5.4.1.4.1  2007/06/13 08:35:52  wimvranken
Removed copyright - is not CCPN

Revision 1.5.4.1  2006/07/21 17:29:01  wimvranken
Now only contain full path imports

Revision 1.5  2005/06/27 16:43:42  wb104
Updated licenses.

Revision 1.4  2004/12/15 17:57:20  tjs23
TJS: Updated licenses.

Revision 1.3  2003/08/08 12:46:08  wfv20
Added comment writing support

Revision 1.2  2003/07/14 09:34:32  wfv20
Modified so list references are not carried through class initialization

Revision 1.1  2003/07/01 12:56:18  wfv20
Jurgen Doreleijers modified Python nmrStar reader. Bug fixes and additions (notably for reading nmrView star files) have been made.

Revision 1.1.1.1  2001/11/02 20:16:40  jurgen
Initial package capable of read/write access to STAR files without nested loops

"""
## BMRB modules
#import ccp.format.nmrStar.bmrb as STAR
##from File           import *
from ccp.format.nmrStar.bmrb.TagTable       import * #@UnusedWildImport
##from Text           import *
from ccp.format.nmrStar.bmrb.Utils          import * #@UnusedWildImport

"""
Saveframe class
"""
class SaveFrame (Lister):
    def __init__( self,
                  title     = 'general_sf_title',
                  tagtables = None,
                  text      = '',
                  verbosity = 2,
                  comment = ''):
        self.title      = title
        
        # Modified tagtables initialization so list references
        # are not carried through (Wim 14/07/2002)
        self.tagtables = tagtables
        
        if self.tagtables == None:
          self.tagtables  = []
          
        self.text       = text
        self.verbosity  = verbosity
        self.comment = comment          # Comment attribute added to node (Wim 2003/08/05)
        
    "Returns the STAR text representation"
    def star_text (self,
                   flavor = 'NMR-STAR'
                   ):
        str = "\n"
        str = str + 'save_%s\n' % self.title
        
        for tagtable in self.tagtables:
            str = str + tagtable.star_text( flavor=flavor )
            
        str = str + '\nsave_\n'
        return str
    
    "Simple checks on integrity"
    def check_integrity( self,  recursive = 1  ):
        if recursive:
            for tagtable in self.tagtables:
                if tagtable.check_integrity():
                    print "ERROR: integrity check failed for tagtable"
                    return 1
        if self.verbosity >= 9:
            print 'Checked integrity of SaveFrame(%2s tagtables, recurs.=%s)  : OK [%s]' % (
                len(self.tagtables), recursive, self.title )

        
###############################################################################
if __name__ == '__main__':
    if 1:
        sf = SaveFrame( verbosity=9 )
        sf.check_integrity()
        print "sf:", sf.star_text()
        print "sf:", sf
        print "FINISHED"
