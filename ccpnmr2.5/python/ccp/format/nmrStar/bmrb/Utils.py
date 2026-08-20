"""
Just a few utilities that can be of more general use.
"""

__author__    = "$Author: wimvranken $"
___revision__ = "$Revision: 1.4 $"
___date__     = "$Date: 2007-07-09 12:22:27 $"

"""
$Log: not supported by cvs2svn $
Revision 1.3.8.1  2007/06/13 08:35:52  wimvranken
Removed copyright - is not CCPN

Revision 1.3  2005/06/27 16:43:42  wb104
Updated licenses.

Revision 1.2  2004/12/15 17:57:20  tjs23
TJS: Updated licenses.

Revision 1.1  2003/07/01 12:56:19  wfv20
Jurgen Doreleijers modified Python nmrStar reader. Bug fixes and additions (notably for reading nmrView star files) have been made.

Revision 1.1.1.1  2001/11/02 20:16:40  jurgen
Initial package capable of read/write access to STAR files without nested loops

"""

class Lister:
    """Example from 'Learning Python from O'Reilly publisher'"""
    def __repr__(self):
        return ("<Instance of %s, address %s:\n%s>" %
           (self.__class__.__name__, id(self), self.attrnames()))

    def attrnames(self):
        result=''
        keys = self.__dict__.keys()
        keys.sort()
        for attr in keys:
            if attr[:2] == "__":
                result = result + "\tname %s=<built-in>\n" % attr
            else:
                result = result + "\tname %s=%s\n" % (attr, self.__dict__[attr])
        return result        


"""
A fast transposing algorithm from the python mailing list
"""
def transpose ( matrix ):
    if len( matrix ) < 1:
        print 'ERROR: trying to transpose an empty matrix'
        return 1
    elif len( matrix ) == 1:
        if len(matrix[0]) == 0:
            print 'ERROR: trying to transpose an empty matrix, shape would be lost'
            print 'ERROR: [[]] would become []'
            return 1
        else:
            return map( lambda y : (y,), matrix[0] )
    else:
        return apply( map, [None,] + list(matrix) )

###############################################################################

if __name__ == '__main__':
    if 1:
        m = [ [1,2], [3,4] ]
        print m, 'is transposed:', transpose(m)
