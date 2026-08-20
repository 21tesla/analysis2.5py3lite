
"""
======================COPYRIGHT/LICENSE START==========================

makeAllDoc.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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

#>>> from testCode.ccpnmr.analysis.testPopup import openPopups
#>>> openPopups(top)

def openAllAnalysisPopupsMacro(argServer):

  from testCode.ccpnmr.analysis.testPopup import openPopups
  openPopups(argServer.parent)
  
if __name__ == '__main__':

  try:
    import sphinx
  except ImportError:
    raise Exception('You must have the Sphinx module installed to run this script')

  import subprocess
  from ccpnmr.analysis.doc.makeAnalysisDocRst import makeAnalysisDoc
  from ccpnmr.analysis.doc.makeCoreLibraryRst import makeCoreDoc
  from ccpnmr.analysis.AnalysisGui import main
  from testCode.ccpnmr.analysis.testPopup import openPopups


  # Open Analysis

  top = main()

  # New, blank project to get things active

  top.newProject('DocExtract')

  # Open all popups

  openPopups(top, False)

  # Extract popup rst

  makeAnalysisDoc(top)

  # Extract code library rst

  makeCoreDoc()

  # run sphinx

  #subprocess.Popen(['make', 'clean']) # Wipes out CVS dirs
  subprocess.Popen(['make', 'html'])

  # Transfer to mammoth

  #PUBLISH_LOCATION = 'ccpn@mammoth:/data/ccpn/www/htdocs/documentation/analysis/'
  #subprocess.Popen(['scp', '-r', './build/html/*', PUBLISH_LOCATION])

  cmd1 = "scp -r ./build/html/* ccpn@mammoth:/data/ccpn/www/htdocs/documentation/analysis/"
  cmd2 = "scp -r ./build/html/* ccpn@mammoth:/data/ccpn/www/htdocs/documentation/analysisEdge/"
  
  print "To update in-program documentation website for stable issue:\n%s\n" % cmd1
  print "To update in-program documentation website for leading-edge issue:\n%s\n" % cmd2
  
