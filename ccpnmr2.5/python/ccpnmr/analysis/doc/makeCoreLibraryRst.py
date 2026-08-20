"""
======================COPYRIGHT/LICENSE START==========================

makeCoreLibraryRst.py: Part of the CcpNmr Analysis program

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
import os

from memops.universal.Io import getPythonDirectory

RST_PATH = 'ccpnmr/analysis/doc/source/core/'
CODE_PATH = 'ccpnmr/analysis/core/'
TOP_DIR = getPythonDirectory()

def makeCoreDoc():

  docDir = os.path.join(TOP_DIR, RST_PATH)
  if not os.path.exists(docDir):
    os.makedirs(docDir)

  files = [f for f in os.listdir(docDir) if f.endswith('.rst')]
  for file in files:
    fileName = os.path.join(docDir, file)
    os.unlink(fileName)

  codeDir = os.path.join(TOP_DIR, CODE_PATH)
  #files = [f for f in os.listdir(codeDir) if f.endswith('.py') and not f.startswith('__')]
  #names = [f.split('.')[0] for f in files]
  names = ['AssignmentBasic',
           'ChemicalShiftBasic',
           'ConstraintBasic',
           'ConstraintBasic',
           'CouplingBasic',
           'DataAnalysisBasic',
           'ExperimentBasic',
           'MarkBasic',
           'MergeObjects',
           'MoleculeBasic',
           'PeakBasic',
           'PrintBasic',
           'QualityControlBasic',
           'StructureBasic',
           'ValidationBasic',
           'WindowBasic',
           'Util']

  for name in names:
    fileName = os.path.join(docDir, name + '.rst')
    fileObj = open(fileName, 'w')

    n = len(name)
    head = '=' * n

    fileObj.write(head+'\n'+name+'\n'+head+'\n\n')
    fileObj.write('.. automodule:: ccpnmr.analysis.core.'+name+'\n   :members:\n\n')
    fileObj.close()

if __name__ == '__main__':

  makeCoreDoc()
