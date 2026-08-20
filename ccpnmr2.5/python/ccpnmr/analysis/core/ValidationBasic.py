LICENSE = """
======================COPYRIGHT/LICENSE START==========================

ValidationBasic.py: Part of the CcpNmr Analysis program

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
ANALYSIS_RMSD_CONTEXT = 'RMSD'

RMSD_KEYWORDS = ['backbone','all','CA','CB','H','O',]

from ccpnmr.analysis.core.Util import getSoftware
from ccp.util import Validation

def getEnsembleValidationStore(ensemble, context, keywords=None,
                               definitions=None, synonyms=None,
                               software=None):
  """Get a CCPN object to store validation results for an ensemble
             in a given program context. Requires a list of keywords which will
             be used in this context. Allows optional lists of definitions and
             user-friendly synonyms for these keywords.
             Optional argument for passing software specification (otherwise
             defaults to current CcpNmr Analysis)
\n.. describe:: Input\n\nMolStructure.StructureEnsemble, Word, List of Words,
             List of Lines, List of Words, Method.Software  
\n.. describe:: Output\n\nValidation.ValidationStore
  """
  
  if not software:
    software = getSoftware(ensemble.root)
    
  getFunc = Validation.getEnsembleValidationStore
 
  validStore = getFunc(ensemble, context, keywords,
                       definitions, synonyms, software)
          
  return validStore




