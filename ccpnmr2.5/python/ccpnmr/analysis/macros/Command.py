
"""
======================COPYRIGHT/LICENSE START==========================

Command.py: Part of the CcpNmr Analysis program

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
from memops.general import Implementation

class Command:
 
  def __init__(self, argumentServer, command_name, module_name, func_name):
 
    self.name = command_name
    self.module_name = module_name
    self.func_name = func_name
    self.argumentServer = argumentServer

    try:
      self.module = __import__(module_name)
    except:
      raise Implementation.ApiError('could not import module "' + self.module_name + '"')
 
    try:
      self.func = getattr(self.module, self.func_name)
    except:
      raise Implementation.ApiError('could not find function "' + self.func_name + '" in module "' + self.module_name + '"')

  def run(self):
 
    self.func(argServer=self.argumentServer)
    #try:
    #  self.func(argServer=self.argumentServer)
    #except Implementation.ApiError:
    #  raise
    #except TypeError, e:
    #  raise Implementation.ApiError(str(e))
    #except:
    #  raise Implementation.ApiError('could not run function "' + self.func_name + '" in module "' + self.module_name + '"')
 
  def reload(self):
 
    try:
      reload(self.module)
    except:
      raise Implementation.ApiError('could not reload module "' + self.module_name + '"')
      return
 
    try:
      self.func = getattr(self.module, self.func_name)
    except:
      raise Implementation.ApiError('could not find function "' + self.func_name + '" in module "' + self.module_name + '"')

