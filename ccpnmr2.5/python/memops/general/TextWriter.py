""" Version for Python version > 2.1
"""

from memops.metamodel import MetaModel
MemopsError = MetaModel.MemopsError
from memops.metamodel import ImpConstants
import memops.universal.Io as uniIo

from memops.general import TextWriter_py_2_1

settings = TextWriter_py_2_1.settings

class TextWriter(TextWriter_py_2_1.TextWriter_py_2_1, object):
  """
  """

  def __init__(self):
    """ parameters are mandatoryInitParams, optionalInitParams
    and the special optional
    """   
    
    for tag in TextWriter_py_2_1.mandatoryAttributes:
      if not hasattr(self, tag):
        raise MemopsError(" TextWriter lacks mandatory %s attribute" % tag)
    
    super(TextWriter, self).__init__()
      
    # special parameters: optional with default values
    if self.rootFileName is None:
      self.rootFileName = ImpConstants.rootPackageDirName
      
    if self.rootDirName is None:
      self.rootDirName = uniIo.getTopDirectory()
      

    self.fp = None
    self.fileName = ''
    self.indent = 0
    self.indents = []
    self.errorMsg = ''

    self.previousLineEmpty = False # used so that do not print out two '\n' in a row
    
  ###########################################################################

  ###########################################################################
