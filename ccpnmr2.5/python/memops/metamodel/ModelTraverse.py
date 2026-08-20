"""  Version for python version > 2.1

Operations on the model in memory, adapting the basic model for a 
particular language and implementation. 

Includes general model querying functions the depend on context.
Includes default version of standard information, such as operations data,
variable names, etc.

Does no output.
"""

from memops.metamodel import ModelTraverse_py_2_1
from memops.metamodel.ModelPortal import ModelPortal
from memops.metamodel import MetaModel
MemopsError = MetaModel.MemopsError

class ModelTraverse(ModelTraverse_py_2_1.ModelTraverse_py_2_1, object):

  def __init__(self):
    """Constructor - sets up for future processing
    """

    # class needed within ObjectDomain (version 2.1)  
    super(ModelTraverse, self).__init__()
    
    for tag in ModelTraverse_py_2_1.mandatoryAttributes:
      if not hasattr(self, tag):
        raise MemopsError(" ModelTraverse lacks mandatory %s attribute" % tag)

    # has to be done this way to allow for different initialisation orders
    if not hasattr(self, 'modelFlavours'):
      self.modelFlavours = {}
    
    # input check
    if not isinstance(self.modelPortal, ModelPortal):
      raise MemopsError("ModelTraverse input %s is not a ModelPortal"
                        % self.modelPortal)
    
    # link varNames for easier access, and check modelPortal has been processed.
    if hasattr(self, 'varNames'):
      # this must have been called from ModelAdapt
      pass
    elif hasattr(self.modelPortal, 'varNames'):
      self.varNames = self.modelPortal.varNames
      self.operationData = self.modelPortal.operationData
      
    
  ###########################################################################
  
  ###########################################################################
