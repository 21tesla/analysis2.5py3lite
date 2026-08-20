from optparse import OptionParser, Option, OptionValueError

"""

Here defining additional types for option parser, so can handle for example a date format.
Should be generally useful! 

"""

def check_date(option, opt, value):

  import re
  dateCheck = re.compile("^\d{4}-\d{2}-\d{2}$")
  
  try:
    if dateCheck.search(value):
      return value
  except:
    raise OptionValueError("Option %s: invalid date value: %r. Please use YYYY-MM-DD" % (opt, value))
  
  raise OptionValueError("Option %s: invalid date format: %r. Please use YYYY-MM-DD" % (opt, value))

def check_namingSystemName(option, opt, value):
  
  from ccpnmr.format.general.Constants import namingSystemNames
  
  try:
    if value in namingSystemNames[1:]:
      return value
  except:
    raise OptionValueError("Option %s: invalid naming system value: %r" % (opt, value))
  
  print "Valid naming system names: %s" % (', '.join(namingSystemNames[1:]))
  raise OptionValueError("Option %s: invalid naming system value: %r" % (opt, value))

def check_yesNo(option, opt, value):
  
  import re
  yesCheck = re.compile("^y(es)?$")
  noCheck = re.compile("^n(o)?$")
  
  try:
    if yesCheck.search(value):
      return True
    elif noCheck.search(value):
      return False
  except:
    raise OptionValueError("Option %s: invalid yes/no value: %r" % (opt, value))
  
  raise OptionValueError("Option %s: invalid yes/no format: %r" % (opt, value))

def check_dataType(option, opt, value):

  # Dummy class
  return None

from copy import copy

class CustomOption(Option):
  TYPES = Option.TYPES + ("date","dataType","namingSystem",'yesNo')
  TYPE_CHECKER = copy(Option.TYPE_CHECKER)
  TYPE_CHECKER["date"] = check_date
  TYPE_CHECKER["dataType"] = check_dataType
  TYPE_CHECKER['yesNo'] = check_yesNo
  TYPE_CHECKER["namingSystem"] = check_namingSystemName


"""
 
Generic script handler class

"""

class ScriptHandler:

  ScriptHandlerError = StandardError
  
  """
  In subclasses, below variables should be redefined!
  """
  
  programDescription = "Set this in subclass"
  programVersion = '1.0'
  programUsage = 'Usage: %prog <options>'

  # Note: can put [%default] in helpText for listing default values!
  optionTuple = (('longArgumentName',   # Long argument name, e.g. inputFile
                  's',                  # Short argument name, e.g. i
                  'ArgName',            # Argument name, e.g. FILE or DATE. Set to None for Boolean arguments.',
                  'isMandatory',        # False/True,
                  'dataType',           # Standard: string (default for None), int, long, choice, float, complex
                                        # User defined: date, dataType, namingSystemName
                  'helpText',           # Help text, e.g.'Input file'
                  ),
                 )

  def __init__(self,scriptInit=True,**keywds):
  
    """

    This initialisation handles the system arguments passed to the script.
    Specifics are set in subclasses of this class (e.g. annotation/convertCoords.py)

    """
    
    if scriptInit: 
    
      # This creates the info for getopt based on the defaultOptionDict and optionDict definitions.
      self.parseOptions()
      
      self.handleOptions()
      
      # This can be used by subclasses for handling sub-class specific information    
      self.initialise(keywds)
  
    else:
      
      self.setDefaultOptions()
  
  def parseOptions(self):
  
    # Start the option parser with a customised set of data types
    self.optionParser = OptionParser(description=self.programDescription,
                                     version='%%prog version %s' % self.programVersion,
                                     option_class=CustomOption,
                                     usage=self.programUsage)
  
    # Loop over the optionTuple (should be defined in subclass!)
    for optionInfo in self.optionTuple:
    
      (longOpt,shortOpt,argName,isMandatory,dataType,helpText) = optionInfo
      
      # Set the keywords for the option to be added
      optKeywords = {'dest': longOpt, 'help': helpText}
      
      if type(argName) != type(''):
        # Expecting boolean
        optKeywords['default'] = argName
        if argName == True:
          optKeywords['action']  = 'store_false'
        elif not argName: # Will work for both None and False..
          optKeywords['action']  = 'store_true'
      else:
        optKeywords['metavar'] = argName
        
      if dataType:
        optKeywords['type'] = dataType
        
      # Add the option
      self.optionParser.add_option("-%s" % shortOpt, "--%s" % longOpt, **optKeywords)

      # Use nargs=2 for multiple arguments with option!
    
    # All set, get the information from sys.argv    
    (self.options, self.args) = self.optionParser.parse_args()    
      
  def handleOptions(self):
  
    # Here the option values are set as self. variables

    for optionInfo in self.optionTuple:
    
      (longOpt,shortOpt,argName,isMandatory,dataType,helpText) = optionInfo
      
      argValue = getattr(self.options,longOpt)
      
      if isMandatory and argValue is None:     
        print "A mandatory option is missing\n"
        self.optionParser.print_help()
        exit(-1)
                    
      setattr(self,longOpt,argValue)
      
      #print "Set %s to %s." % (longOpt, str(argValue))

  def setDefaultOptions(self):
  
    #
    # Make sure self.attrs are set even if script init ignored
    #

    for optionInfo in self.optionTuple:
    
      (longOpt,shortOpt,argName,isMandatory,dataType,helpText) = optionInfo
      
      if type(argName) != type(''):
        if argName == True:
          argValue = True
        elif not argName:
          argValue = False
      else:
        argValue = argName
        
      setattr(self,longOpt,argValue)

  def initialise(self,keywds):
  
    # Subclass specific initialisation options
    pass
  
  #
  # Functions that can be used in self.initialise()
  #
  
  def getExecutionInfo(self):
     
    import time, getpass
    
    self.currentTimeStamp = time.strftime("%Y-%m-%d.%H:%M")
    self.userName = getpass.getuser()

  def getInputFile(self):
    
    #
    # Define input file from argument list
    #

    if len(self.args) != 1:
      self.optionParser.print_help()
      raise self.ScriptHandlerError("This script needs an inputFile argument! See usage above.")
    
    self.inputFile = self.args[0]
  
  def getInputOutputFile(self):
    
    #
    # Define input and output file from argument list
    #

    if len(self.args) != 2:
      self.optionParser.print_help()
      raise self.ScriptHandlerError("This script needs inputFile and outputFile arguments! See usage above.")
    
    self.inputFile = self.args[0]
    self.outputFile = self.args[1]
  
    if self.inputFile == self.outputFile:
      raise self.ScriptHandlerError("Input and output file have to be different.")

    
if __name__ == '__main__':

  scriptHandler = ScriptHandler()
