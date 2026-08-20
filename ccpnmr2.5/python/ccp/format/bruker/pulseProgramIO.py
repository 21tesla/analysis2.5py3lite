from ccp.format.bruker.generalIO import BrukerGenericFile

#
# Need experiment to map to CCPN!
#
# TODO: Problem with putting info in pulse program file: what if it changes? Make checks to see if 'original' one
# or not... .
#
# ;CCPN HNcoCACB F1:3,F2:2,F3:1 for mapping to EXPERIMENT side (can also use indexes, probably better to use dims!).

class PulseProgramFile(BrukerGenericFile):

  def initialize(self):
  
    self.fileProgramName = None
    self.comments = []
    self.dimIncrements = {}
    self.increments = []
    self.fnModes = []
    self.programName = None
    
    self.phaseSensitive = []
    self.semiconstantTime = []
    self.constantTime = []
    
    self.isEchoAntiEcho = 0
    
    self.pathwayInfo = []
    
    self.pathwayPointer = '->'

  def read(self):

    #
    # TODO TODO HERE have to use dictionary information, if available, to select
    # the correct parameters. If they were found, use those in code below!!!
    #

    fin = open(self.name)

    #line = fin.readline()

    #pulseLine = ""
    avanceVersion = False

    for line in fin:
      
      #
      # Look for hash comments
      #
      # RHFogh 25/7/08 none are useful
      
      #hashSearch = self.patt[self.format + 'HashComment'].search(line)
      
      # RHFogh 25/7/08 not needed
      #if not pulseLine and hashSearch:
      #  pulseLine = hashSearch.group(1).strip()
        
      #
      # Semicolon comment lines
      #
      
      semiColonSearch = self.patt[self.format + 'SemiColonComment'].search(line)
 
      if semiColonSearch:
        lineString = semiColonSearch.group(1).strip()
        # First line should be pulse program name...
        if not self.fileProgramName:
          self.fileProgramName = lineString
          
          # RHFogh 25/7/08 start addition
          continue
        
        elif lineString .startswith('avance-version'):
          # needed for program name
          avanceVersion = True
          continue
        
        elif avanceVersion and not self.programName:
          if lineString:
            self.programName = lineString
          else:
            self.programName = 'generic'
          
        
        # RHFogh 25/7/08 end of addition
        
        
        else:
          if lineString.count(self.pathwayPointer):
            while (1):
              pathwayElementSearch = self.patt[self.format + 'PathwayElements'].search(lineString)
              
              if pathwayElementSearch:
                
                pathwayElementInfo = (pathwayElementSearch.group(1),pathwayElementSearch.group(2))
                lineString = lineString[pathwayElementSearch.end():]            
                self.pathwayInfo.append(pathwayElementInfo)
               
              else:
                break
                
          else:
            self.comments.append(lineString)
      
      #
      # Try to get program name (????)
      #
      # RHFogh 35/7/08 not correct. Done differently
      
      #if not self.programName and pulseLine and self.comments:
      #  pulseNameSearch = self.patt[self.format + 'PulseName'].search(pulseLine)
      #  if pulseNameSearch:
      #    self.programName = pulseNameSearch.group(1)
      #  else:
      #    self.programName = 'generic'
          
      #
      # Get fnMode info from comment lines
      #

      fnModeSearch = self.patt[self.format + 'FnMode'].search(line)
      if fnModeSearch:
          self.fnModes.append(fnModeSearch.group(1).strip())


        # TODO look for dictionary match to set parameters correctly!!
        # self.pulseProgramInfo = ...
        
      #
      # Dimension increment search - in newer files only!
      #
      
      dimIncrSearch = self.patt[self.format + "DimIncrement"].search(line)
      
      if dimIncrSearch:
        incNum = dimIncrSearch.group(1)
        dimIncNum = dimIncrSearch.group(2)
        factor = dimIncrSearch.group(3)
        
        self.dimIncrements[dimIncNum] = (incNum, factor)
        
        # NOTE: incNum = dimIncNum / factor!!!

      #
      # Increment/Decrement search here - could be obsolete. Also have to add other searches
      # based on dictionary if match found!!
      #

      incrSearch = self.patt[self.format + 'InOrDecrement'].search(line)

      # TODO LET 0,10 take preference over everything else?!?!
      # DO bit of after-processing in case there are non 10/0 ones!
      # ALSO: try to use central dictionary to get this information out...
      if incrSearch:
        incNum = incrSearch.group(2)
        self.increments.append(incNum)
        
      #
      # Phase sensitive comment info
      #
      
      if self.patt[self.format + 'phaseSensitive'].search(line):
        
        tValues = self.patt[self.format + 'tSearch'].findall(line)

        if tValues:
          self.phaseSensitive.extend(tValues)
      
      #
      # Constant time comment info
      #
      for tag in ('semiconstantTime', 'constantTime'):
        if self.patt[self.format + tag].search(line):
          tValues = self.patt[self.format + 'tSearch'].findall(line)
          if tValues:
            getattr(self,tag).extend(tValues)
          break

      #
      # See if acquired with Echo-AntiEcho
      # TODO: this is not foolproof, for newer pulse programs try to use name or 
      #       internal information!!
      #

      echoAntiEchoSearch = self.patt[self.format + 'EchoAntiEcho'].search(line)

      if echoAntiEchoSearch:
        self.isEchoAntiEcho += 1

      #line = fin.readline()
      
    return True
