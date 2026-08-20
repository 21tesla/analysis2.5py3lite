import re

braceRe = re.compile("\{([^{}]*)\}")

ANSIG_INT = 0
ANSIG_FLOAT = 1

class AnsigParseException(Exception):
  pass

# keyword dicts: keyword --> (block start?, number of args (tuple if several), types)
def parseAnsigFile(fileName, callback, keywordDict,
                   directoryDict, canDefineDirectory=False):

  fp = open(fileName, 'rU')
  lines = fp.readlines()
  fp.close()

  content = []
  for (lineNumber, line) in enumerate(lines):
    n = line.find('!')
    if n >= 0:
      comment = line[n+1:].strip()
      text = line[:n].strip()
    else:
      comment = None
      text = line.strip()
    if text:
      content.append((lineNumber+1, line, text, comment))

  ncontent = len(content)
  blockStack = []

  endBlock = None
  remValues = None
  contentIndex = 0
  try:
    while contentIndex < ncontent:
      keywords = keywordDict.keys()

      if remValues:
        fields = remValues
      else:
        (lineNumber, line, text, comment) = content[contentIndex]
        fields = text.split()
      key = fields[0]
      values = fields[1:]

      # check if key is really just a directory definition
      if canDefineDirectory:
        obj = braceRe.match(key)
        if obj:
          # just a directory definition
          dirKey = obj.group(1)
          if len(values) != 1:
            raise AnsigParseException('dir_symbol "%s" has values = %s, must be of length 1' % (key, values))
          directoryDict[dirKey] = substituteDirectory(values, directoryDict)[0]
          contentIndex += 1
          continue

      if key in keywords:
        (blockStart, nreqdValues, types) = keywordDict[key]
      else:
        blockStart = False

      # check for start of block
      if blockStart:
        # start of block is special: assume know exactly number of required values
        # (because don't know what keywords will follow since not yet specified)
        while len(values) < nreqdValues:
          contentIndex += 1
          if contentIndex >= ncontent:
            raise AnsigParseException('key "%s" only has %d arguments, requires %d' % (len(values), nreqdValues))
          text = content[contentIndex][2]
          if comment is None:
            comment = content[contentIndex][3]
          fields = text.split()
          values.extend(fields)
        remValues = values[nreqdValues:]
        values = values[:nreqdValues]
        values = substituteDirectory(values, directoryDict)
        values = convertTypes(values, types)
        blockStack.append((keywordDict, endBlock))
        keywordDict = callback(key, values, comment)
        endBlock = 'end_' + key

      # check for end of block
      elif key == endBlock:
        remValues = values
        callback(key, None, comment)
        (keywordDict, endBlock) = blockStack.pop()

      # ought to be keyword
      elif key in keywords:
        if type(nreqdValues) == type(()) or nreqdValues == -1:
          # look for next keyword
          # assume everything between this keyword and next keyword is a value
          n = keywordIndex(values, keywords, endBlock)
          while n < 0:
            contentIndex += 1
            if contentIndex >= ncontent:
              break
            text = content[contentIndex][2]
            if comment is None:
              comment = content[contentIndex][3]
            fields = text.split()
            values.extend(fields)
            n = keywordIndex(values, keywords, endBlock)
          if n >= 0:
            remValues = values[n:]
            values = values[:n]
          else:
            remValues = []
          if nreqdValues != -1:
            if len(values) not in nreqdValues:
              raise AnsigParseException('key "%s" has %d arguments, requires one of %s' % (key, len(values), nreqdValues))
        else:
          while len(values) < nreqdValues:
            contentIndex += 1
            if contentIndex >= ncontent:
              raise AnsigParseException('key "%s" only has %d arguments, requires %d' % (key, len(values), nreqdValues))
            text = content[contentIndex][2]
            if comment is None:
              comment = content[contentIndex][3]
            fields = text.split()
            values.extend(fields)
          remValues = values[nreqdValues:]
          values = values[:nreqdValues]
        values = substituteDirectory(values, directoryDict)
        if nreqdValues == -1:
          types = len(values) * types
        values = convertTypes(values, types)
        callback(key, values, comment)

      # problem, assume it's an unknown keyword
      else:
        raise AnsigParseException('unknown keyword "%s"' % key)

      # increment contentIndex if need be
      if not remValues:
        contentIndex += 1

  except AnsigParseException, e:
    raise IOError('In file "%s", line number %d: %s:\n  %s' % (fileName, lineNumber, str(e), line))

def substituteDirectory(values, directoryDict):

  newvalues = []
  for value in values:
    obj = braceRe.match(value)
    if obj:
      key = obj.group(1)
      if not directoryDict.has_key(key):
        raise AnsigParseException('key "%s" missing from directory dict' % key)
      value = value.replace('{%s}' % key, directoryDict[key])
    newvalues.append(value)

  return newvalues

def convertTypes(values, types):

  newvalues = []
  for (n, value) in enumerate(values):
    typ = types[n]
    if typ == ANSIG_INT:
      try:
        value = int(value.strip())
      except ValueError, e:
        raise AnsigParseException('value "%s" not an integer' % value)
    elif typ == ANSIG_FLOAT:
      try:
        value = float(value.strip())
      except ValueError, e:
        raise AnsigParseException('value "%s" not a real' % value)
    newvalues.append(value)

  return newvalues

def keywordIndex(values, keywords, endBlock):

  for (n, value) in enumerate(values):
    if value == endBlock or value in keywords:
      return n

  return -1

if __name__ == '__main__':

  import sys

  if len(sys.argv) != 2:
    print 'need argument: ansig ctr file'
    sys.exit()

  fileName = sys.argv[1]

  keywords = ()
  blockKeywordDict = {'control': 0}
  directoryDict = {}

  def callback(key, values, comment):
    print 'callback: key = %s: values = %s: comment = %s' % (key, values, comment)
    if key == 'control':
      controlKeywords = ('help_dir', 'html_viewer', 'dictionary', 'gr_translate', 'sequence',
              'spectra', 'crosspeaks', 'backup', 'initialize', 'exit', 'binaryformat', 'font')
      controlBlockKeywordDict = {}
      return (controlKeywords, controlBlockKeywordDict)

  parseAnsigFile(fileName, callback, keywords, blockKeywordDict, directoryDict, canDefineDirectory=True)

