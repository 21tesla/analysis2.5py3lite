import sys

def parseStarDict():

  global starSfNameDict, starSfTableDict

  from msd.nmrStar.IO import nmrStarDict

  starDict = nmrStarDict.sfDict

  starSfNameDict = {}
  starSfTableDict = {}

  for sfTitle in sorted(starDict.keys() ):
    sfName = starDict[sfTitle]['name']

    starSfNameDict[sfName] = sfTitle

    if 'tableNames' in starDict[sfTitle]:
      for tableName in starDict[sfTitle]['tableNames']:
        starSfTableDict[tableName] = sfTitle

def parseFile(fileName):

  global aditMandDict

  file = open(fileName)

  lines = file.readlines()

  aditMandDict = {}

  for line in lines[1:]:
    if line[0] == '#':
      continue

    row = line.strip('\n').split('\t')

    #if len(row) == 8:
    #  print 'ROW: [%s]' % row

    if len(row) > 5:
      (sfTable, tagName) = row[4].split('.')

      sfTable = sfTable[1:]

      if sfTable in starSfNameDict:
        sfName = starSfNameDict[sfTable]

        #print 'SF: [%s] [%s]' % (sfName, tagName)

        if sfName not in aditMandDict:
          aditMandDict[sfName] = {}

          aditMandDict[sfName]['tags'] = {}
          aditMandDict[sfName]['tagNames'] = set()

          aditMandDict[sfName]['tables'] = {}
          aditMandDict[sfName]['tableNames'] = set()

        if tagName not in aditMandDict[sfName]['tags']:
          aditMandDict[sfName]['tags'][tagName] = row[5:] + [row[3]]
          aditMandDict[sfName]['tagNames'].add(tagName)

      else:
        tableName = sfTable
        sfName = starSfTableDict[tableName]

        if sfName not in aditMandDict:
          aditMandDict[sfName] = {}
          aditMandDict[sfName]['tags'] = {}
          aditMandDict[sfName]['tagNames'] = set()
          aditMandDict[sfName]['tables'] = {}
          aditMandDict[sfName]['tableNames'] = set()

        if tableName not in aditMandDict[sfName]['tables']:
          aditMandDict[sfName]['tables'][tableName] = {}

          aditMandDict[sfName]['tables'][tableName]['tags'] = {}
          aditMandDict[sfName]['tables'][tableName]['tagNames'] = set()

          aditMandDict[sfName]['tableNames'].add(tableName)

        if tagName not in aditMandDict[sfName]['tables'][tableName]['tags']:
          aditMandDict[sfName]['tables'][tableName]['tags'][tagName] = row[5:] + [row[3]]
          aditMandDict[sfName]['tables'][tableName]['tagNames'].add(tagName)

        #print 'TABLE: [%s] [%s] [%s]' % (sfName, tableName, tagName)

  print 'aditMandDict = {\n'

  for sfName in aditMandDict.keys():

    if aditMandDict[sfName]['tagNames'] or aditMandDict[sfName]['tableNames']:
      print "  '%s': {\n" % sfName

      if aditMandDict[sfName]['tagNames']:
        print  "    'tags': {\n"

        for tagName in sorted(aditMandDict[sfName]['tagNames']):
          print "      '%s': %s," % (tagName, aditMandDict[sfName]['tags'][tagName])

        print '\n      },\n'

        print "    'tagNames': %s,\n" % sorted(list(aditMandDict[sfName]['tagNames']) )

      if aditMandDict[sfName]['tableNames']:
        print "    'tables': {"

        for tableName in aditMandDict[sfName]['tableNames']:
          print "\n      '%s': {\n" % tableName

          if aditMandDict[sfName]['tables'][tableName]['tagNames']:
            print "        'tags': {\n"

            for tagName in aditMandDict[sfName]['tables'][tableName]['tagNames']:
              print "          '%s': %s," % (tagName, aditMandDict[sfName]['tables'][tableName]['tags'][tagName])

            print '\n          },\n'

            print "        'tagNames': %s,\n" % sorted(list(aditMandDict[sfName]['tables'][tableName]['tagNames']) )

          print '        },'

        print '      },\n'

        print "    'tableNames': %s,\n" % sorted(list(aditMandDict[sfName]['tableNames']) )

      print '    },\n'

  #print aditMandDict

  print '  }'

if __name__ == '__main__':

  parseStarDict()

  argv = sys.argv
  argc = len(argv)

  if argc > 1:
    fileName = argv[1]

    parseFile(fileName)
