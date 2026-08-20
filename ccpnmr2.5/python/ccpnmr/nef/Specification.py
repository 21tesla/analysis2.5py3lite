"""Code for handling NEF specification and metadata

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2019"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2019-12-13 19:01:25 +0000 (Fri, December 13, 2019) $"
__version__ = "$Revision: 3.0.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import sys

from . import GenericStarParser
from . import StarIo


INFOPREFIX = 'INFO: '


# TODO, This is a DRAFT only - not used and not currently functional.
# May be upgraded later, for specification-aware NEF I/O


def getCcpnSpecification(filePath):
    """Get NEF specification summary with ccpn-specific additions"""

    converter = CifDicConverter(open(filePath).read(),
                                additionalBlocks=('ccpn_additions',))
    #
    return converter.convertToNef()


class CifDicConverter(object):
    """Converts mmcif .dic file, with program-specific additions datablocks
    into a single NEF data structure, containing:

    1) a nef_specification saveframe, containing a dictionary_history loop and a
    item_type_list loop

    2) A saveframe for each saveframe_dategory in the specification. Each saveframe contains

        items:
        _nef_saveframe.sf_framecode
        _nef_saveframe.sf_category
        _nef_saveframe.is_mandatory
        _nef_saveframe.description
        _nef_saveframe.example

        A table for contained loops:

        loop_
           _nef_loop.category
           _nef_loop.is_mandatory
           _nef_loop.description
           _nef_loop.example

          And a table for contained items and loop columns:

        loop_
           _nef_item.name
           _nef_item.loop_category
           _nef_item.type_code
           _nef_item.is_mandatory
           _nef_item.is_key
           _nef_item.example_1
           _nef_item.example_2
           _nef_item.description

           The loop_category defines which loop the item belongs to
           (if empty it belongs directly inside the saveframe)
        """

    def __init__(self, inputText, skipExamples=True, additionalBlocks=()):
        self.specification = GenericStarParser.parse(inputText)
        self.additionalBlocks = additionalBlocks
        self.keyTags = {}
        self.result = None
        self.skipExamples = skipExamples
        self._category2SaveFrame = {}

    def _logging(self, *args):
        """Log messages as required
        """
        try:
            print('{}{}'.format(INFOPREFIX, ' '.join([str(arg) for arg in args])))
        except Exception as es:
            print('{}>>> Error during logging: {}'.format(INFOPREFIX, str(es)))

    def convertToNef(self):
        """Convert RCSB .cif file into a nef specification summary file
        """
        # NOTE - this assumes a single datablock.

        # set up
        rcsbDataBlock = list(self.specification.values())[0]
        result = self.result = StarIo.NmrDataBlock(name='specification')

        # make specific content saveframes
        self.extractGeneralDataFrame(rcsbDataBlock)

        dataBlocks = [rcsbDataBlock] + list(self.specification.get(tag)
                                            for tag in self.additionalBlocks)
        if None in dataBlocks:
            ii = dataBlocks.index(None)
            raise ValueError("Specification file has no data block matching %s"
                             % repr(self.additionalBlocks[ii - 1]))

        for dataBlock in dataBlocks:
            toSaveFrames, toLoops, toItems = extractByCategories(dataBlock)
            self._logging('%s SAVEFRAMES' % len(toSaveFrames))
            for xx in toSaveFrames:
                self.extractSaveFrameDescription(xx)
            self._logging('saveframes are', list(self.result.keys()))
            self._logging('%s LOOPS' % len(toLoops))
            for xx in toLoops:
                self.extractLoopDescription(xx)
            self._logging('%s ITEMS' % len(toItems))
            for xx in toItems:
                self.extractItemDescription(xx)

        # error check:
        if self.keyTags:
            self._logging("Error. unused keys:")
            for tt in self.keyTags:
                self._logging(tt)
            self._logging()

        #
        return result

    def extractGeneralDataFrame(self, rcsbDataBlock):
        """Extract general data saveframe
        """
        saveFrame = self.result.newSaveFrame('nef_specification', category='nef_specification')
        saveFrame.addItem('version', rcsbDataBlock.get('_dictionary.version'))

        # VersionHistory loop
        transferLoop(rcsbDataBlock, saveFrame, ('_dictionary_history.version',
                                                '_dictionary_history.update',
                                                '_dictionary_history.revision'))
        # ItemType loop
        typeLoop = transferLoop(rcsbDataBlock, saveFrame, ('_item_type_list.code',
                                                           '_item_type_list.primitive_code',
                                                           '_item_type_list.construct',
                                                           '_item_type_list.detail'))
        # Strip off spaces
        for row in typeLoop.data:
            detail = row['detail'].strip()
            if detail:
                ll = detail.splitlines()
                row['detail'] = '\n'.join(x.strip() for x in ll) + '\n'
        #
        return saveFrame

    def extractSaveFrameDescription(self, inputSaveFrame):
        """Extract saveframe description
        """
        expectedTags = ('_category.description', '_category.id', '_category.mandatory_code',
                        '_category_group.id', '_category_key.name', '_category_examples.case',
                        '_category_examples.detail',)

        metaCategory = 'nef_saveframe'
        category = inputSaveFrame['_category.id']
        name = '%s_%s' % (metaCategory, category)
        saveFrame = self.result.newSaveFrame(name, category=metaCategory)
        self._category2SaveFrame[category] = saveFrame
        saveFrame.addItem('is_mandatory',
                          inputSaveFrame.get('_category.mandatory_code') == 'yes')

        saveFrame.addItem('description', inputSaveFrame.get('_category.description'))

        data = inputSaveFrame.multiColumnValues(('_category_examples.detail',
                                                 '_category_examples.case',))

        examples = [x['_category_examples.case'] for x in data]
        if len(examples) == 1:
            if self.skipExamples:
                example = 'omitted'
            else:
                example = examples[0]
            saveFrame.addItem('example', example)
        elif examples:
            self._logging("Multiple examples for %s" % name)
            for dd in data:
                self._logging(dd['_category_examples.detail'], dd['_category_examples.case'])

        # Get keytags for later use
        keyNamesData = inputSaveFrame.multiColumnValues(('_category_key.name',))
        for dd in keyNamesData:
            tt = list(dd.values())[0].split('.', 1)
            self.keyTags[(tt[0][1:], tt[1])] = name

        # Check for untreated tags
        for tag in inputSaveFrame:
            if tag not in expectedTags:
                self._logging("Unexpected item in %s:" % inputSaveFrame['_category.id'], tag,
                              inputSaveFrame.get(tag))

    def extractLoopDescription(self, inputSaveFrame):
        """Extract loop description
        """
        expectedTags = ('_category.description', '_category.id', '_category.parent_category_id',
                        '_category.mandatory_code', '_category_group.id', '_category_key.name',
                        '_category_examples.case', '_category_examples.detail',)
        name = inputSaveFrame['_category.id']
        parentCategory = inputSaveFrame.get('_category.parent_category_id')
        if parentCategory is None:
            self._logging("loop is missing _category.parent_category_id:", name)
        else:
            parent = self._category2SaveFrame.get(parentCategory)
            if parent is None:
                self._logging("loop is missing parent saveFrame:", name, parentCategory,
                              list(self._category2SaveFrame.keys()))
            else:
                self._category2SaveFrame[name] = parent

                # get example
                data = inputSaveFrame.multiColumnValues(('_category_examples.detail',
                                                         '_category_examples.case',))
                examples = [x['_category_examples.case'] for x in data]
                if examples:
                    example = 'omitted'
                    if not self.skipExamples:
                        if len(examples) == 1:
                            example = examples[0]
                        else:
                            self._logging("Multiple examples for %s" % name)
                            for dd in data:
                                self._logging(dd['_category_examples.detail'], dd['_category_examples.case'])
                else:
                    example = None

                # make loop
                loop = parent.get('nef_loop')
                if loop is None:
                    loop = parent.newLoop('nef_loop', ('category', 'is_mandatory', 'description', 'example'))
                loop.newRow(dict(
                        category=name,
                        is_mandatory=inputSaveFrame.get('_category.mandatory_code') == 'yes',
                        description=inputSaveFrame.get('_category.description'),
                        example=example)
                        )

        # Get keytags for later use
        keyNamesData = inputSaveFrame.multiColumnValues(('_category_key.name',))
        for dd in keyNamesData:
            tt = list(dd.values())[0].split('.', 1)
            if len(tt) != 2:
                self._logging("key lacks internal '.'", parentCategory, name, tt)
            self.keyTags[(tt[0][1:], tt[1])] = name

        # Check for untreated tags
        for tag in inputSaveFrame:
            if tag not in expectedTags:
                self._logging("Unexpected item in %s:" % inputSaveFrame['_category.id'], tag,
                              inputSaveFrame.get(tag))

    def extractItemDescription(self, inputSaveFrame):
        """Extract item description
        """
        expectedTags = ('_item_description.description', '_item.name', '_item.mandatory_code',
                        '_item.category_id', '_item_type.code', '_item_examples.case',
                        '_item_examples.detail',)

        # get data
        name = inputSaveFrame.get('_item.name').split('.', 1)[1]
        category = inputSaveFrame.get('_item.category_id')
        isKey = (category, name) in self.keyTags
        if isKey:
            del self.keyTags[(category, name)]
        saveFrame = self._category2SaveFrame.get(category)
        if saveFrame is None:
            raise ValueError("SAvFrame named %s not found in list: %s"
                             % (category, list(self._category2SaveFrame.keys())))
        if saveFrame.name == 'nef_saveframe_' + category:
            # item lives in a saveframe, not a loop
            category = None
        isMandatory = inputSaveFrame.get('_item.mandatory_code') == 'yes'
        description = inputSaveFrame.get('_item_description.description')
        typeCode = inputSaveFrame.get('_item_type.code')

        # get examples
        data = inputSaveFrame.multiColumnValues(('_item_examples.detail',
                                                 '_item_examples.case',))
        examples = [x['_item_examples.case'] for x in data]
        if len(examples) > 2:
            self._logging("More than two examples for %s" % name)
            # for dd in data:
            #     self._logging(dd['_item_examples.detail'], dd['_item_examples.case'])
        while len(examples) < 2:
            examples.append(None)

        # Add iten to loop, making it if necessary
        specificationLoop = saveFrame.get('nef_item')
        if specificationLoop is None:
            specificationLoop = saveFrame.newLoop('nef_item', ('name', 'loop_category', 'type_code',
                                                               'is_mandatory', 'is_key', 'example_1', 'example_2',
                                                               'description'))

        specificationLoop.newRow(dict(name=name, loop_category=category,
                                      type_code=typeCode, is_mandatory=isMandatory, is_key=isKey,
                                      example_1=examples[0], example_2=examples[1],
                                      description=description))


def extractByCategories(rcsbDataBlock):
    """Get saveFrames describing SaveFrames, Loops, and items, respectively
    """
    toSaveFrames = []
    toLoops = []
    toItems = []
    for tag, saveFrame in rcsbDataBlock.items():

        if tag.startswith('save_'):
            category = saveFrame.get('_category.id')
            if category is None:
                toItems.append(saveFrame)
            else:
                saveFrameCodeTag = '_%s.sf_framecode' % category
                keyNamesData = saveFrame.multiColumnValues(('_category_key.name',))
                if any(x for x in keyNamesData if list(x.values())[0] == saveFrameCodeTag):
                    toSaveFrames.append(saveFrame)
                else:
                    toLoops.append(saveFrame)
    #
    return (toSaveFrames, toLoops, toItems)


def transferLoop(genericContainer, saveFrame, inputTags):
    """Transfer _category.tag_x, ... to loop named category with tags tag_x etc.
    """
    set1 = set()
    columns = []
    for tag in inputTags:
        tt = tag.split('.')
        if len(tt) == 2 and tt[0][0] == '_':
            columns.append(tt[1])
            set1.add(tt[0][1:])
        else:
            raise ValueError("Tag %s is not of form _xyz.abc")

    if len(set1) == 1:
        category = set1.pop()
        data = genericContainer.multiColumnValues(inputTags)
        if data:
            loop = saveFrame.newLoop(category, columns=columns)
            for row in data:
                loop.newRow(list(row.get(tag) for tag in inputTags))

            return loop
    else:
        raise ValueError("tags have more than on prefix: %s" % sorted((set1)))
    #
    return None


if __name__ == '__main__':

    args = sys.argv
    if len(args) < 2:
        print("Error, input file name is mandatory")
    else:
        infile = sys.argv[1]

        with open(infile) as fp:
            data = fp.read()
        converter = CifDicConverter(data)
        converter.convertToNef()
        print(converter.result.toString())
