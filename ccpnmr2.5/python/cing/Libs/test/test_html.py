"""
Unit test execute as:
python $CINGROOT/python/cing/PluginCode/test/test_html.py
"""
import unittest
import os
from unittest import TestCase

from cing import cingDirTmp
from cing.core.classes import Project
from cing.core.molecule import Ensemble, Molecule
from cing.Libs.html import HTML_TAG_PRE, HTML_TAG_PRE2, HTMLfile, MakeHtmlTable, removePreTagLines
from cing.Libs.NTutils import *  #@UnusedWildImport


class AllChecks(TestCase):

    def testRemovePreTagLines(self):
        spuriousSpaceMsg = 'something     with     many spaces'
        msg = '\n'.join([HTML_TAG_PRE, spuriousSpaceMsg, HTML_TAG_PRE2 ])
        self.assertNotEqual(msg, spuriousSpaceMsg)
        self.assertEqual(removePreTagLines(msg), spuriousSpaceMsg)

    def setupSimplestProject(self):
        cingDirTmpTest = os.path.join( cingDirTmp, 'test_html' )
        os.makedirs( cingDirTmpTest , exist_ok=True)
        os.chdir(cingDirTmpTest)
        entryId = 'test'
        project = Project(entryId)
        self.assertFalse(project.removeFromDisk())
        project = Project.open(entryId, status='new')
        molecule = Molecule(name='moleculeName')
        molecule.ensemble = Ensemble(molecule) # Needed for html.
        project.appendMolecule(molecule) # Needed for html.
        c = molecule.addChain('A')
        r1 = c.addResidue('ALA', 1, Nterminal=True)
        if r1:
            r1.addAllAtoms()

        molecule.updateAll()
        project.setupHtml() # Needed for creating the sub dirs.
        return project

    def testMakeHtmlTableWithJS(self):
#        CSS and Javascript is going to determine much of the formatting.
        project = self.setupSimplestProject()
        h = HTMLfile(project.htmlPath('test.html'), project)
        columnFormats = [ ('col1', {}),
                          ('col2', {}),
                          ('col3', {})
                                 ]

        t = MakeHtmlTable(h.main, showHeader=False, classId="testJsTable", id="testJsTableId", columnFormats=columnFormats,
                          bla="0")
        for row in t.rows(range(2)):
            rStr = str(row)
            t.nextColumn()
            t(None, rStr)
            t.nextColumn()
            t('a', rStr + "." + str(2), href='someRefHere')
            t.nextColumn()  # empty one
            t.nextColumn()
            t(None, rStr + "." + str(4))

        #end for
        h.main('h3', 'Text to start below the table.')
        h.render()


if __name__ == "__main__":
    cing.verbosity = verbosityDebug
#    nTdebug("Starting...")
    unittest.main()
