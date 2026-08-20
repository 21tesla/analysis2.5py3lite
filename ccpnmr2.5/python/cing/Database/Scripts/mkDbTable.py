#@PydevCodeAnalysisIgnore # pylint: disable-all
from cing import NTdb
from cing.Libs.NTutils import *

stream = open('dbTable.new', 'w')
NTdb.exportDef(stream=stream, convention='INTERNAL_1')
stream.close()
