"""
Unit test execute as:
python $CINGROOT/python/cing/Scripts/cingProfile.py
"""
import profile
import pstats

import cing  # py3: `from cing.X import Y` no longer binds the top-level name

from cing.Libs.NTutils import *  #@UnusedWildImport


def run():
    print("hello")
    cing.verbosity = verbosityDebug
    nTdebug( "hello again" )

if True:
    # Commented out because profiling isn't part of unit testing.
    fn = 'fooprof'
    profile.runctx('run()', globals(), locals(), fn)  # py3: profile.run() runs in a bare namespace; runctx sees module globals
    p = pstats.Stats(fn)
    #p.sort_stats('time').print_stats(100)
    p.sort_stats('cumulative').print_stats(20)
else:
    run()
