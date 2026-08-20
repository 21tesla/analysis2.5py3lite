
"""
======================COPYRIGHT/LICENSE START==========================

AnalysisGui.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
import os
import sys
import traceback

import Tkinter

# IsWindowsOS code - (not imported yet due to dependencies
# 1. Win32 binary over-writes sys.path and needs re-setting
# 2. We don't want to add dlls to system dir so add CCPN dir to $PATH 
if sys.platform[:3].lower() == 'win':
  if sys.path[0].count('library.zip') == 1:
    py_dir = os.path.join( os.path.dirname( sys.path[0] ), 'python' )
    sys.path.append( py_dir )
    del py_dir
  os.environ['PATH'] += ';' + os.path.dirname( sys.path[0] )

try:
  from memops.universal.Io import normalisePath
except ImportError:
  #print
  #print os.environ['PYTHONPATH']
  #print

  print 'Error, cannot import core CCPN Python modules:'
  print 'Maybe your PYTHONPATH environment variable is not set or'
  print 'does not contain the current CCPN installation directory.'
  raise

from memops.general.Implementation import ApiError

from memops.gui.MessageReporter import showError

from memops.universal.Util import isWindowsOS, isMacOS

from ccp.gui.Io import loadProject

from ccpnmr.analysis.AnalysisPopup import AnalysisPopup

if isWindowsOS():
  # create interactive session when using MS Windows
  os.environ["PYTHONINSPECT"]="x"

top = None

import atexit


def main(projectDir=None, cache_size=64, glDirect=None):

  global top

  #print 'cache_size =', cache_size

  root = Tkinter.Tk()
  root.withdraw()
  # root.option_add("*Background", "grey90")
  root.option_add("*Font", "Helvetica -12")
  root.option_add("*Dialog.msg.wrapLength", '6i')

  top = AnalysisPopup(root, cache_size=cache_size, glDirect=glDirect)
  #top.option_add("*Cursor", "crosshair")
  #top.option_add("*Cursor", "watch")
  #top.configure(cursor="crosshair")

  def _checkPythonQuit():
    if top is not None and not hasattr(top, '_hasQuitNormally'):
      # # NOTE:ED - can remove the threading Thread.__delete keyError
      # #           but not sure want to remove threading module
      # #           can monkey patch Thread before module loaded?
      # #           error caused because trying to exit twice and top already deleted...
      # if 'threading' in sys.modules:
      #   del sys.modules['threading']
      top.quit()

  # add handler to trap quiting from the Python menu in MacOS - shouldn't occur in Windows and Linux
  atexit.register(_checkPythonQuit)

  project = None
  if projectDir:
    try:
      project = loadProject(top, path=projectDir)
    except ApiError, e:
      showError('Loading project', e.error_msg)
      raise

  # lift(), called from open(), is required for some reason in order to get
  # the top widget first on the toolbar and lift() is slow for some reason
  #top.open() # very slow for some reason
  top.update_idletasks() # much faster
  top.initProject(project)

  # return top

  if isWindowsOS():
    root.mainloop()

  return top

def usage():

  print 'Allowed arguments:'
  print '  [ -m memory_size_in_megabytes ] [ -glDirect gl_rendering_direct (0 or 1) ] [ project_directory ]'
  sys.exit()

def getOptArg(argv, flag, defaultValue, conversionFunc = None, validArg = ''):

  n = len(argv)
  k = [i for i in range(n) if argv[i] == flag]

  if len(k) > 1:
    print 'Multiple occurrences of flag "%s"' % flag
    usage()

  if k:
    k = k[0]
    if k == (n-1):
      print 'Flag "%s" requires argument' % flag
      usage()

    value = argv[k+1]
    if conversionFunc:
      try:
        value = conversionFunc(value)
      except Exception:
        if validArg:
          validArg = validArg + ' '
        print 'Flag "%s" requires valid %sargument' % (flag, validArg)
        usage()

    del argv[k:k+2]

  else:
    value = defaultValue

  return value

if (__name__ == '__main__'):

  # startup error messages courtesy of Gary Thompson, Leeds
  startupExecError = """
*************************************************************************
ERROR: an exception occurred while executing the analysis startup script:

Exception:
%s

The startup file is defined by the environment variable ANALYSIS_STARTUP
and has the following path:
%s

*************************************************************************
"""

  missingStartupFileError = """
*************************************************************************
WARNING: couldn't find the Analysis python startup file

Startup file path:
'%s',

Either:

   1. check that the environment variable ANALYSIS_STARTUP
      points to a valid readable python file or
   2. unset the environment variable ANALYSIS_STARTUP

Continuing...
*************************************************************************
"""

  startupFile = os.environ.get('ANALYSIS_STARTUP')
  if startupFile:
    if os.path.isfile(startupFile):
      try:
        execfile(startupFile)
      except Exception, e:
        print startupExecError % (traceback.format_exc(), startupFile)
        print e
        sys.exit()
    else:
      print missingStartupFileError % startupFile

  argv = sys.argv[:]

  max_size = getOptArg(argv, flag='-m', defaultValue=128, conversionFunc=int, validArg='integer')
  glDirect = getOptArg(argv, flag='-glDirect', defaultValue=None, conversionFunc=int, validArg='integer')

  n = len(argv)

  projectDir = None
  if (n > 2):
    if n > 3:
      s = 's'
    else:
      s = ''
    print 'Have extra arg%s: "%s"' % (s, ', '.join(argv[1:-1]))
    usage()
  elif (n == 2):
    projectDir = argv[1]
    if not os.path.isdir(projectDir):
      print 'Path "%s" does not exist' % projectDir
      usage()
    elif not os.path.isdir(projectDir):
      print 'Path "%s" is not a directory' % projectDir
      usage()

  try:
    main(projectDir, max_size, glDirect)
  except Exception as es:
    # ignore any errors on exit
    pass
  finally:
    # NOTE:ED - this exits correctly now
    if isWindowsOS():
      os._exit(0)