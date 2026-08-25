"""
Module Documentation here
"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (https://www.ccpn.ac.uk) 2014 - 2022"
__credits__ = ("Ed Brooksbank, Joanna Fox, Victoria A Higman, Luca Mureddu, Eliza Płoskoń",
               "Timothy J Ragan, Brian O Smith, Gary S Thompson & Geerten W Vuister")
__licence__ = ("CCPN licence. See https://ccpn.ac.uk/software/licensing/")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Geerten Vuister $"
__dateModified__ = "$dateModified: 2022-02-07 17:13:54 +0000 (Mon, February 07, 2022) $"
__version__ = "$Revision: 3.1.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Ed Brooksbank $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import sys
from functools import wraps

from .GenericStarParser import StarSyntaxError

NEF_STANDARD = 'standard'
NEF_SILENT = 'silent'
NEF_STRICT = 'strict'

NEFVALID = 0
NEFERROR_GENERICGETTABLEERROR = -1
NEFERROR_TABLEDOESNOTEXIST = -2
NEFERROR_SAVEFRAMEDOESNOTEXIST = -3
NEFERROR_ERRORLOADINGFILE = -4
NEFERROR_ERRORSAVINGFILE = -5
NEFERROR_BADTOSTRING = -6
NEFERROR_BADFROMSTRING = -7
NEFERROR_BADMULTICOLUMNVALUES = -8
NEFERROR_BADTABLENAMES = -9
NEFERROR_LISTTYPEERROR = -10
NEFERROR_BADLISTTYPE = -11
NEFERROR_BADFUNCTION = -12
NEFERROR_BADCATEGORIES = -13
NEFERROR_BADADDSAVEFRAME = -14
NEFERROR_READATTRIBUTENAMES = -15
NEFERROR_READATTRIBUTE = -16
NEFERROR_BADKEYS = -17


class ErrorLog():
    """
    A class to facilitate Logging of errors to stderr.

    example:
      errorLogging.logError('Error: %s' % errorMessage)

    functions available are:

      logError              write message to the current output
      func = logger         return the current logger
      logger = func         set the current logger
      loggingMode = mode    set the logging mode where mode is:
                            'standard', 'silent', 'strict'

            current modes are:
                'standard'      errors are written to the stderr, no errors are raised
                'silent'        no errors are raised, no output to the stderr
                'strict'        errors are logged to stderr and errors are raised
                                to be handled by the calling functions

      mode = loggingMode     return the current mode.
    """
    _availableModes = (NEF_STANDARD, NEF_SILENT, NEF_STRICT)
    NEFERRORS = {NEFERROR_BADADDSAVEFRAME      : 'bad add saveFrame',
                 NEFERROR_BADCATEGORIES        : 'bad categories',
                 NEFERROR_BADFUNCTION          : '',
                 NEFERROR_BADLISTTYPE          : 'bad listType',
                 NEFERROR_BADTABLENAMES        : 'bad table names',
                 NEFERROR_LISTTYPEERROR        : 'list type error',
                 NEFERROR_BADMULTICOLUMNVALUES : 'bad multiColumnValues',
                 NEFERROR_BADFROMSTRING        : 'bad convert from string',
                 NEFERROR_BADTOSTRING          : 'bad convert to string',
                 NEFERROR_ERRORSAVINGFILE      : 'error saving file',
                 NEFERROR_ERRORLOADINGFILE     : 'error loading file',
                 NEFERROR_SAVEFRAMEDOESNOTEXIST: 'saveFrame does not exist',
                 NEFERROR_TABLEDOESNOTEXIST    : 'table does not exist',
                 NEFERROR_GENERICGETTABLEERROR : 'table error',
                 NEFERROR_READATTRIBUTENAMES   : 'error reading attribute names',
                 NEFERROR_READATTRIBUTE        : 'error reading attribute',
                 NEFERROR_BADKEYS              : 'error reading keys'}

    def __init__(self, logOutput=sys.stderr.write, loggingMode=NEF_STANDARD, errorCode=NEFVALID):
        """
        Initialise a new eror logging object
        :param logOutput:
        :param loggingMode:
        :param errorCode:
        """
        self._logOutput = logOutput
        self._loggingMode = loggingMode
        self._lastError = errorCode

    def __call__(self, func):
        """
        Method that is called when used as a decorator for functions
        Facilitates the logging or raising of errors depending on the error logging mode
        """

        @wraps(func)
        def errortesting(obj, *args, **kwargs):
            try:
                obj._logError(errorCode=NEFVALID)
                return func(obj, *args, **kwargs)

            except (RuntimeError, StarSyntaxError) as es:
                _type, _value, _traceback = sys.exc_info()
                _errString = '%s.%s: %s' % (obj.__class__.__name__,
                                            func.__name__,
                                            str(es)
                                            )
                obj._logError(errorCode=self._lastError, errorString=_errString)
                return None

        return errortesting

    @property
    def logger(self):
        """
        Return the current logging function
        :return func; defaults to sys.stderr.write
                profile of func:
                func(value:str)
        """
        # return the current logger
        return self._logOutput

    @logger.setter
    def logger(self, func):
        """
        Set the current logging function
        :param func:
               profile of func:
               func(value:str)
        """
        self._logOutput = func

    @property
    def loggingMode(self):
        """
        Return the current logging Mode
        current modes are:
            'standard'      errors are written to the stderr, no errors are raised
            'silent'        no errors are raised, no output to the stderr
            'strict'        errors are logged to stderr and errors are raised
                            to be handled by the calling functions
        :return string:
        """
        return self._loggingMode

    @loggingMode.setter
    def loggingMode(self, loggingMode):
        """
        Set the current logger to a valid mode, invalid modes are ignored
        :param loggingMode:
        """
        if loggingMode in self._availableModes:
            self._loggingMode = loggingMode

    @property
    def lastErrorString(self):
        """
        Return the error string of the last action
        :return string:
        """
        return self._lastErrorString

    @property
    def lastError(self):
        """
        Return the error code of the last action
        :return int:
        """
        return self._lastError

    def _clearError(self):
        """
        Clear the last error code and error string
        """
        self._lastError = NEFVALID
        self._lastErrorString = ''

    def _logError(self, errorCode=NEFVALID, errorString=''):
        """
        Log errorCode to the current logger
        :param errorCode:
        :param errorString:
        """
        self._clearError()
        if errorCode != NEFVALID:
            self._lastError = errorCode
            if not errorString:
                errorString = self.NEFERRORS[errorCode]
            self._lastErrorString = 'NefRunTimeError[%r]:\n%s' % \
                                    (self.NEFERRORS[errorCode], errorString)

            if self._loggingMode != NEF_SILENT:
               self._logOutput(self._lastErrorString)

            if self._loggingMode == NEF_STRICT:
                raise RuntimeError(str(self._lastErrorString))
