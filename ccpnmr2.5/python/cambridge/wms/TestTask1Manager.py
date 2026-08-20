"""

 Example slave process. For now this has been set up as a separate
 thread although this is to make things clearer form an admin point
 of view as much as anything else. 

 The details of the architecture need to be settled. In the long run
 this needs to tie in with ExtendNMR.

 The main question is whether we need this specific manager layer and
 how we link this with the task code layer underneath.

"""


__author__ = "John Ionides"
__copyright__ = "Copyright (C) 2009 CCPN"

__revision__ = "$Id: TestTask1Manager.py,v 1.2 2010-09-09 16:04:50 wb104 Exp $"
__version__ = "0.1"



import os, sys, time

from threading import Thread

from SharedBeanService_services import *
from WSString import *

from Test1 import Test1

# TODO

# we currently have this set up so that the task MUST be
# registered in the WMS. Is this really correct? Or do we
# want to be able to run tasks from a local file location?
# Ignore this for now but consider later


class TestTask1Manager(Thread):

    # TODO

    # we have to decide how to pass the details about the task across. It might be best
    # just to pass the task ID and the repository (null for now) across and then pick
    # up the rest from a query.

    def __init__(self, id, status, repository=None):

        print 'creating new TestTask1Manager ', id, ', ', status

        # we need to know the repository because client side
        # tasks could be running off any repository (or none
        # at all!). For now this will be none and the repository
        # is implicit within the imported ws stubs
        self._repository = repository
        
        # task serial number
        self._taskId = id

        # get the info
        loc = SharedBeanServiceLocator()
        port = loc.getSharedBean()
        request1 = getFieldsWithMethod();
        request1._arg0 = 'org.pimslims.applet.server.TaskBean';
        request1._arg1 = 'getTest1Fields'
        h2 = {'serial': self._taskId }
        wsstr_in = WSString(h2)
        request1._arg2 = wsstr_in.str
        response1 = port.getFieldsWithMethod(request1)
        wsstr_out = WSString(response1._return)
        ss = wsstr_out.getStructHM()
        
        print 'got task details ', ss

        # input project Name. Needs to be
        if ss.has_key('inputVersion'):
            self._inputProjectName = ss['inputVersion']['project']['name'].__str__()
            self._inputVersionTag = ss['inputVersion']['versionTag'].__str__()
        else:
            # need to throw an exception. We don't have any input information!
            self.abort()

        # task serial number
        self._taskStatus = status 

        #checking
        self.isComplete = False

        #checking
        self.isRunning = True

        # for admin purposes. Check again precisly what we want here
        self._startTime = time.time()

        ##############################################
        #
        # EXPECT A LOT OF TASK-SPECIFIC PARAMETERS ETC

        # simple mimic for now
        self.status = None

        Thread.__init__(self)


    def run(self):

        # check the current status. Could use a dedicated tag if need be
        if self.status == 'FAILED':
            return

        # for now we don't need to do anything. We are just going
        # to set flags to show that a run has finished. This method
        # will call the actual exec (be it through WS or directly)
        # and will essentially tie in to the extend NMR code base

        print 'running ', self


        # important to set the status to RUNNING so we don't try
        # to start another process with the same input params!
        
        h2 = {'serial': self._taskId,
              'status': 'RUNNING' }

        loc = SharedBeanServiceLocator()
        port = loc.getSharedBean()

        request1 = record();

        request1._arg0 = 'org.pimslims.applet.server.TaskBean';
        request1._arg1 = 'updateStatus';

        wsstr_in = WSString(h2)

        print 'trying to update with ', h2, ', ', wsstr_in.str
        request1._arg2 = wsstr_in.str

        port.record(request1)


        # first get the parameters
        # NO: should have got these during instantiation
        # should be all ready to export the file

        # This needs to be set off the data model (or could be prompted)
        expDir = '/home/jionides/work/CCPN/test_EXP_data'
        print 'run: STARTING EXPORT'
        self._repository.export_project(self._inputProjectName, self._inputVersionTag,expDir)
        print 'run: FINISHED EXPORT'

        # For now run directly off here
        print 'run: STARTING TASK'
        projDir = expDir + '/' + self._inputProjectName
        task = Test1(projDir)
        task.process()
        print 'run: FINISHED TASK'

        self.isComplete = True
        self.isRunning = False
        print 'run: DONE'


    def abort(self):

        # need to set a validate flag that the run() method can pick up
        self.isComplete = True
        self.isRunning = False
        self.status = 'FAILED'

    """

    To begin with just update the status for the Task. Next step will be
    to make some trivial changes to the file and then check this in as a
    revised version. At that stage we have the architecture for the entire
    client-side daemon system

    """


    def recordResults(self):


        # re-import the results
        newVersionTag = self._repository.increment_version(self._inputVersionTag)
        expDir = '/home/jionides/work/CCPN/test_EXP_data/' +  self._inputProjectName
        print 'run: STARTING IMPORT'
        self._repository.import_project(self._inputProjectName, newVersionTag,expDir)
        print 'run: FINISHED IMPORT'

        # then set the output version
        h2 = {'serial': self._taskId,
              'generatedVersion': {'versionTag': newVersionTag,
                                   'project': { 'name' : self._inputProjectName }}}

        loc = SharedBeanServiceLocator()
        port = loc.getSharedBean()

        request1 = record();

        request1._arg0 = 'org.pimslims.applet.server.TaskBean';
        request1._arg1 = 'updateOutput';

        wsstr_in = WSString(h2)

        print 'trying to update output with ', h2, ', ', wsstr_in.str
        request1._arg2 = wsstr_in.str

        port.record(request1)


        # finally set status to OK
        h2 = {'serial': self._taskId,
              'status': 'OK' }

        loc = SharedBeanServiceLocator()
        port = loc.getSharedBean()

        request2 = record();
        
        request2._arg0 = 'org.pimslims.applet.server.TaskBean';
        request2._arg1 = 'updateStatus';

        wsstr_in = WSString(h2)

        print 'trying to update status with ', h2, ', ', wsstr_in.str
        request2._arg2 = wsstr_in.str

        port.record(request2)
        

    """

    Method that should validate the parameters. This means reading the input files
    and checking that everything is OK to go

    For now just return true

    """

    def exit(self):

        super(Thread)


    def validate(self):

        return True


    def getTimeStarted(self):

        return self._startTime


    """

    At some point we will want to build in a timeout. Will need this

    """


    def getTimeRunning(self):

        return (time.time() - self._startTime)
