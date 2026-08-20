"""

 Master daemon process for monitoring the state of tasks
 and farming them out as necessary. This will start as
 a very simple process and, if everything goes to plan,
 should end up as a very sophisticated piece of code

 This needs to run in the baskground and periodically monitor
 the task table.

 Note: we could run this at the server end. Indeed, this is a major
 design question that needs to be put to the design panel. The
 advantage of having it at the client end is that it makes it
 easier for client-side software to be plugged in.

 If we want to be clever, we should probably support both client
 and server side code and then there needs to be a flag in the
 Task table to indicate which side should pick up a given task.

 Also note that in this context we have two different types of
 daemon. The server-side daemon will run in the background and
 pick up tasks (for that repository only) ALL the time, whereas
 the client side daemon is essentially just a background process
 sitting behind the GUI

"""


#__author__ = "John Ionides"
#__copyright__ = "Copyright (C) 2009 CCPN"

#__revision__ = "$Id: TaskManager.py,v 1.2 2010-09-09 16:04:50 wb104 Exp $"
#__version__ = "0.1"


import os, sys, time

from threading import Thread

from cambridge.wms.TestTask1Manager import TestTask1Manager

from SharedBeanService_services import *
from WSString import *

"""

Not only does this have to run as a daemon, but it also has to
keep track of the slave processes

At the moment we are tracking just a single repository. For client
side tasks we need to track ALL repositories. This means that we
need to maintain a link between each task and the repository from
which it comes. At the moment we don't need to worry too much about
this because the repository is implicit with the pyton ws stubs. Once
we fix ZSI to allow direct calls then we will need to look at this a
lot more carefully

"""

class TaskManager(Thread):

    def __init__(self, parent):

        # owner of the daemon. This will allow us to track back
        # simply into the Wms structure
        self.parent = parent

        # read from config. We need to be a lot cleverer
        # about how we set these up. There needs to be a ranking
        # and all sort of things        
        self._slaveTaskTypes = ['TestTask1']

        # max number of slave threads at once
        self._MAX_SLAVE_TASKS = 6

        # frequency of scan (in seconds)
        self._SCAN_FREQ = 5

        # list of slave processes
        self._slaves = []

        # checking. Set to false to stop any processing
        self.active = True

        # checking. Set to false to automatically redraw TasksFrame
        # This should maybe be run through a different daemon but run
        # it through this for now
        self.redraw = True

        # checking. Set to false to let the Thread exit
        self.alive = True

        Thread.__init__(self)


    """

    Main loop. Check to see whether any existing processes have finished. If
    so, the results will be recorded and task removed from the slave list.

    Then new tasks will be added up to a specified number.

    """

    def run(self):

        i = 1

        while self.alive:

            #print 'looping in main loop: iteration ', i, ', ', self

            time.sleep(self._SCAN_FREQ)

            if self.active:
                # first check the existing processes to see how they are doing
                print 'checking slave processes', self._slaves
                self._check_slaves()

                print 'topping up slave array ', self._slaves
                # then use up any spare resources on new one
                if (self.active and len(self._slaves) < self._MAX_SLAVE_TASKS):
                    self._allocate_slaves()

            i += 1



    """

    Method to check the progress of existing slave processes. This
    will pick up any tasks that have either failed and/or completed
    and register them. We could also consider adding a timeout value
    for each task and picking up stalled tasks too

    FIXME JMCI

    The Task class should have some details block to put failure
    comments in

    """


    def _check_slaves(self):

        tmp = []

        for slave in self._slaves:

            print 'checking slave ', slave
            print 'DICT ', slave.__dict__
            

            if slave.isComplete:

                slave.recordResults()
                slave.exit()

            else:

                tmp.append(slave)

            # FIXME JMCI

            # add timeout block here too?
                    
        self._slaves = tmp

    """

    This is the brains of the operation. We need to identify the tasks that
    need running and then set up jobs. For this we first need to get the
    complete list of tasks outstanding, along with their status and task
    type.

    """
        
    def _allocate_slaves(self):

        print "allocating tasks ", time.time()

        # get the tasks

        # This should end up as a direct call based on a repository
        # connection. Note that we will have to cope with multiple
        # repositories. For this (python) code we are only interested
        # in the client-side daemon so this should be linked to a WmsGui
        # and thus have access to all the connections etc
        loc = SharedBeanServiceLocator()
        port = loc.getSharedBean()

        request1 = getList();

        # need to decide precisely how to handle these. Should we
        # be able to choose which status we want to pick up?

        a2 = ['PROVISIONAL','PENDING']
        h2 = {'status' : a2 }
        wsstr_in = WSString(h2)
        

        request1._arg0 = 'org.pimslims.applet.server.TaskBean';
        request1._arg1 = 'getListWithFields';
        request1._arg2 = wsstr_in.str

        response1 = port.getList(request1)
        
        wsstr1 = WSString(response1._return)
        ss1 = wsstr1.getStruct()

        print 'got raw array ', ss1
        print 'setting up jobs ', self._slaves, ', ', self._MAX_SLAVE_TASKS

        # FIXME
        # this is a HORRIBLE hack. We clearly need to track the repository for
        # each client side process
        repository = self.parent.repList.currentRepository
        i = 0
        while len(self._slaves) < self._MAX_SLAVE_TASKS and i < len(ss1):
            
            print 'calling ', ss1[i]['status'], ' on task ',  ss1[i]['serial'],'(', ss1[i]['type'], ')'
            print 'monitoring length ', len(self._slaves)

            # Constructor needs to depend on type
            slave = TestTask1Manager(int(ss1[i]['serial'].__str__()), ss1[i]['status'].__str__(), repository)
            print 'created new slave process ', slave
            slave.start()
            self._slaves.append(slave)

            print 'monitoring final length ', len(self._slaves), ', ', len(ss1)
            # finally, increment i 
            i += 1

            print 'monitoring final contents ', self._slaves


    
