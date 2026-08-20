

# class that represents a view of a repository from
# the perspective of a given user (self.user).

# FIXME JMCI

# consider the model carefully. we could have two systems,
# one which represents a repository and one which represents
# a repository log in. There might be benefits to doing things
# this way.

# UPDATE JMCI

# The way things are currently going, it looks as though we
# are going to have a VERY thin client -- essentially only
# handling the display information. This means that we probably
# only want data in the simplest possible form here.

# For now, we can get away with the natural keys of the current
# information. All data should be retrieved from the DB every
# time anyway -- otherwise we start to get into the whole realm
# of synchronisation.

# What we may well need are other classes with static methods
# (essentially Beans) to abstract out common functionality.

import os

from memops.general.Io import loadProject

from cambridge.wms.Project import Project
from cambridge.wms.Task import Task

# required for WS layer
from SharedBeanService_services import *
from WSString import *


class Repository:

    def __init__(self, repList, name, connect, user, type=None):

        self.repList = repList
        self.name = name
        self.connect = connect
        self.type = type
        
        # it is not clear whether this belongs here or not.
        # yes, I think that we do somehow need to associate all
        # the connection and password information for the repository
        # connection with the repository somehow. I might be sensible
        # not to load the passwd info:: need to consider the security
        # implications
        self.user = user

        self.projects=[]
        self.basket=[]
        self.tasks=[]

        # If we are going to have a very thin client then we need
        # something more like this. These should be read in from a
        # local config file and saved on exit
        
        self.currentVersionTag = None
        self.currentProjectName = None
        self.currentTask = None
        self.currentBasket=None

        # for display purposes
        self.fullDisplay = True

        # for knitting together reference data. This should be read in
        # from a local data file and saved on exit        
        self.refDataLocations = []


    def import_project(self, proj_name, proj_ver, proj_dir, size = 100000):

        print 'in import project with params ', proj_name, ', ', proj_dir

        # These should be stored within this layer. They
        # should not come from an import but should be
        # built during the __init__() call
        loc = SharedBeanServiceLocator()
        port = loc.getSharedBean()

        # first check to see whether the project exists

        hm1 = {'name': proj_name}
        wsstr1 = WSString(hm1)

        req1 = isExist()
        req1._arg0 = 'org.pimslims.applet.server.ProjectBean'
        req1._arg1 = 'isExist'
        req1._arg2 = wsstr1.str
        projectExists = reponse1 = port.isExist(req1)

        print 'in import project with initial response', reponse1._return

        # FIXME
        # need lots more logic here once we have decided on a standard        
        if proj_ver == None:
            proj_ver = '1.1'

        # now actually store the file
        req2 = storeDocWithParams();

        # it seems as though the best practice is to just find all the files
        # in the directory and copy them across. We are not bothering with
        # the reference data for now; that will be handled later.
        file_list = []
        for root, dirs, files in os.walk(proj_dir):
            for filename in files:
                file_list.append(os.path.join(root, filename))


        print 'IMPORT PROJECT: FILES: ', file_list

        # going to upload in chunks

        for file in file_list:

            # need to handle possible tailing '/'
            # should really do this be regex and catch any odd path mismatches
            local_path = file[len(proj_dir)+1:]
            print 'copying ', file, ' to local path ', local_path

            f = open(file, 'r')
            f_str = f.read()
    
            size = 100000
            offset = 0
            send = True
            
            while (send):
                start = (offset * size)
                end = ((offset + 1) * size ) 
                if (end >= len(f_str)):
                    end = len(f_str)
                    send = False

                h2 = {'fileName': local_path,
                      'size': size,
                      'offset': offset,
                      'encoded': False,
                      'check': 'none_yet',
                      'projectName': proj_name,
                      'versionTag': proj_ver,
                      'segment': 'jionides',             
                      'raw': False}
                
                wsstr2 = WSString(h2)

                print 'sending hash ', h2
                print 'sending hash text ', wsstr2.str

                req2._arg0 = 'org.pimslims.applet.server.CcpnFileBean'
                req2._arg1 = 'appendChunk'
                req2._arg2 = f_str[start:end]
                req2._arg3 = wsstr2.str
                resp2 = port.storeDocWithParams(req2)

                offset += 1
        
            f.close()


            # register the project and projectVersion. This should
            # probably be done within the bean. should maybe look at
            # whether there is a better way of organising the
            # interface (e.g. have a loadProject bean that does all
            # the hard work and then registers the project)

            # need to know whether the project should already exist or
            # not !!
    
            hm3 = {'name': proj_name,
                   'location': '/home/jionides/work/CCPN/test_WMS_archive/jionides'}
            wsstr3 = WSString(hm3)

            print "going to create new record ", hm3

            req3 = record()
            req3._arg0 = 'org.pimslims.applet.server.ProjectBean'
            # need to get this name from either popup or directory name
            req3._arg1 = 'record'
            req3._arg2 = wsstr3.str
            resp3 = port.record(req3)
      
            
            hm4 = {'parent': {'name' : proj_name},
                   'versionTag': proj_ver,
                   'status': 'AVAILABLE',
                   'creationTime' : '9 FEB 2008:19:09'}
            wsstr4 = WSString(hm4)

            print "going to create new pv record ", hm4

            req4 = record()
            req4._arg0 = 'org.pimslims.applet.server.ProjectVersionBean'
            # need to get this name from either popup or directory name
            req4._arg1 = 'record'
            req4._arg2 = wsstr4.str
            resp4 = port.record(req4)        


    def export_project(self, projectName, versionTag, expDir):

        print 'EXPORT_PROJECT ', projectName, ', ' ,versionTag

        loc = SharedBeanServiceLocator()
        port = loc.getSharedBean()
        
        # first get a list of all the files required
        
        req1 = getList()
        
        hm1 = {'projectName': projectName,
               'versionTag': versionTag,
               'segment': 'jionides' }
        wsstr1 = WSString(hm1)
        
        # these are actually static
        req1._arg0 = 'org.pimslims.applet.server.CcpnFileBean';
        req1._arg1 = 'getList';
        req1._arg2 = wsstr1.str
        
        # get the response
        resp1 = port.getList(req1)
        
        # this is a hack at present. It needs to be written properly
        wsstr = WSString(resp1._return)
        ss = wsstr.getStruct()
        
        print 'got array ', ss
        
        for strg in ss:
            
            print 'trying to handle file ', strg
            
            hm2 = {'projectName': projectName,
                   'versionTag': versionTag,
                   'segment': 'jionides',
                   'fileName' : '/' + strg.__str__(),
                   'random' : 'test',
                   'raw' : False}
                
            print 'trying to handle file 2 ', strg
            print 'trying to handle file 2 ', hm2
                
            req2 = getDocWithParams()
                
            wsstr2 = WSString(hm2)
            
            # these are actually static
            req2._arg0 = 'org.pimslims.applet.server.CcpnFileBean';
            req2._arg1 = 'getDocWithParams'
            req2._arg2 = wsstr2.str
            
            # get the response
            resp2 = port.getDocWithParams(req2)
            
            local_file_path = expDir + '/' + projectName + '/' + strg
            
            idx = local_file_path.rfind('/')
            local_dir = local_file_path[:idx]
            print 'trying to create ', local_dir

            if (os.path.exists(local_dir) == 0):
                os.makedirs(local_dir)
            
            print 'writing to local area; ', local_file_path

            f = open(local_file_path, 'w')
            f.write(resp2._return)

        # now set the local repository paths correctly by
        # loading in and saving. We may need some much more
        # clever stuff added here to handle errors and to
        # check that all the reference data is in place

        # need to think about the requirements for exporting
        # projects multiple times. Do we need to provide the
        # export project name from the export dialogue box?
        projectDirectory = expDir + '/' + projectName
        r = loadProject(projectDirectory)
        r.saveModified()

    # FIXME

    # This is a quick hack. Needs to be done properly once we
    # decide on a proper versioning system
 
    def increment_version(self, versionTag):

        r = versionTag[:-1]
        l = int(versionTag[-1:]) + 1
        return r + l.__str__()
           
