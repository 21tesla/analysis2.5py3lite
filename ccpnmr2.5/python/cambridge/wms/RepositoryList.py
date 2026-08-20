
"""

  This is effectively the object that represents the GUI level. If
  contains a set of repositories and any other information that
  need to be stored concerning the GUI layer should probably go
  in here. This includes the TaskManager information



"""



from cambridge.wms.Repository import Repository
from cambridge.wms.User import User


class RepositoryList:

    def __init__(self, user='guest'):

        self.repositories=[]
        self.user = User(user)
        self.currentRepository = None

        self.current_export_dir = None
        self.current_import_dir = None

        # FIXME JMCI

        # we also want some other daemons that sit and refresh screen periodically

    def addRepository(self, rep):

        self.repositories.append(rep)

            
    def newRepository(self, repName, repConnect, user, type=None):

        rep = Repository(self, repName, repConnect, user, type)
        self.repositories.append(rep)

        return rep

    def getRepositories(self):

        return self.repositories
            
