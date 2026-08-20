
# class that represents a user in the context of the GUI. The idea is
# that the GUI will have a lot of user parameters stored locally
# (there are security issues here). On login these are read and this
# sets up the GUI. We then have to decide whether it is meaningful to
# run the GUI without a user (i.e. should this be a mandatory first
# step?)



class User:

    def __init__(self, name, statusFile=None):

        self.name = name

        # do we need this? no, probably not. let us
        # assume that all the status info is stored
        # together in an encrypted file
        self.statusFile = statusFile
