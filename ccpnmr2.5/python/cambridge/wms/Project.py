

from cambridge.wms.ProjectVersion import ProjectVersion


class Project:

    def __init__(self, rep, name):

        self.rep = rep
        self.name = name
        self.children = []

        # need to have status for
        # READ
        # READ_WRITE
        self.status=None

    def addVersion(self, ver):

        self.children.append(ver)

    def newVersion(self, verNum):

        ver = ProjectVersion(self, verNum)
        self.children.append(ver)

        return ver

