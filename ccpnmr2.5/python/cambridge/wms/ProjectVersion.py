

class ProjectVersion:

    def __init__(self, projRoot, verNum, parent=None, comment=None):

        self.proj = projRoot
        self.verNum = verNum
        self.children = []
        self.parent = parent
        self.comment = comment

        # need to have status for
        # AVAILABLE
        # CHECKED OUT
        # LOCAL CHANGES
        # UNAVAILABLE

        self.status = None

        # for task management

        self.createdBy = None
        self.usedFor = []

    def addChildVersion(self, verNum):

        newVer = ProjectVersion(self.proj, verNum, self)
        self.children.append(newVer)

        return newVer

    def getChildren(self):

        return self.children

