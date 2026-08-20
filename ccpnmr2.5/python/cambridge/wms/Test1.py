
# The first bit (resetting the repository paths) should really
# be part of the export_project call, so that this sets up the data
# in a state that is know to be good for reading.

from memops.general.Io import loadProject


class Test1:

    def __init__(self, projectDir):

        self._projectDir = projectDir

    def process(self):
        
        rootProject = loadProject(self._projectDir)

        rootProject.saveModified()



if __name__ == '__main__':

    test = Test1('/home/jionides/work/CCPN/test_EXP_data/1ai0')

    test.process()

      


