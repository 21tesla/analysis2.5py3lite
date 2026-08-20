
class Task:

    def __init__(self, rep, name, status=None):

        self.rep = rep
        self.name = name

        # need status

        # OK
        # PENDING
        # FAILED
        # RUNNING
        
        self.status = status

        # FIXME JMCI

        # need to look at the cardinality here

        self.input = input
        self.output= output

    
