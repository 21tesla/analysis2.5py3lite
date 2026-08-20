
# do we actually need this on the client side? Can probably
# handle with simple structures

class Test1Task(Task):

    def __init__(self, rep, name, status=None):

        super()

        # does it make sense to encode the various i/o ports in
        # specific role associations to the task subtypes? what are
        # the disadvantages of this approach. Presumably we could
        # write a small utility that queries the model and produces
        # GUI widgets?
        

        self.role1InputFile = None
        self.role2InputFile = None

        self.role1OutputFile = None
        self.role2OutputFile = None
        

    
