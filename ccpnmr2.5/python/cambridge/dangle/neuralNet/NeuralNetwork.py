# Back-Propagation Neural Networks
# 
# Written in Python.  See http://www.python.org/
# Placed in the public domain.
# Neil Schemenauer <nas@arctrix.com>

import math
import random
import string
import time

from numpy import array

random.seed(time.time())

# calculate a random number where:  a <= rand < b
def rand(a, b):
    return (b-a)*random.random() + a

# Make a matrix (we could use NumPy to speed this up)
def makeMatrix(I, J, fill=0.0):
    m = []
    for i in range(I):
        m.append([fill]*J)
    return m

#def makeMatrix(I, J):
#    m = numpy.zeros([I,J])
#    return m

def triggerFuncTanh(x):
    return math.tanh(x)

# derivative of our triggerFunc function
def triggerGradTanh(y):
    return 1.0-y*y


# our self.weightsOutput function, tanh is a little nicer than the standard 1/(1+e^-x)
def triggerFunc(x):
    
    v =  0.8920620580763845 * math.exp(-2.5*x*x)
    return v
    #return math.tanh(x)

# derivative of our triggerFunc function
def triggerGrad(y):
    return -4.4603102903819227 * y * math.exp(-2.5*y*y)

class NN:
    def __init__(self, ni, nh, no, testFunction=None):
        # number of input, hidden, and output nodes
        self.nInput  = ni + 1 # +1 for bias node
        self.nHidden = nh
        self.nOutput = no

        self.testFunction = testFunction or self.test

        self.missingInput = {}

        # activations for nodes
        self.signalInput  = [1.0]*self.nInput
        self.signalHidden = [1.0]*self.nHidden
        self.signalOutput = [1.0]*self.nOutput
        
        # create weights
        self.weightsInput = makeMatrix(self.nInput, self.nHidden)
        self.weightsOutput = makeMatrix(self.nHidden, self.nOutput)
        
        # set them to random vaules
        for i in range(self.nInput):
            for j in range(self.nHidden):
                self.weightsInput[i][j] = rand(-1, 1)
                
        for j in range(self.nHidden):
            for k in range(self.nOutput):
                self.weightsOutput[j][k] = rand(-1, 1)

        # last change in weights for momentum   
        self.ci = makeMatrix(self.nInput, self.nHidden)
        self.co = makeMatrix(self.nHidden, self.nOutput)
        self.cm = [1.0]*self.nInput

        # Record best weights
        self.bestScore = None
        self.weightsBestIn = None
        self.weightsBestOut = None

    def update(self, inputs):
        wo = self.weightsOutput
        wi = self.weightsInput
        so = self.signalOutput
        si = self.signalInput
        sh = self.signalHidden
        nOut = range(self.nOutput)
        nHid = range(self.nHidden)
        nInp = range(self.nInput)
        
        self.missingInput = {}
        
        if len(inputs) != self.nInput-1:
            print len(inputs), self.nInput-1
            raise ValueError, 'wrong number of inputs'

        # input activations
        for i in range(self.nInput-1):
            #self.signalInput[i] = triggerFunc(inputs[i])
            v = inputs[i]
            if v is None:
              si[i] = 0.0
              self.missingInput[i] = True
            else:
              si[i] = inputs[i]

        # hidden activations
        for j in nHid:
            totalSignal = 0.0
            for i in nInp:
                totalSignal += si[i] * wi[i][j]
            sh[j] = triggerFuncTanh(totalSignal)

        # output activations
        for k in nOut:
            totalSignal = 0.0
            for j in nHid:
                totalSignal += sh[j] * wo[j][k]
            so[k] = triggerFuncTanh(totalSignal)

        return so[:]


    def backPropagate(self, targets, N, M):
        wo = self.weightsOutput
        wi = self.weightsInput
        co = self.co
        ci = self.ci
        cm = self.cm
        so = self.signalOutput
        si = self.signalInput
        sh = self.signalHidden
        nOut = range(self.nOutput)
        nHid = range(self.nHidden)
        nInp = range(self.nInput)
    
        if len(targets) != self.nOutput:
            raise ValueError, 'wrong number of target values'

        # calculate error terms for output
        output_deltas = [0.0] * self.nOutput
        for k in nOut:
            error = targets[k]-so[k]
            output_deltas[k] = triggerGradTanh(so[k]) * error

        # calculate error terms for hidden
        hidden_deltas = [0.0] * self.nHidden
        for j in nHid:
            error = 0.0
            for k in nOut:
                error += output_deltas[k]*wo[j][k]
            hidden_deltas[j] = triggerGradTanh(sh[j]) * error

        # update output weights
        for j in nHid:
            for k in nOut:
                change = output_deltas[k]*sh[j]
                wo[j][k] += N * change + M*co[j][k]
                co[j][k] = change
                #print N*change, M*self.co[j][k]

        # update input weights
        for i in nInp:
            for j in nHid:
                change = hidden_deltas[j]*si[i]
                wi[i][j] += N*change + M*ci[i][j]
                ci[i][j] = change


        # calculate error
        #error = 0.0
        #for k in range(len(targets)):
        #    error += 0.5*(targets[k]-so[k])**2
            
            
        #return error


    def test(self, nn, patterns, verbose=True):
    
        n = 0 
        m = len(patterns)
        err = 0.0
        for data, known in patterns:
          #print p[0], '->', self.update(p[0])
         
          predict =  nn.update(data)
          
          i = known.index(max(known))
          j = predict.index(max(predict))
         
          for k in range(len(known)):
            err += 0.5*(known[k]-predict[k])**2

          if i == j:
             n += 1

        pc = (100.0*n)/m
        
        if verbose:
          print "Success rate %d from %d : %.2f%% error %.2f" % (n,m,pc,err)
        return pc

    def weights(self):
        print 'Input weights:'
        for i in range(self.nInput):
            print self.weightsInput[i]
        print
        print 'Output weights:'
        for j in range(self.nHidden):
            print self.weightsOutput[j]

    def train(self, patterns, iterations=10, N=0.5, M=0.1):

      #self.monteCarloTrain(patterns)

      t0 = time.time()

      n = len(patterns)

      """
      """
      print "Phase 1"
      nSplit = 10
      for i in range(nSplit):
        print 'Subset: %d' % i
        random.shuffle(patterns)
        self.backPropagateTrain(patterns[:int(n/nSplit)], 7, 0.5,  0.2)
      
      print "Phase 2"
      self.backPropagateTrain(patterns, iterations, 0.5,  0.2)
      print "Phase 3"
      self.backPropagateTrain(patterns, iterations, 0.5,  0.1)
      print "Phase 4"
      self.backPropagateTrain(patterns, iterations, 0.5, 0.05)
      print "Phase 5"
      self.bestScore = None
      self.backPropagateTrain(patterns, iterations, 0.5, 0.025)

      self.weightsInput  = [x[:] for x in self.weightsBestIn]
      self.weightsOutput = [x[:] for x in self.weightsBestOut]

      print "Time taken:", time.time()-t0


    def backPropagateTrain(self, patterns, iterations=10, N=0.5, M=0.1):
      
      
      # N: learning rate
      # M: momentum factor
      for i in xrange(iterations):
        #error = 0.0
        
        #print i
        
        random.shuffle(patterns)

        for inputs, targets in patterns:
          self.update(inputs)
           #error += self.backPropagate(targets, N, M)
          self.backPropagate(targets, N, M)
          
        #if i % 1 == 0:
        score = self.testFunction(self, patterns)

        if self.bestScore is None:
          self.bestScore = score
         
        if score >= self.bestScore:
          self.bestScore = score
          self.weightsBestIn = [x[:] for x in self.weightsInput] 
          self.weightsBestOut = [x[:] for x in self.weightsOutput]
   
          #print 'Back Propagate Iteration %d error %-14f' % (i+1,error)


def demo():
    # Teach network XOR function
    pat = [
          [[0,0], [0]],
          [[0,1], [1]],
          [[1,0], [1]],
          [[1,1], [0]]
          ]

    # create a network with two input, two hidden, and one output nodes
    n = NN(2, 2, 1)
    # train it with some patterns
    n.train(pat, iterations=100, N=0.5, M=0.1)
    # test it
    n.test(n, pat)



if __name__ == '__main__':
    t1 = time.time()
    demo()
    print time.time() - t1
