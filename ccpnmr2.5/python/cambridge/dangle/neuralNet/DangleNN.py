from random import random, randint

from math import atan2, tanh, cos, sin, sqrt

from NeuralNetwork import NN

aminoAcidLetters = ['I','L','V','F','W',
                    'Y','C','A','G','P',
                    'T','S','M','H','Q',
                    'N','E','D','K','R']

secStructTypes = ['H','E','C']

pi = 3.14159265358979
degrees45 = pi/4.0
degrees30 = pi/6.0
degrees20 = pi/9.0
degrees10 = pi/18.0

def readDatabase(fileName):

  data = []
  fileObj = open(fileName, 'r')
  line  = fileObj.readline()

  j = 0
  while line:
    line.strip()
     
    if line[0] == '!':
      line = fileObj.readline()
      continue
    
    datum = line.split()
    
    for i in range(1,len(datum)):
      if datum[i] == 'None':
        datum[i] = None
      elif datum[i] in secStructTypes:
        pass
      else:
        datum[i] = float(datum[i])
    
    data.append(datum)
     
    line = fileObj.readline()
    j += 1
  
  #print "Read %d lines" % j
    
  return data


def makeTrainingDataAngles(data):

  aaDict = {}
  i = 0
  for aa in aminoAcidLetters:
    aaDict[aa] = i
    i += 1

  missing = 0

  blankVec = [0.0] * 20


  trainingData = []

  for datum in data:
    seq  = datum[0]
    
    vec = []
    for i in range(5):
      aaVec = blankVec[:]
      aaVec[aaDict[seq[i]]] = 0.1
      vec.extend(aaVec)
            
    # Ha,Ca,Cb,C,N
    shifts = []
    for i in range(1,26):
      
      v = datum[i]
    
    
      if v is None:
        #shifts.append(0.0)
        shifts.extend([0.0, -0.2])
        
      else:
        #shifts.append(v/25.0)
        shifts.extend([v/25.0, 0.2])

    phi, psi = datum[26:28]
    secStruc = datum[28]
    
    secStructs = [0.0] * len(secStructTypes)
    
    index = secStructTypes.index(secStruc)
    secStructs[index] = 0.1
    
    phi *= pi/180.0
    psi *= pi/180.0
    
    angles = [cos(phi), sin(phi), cos(psi), sin(psi)]
    
    if None not in shifts:
      vec.extend(shifts)
  
      trainingData.append([vec, secStructs])
    else:
      missing += 1

  print "Excluded inputs:", missing
  return trainingData

def makeTrainingDataMissingShift(data):

  aaDict = {}
  i = 0
  for aa in aminoAcidLetters:
    aaDict[aa] = i
    i += 1

  missing = 0

  blankVec = [0.0] * 20


  trainingData = []

  for datum in data:
    seq  = datum[0]
    phi, psi = datum[26:28]
    phi *= pi/180.0
    psi *= pi/180.0
    
    angles = [cos(phi), sin(phi), cos(psi), sin(psi)]

    secStruc = datum[28]
    secStructs = [0.0] * len(secStructTypes)
    
    index = secStructTypes.index(secStruc)
    secStructs[index] = 0.1
    
    vec = []
    for i in range(2,3):
      aaVec = blankVec[:]
      aaVec[aaDict[seq[i]]] = 0.1
      vec.extend(aaVec)
      
    # 11 12 13 14 15   
    # Ha,Ca,Cb,C,N
    shifts = []
    indices = range(1,26)
    indices.remove(11)
    
    for i in indices:
      v = datum[i]
      if v is None:
        shifts.append(None)
      else:
        shifts.append(v/10.0)

    missing = datum[11]
    
    if (None not in shifts) and (missing is not None):
      missing /= 10.0
      
      vec = secStructs + angles #+ vec # + shifts
      datum = [vec, [missing,]]
      #print datum
      trainingData.append(datum)
 

  return trainingData

def angleDifference(a1, a2, oneTurn=6.2831853071796):

  delta = abs(a1- a2) % oneTurn
  
  if delta > (oneTurn/2.0):
    delta = oneTurn - delta

  return delta


def testFuncSS(nn, inputs, verbose=True):


  numSSGood = 0
  for data, known in inputs:
    
    predict = nn.update(data)
    
    helix1, sheet1, coil1  = known
    helix2, sheet2, coil2  = predict
        
    ssVec1 = [helix1, sheet1, coil1]
    ssVec2 = [helix2, sheet2, coil2]
    
    ss1 = ssVec1.index(max(ssVec1))
    ss2 = ssVec2.index(max(ssVec2))
    
    if ss1 == ss2:
      numSSGood += 1
       
  N = float(len(inputs))   

   
  pcSS =  (100.0*numSSGood)/N

  data = (pcSS,)
  
  if verbose:
    print "Success rate : %5.2f%% " % data
  
  return pcSS
 
def testFuncAngles(nn, inputs, verbose=True):

  meanPsi = 0.0
  meanPhi = 0.0
  numGood10 = 0
  numGood20 = 0
  numGood30 = 0
  numGood45 = 0

  numGoodPhi30 = 0
  numGoodPsi30 = 0
  
  numSSGood = 0
  
  err = 0.0
  for data, known in inputs:
    
    predict = nn.update(data)
    
    cosPhi1,  sinPhi1,  cosPsi1,  sinPsi1, helix1, sheet1, coil1  = known
    cosPhi2,  sinPhi2,  cosPsi2,  sinPsi2, helix2, sheet2, coil2  = predict
    
    phiKnown = atan2(sinPhi1,cosPhi1)
    psiKnown = atan2(sinPsi1,cosPsi1)
    
    phiPred = atan2(sinPhi2,cosPhi2)
    psiPred = atan2(sinPsi2,cosPsi2)
    
    dPhi = abs(angleDifference(phiPred,phiKnown))
    dPsi = abs(angleDifference(psiPred,psiKnown))
    
    meanPsi += dPhi
    meanPhi += dPsi
    
    ssVec1 = [helix1, sheet1, coil1]
    ssVec2 = [helix2, sheet2, coil2]
    
    ss1 = ssVec1.index(max(ssVec1))
    ss2 = ssVec2.index(max(ssVec2))
    
    if ss1 == ss2:
      numSSGood += 1
    
    for k in range(len(known)):
      err += 0.5*(known[k]-predict[k])**2

    if (dPhi<degrees45) and (dPsi<degrees45):
      if (dPhi<degrees30) and (dPsi<degrees30):
        if (dPhi<degrees20) and (dPsi<degrees20):
          if (dPhi<degrees10) and (dPsi<degrees10):
            numGood10 += 1
          numGood20 += 1
        numGood30 += 1
      numGood45 += 1
    
    if dPhi<degrees30:
      numGoodPhi30 += 1
      
    if dPsi<degrees30:
      numGoodPsi30 += 1
      
  N = float(len(inputs))   

  meanPsi /= N 
  meanPhi /= N
  meanPsi *= 180.0/pi
  meanPhi *= 180.0/pi
  
  pcSS =  (100.0*numSSGood)/N
  
  pcPhi = (100.0*numGoodPhi30)/N
  pcPsi = (100.0*numGoodPsi30)/N
 
  pc = (100.0*numGood30)/N
  data = (pc,pcSS, pcPhi,pcPsi,numGood10,numGood20,numGood30,numGood45,N,err,meanPhi, meanPsi)
  
  if verbose:
    print "Success rate : %5.2f%% :: SS %5.2f%% Phi %5.2f%% Psi %5.2f%%  %5d|%5d|%5d|%5d of %5d error %5.2f mean %5.2f, %5.2f" % data
  
  return pc
    

def testFuncMissingShift(nn, inputs):

  meanDiff = 0.0
  meanDiff2 = 0.0
  numGood  = 0
  err      = 0.0
  
  for data, known in inputs:
    
    predict = nn.update(data)
    delta = abs(known[0]-predict[0])*10.0
    meanDiff2 += delta*delta
    meanDiff += delta
    
    for k in range(len(known)):
      err += 0.5*(known[k]-predict[k])**2
    
    if delta < 0.25:
      numGood += 1
      
      
  N = float(len(inputs))   

  meanDiff2 /= N 
  meanDiff /= N 
  
  
  pc = (100.0*numGood)/N
  data = (pc,numGood,N,err,sqrt(meanDiff2),meanDiff)
  print "Success rate : %5.2f%% - %5d of %5d error %5.2f RMS %5.2f mean %5.2f" % data
  return pc

if __name__ == '__main__':
  
  # To-Do
  
  # Missing 1 or 2 of Ha, Ca, Cb, Co, N 
  
  # Train missing shifts
  
  print 'Reading database'
  #data = readDatabase('networkVarInput.txt')
  data = readDatabase('DB_TALOS_186_secStruct_nn')
  
  print 'Making training set'
  trainingSet = makeTrainingDataMissingShift(data)
  
  nIn  = len(trainingSet[0][0])
  nOut = len(trainingSet[0][1])
  nHid = 3
  
  print "Inputs", nIn
  print "Hidden", nHid
  print "Output", nOut

  nn = NN(nIn, nHid, nOut, testFuncMissingShift)

  nn.train(trainingSet, iterations=50, N=0.5, M=0.2)
 
  #testFuncMissingShift(nn, trainingSet)
  """

  for iter in range(1):

    print iter

    # Train predict angles
    print 'Reading database'
    data = readDatabase('DB_TALOS_186_secStruct_nn')
    print 'Read %d entries' % len(data)

    nTest = int(len(data)/2.0)

    testData = data[:nTest]
    data = data[nTest:]
    
    #for i in xrange():
    #  n = len(data)
    #  j = randint(0, n-1)
    #  testData.append(data.pop(j))



    print 'Making training set'
    trainingSet = makeTrainingDataAngles(data)
    print 'Training Database size: %d' % len(trainingSet)

    print 'Making test set'
    testSet = makeTrainingDataAngles(testData)
    print 'Test Database size: %d' % len(testSet)
 
    nIn  = len(trainingSet[0][0])
    nOut = len(trainingSet[0][1])
    nHid = 4
 
    print "Inputs", nIn
    print "Hidden", nHid
    print "Output", nOut

    nn = NN(nIn, nHid, nOut, testFuncSS)

    nn.train(trainingSet, iterations=20, N=0.5, M=0.2)
 
    print 'BEST OF TRAINING:'
 
    testFuncSS(nn, trainingSet, verbose=True)
 
    print 'TEST RESULT:'
 
    testFuncSS(nn, testSet, True)

    nn.weights
  """
