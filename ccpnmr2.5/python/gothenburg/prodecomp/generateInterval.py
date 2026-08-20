#======================================================================
# FUNCTION for generating list of intervals 
# and components from Peak list
# 
# fin -> file with a peak list -> string
#    It assumes that 1st and 2nd column are No of peak and 
#    HN peak position (x0) in points, respectively
# col -> No of column in the list with HN (ppm) -> integer 
# colN -> No of column in the list with N (ppm) -> integer
# colI -> No of column in the list with Peak Intensity -> integer
# m -> number of points for an interval (half range)-> float
# mpar -> shifts in an interval of x0+/-mpar*m will be taken -> float
# PeakFunc -> type of function: Gaussian or Lorentzian -> string
# n -> half-width at half-maximum (HWHM) in points -> float
# thresh -> threshold for the ratio: Lorentzian(x)/Intensity(x0) -> float 
#    such that a neighbour peak would give a component for x0
#======================================================================

# NBNB renamed. Used to be 'peaksToInterval.py'

def peak2comp(fin, col, colN, colI, m, mpar, PeakFunc, n, thresh):
  """
  Wrapper to generate intervals from file input
  """

  # Reading the data
  f = open(fin)
  buf = f.readlines()
  f.close()
  
  # Pini: HN peak position in points; IPini: Peak Intensity 
  Pini, IPini = [], []
  # HN, N - positions in ppm; pN - No of the peak
  csHN, csN, pN = [],[],[]
  #Finding the first line with peaks
  
  for line in buf:
    lineData = line.split()
  
    if lineData:

      #looks only into lines starting with number     
      if lineData[0][0].isdigit():
        pN.append(int(lineData[0]))
        Pini.append(float(lineData[1])-1)
        IPini.append(float(lineData[colI]))    
        csHN.append(float(lineData[col]))
        csN.append(float(lineData[colN]))
       
  
  
  peakData = (pN, Pini, IPini, csHN, csN)
  
  pNlst, stI, endI, comp, csHNlst, csNlst = peak2compProcess(peakData, m, mpar, PeakFunc, n, thresh)

  # Preparing the output
  IntlRes = ['NoPeak'+' '+'StartIntl'+' '+'EndIntl'+' '+'NoComp'+' '+'HN-shift'+' '+'N-shift'+'\n']
  [IntlRes.append(str(pNlst[x])+' \t'+str(stI[x])+' \t'+str(endI[x])+' \t'+str(comp[x])+' \t'+str(csHNlst[x])+' \t'+str(csNlst[x])+'\n') \
     for x in range(len(pNlst))]
  
  return IntlRes 
  

def peak2compProcess(peakData, m, mpar, PeakFunc, n, thresh):

  #**********************************************************
  # Lorentzian/Gaussian function
  # x0: peak position; Lx0: Intensity at x0; 
  # hwhm: half-width at half-maximum; x: position of interest
  #**********************************************************
  def LorF(Lx0,x0,hwhm,x):        
     return Lx0/(1+((x-x0)/hwhm)**2)
  
  def GausF(Lx0,x0,hwhm,x):
     return Lx0* math.exp(-LN2*((x-x0)/hwhm)**2)
  
  pN, Pini, IPini, csHN, csN = peakData
  
  # alternative, numpy-free sorting:
  # make list of tuples, starting with csHN
  ll = zip(csHN, csN, pN, Pini, IPini)
  # sort, csHN will be major sorting key
  ll.sort()
  # transform back to list of tuples (now sorted): csHN, csN, pN, Pini, IPini
  ll = zip (*ll)
  # extract individual sorted lists
  csHNlst = list(ll[0])
  csNlst = list(ll[1])
  pNlst = list(ll[2])
  Plst = list(ll[3])
  IPlst = list(ll[4])
  # remove ll, just in cas we have memory problems
  del ll
  
  
  stI  = [int(round(x-m)) for x in Plst]
  endI = [int(round(x+m)) for x in Plst]

  # Shifts in a range of +/-dmin around x0 will be considered 
  dmin = mpar*m
  k = range(len(Plst))
  comp = []
  #counts how many components are found in an interval using intensities of neighbouring resonances
  for j in range(len(Plst)):   
    #indcomp = [x for x in k if logical_and(abs(Plst[x]-Plst[j])<=dmin,x!=j)]
    indcomp = [x for x in k if abs(Plst[x]-Plst[j])<=dmin and x!=j]
    Prest   = [Plst[x] for x in indcomp]
    IPrest  = [IPlst[x] for x in indcomp]    
        
    if Prest:
        # split the neighbours in 2: right and left from the j-one 
        RPrest  = [Prest[x] for x in range(len(Prest)) if Prest[x] > Plst[j]]
        IRPrest = [IPrest[x] for x in range(len(Prest)) if Prest[x] > Plst[j]]
        LPrest  = [Prest[x] for x in range(len(Prest)) if Prest[x] < Plst[j]]
        ILPrest = [IPrest[x] for x in range(len(Prest)) if Prest[x] < Plst[j]]
     
        Rcomp, Lcomp = 0, 0
        if RPrest:
           if PeakFunc=='G':
              Icalc = [GausF(IRPrest[x],RPrest[x],n,endI[j]) for x in range(len(RPrest))]
           elif PeakFunc=='L':
              Icalc = [LorF(IRPrest[x],RPrest[x],n,endI[j]) for x in range(len(RPrest))]
           res = [x for x in Icalc if x > thresh*IPlst[j]]
           Rcomp = len(res)
           
        if LPrest:
           if PeakFunc=='G':
              Icalc = [GausF(ILPrest[x],LPrest[x],n,stI[j]) for x in range(len(LPrest))]
           elif PeakFunc=='L':
              Icalc = [LorF(ILPrest[x],LPrest[x],n,stI[j]) for x in range(len(LPrest))]
           res = [x for x in Icalc if x > thresh*IPlst[j]]
           Lcomp = len(res)
           
        comp.append(Rcomp+Lcomp+1)              
    else:
        # no neighbouring resonances in a range
        comp.append(1)

  return pNlst, stI, endI, comp, csHNlst, csNlst
#======
# MAIN
#======
import math
LN2 = math.log(2)

if __name__ == "__main__":
        # for test    fin, col, colN, colI, m, mpar, PeakFunc, n, thresh    
        result = peak2comp('pk1assign.tab', 5, 6, 17, 2.5, 10., 'G', 2., 0.1)
        for line in result:
          print line
        #f = open('ftest','w')
        #f.writelines(result)
        #f.close()


