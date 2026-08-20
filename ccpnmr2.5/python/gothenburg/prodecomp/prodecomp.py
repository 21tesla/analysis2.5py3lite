"""
======================COPYRIGHT/LICENSE START==========================

prodecomp.py: Part of PRODECOMP program

Copyright (C) 2008 Doroteya K. Staykova, Daniel Malmodin, Martin Billeter 
(University of Gothenburg, Sweden)

=======================================================================

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
======================COPYRIGHT/LICENSE END============================

for further information, please contact :
- our website (http://www.lundberg.gu.se/nmr/)
- dkstaykova@gmail.com, martin.billeter@chem.gu.se 

=======================================================================

If you are using this software for academic purposes, we suggest
quoting (one of) the following references:

===========================REFERENCE START=============================
Malmodin, D. and Billeter, M. (2005) Multiway Decomposition of NMR Spectra 
with Coupled Evolution Periods. J. Am. Chem. Soc., 127, 13486-13487. 

Malmodin, D. and Billeter, M. (2006) Robust and versatile interpretation 
of spectra with coupled evolution periods using multi-way decomposition. 
Magn. Reson. Chem., 44, S185-S195.

Staykova, D. K.; Fredriksson, J.; Bermel, W. and Billeter, M. (2008) 
Assignment of protein NMR spectra based on projections, multi-way decomposition 
and a fast correlation approach. J. Biomol. NMR., 42, 87-97.

===========================REFERENCE END===============================
"""

import time

from numpy import array, asmatrix, concatenate, diag, dot, finfo, flipud, \
                  linalg, logical_and, matrix, mod, ndarray, ones, outer, \
                  put, r_, random, roll, sometrue, squeeze, where, zeros
                  
# "sum" is a python built-in
from numpy import sum as nSum

linalgInv    = linalg.inv
randomSeed   = random.seed
randomSample = random.random_sample

publicDocumentation = """PRODECOMP Version 3.0 (07 Nov 2008 05:24:36)
Decomposition of 2D projections, representing a high-dimensional spectrum,
to a set of components (defined by one-dimensional "shapes".)

Copyright (C) 2008 Doroteya K. Staykova, Daniel Malmodin, Martin Billeter
(University of Gothenburg, Sweden)

For further help, please visit: http://www2.chem.gu.se/bcbp/nmr/
or contact: dkstaykova@gmail.com, martin.billeter@chem.gu.se 
        """

def prodecomp(Pmx, defs, intl, cmps, rglf, itrs): 
    # Input arguments/Type
    # Pmx ('numpy.ndarray'); defs ('numpy.ndarray'); intl ('list'); cmps ('list'); rglf('float'); itrs('int') 
    def fastnnls(XtX,Xty):
        # FASTNNLS function
    	# Translated from Matlab version (fastnnls.m)
        XtX = asmatrix(XtX)
        Xty = asmatrix(Xty)
        nrm = nSum(abs(XtX),1).max(0)   
        sz  = array(XtX.shape)
        tol = 10*finfo(float).eps*nrm*sz.max(0)    
        n = XtX.shape[1]
        P = zeros((1,n))
        Z = asmatrix(range(1,n+1))    
        x = P.copy()
        # range for indecies
        ZZ = Z-1
        w = Xty - dot(XtX,x.T)    
        # set up iteration criterion        
        it = 0
        itmax = 30*n
        
        # outer loop to put variables into set to hold positive coefficients
        z = zeros((n))
        while logical_and(sometrue(Z),sometrue(w[ZZ.tolist()] > tol)):       
            wt = w[ZZ.tolist()].max(0)        
            t = where(w==wt)[0]        
            put(P, t, t+1)
            put(Z, t, 0)
            PP = where(P!=0)[1] 
            ZZ = where(Z!=0)[1]       
            nzz = ZZ.shape
            var1 = Xty[PP].T
            var2 = asmatrix(XtX[array([PP]).T,PP].T)       
            var3 = dot(var1,linalgInv(var2))           
            
            in1 = asmatrix(PP).T        
            # put considers matrix as flatten        
            put(z,in1,var3)
            put(z,ZZ,zeros((nzz[0],nzz[1])))        
            # z comes out as row         
            #----------------------------------------------------
            # inner loop to remove elements from the positive set
            while logical_and(sometrue(z[PP] <= tol), it < itmax):
                it = it + 1           
                QQ = where(logical_and((z <= tol),P))[1]
                f = x.flatten(1)[QQ]/(x.flatten(1)[QQ] - z[QQ])            
                alpha = f.min()            
                x = x+alpha*(z-x)            
                ij = where(logical_and((abs(x) < tol),(P!=0)))[1]            
                put (Z, ij, ij.T+1)            
                # P[ij]; ij has to be column vector to make the choice of elements in a column of P
                put(P, ij, zeros((len(ij))).reshape(-1,1))                
                PP = where(P != 0)[1] 
                ZZ = where(Z!=0)[1] 
                nzz = ZZ.shape    
                var1 = Xty[PP].T
                var2 = asmatrix(XtX[array([PP]).T,PP].T)            
                if var2.shape == (0,0): var3 = []
                else: var3 = dot(var1,linalgInv(var2))
                                          
                put(z,PP,var3)
                put(z,ZZ,zeros((nzz[1],nzz[0]))) 
                
            x = z.copy()        
            w = Xty-dot(XtX,asmatrix(x).T)
        
        # comes as row
        return asmatrix(x) 

    #**************
    def midmatrix(Fin,addorsub): 
        # [m n p]=size(Fin); n is always 1 since Fin comes as (m,1,p)   
        Fin = Fin[:,0,:]
        [m,p] = Fin.shape
        # Folding of the spectra
        if mod(m,2)==0: fold = (m+2)/2
        else: fold = (m+1)/2
        
        if addorsub[0,0] == -1: Fin = flipud(Fin)

        Fout = zeros((m,m,p))
        Fin_sh = roll(Fin,fold,0)
        for j in range(m): Fout[:,j,:] = roll(Fin_sh,j,0) 
        
        if addorsub[0,1] == -1:
            [a,b,c] = Fout.shape
            # Fout=flipdim(Fout,2);
            Fout = Fout[0:a,r_[b-1:-1:-1],0:c]

        return Fout

    #**************
    def signMx(a):
        b = a.copy()
        i1 = where(b > 0)
        b[i1] = 1
        i2 = where(b == 0)
        b[i2] = 0
        i3 = where(b < 0)
        b[i3] = -1

        return b
 
    
    #==MAIN PRODECOMP FUNCTION==
    [x1,x2,x3] = Pmx.shape
    # NORMALIZATION of the input data to the MAX value for all spectra
    
    mtot = Pmx.max()
    for i in range(x3):
        Pmx[:,:,i] = Pmx[:,:,i]/mtot    
    
    interval = matrix(intl)
    ncomps = matrix(cmps)    
    Nb = defs.shape[1]
    branches = array([ones(Nb+1)])
    regfac = rglf
    iterations = itrs
    # PRODECOMP
    rng = range(int(interval[0,0]),int(interval[0,1])+1)
    projs = Pmx[:,rng,:].copy()
    [mp,np,op] = projs.shape
    [op,nd] = defs.shape
        
    # Cleaning the memory
    Ldata, Ldefs, Pmx = [[],[],[]]
    # Positioning the a(0) of a vector, i.e. if len(a)=5 points
    # midpos = 3: a(-2) a(-1) a(0) a(1) a(2)
    # if len(a) = 4; midpos = 2; a(-1) a(0) a(1) a(2)
    if mod(mp,2)==0: midpos = mp/2    
    else: midpos = (mp+1)/2

    for comps in range(int(ncomps[0,0]),int(ncomps[0,1])+1):
        comps = array([comps])
        branches = branches.flatten(1)

        randomSeed([0])
        f = randomSample((mp,nd,comps)) # initialization of f 
        randomSeed([0])
        fdir = randomSample((np,comps)) # initialization of fdir             
        n = 0
        it = 0    
        while it < iterations:       
            if n < nd: n = n+1            
            else: it = it+1
               
            start_time = time.time()
            for i in range(n,0,-1):        
                P = []
                M1 = []
                for j in range(1,op+1):
                    # defs(j,n+1:nd); (n-1)+1=n; (nd-1)+1=nd - range to nd-1
                    if logical_and(nSum(abs(defs[j-1,r_[n:nd]])) == 0, abs(defs[j-1,i-1]) == 1):
                            if P==[]: P = projs[:,:,j-1]
                            else: P = concatenate((P,projs[:,:,j-1]))
                            ff = zeros((mp,1,comps))
                            ff[midpos-1,0,:] = 1
                            # [1:i-1 i+1:n]
                            rn1 = range(1,i)
                            rn2 = range(i+1,n+1)
                            
                            if logical_and(rn1!=[],rn2!=[]): rng = rn1 + rn2
                            elif rn1!=[]: rng = rn1
                            elif rn2!=[]: rng = rn2
                            else: rng = None                                           

                            if rng != None:
                                for jj in rng:
                                    if abs(defs[j-1,jj-1]) > 0:                                    
                                        M2 = midmatrix(ff,asmatrix([1,defs[j-1,jj-1]]))                                    
                                        for jjj in range(1,comps+1):
                                            inff = dot(M2[:,:,jjj-1],f[:,jj-1,jjj-1])
                                            ff[:,0,jjj-1] = inff.copy()                                    
                                            
                            rn3 = range(1,int(comps/branches[i-1]+1))
                            [m1,m2,m3] = ff.shape
                            fff = zeros((m1,1,len(rn3)))                    
                            for jjj in rn3:                            
                                dim3 = range(int(branches[i-1]*jjj+1-branches[i-1]),int(branches[i-1]*jjj+1))
                                d3 = [d-1 for d in dim3]
                                # ff[:,:,:] = ff[:,0,:]; sum along axes 3 gives 2D matrix
                                fff[:,:,jjj-1] = ff[:,:,d3].sum(axis=2)                        
                        
                            M2 = midmatrix(fff,asmatrix([1,defs[j-1,i-1]]))
                            # reshape in Python takes row wise, reshape in Matlab: columnwise
                            # used another approach to simulate matlab result of reshape                       
                            n3 = M2.shape[2]
                            rM2 = []
                            for c in range(0,n3):
                                if rM2 == []: rM2 = M2[:,:,0]
                                else: rM2 = concatenate((rM2,M2[:,:,c]), axis=1) 
                            
                            if M1==[]: M1 = asmatrix(rM2)
                            else: M1 = concatenate((asmatrix(M1),asmatrix(rM2)))
                            
                rn4 = range(1,int(comps/branches[i-1]+1))
                # calculation of possible size of M3 for its initialization
                lim2 = int(comps/branches[i-1]*mp)
                [l1,l2] = fdir.shape
                [l3,l4] = M1.shape

                # alternative approach for reducing the calculation time for MTM and M3.T*x
                MTM = zeros((lim2,lim2))
                m3tx = zeros((lim2,1))
                Lavec = [[] for x in range(len(rn4))]
                Lbmx = [[] for x in range(len(rn4))]
                Lc = -1
                for jjj in rn4:                 
                        dim2 = range(int((jjj-1)*mp+1),int(jjj*mp+1))
                        d2 = [k-1 for k in dim2] 
                        Lc = Lc+1
                        Lavec[Lc] = fdir[:,jjj-1].reshape(-1,1).copy()
                        Lbmx[Lc] = M1[:,d2].copy()
                        pmx = zeros((lim2/len(rn4),1))
                        for count in range(Lavec[Lc].shape[0]):
                                if Lavec[Lc][count,0]!=0:                                       
                                        pmx = pmx + dot(Lbmx[Lc].T,P[:,count].reshape(-1,1))*Lavec[Lc][count,0]
                        m3tx[Lc*lim2/len(rn4):(Lc+1)*lim2/len(rn4)] = pmx.copy()
                        pmx = []

                for iL in range(len(Lavec)):
                        for jL in range(len(Lavec)):
                                # unit element of MTM
                                smx = zeros((lim2/len(rn4),lim2/len(rn4)))                                               
                                for count in range(Lavec[0].shape[0]):
                                        if logical_and(Lavec[iL][count,0]!=0,Lavec[jL][count,0]!=0):                                    
                                                smx = smx + dot(Lbmx[iL].T,Lbmx[jL])*Lavec[iL][count,0]*Lavec[jL][count,0]
                                MTM[iL*lim2/len(rn4):(iL+1)*lim2/len(rn4),jL*lim2/len(rn4):(jL+1)*lim2/len(rn4)] = smx.copy() 
                                smx = []
   
                MTM=MTM+regfac*diag(diag(MTM, k=0)*((signMx(diag(MTM, k=0))+1)/2), k=0);          
                ffff = fastnnls(MTM,m3tx).reshape(int(comps/branches[i-1]),mp)                           
                ffff = ffff.T
                # clean the memory
                MTM = []
                            
                for jjj in range(1,int(branches[i-1]+1)):
                    slc = range(jjj,int(comps+jjj-branches[i-1]+1),int(branches[i-1]))
                    islc = [s-1 for s in slc]                 
                    f[:,i-1,islc] = ffff.copy()            
                
            P = []
            M1 = []
            for j in range(1,op+1):
                if nSum(abs(defs[j-1,r_[n:nd]])) == 0:
                    if P==[]: P = projs[:,:,j-1]
                    else: P = concatenate((P,projs[:,:,j-1]))
                    ff = zeros((mp,1,comps))
                    ff[midpos-1,0,:] = 1
                    for jj in range(1,n+1):
                        if abs(defs[j-1,jj-1]) > 0:
                            M2 = midmatrix(ff,asmatrix([1,defs[j-1,jj-1]]))
                            for jjj in range(1,comps+1):
                                inff2 = dot(M2[:,:,jjj-1],f[:,jj-1,jjj-1])
                                ff[:,0,jjj-1] = inff2.copy()                        
                

                    rn3 = range(1,int(comps/branches[n]+1))
                    [m1,m2,m3] = ff.shape
                    fff = zeros((m1,1,len(rn3)))                    
                    for jjj in rn3:                        
                        dim3 = range(int(branches[n]*jjj+1-branches[n]),int(branches[n]*jjj+1))
                        d3 = [d-1 for d in dim3]
                        # ff[:,:,:] = ff[:,0,:]; sum along axes 3 gives 2D matrix
                        fff[:,:,jjj-1] = ff[:,:,d3].sum(axis=2)
                    
                    if M1==[]: M1 = fff
                    else: M1 = concatenate((M1,fff))
                    
            M3 = squeeze(M1).copy()
            M3 = asmatrix(M3)    
            [rM3,cM3] = M3.shape
            if rM3 == 1: M3 = M3.reshape(-1,1)

            MTM = outer(zeros(M3.shape[1]),zeros(M3.shape[1]))               
            MTM = dot(M3.T,M3)           
            MTP = dot(M3.T,P)    
            mx = MTM.shape[0]
            ffff = zeros((np,mx))
            for jjj in range(1,np+1):        
                ffff[jjj-1,:] = fastnnls(MTM,MTP[:,jjj-1]) 
            
            for jjj in range(1,int(branches[n]+1)):
                slc = range(jjj,int(comps+jjj-branches[n]+1),int(branches[n]))
                islc = [s-1 for s in slc]            
                fdir[:,islc] = ffff.copy()
            
            # output time for each iteration
            elapsed_time = time.time()-start_time
            print 'iteration', it, '->', elapsed_time, 'seconds'
            # clean the memory
            M3, MTM, MTP = [[],[],[]]      
    
    return fdir, f

