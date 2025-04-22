import numpy as np
import matplotlib.pyplot as plt

# material parameters
K = [0.8, 2.2, 2.5, 2.5, 2.4]  # thermal conductivity per domain ... W / (m K)
C = [310.0 * 2400, 300.0 * 2400, 810.0 * 2400, 810.0 * 2400, 970.0 * 2400]  # heat capacity per domain 

# generate conduction data
lowbound = [-1e2, -1e2, -1e2, -1e2, -1e2]  # lower bound for lin space for each material
upbound  = [1e2, 1e2, 1e2, 1e2, 1e2]  # upper bound for lin space for each material
nDB = [int(1e5+1), int(1E5+1), int(1E5+1), int(1E5+1), int(1E5+1)]
# the above code will be later replaced by real data
# then nDB needs to be the length of each data set

# geometrical parameters
L = [0.04, 0.14, 0.22, 0.16, 11.44]  # domain widths in m
NN = [2, 4, 5, 3, 16]  # number of elements per domain

N = 0  # total number of elements used, is calculated
NNN = np.zeros(len(NN))  # assignment which element numbers to which domain
count = 0
for i in NN:
    N += i
    NNN[count] = N
    count += 1

LL = 0.0  # total length
for i in L:
    LL += i

mats = len(K)
mnDB = int(max(nDB))
gData = np.zeros((mats, mnDB))
qData = np.zeros((mats, mnDB))
for i in range(mats):
    a = np.linspace(lowbound[i], upbound[i], nDB[i])
    gData[i, 0:nDB[i]] = a[0:nDB[i]]
    qData[i, 0:nDB[i]] = K[i] * gData[i, 0:nDB[i]]

# restrain bottom
bflag = True
tbottom = 10.0  # 

# define a Data Driven metric
Kddfac = 0.4
Kdd = np.zeros(mats)
mmm = np.ones(mats) * Kddfac
Kdd = mmm * K  # this matrix is chosen

# result file
resfile = "00_DD_res.dat"
    
# initial temp
inittemp = True
itemps = [12.17112, 12.75, 8.188, 9.625, 9.625]

# load function
loadtype = 1  # 0 is by linear ramp, 1 is by file
loadfile = "load.dat"
    
# material parameters
if loadtype == 0:
    tau0 = 1.0*LL**2/20.0 # characteristic time
    
    # loading parameters
    Tmax = 10.0        # maximum temperature
    tRamp = 0.1*tau0   # ramp time
    tFinal = 2.0*tau0  # final time
    nSteps = 100       # number of time steps
    dt = tFinal/nSteps # time step
    
elif loadtype == 1:
    loaddata = np.loadtxt(loadfile)
    tFinal = loaddata[-1, 0]  # second number is column
    dt = 60.0  # time increment
    nSteps = tFinal / dt
    print(nSteps)
    if abs(nSteps - int(nSteps)) > 0:
        print("tFinal not divisable by dt")
        # exit()
    nSteps = int(nSteps)
    xp = loaddata[:, 0]
    fp = loaddata[:, 1]
    
else:
    print("No such load type")

# assemble FE matrices
Bmat = np.zeros((N,N+1))
Cmat = np.zeros((N+1,N+1))
Kref = np.zeros((N+1,N+1))
Mref = np.zeros((N+1,N+1))
Dmat = np.zeros((N,N+1))

# Le = L/N
for e in range(N):
    for i in range(len(NN)):
        if e < NNN[i]:
            ei = i
            break
    Le = L[ei] / NN[ei]
    Ce = C[ei]
    Ke = K[ei]
    
    # gradient matrix
    Bmat[e,e]  = -1.0 / Le
    Bmat[e,e+1] = 1.0 / Le
    
    # Dmat
    Dmat[e,e]  = -1.0
    Dmat[e,e+1] = 1.0
    
    # conductivity matrix
    Cmat[e,e] += 1.0/3.0 * Ce * Le
    Cmat[e+1,e+1] += 1/3.0 * Ce * Le
    Cmat[e,e+1] += 1.0/6.0 * Ce * Le
    Cmat[e+1,e] += 1.0/6.0 * Ce * Le
    
    # Kref
    Kref[e,e] += 1.0
    Kref[e+1,e+1] += 1.0
    Kref[e,e+1] += -1.0
    Kref[e+1,e] += -1.0
    
    # Mref derived from Kref for DD
    Mref[e,e] += 1.0 / Le * Kdd[ei]
    Mref[e+1,e+1] += 1.0 / Le * Kdd[ei]
    Mref[e,e+1] += -1.0 / Le * Kdd[ei]
    Mref[e+1,e] += -1.0 / Le * Kdd[ei]

# initialize arrays
T_DD = np.zeros((nSteps+1,N+1))  # temperature data for all nodes and all time steps ... also 0th time step
gDD = np.zeros((nSteps+1,N))  # solutions for -grad(T) for all elements
qDD = np.zeros((nSteps+1,N))  # solutions for temperature fluxes
idx = np.zeros((nSteps+1,N),dtype=int)  

if inittemp:
    # always set the whole temperature of the beam
    for i in range(len(NN)):
        numpoints = int(NN[i])  # number of points for each segment
        lowfp = itemps[i]
        if i == 0:
            lowxp = 0  # catch the case of first segment
        else:
            lowxp = NNN[i - 1]  # either way the first entry in T is the number in NNN of i - 1
            
        highxp = NNN[i]  # upper interp for x directly as number of elements
        if i == len(NN) - 1:  # if last i, take the bottom temp to interpolate, else take the next temp
            if not bflag:
                tbottom = 0.0  # safety measure
            highfp = tbottom
        else:
            highfp = itemps[i + 1]
            
        # interpolate
        xxp = [lowxp, highxp]
        ffp = [lowfp, highfp]
        for j in range(int(lowxp), int(highxp)):
            T_DD[0, j] = np.interp(j, xxp, ffp)
            
    if bflag:  # if the bottom is fixed, set also the previous bottom temp
        T_DD[0, -1] = tbottom
    else:
        T_DD[0, -1] = 0.0

# get average Kdd
Kdda = np.average(Kdd)

# initialize pairing (zero gradient and flux), here average works, since init only, and just first data set should work since init
dist = 0.5*(Kdda*(gData[0, :]*gData[0, :])+(1.0/Kdda)*(qData[0, :]*qData[0, :]))
idx[0,:] = np.argmin(dist)


# Data Driven iteration matrix ... first node not, since fixed?
if bflag:
    Add = np.block([[Cmat[1:-1,1:-1],             dt*Mref[1:-1,1:-1]],
                    [dt*Mref[1:-1,1:-1], -Cmat[1:-1,1:-1]]])
else:
    Add = np.block([[Cmat[1:,1:],             dt*Mref[1:,1:]],
                    [dt*Mref[1:,1:], -Cmat[1:,1:]]])
if bflag:
    # time integration loop
    for n in range(nSteps):
    
        # init
        idx[n+1,:] = idx[n,:]
        dist_ref = 1.0
        
        t = (n+1)*dt
        
        # load function 
        if loadtype == 0:
            if t < tRamp:
                T_DD[n+1, 0] = Tmax*t/tRamp
                T_DD[n+1,-1] = tbottom
            else:
                T_DD[n+1,0 ] = Tmax
                T_DD[n+1,-1] = tbottom
        elif loadtype == 1:
            T_DD[n+1, 0 ] = np.interp(t, xp, fp)
            T_DD[n+1, -1] = tbottom
            
        print("Step %d - Time = %f:"%(n+1,t))
        
        # DD iterative loop
        itMax = 10
        iter = 0
        New = N-1
        while (iter < itMax):
            # compute RHS
            Rdd = np.zeros(2*New)  # reduce since one more node restrained
            for e in range(New):  
                for i in range(len(NN)):
                    if e < NNN[i]:
                        ei = i
                        break
                Rdd[0:New]     += dt*Dmat[e,1:-1]*qData[ei, idx[n+1,e]] 
                Rdd[New:2*New] -= dt*Dmat[e,1:-1]*Kdd[ei]*gData[ei, idx[n+1,e]]  
                
            adebug = (np.dot(Kref[1:-1,1:-1],T_DD[n,1:-1])+Kref[1:-1,0]*T_DD[n+1,0]+Kref[1:-1,-1]*T_DD[n+1,-1])  # exclude last node
            for e in range(New):  # loop over elements stays
                for i in range(len(NN)):
                    if e < NNN[i]:
                        ei = i
                        break
                Le = L[ei] / NN[ei]  # get element length for this material back
                Rdd[e]     -= Cmat[e + 1,0]*(T_DD[n+1,0]-T_DD[n,0])
                Rdd[e + New] -= dt*(Kdd[ei]/Le)* adebug[e]
            
            dZ = np.linalg.solve(Add,Rdd)
            T_DD[n+1,1:-1] = T_DD[n,1:-1]+dZ[0:New]  # here just update free nodal temperatures
            gDD[n+1,:] = -np.dot(Bmat,T_DD[n+1,:])  # gradients are always calculated over all existing nodes
            for e in range(N):
                for i in range(len(NN)):
                    if e < NNN[i]:
                        ei = i
                        break
                qDD[n+1,e] = qData[ei, idx[n+1,e]] - Kdd[ei]*np.dot(Bmat[e,1:-1],dZ[New:2*New]) 
            
            # update pairing
            dist_tot = 0.0
            for e in range(N):
                for i in range(len(NN)):
                    if e < NNN[i]:
                        ei = i
                        break
                Le = L[ei] / NN[ei]
                dist = 0.5*(Kdd[ei]*(gData[ei, :]-gDD[n+1,e])*(gData[ei, :]-gDD[n+1,e]) 
                           +(1.0/Kdd[ei])*(qData[ei, :]-qDD[n+1,e])*(qData[ei, :]-qDD[n+1,e]))
                idx[n+1,e] = np.argmin(dist)
                if idx[n+1, e] == 0:
                    print("W A R N I N G: M I N  F O U N D")
                elif idx[n+1, e] == len(gData[ei,:]):
                    print("W A R N I N G: M A X  F O U N D")
                dist_tot += Le*dist[idx[n+1,e]]
    
            # convergence test
            print("\titer %d: dist=%f"%(iter,dist_tot))
            if (iter == 0):
                dist_ref = dist_tot
            else:
                if (dist_tot < dist_ref):
                    if ((dist_ref-dist_tot) < 1.0e-3*dist_ref): break
                    dist_ref = dist_tot
    
            iter = iter+1
            
else:
    # time integration loop
    for n in range(nSteps):
    
        # init
        idx[n+1,:] = idx[n,:]
        dist_ref = 1.0
        
        t = (n+1)*dt
        
        # load function
        if loadtype == 0:
            if t < tRamp:
                T_DD[n+1,0] = Tmax*t/tRamp
            else:
                T_DD[n+1,0] = Tmax
        elif loadtype == 1:
            T_DD[n+1, 0 ] = np.interp(t, xp, fp)
        print("Step %d - Time = %f:"%(n+1,t))
        
        # DD iterative loop
        itMax = 10
        iter = 0
        while (iter < itMax):
            # compute RHS
            Rdd = np.zeros(2*N)
            for e in range(N):
                for i in range(len(NN)):
                    if e < NNN[i]:
                        ei = i
                        break
                
                Rdd[0:N]   += dt*Dmat[e,1:]*qData[ei, idx[n+1,e]]  
                Rdd[N:2*N] -= dt*Dmat[e,1:]*Kdd[ei]*gData[ei, idx[n+1,e]]
                
            adebug = (np.dot(Kref[1:,1:],T_DD[n,1:])+Kref[1:,0]*T_DD[n+1,0])
            for e in range(N):
                for i in range(len(NN)):
                    if e < NNN[i]:
                        ei = i
                        break
                Le = L[ei] / NN[ei]  # get element length for this material back
                Rdd[e]     -= Cmat[e + 1,0]*(T_DD[n+1,0]-T_DD[n,0])
                Rdd[e + N] -= dt*(Kdd[ei]/Le)* adebug[e]
            
            dZ = np.linalg.solve(Add,Rdd)
            T_DD[n+1,1:] = T_DD[n,1:]+dZ[0:N]
            gDD[n+1,:] = -np.dot(Bmat,T_DD[n+1,:])
            for e in range(N):
                for i in range(len(NN)):
                    if e < NNN[i]:
                        ei = i
                        break
                qDD[n+1,e] = qData[ei, idx[n+1,e]] - Kdd[ei]*np.dot(Bmat[e,1:],dZ[N:2*N])
            
            # update pairing
            dist_tot = 0.0
            for e in range(N):
                for i in range(len(NN)):
                    if e < NNN[i]:
                        ei = i
                        break
                Le = L[ei] / NN[ei]
                dist = 0.5*(Kdd[ei]*(gData[ei, :]-gDD[n+1,e])*(gData[ei, :]-gDD[n+1,e]) 
                           +(1.0/Kdd[ei])*(qData[ei, :]-qDD[n+1,e])*(qData[ei, :]-qDD[n+1,e]))
                idx[n+1,e] = np.argmin(dist)
                dist_tot += Le*dist[idx[n+1,e]]
    
            # convergence test
            print("\titer %d: dist=%f"%(iter,dist_tot))
            if (iter == 0):
                dist_ref = dist_tot
            else:
                if (dist_tot < dist_ref):
                    if ((dist_ref-dist_tot) < 1.0e-3*dist_ref): break
                    dist_ref = dist_tot
    
            iter = iter+1

# save important results
Tout = np.zeros((nSteps + 1, len(NN) + 3))
times = np.linspace(0,tFinal,nSteps+1)
Tout[:, 0 ] = times  # Time for temperature time plot versus sensors
Tout[:, 1 ] = T_DD[:, 0]  # temperature at the top, prescribed
Tout[:, -1] = T_DD[:, -1]  # temperature at the bottom, prescribed
for i in range(len(NN)):
    Tout[:, i + 2] = T_DD[:, int(NNN[i])]
    
np.savetxt(fname="00_DD_dat", X=Tout)

# plot results
plt.plot(times,T_DD[:,0],label='x=0')
for i in NNN[:-1]:
     lab = str("Number " + str(int(i)))
     plt.plot(times, T_DD[:, int(i)], label=lab)
plt.plot(times,T_DD[:,-1],label='x=L')
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.legend()
plt.savefig("test.png", dpi=600)
plt.show()

# save all results
out = np.zeros((nSteps+1,N+2))
out[:,1:] = T_DD
out[:,0 ] = times
np.savetxt(fname=resfile, X=out)




