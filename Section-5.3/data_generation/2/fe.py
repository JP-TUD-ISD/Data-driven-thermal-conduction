import numpy as np
import matplotlib.pyplot as plt

# material parameters
K = [1.3, 0.5, 2.3, 2.3, 2.5]  # thermal conductivity per domain ... W / (m K)
# order is sma, concrete/asphalt, FSS, subsoil
C = [650.0 * 2400, 310.0 * 2400, 300.0 * 2400, 300.0 * 2400, 1000.0 * 2400]  # heat capacity per domain 
# values are specific in J/(kg K), multiplied with density in kg /m^3

# geometrical parameters
L = [0.04, 0.14, 0.22, 0.16, 11.44]  # domain widths in m
NN = [10, 20, 30, 20, 150]  # number of elements per domain

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

print(N)
print(NN)
print(NNN)

# restrained bottom side
bflag = True
tbottom = 10.0  # source: Bodenzustand Deutschland of Umweltbundesamt

# load function
loadtype = 1  # 0 is by linear ramp, 1 is by file
loadfile = "load.dat"

# initial temp
inittemp = True
itemps = [14.26218, 15.26239, 11.5, 9.5, 9.5]

# result file containing all temperatures and times
resfile = "T.dat"

if loadtype == 0:
    # time
    tau0 = C[0]*LL**2/K[0] # characteristic time
    
    # loading parameters
    Tmax = 10.0        # maximum temperature
    tRamp = 0.1*tau0   # ramp time
    tFinal = 2.0*tau0  # final time
    nSteps = 100       # number of time steps
    dt = tFinal/nSteps # time step
elif loadtype == 1:
    loaddata = np.loadtxt(loadfile)
    tFinal = loaddata[-1, 0]  # second number is column
    dt = 120  # time increment
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
BforQ = np.zeros((N,N+1))
Cmat = np.zeros((N+1,N+1))
Kref = np.zeros((N+1,N+1))
for e in range(N):
    for i in range(len(NN)):
        if e < NNN[i]:
            ei = i
            break
    Le = L[ei] / NN[ei]
    Ce = C[ei]
    Ke = K[ei]
    Bmat[e,e]  = -1.0 / Le
    Bmat[e,e+1] = 1.0 / Le
    BforQ[e,e]  = -1.0 * Ke / Le
    BforQ[e,e+1] = 1.0 * Ke / Le
    Cmat[e,e] += 1.0/3.0 * Ce * Le
    Cmat[e+1,e+1] += 1/3.0 * Ce * Le
    Cmat[e,e+1] += 1.0/6.0 * Ce * Le
    Cmat[e+1,e] += 1.0/6.0 * Ce * Le
    Kref[e,e] += 1.0 * Ke / Le
    Kref[e+1,e+1] += 1.0 * Ke / Le
    Kref[e,e+1] += -1.0 * Ke / Le
    Kref[e+1,e] += -1.0 * Ke / Le
Kmat = Kref

T = np.zeros((nSteps+1,N+1))
QQ = np.zeros((nSteps+1,N))

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
            T[0, j] = np.interp(j, xxp, ffp)
            
    if bflag:  # if the bottom is fixed, set also the previous bottom temp
        T[0, -1] = tbottom
    else:
        T[0, -1] = 0.0
    QQ[0, :] = np.matmul(BforQ[:,:], T[0, :])
    
# time integration
if bflag:  # restrain second node as well
    if loadtype == 0:
        beta = 0.5
        for n in range(nSteps):
            t = (n+1)*dt
            if t < tRamp:
                T[n+1,0] = Tmax*t/tRamp
                T[n+1,-1]= tbottom
            else:
                T[n+1,0] = Tmax
                T[n+1,-1]= tbottom
            T[n+1,1:-1] = np.linalg.solve(Cmat[1:-1,1:-1]+beta*dt*Kmat[1:-1,1:-1],  # matrix to be solved
                                          (Cmat[1:-1,:]-(1-beta)*dt*Kmat[1:-1,:])@T[n]  # initial load ...
                                           -(Cmat[1:-1,0 ]+beta*dt*Kmat[1:-1,0 ])*T[n+1,0]  # condensing first node
                                           -(Cmat[1:-1,-1]+beta*dt*Kmat[1:-1,-1])*T[n+1,-1])  # condensing last node
    elif loadtype == 1:
        beta = 0.5
        for n in range(nSteps):
            t = (n+1)*dt
            
            # get load function value
            T[n + 1, 0] = np.interp(t, xp, fp)
            T[n + 1,-1] = tbottom
            
            # get solution
            T[n+1,1:-1] = np.linalg.solve(Cmat[1:-1,1:-1]+beta*dt*Kmat[1:-1,1:-1],  # matrix to be solved
                                          (Cmat[1:-1,:]-(1-beta)*dt*Kmat[1:-1,:])@T[n]  # initial load ...
                                           -(Cmat[1:-1,0 ]+beta*dt*Kmat[1:-1,0 ])*T[n+1,0]  # condensing first node
                                           -(Cmat[1:-1,-1]+beta*dt*Kmat[1:-1,-1])*T[n+1,-1])  # condensing last node
            QQ[n+1, :] = np.matmul(BforQ[:,:], T[n+1, :])
    
else:
    if loadtype == 0:
        beta = 0.5
        for n in range(nSteps):
            t = (n+1)*dt
            if t < tRamp:
                T[n+1,0] = Tmax*t/tRamp
            else:
                T[n+1,0] = Tmax
            T[n+1,1:] = np.linalg.solve(Cmat[1:,1:]+beta*dt*Kmat[1:,1:],
                                        (Cmat[1:,:]-(1-beta)*dt*Kmat[1:,:])@T[n]-(Cmat[1:,0]+beta*dt*Kmat[1:,0])*T[n+1,0])
    elif loadtype == 1:
        beta = 0.5
        for n in range(nSteps):
            t = (n+1)*dt
            
            # get load function value
            T[n + 1, 0] = np.interp(t, xp, fp)
            T[n + 1,-1] = tbottom
            
            # get solution
            T[n+1,1:] = np.linalg.solve(Cmat[1:,1:]+beta*dt*Kmat[1:,1:],
                                        (Cmat[1:,:]-(1-beta)*dt*Kmat[1:,:])@T[n]-(Cmat[1:,0]+beta*dt*Kmat[1:,0])*T[n+1,0])
                                
# save important results
Tout = np.zeros((nSteps + 1, len(NN) + 3))
times = np.linspace(0,tFinal,nSteps+1)
Tout[:, 0 ] = times  # Time for temperature time plot versus sensors
Tout[:, 1 ] = T[:, 0]  # temperature at the top, prescribed
Tout[:, -1] = T[:, -1]  # temperature at the bottom, prescribed
for i in range(len(NN)):
    Tout[:, i + 2] = T[:, int(NNN[i])]
    
np.savetxt(fname="00_dat", X=Tout)

# plot results
plt.plot(times,T[:,0],label='x=0')
for i in NNN[:-1]:
     lab = str("Number " + str(int(i)))
     plt.plot(times, T[:, int(i)], label=lab)
plt.plot(times,T[:,-1],label='x=L')
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.legend()
plt.savefig("test.png", dpi=600)
out = np.zeros((nSteps+1,N+2))
out[:,1:] = T
out[:,0 ] = times
np.savetxt(fname=resfile, X=out)
out2 = np.zeros((nSteps+1,N+1))
out2[:,0 ] = times
out2[:,1:] = QQ
np.savetxt(fname="Q.dat", X=out2)



