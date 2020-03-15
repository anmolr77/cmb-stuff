import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.fftpack as fourrier

colors = [(1,0,0,c) for c in np.linspace(0,1,100)]
cmapred = mcolors.LinearSegmentedColormap.from_list('mycmap',colors,N=5)
#colors = [(0,0,1,c) for c in np.linspace(0,1,100)]

N = 11

x = np.linspace(0,1,N)
y = np.linspace(0,1,N)

xx, yy = np.meshgrid(x,y)


##################################

muT    = 2.7260
sigmaT = 6e-4
Tval = np.random.normal(muT,sigmaT,(N,N)) #creates a NxN matrix of grv's that are uncorrelated
print(Tval)

theta = (Tval - muT)/Tval

plt.scatter(xx,yy,s=1,color='black')
#plt.set_aspect('equal')
plt.contourf(xx,yy,Tval,cmap=cmapred)
#plt.colormesh(X,Y,Tval,cmap=cmapblue)
#plt.contourf(X,Y,uval,cmap='jet')
plt.title('Absolute Temperature')
plt.colorbar()
plt.show()

plt.scatter(xx,yy,s=1,color='black')
#plt.set_aspect('equal')
plt.contourf(xx,yy,theta,cmap=cmapred)
#plt.colormesh(X,Y,Tval,cmap=cmapblue)
#plt.contourf(X,Y,uval,cmap='jet')
plt.title('Absolute Temperature')
plt.colorbar()
plt.show()


###################################



def dist(i,j):
    return np.sqrt((xx-x[i])**2 + (yy-y[j])**2)

def tpCF(r):
    Ttmp=np.zeros((N,N))
    for i in range(N):
        for j in range(N):

            ctr=0
            for l in range(N):
                for m in range(N):

                    if r <= dist(i,j)[l][m] and dist(i,j)[l][m] <= r+dr:
                        Ttmp[i][j] = Ttmp[i][j] + theta[i][j]*Tval[l][m]
                        ctr = ctr + 1

            Ttmp[i][j] = Ttmp[i][j]/ctr

    return np.mean(Ttmp.flatten())

###################################

L = 50
dr = 0.08
rmax = np.sqrt(2)-dr
rvals = np.linspace(0,rmax,L)

tpCF_vals = np.zeros(L)

itr = 0
for r in rvals:
    tpCF_vals[itr] = tpCF(r)
    itr = itr+1


#tpCF_vals = tpCF_vals[~np.isfinite(tpCF_vals)]
print(tpCF_vals)

plt.plot(rvals,tpCF_vals)
plt.show()


#val = fourrier.ifft(tpCF_vals)
#plt.plot(rvals,np.imag(val))
#plt.show()
