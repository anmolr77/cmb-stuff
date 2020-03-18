import sys
import camb
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import astropy.io.fits as fits

## creates a parameter class
pars = camb.CAMBparams()

pars.set_cosmology(H0=70,TCMB=2.7255,ombh2=0.0226,omch2=0.112,omk=0.0,YHe=0.24,tau=0.09)
pars.set_dark_energy(w=-1,cs2=1,wa=0,dark_energy_model='fluid')
pars.set_for_lmax(2500,lens_potential_accuracy=2)
pars.InitPower.set_params(As=2e-9,ns=1,r=0)

## results gives the class structure that we just created, unmark print() to see how it looks like
results = camb.get_results(pars)
#print(results)

## powers is the dictionary containing the various types of pwerspectrums
## i.e total, unlensed scalar, unlensed tensor etc, unmark print() to see how it looks like
powers = results.get_cmb_power_spectra(pars, CMB_unit='muK')
#print(powers)

totCL = powers['total']
unlensedCL = powers['unlensed_scalar']
#print(unlensedCL)

DlTT = totCL[:,0]
#print(ClTT)
ls = np.arange(totCL.shape[0])
#print(ls)

'''
## me being a degenerate (pt 1)

lsize = np.size(ls)
DlTT = np.zeros(lsize)

for i in range(lsize):
    DlTT[i] = 2*np.pi*ClTT[i] / (ls[i]*(ls[i]+1))
'''

'''
ClTT = (( 2*np.pi*DlTT ) / (ls*(ls+1)))

plt.plot(ls,ClTT)
#plt.xscale('log')  ## to make it similar to Wayne Hu's plots
plt.title('TT')
plt.xlabel(r'$l$')
plt.ylabel(r'$C_{l}(\mu K)$')
plt.show()
'''

plt.plot(ls,DlTT)
plt.title('TT')
plt.xlabel(r'$l$')
plt.ylabel(r'$D_{l}(\mu K)$')
plt.show()


## variables to set up the size of the map
N = 2**10  # this is the number of pixels in a linear dimension
            ## since we are using lots of FFTs this should be a factor of 2^N
pix_size  = 0.5 # size of a pixel in arcminutes

## variables to set up the map plots
c_min = -400  # minimum for color bar
c_max = 400   # maximum for color bar
X_width = N*pix_size/60.  # horizontal map width in degrees
Y_width = N*pix_size/60.  # vertical map width in degrees


def make_CMB_T_map(N,pix_size,ls,DlTT):
    "makes a realization of a simulated CMB sky map given an input DlTT as a function of ls,"
    "the pixel size (pix_size) required and the number N of pixels in the linear dimension."
    #np.random.seed(100)
    # convert Dl to Cl
    ClTT = DlTT * 2 * np.pi / (ls*(ls+1.))
    ClTT[0] = 0. # set the monopole and the dipole of the Cl spectrum to zero
    ClTT[1] = 0.

    # make a 2D real space coordinate system
    onesvec = np.ones(N)
    inds  = (np.arange(N)+.5 - N/2.) /(N-1.) # create an array of size N between -0.5 and +0.5
    # compute the outer product matrix: X[i, j] = onesvec[i] * inds[j] for i,j
    # in range(N), which is just N rows copies of inds - for the x dimension
    X = np.outer(onesvec,inds)
    # compute the transpose for the y dimension
    Y = np.transpose(X)
    # radial component R
    R = np.sqrt(X**2. + Y**2.)

    # now make a 2D CMB power spectrum
    pix_to_rad = (pix_size/60. * np.pi/180.) # going from pix_size in arcmins to degrees and then degrees to radians
    ls_scale_factor = 2. * np.pi /pix_to_rad  # now relating the angular size in radians to multipoles
    ls2d = R * ls_scale_factor # making a fourier space analogue to the real space R vector
    ClTT_expanded = np.zeros(int(ls2d.max())+1)
    # making an expanded Cl spectrum (of zeros) that goes all the way to the size of the 2D ls vector
    ClTT_expanded[0:(ClTT.size)] = ClTT # fill in the Cls until the max of the ClTT vector

    # the 2D Cl spectrum is defined on the multiple vector set by the pixel scale
    CLTT2d = ClTT_expanded[ls2d.astype(int)]
    #plt.imshow(np.log(CLTT2d))


    # now make a realization of the CMB with the given power spectrum in real space
    random_array_for_T = np.random.normal(0,1,(N,N))
    FT_random_array_for_T = np.fft.fft2(random_array_for_T)   # take FFT since we are in Fourier space

    FT_2d = np.sqrt(CLTT2d) * FT_random_array_for_T # we take the sqrt since the power spectrum is T^2
    plt.imshow(np.real(FT_2d))


    ## make a plot of the 2D cmb simulated map in Fourier space, note the x and y axis labels need to be fixed
    #Plot_CMB_Map(np.real(np.conj(FT_2d)*FT_2d*ls2d * (ls2d+1)/2/np.pi),0,np.max(np.conj(FT_2d)*FT_2d*ls2d * (ls2d+1)/2/np.pi),ls2d.max(),ls2d.max())  ###

    # move back from ls space to real space
    CMB_T = np.fft.ifft2(np.fft.fftshift(FT_2d))
    # move back to pixel space for the map
    CMB_T = CMB_T/(pix_size /60.* np.pi/180.)
    # we only want to plot the real component
    CMB_T = np.real(CMB_T)

    ## return the map
    return(CMB_T)
  ###############################

def Plot_CMB_Map(Map_to_Plot,c_min,c_max,X_width,Y_width):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    print("map mean:",np.mean(Map_to_Plot),"map rms:",np.std(Map_to_Plot))
    plt.gcf().set_size_inches(10, 10)
    im = plt.imshow(Map_to_Plot, interpolation='bilinear', origin='lower',cmap=cm.RdBu_r)
    im.set_clim(c_min,c_max)
    ax=plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = plt.colorbar(im, cax=cax)
    #cbar = plt.colorbar()
    im.set_extent([0,X_width,0,Y_width])
    plt.ylabel('angle $[^\circ]$')
    plt.xlabel('angle $[^\circ]$')
    cbar.set_label('tempearture [uK]', rotation=270)

    plt.show()
    return(0)
  ###############################

## make a CMB T map
CMB_T = make_CMB_T_map(N,pix_size,ls,DlTT)
Plot_CMB_Map(CMB_T,c_min,c_max,X_width,Y_width)
plt.savefig('cmb1.png')
plt.clf()
