import astropy.cosmology
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import camb
from camb import model, initialpower

# setting up a set of parameters for CAMB

pars = camb.CAMBparams()

# .CAMBparams sets up CosmoMC (Monte-Carlo) like setting, with one
#  massive neutrino and helium set using BBN consistency

pars.set_cosmology(H0=70,ombh2=0.06,omch2=0.122,omk=0.0)
pars.set_dark_energy(w=-1,cs2=1.0,wa=0,dark_energy_model='fluid')
pars.InitPower.set_params(As=2e-9,ns=1,r=0)
pars.set_for_lmax(2500,lens_potential_accuracy=2)

# calculate the results for the model

results = camb.get_results(pars)
powers = results.get_cmb_power_spectra(pars, CMB_unit = 'muK')

#plot the total lensed CMB power spectra versus unlensed, and fractional difference

totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print(totCL.shape)

#Python CL arrays are all zero based (starting at L=0), Note L=0,1 entries will be zero by default.
#The different CL are always in the order TT, EE, BB, TE (with BB=0 for unlensed scalar results).

ls = np.arange(totCL.shape[0])
#fig, ax = plt.subplots(2,2, figsize = (12,12))

'''
ax[0,0].plot(ls,totCL[:,0], color='k')
ax[0,0].plot(ls,unlensedCL[:,0], color='r')
ax[0,0].set_title('TT')
ax[0,0].set_xscale('log')
'''

plt.plot(ls,totCL[:,0], color='k')
plt.plot(ls,unlensedCL[:,0], color='r')
plt.title('TT')
plt.xscale('log')

'''
ax[0,1].plot(ls[2:], 1-unlensedCL[2:,0]/totCL[2:,0]);
ax[0,1].set_title(r'$\Delta TT$')
ax[1,0].plot(ls,totCL[:,1], color='k')
ax[1,0].plot(ls,unlensedCL[:,1], color='r')
ax[1,0].set_title(r'$EE$')
ax[1,1].plot(ls,totCL[:,3], color='k')
ax[1,1].plot(ls,unlensedCL[:,3], color='r')
ax[1,1].set_title(r'$TE$');
for ax in ax.reshape(-1): ax.set_xlim([2,2500])
'''

plt.show()


'''

for val in omega_range
plt diff graphs for the same set-up

compare the dtt values for different lensing potential lens_potential_accuracy (0,1,2)
plot (dtt(lpa1)-dtt(lpa2) 01 02 11 12)

'''
