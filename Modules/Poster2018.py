import functions2
import functions
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# filepath = '/Users/ballanr/Desktop/Research/DR15/Wave_and_Flux/6223-56283-099.csv'
# # filepath2 = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/6223/56283/apVisit-r8-6223-56283-099.fits'

# # openfile2 = fits.open(filepath2)
# # wave2 = openfile2[4].data
# # flux2 = openfile2[1].data
# # vbc = openfile2[0].header['BC']

# # print(vbc)
# openfile = pd.read_csv(filepath)

# wave = openfile['Wavelength']
# flux = openfile['Flux']

# a,b,rest = functions.Barycentric_Correction(11,0)

# rest = rest * (10**10)

# centerline = functions.find_nearest(wave,a)

# x,continuum,z = functions2.Br_EqW_Updated(wave,flux,11,centerline)

# L1 = wave[centerline - 301] # ~ 27.42 Angstroms
# L2 = wave[centerline - 150] # ~ 17.21 Angstroms
# R1 = wave[centerline + 150]
# R2 = wave[centerline + 301]


# plt.figure(figsize=(13,10))
# # plt.plot(wave2[0],flux2[0],color='red',linewidth=1.25,label='Original Spectra')
# # plt.plot(wave2[1],flux2[1],color='red',linewidth=1.25)
# # plt.plot(wave2[2],flux2[2],color='red',linewidth=1.25)
# plt.plot(wave,flux,linewidth=1)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.ylabel('Flux (10^-17 erg/s/cm^2/Ang)',fontsize=20)
# plt.xlabel('Wavelength (Ang)',fontsize=20)
# #plt.title('2M03434449+3143092 Spectra',fontsize=24)
# plt.ylim(400,1300)
# #plt.legend(fontsize=18)
# #plt.xlim(wave[centerline]-50,wave[centerline]+50)

# plt.savefig('/Users/ballanr/Desktop/Spectra.pdf',bbox_inches='tight',dpi=300)

#functions2.Updated_Spectra_Plotter('/Users/ballanr/Desktop/Research/DR15/Wave_and_Flux/6223-56283-099.csv')


functions2.multichi('6223-56283-099.csv',True)