import functions
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np

fitsfile = '/Users/ballanr/Desktop/SummerResearch/dr13/apogee/spectro/redux/r6/apo25m/6218/56172/apVisit-r6-6218-56172-094.fits'
openfile = fits.open(fitsfile)
vbc = openfile[0].header['BC']
c = 299792.458
shift = 1 + (vbc/c)
        
fspec = openfile[1]
fwave = openfile[4]
wave = []
flux = []
        
for i in range(len(fwave.data[2])):
    wave.append(fwave.data[2][-i-1])
    flux.append(fspec.data[2][-i-1])
for i in range(len(fwave.data[1])):
    wave.append(fwave.data[1][-i-1])
    flux.append(fspec.data[1][-i-1])
for i in range(len(fwave.data[0])):
    wave.append(fwave.data[0][-i-1])
    flux.append(fspec.data[0][-i-1])
#plt.plot(wave,flux,linewidth=1,color='red',label='Original Spectra')


y = functions.OH_Skylines_Remover(wave,flux)
g = functions.Br_EqW(wave,y,11,vbc)
wave = np.asarray(wave) * shift
gg = g[3]*(10**10)
plt.plot(wave,y,linewidth=1,color = 'blue',label='Removed Skylines')
plt.legend()
plt.xlim(16770,16850)
plt.axhline(618.3127,ls = 'dashed',color='black')
plt.axvline(gg,ls='dashed',color='red')
#plt.ylim(200,600)
plt.show()
'''
for i in range(10):

    g = functions.Br_EqW(wave,y,11+i,vbc)
    print(g)'''