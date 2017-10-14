import functions
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits


filepath = '/Users/ballanr/Desktop/JessicasDwarves/6531-56254-233.csv'
openfile = pd.read_csv(filepath)

wave = openfile['Wavelength']
flux = openfile['Flux']


plt.figure(figsize=(20,10))
plt.plot(wave,flux)

#plt.xlim(16000,16100)
plt.show()


'''filepath = '/Volumes/CoveyData/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/apo25m/6531/56254/apVisit-r8-6531-56254-233.fits'
openfile = fits.open(filepath)

flux = openfile[1].data
wave = openfile[4].data

plt.figure(figsize=(20,10))
plt.plot(wave[0],flux[0])
plt.plot(wave[1],flux[1])
plt.plot(wave[2],flux[2])
plt.show()
'''
#plt.savefig('/Users/ballanr/Desktop/Server.pdf',bbox_inches='tight',dpi=300)

