import functions
import functions2
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd



#functions2.Updated_Spectra_Plotter('/Users/ballanr/Desktop/Testgrid.csv')
#functions2.Updated_Spectra_Plotter()

filepath = '/Users/ballanr/Desktop/Testgrid.csv'
openfile = pd.read_csv(filepath)

wave = openfile['Wavelength']
flux = openfile['Flux']


blah = np.arange(4091,4110)
wave = np.asarray(wave)
flux = np.asarray(flux)
#print(blah,len(wave))

wave1 = np.delete(wave,blah,axis=0)
flux1 = np.delete(flux,blah,axis=0)
#print(wave1[4090:4110])

avg = 0.13207440467
waveright = wave1[4091]
waveleft = wave1[4090]
fluxright = flux1[4091]
fluxleft = flux1[4090]

slope = (fluxright-fluxleft)/(waveright-waveleft)


for i in range(377):

    i += 1
    value = waveleft + i*avg

    fluxval = slope*(i*avg) + fluxleft

    wave1 = np.insert(wave1,4090+i,value)

    flux1 = np.insert(flux1,4090+i,fluxval)

df = pd.DataFrame(columns=['Wavelength','Flux'])
df['Wavelength'] = wave1
df['Flux'] = flux1
df.to_csv('/Users/ballanr/Desktop/Testgrid2.csv',index=False)


