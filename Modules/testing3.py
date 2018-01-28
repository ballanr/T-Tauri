import pandas as pd
from astropy.io import fits
import itertools
from PyAstronomy.pyasl import helcorr as helcorr
import numpy as np
import functions
import functions2
#from rpy2.robjects.packages import importr
#import rpy2.robjects as ro
#import pandas.rpy.common as com

filepathh = '/Users/ballanr/Desktop/File Outputs/DR15/Chip Gap Fixed/List.csv'
totalfile = pd.read_csv(filepathh)

#names = totalfile['']
i = 0
totals = len(totalfile['Filename'])

for index,row in itertools.islice(totalfile.iterrows(),0,None):
    i += 1
    print(str(i) + '/' + str(totals),row['Filename'])
    
    wave,flux = functions2.Chipgap_Fix(row['Filename'])
    
    df = pd.DataFrame(columns=['Wavelength','Flux'])
    df['Wavelength'] = wave
    df['Flux'] = flux
    df.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Chip Gap Fixed/' + str(row['Filename']),index=False)
    