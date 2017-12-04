import pandas as pd
from astropy.io import fits
import itertools
from PyAstronomy.pyasl import helcorr as helcorr
import numpy as np
import functions

'''
filepath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/5534/55847/apVisit-r8-5534-55847-002.fits'
openfile = fits.open(filepath)

header = openfile[0].header

print(header['RVTEFF'])
'''

cols = ['Location ID','2Mass ID','Plate','MJD','Fiber','Confidence','Br11','TEFF','J MAG','H MAG','K MAG','J-H','H-K']
df = pd.DataFrame(columns=cols)

filepath = '/Users/ballanr/Desktop/File Outputs/DR15/forMdot2.csv'
openfile = pd.read_csv(filepath)

g = 0 

for index,row in itertools.islice(openfile.iterrows(),0,None):

    g += 1
    print(g)

    loc = row['Location ID']
    twomass = row['2Mass ID']
    plate = row['Plate']
    mjd = row['MJD']
    fiber = row['Fiber']

    if len(str(fiber)) == 1:
        fiber = '00' + str(fiber) 
    if len(str(fiber)) == 2:
        fiber = '0' + str(fiber)
    if len(str(fiber)) == 3:
        fiber = str(fiber)

    if int(plate) < 8870 :

        fitspath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'+str(plate)+'/'+str(mjd)+'/apVisit-r8-'+str(plate)+'-'+str(mjd)+'-'+str(fiber)+'.fits'
            
    elif int(plate) > 8870 and int(plate) < 9700:

        fitspath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'+str(plate)+'/'+str(mjd)+'/apVisit-apogee2-'+str(plate)+'-'+str(mjd)+'-'+str(fiber)+'.fits'
            
    else:
        fitspath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/lco25m/'+str(plate)+'/'+str(mjd)+'/asVisit-apogee2-'+str(plate)+'-'+str(mjd)+'-'+str(fiber)+'.fits'
                
    fitsvisit = fits.open(fitspath)    
                
    header = fitsvisit[0].header

    J = header['J']
    H = header['H']
    K = header['K']

    J_H = float(J) - float(H)
    H_K = float(H) - float(K)

    try:
        
        teff = header['RVTEFF']

    except:

        teff = 0

    try:

        vbc = header['VHELIO']

    except:

        ra = header['RA']
        dec = header['DEC']
        jd = header['JD-MID']
        telescope = header['TELESCOP']
                    
        if telescope == 'apo25m':
            height = 2788
            longg = -105.4913
            lat = 36.4649
                
        else:
            height = 2380
            longg = -70.413336
            lat = -29.05256

        vbc,hjd = helcorr(longg,lat,height,ra,dec,jd)

    c = 299792.458
    lamshift = 1 + (vbc/c)
                
    fspec = fitsvisit[1]
    ferr = fitsvisit[2]
    fwave = fitsvisit[4]
    wave = []
    flux = []
    error = []
                
    for i in range(len(fwave.data[2])):
        wave.append(fwave.data[2][-i-1])
        flux.append(fspec.data[2][-i-1])
        error.append(ferr.data[2][-i-1])
    for j in range(len(fwave.data[1])):
        wave.append(fwave.data[1][-j-1])
        flux.append(fspec.data[1][-j-1])
        error.append(ferr.data[2][-j-1])
    for k in range(len(fwave.data[0])):
        wave.append(fwave.data[0][-k-1])
        flux.append(fspec.data[0][-k-1])
        error.append(ferr.data[2][-k-1])
                
    fitsvisit.close()

    newflux = functions.skylines_cleaner(wave,flux)
                
    #now we run equiv width calc

    lines = [11,12,13,14,15,16,17,18,19,20]
       
    wave = np.asarray(wave) * lamshift #check which direction shift is going            

    equiv_width,fcontinuum,shift,rest_wavelength,centers = functions.Br_EqW(wave,newflux,lines[0],vbc)
    rest = rest_wavelength*(10**10)

    Confidence = functions.Confidence_Level(wave,newflux,rest)

    data = [int(loc),twomass,int(plate),int(mjd),fiber,Confidence,equiv_width,int(teff),J,H,K,J_H,H_K]

    df.loc[len(df)+1] = data

df.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/somelist1.csv',index=False)
df = df.iloc[0:0]