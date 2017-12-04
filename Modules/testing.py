import functions
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from PyAstronomy.pyasl import helcorr
import numpy as np
import itertools
import functions

filepath = '/Users/ballanr/Desktop/Confidence Plots/Confidence2.csv'
openfile = pd.read_csv(filepath)
counter = 0
for index,row in itertools.islice(openfile.iterrows(),0,None):
    counter += 1
    plate = row['Plate']
    mjd = row['MJD']
    if len(str(row['Fiber'])) == 3:
        fiber = str(row['Fiber'])
    elif len(str(row['Fiber'])) == 2:
        fiber = '0' + str(row['Fiber']) 
    else:
        fiber = '00' + str(row['Fiber'])

    if int(plate) < 9700:
        try:
            serverpath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/' 
            fitspath = serverpath + str(plate) + '/' + str(mjd) + '/apVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
            openfile = fits.open(fitspath)
            header = openfile[0].header
        except:
            serverpath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/' 
            fitspath = serverpath + str(plate) + '/' + str(mjd) + '/apVisit-r8-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
            openfile = fits.open(fitspath)
            header = openfile[0].header

        try:
            vbc = header['VHELIO']
        except:
            ra = header['RA']
            dec = header['DEC'] 
            jd = header['JD-MID']
                            
            height = 2788
            longg = -105.4913
            lat = 36.4649
                            
            vbc,hjd = helcorr(longg,lat,height,ra,dec,jd)
    
    else: 
        serverpath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/lco25m/' 
        fitspath = serverpath + str(plate) + '/' + str(mjd) + '/asVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
        openfile = fits.open(fitspath)
        header = openfile[0].header
        
        try:

            vbc = header['VHELIO']

        except:

            ra = header['RA']
            dec = header['DEC'] 
            jd = header['JD-MID']
                    
            height = 2380
            longg = -70.413336
            lat = -29.05256
                    
            vbc,hjd = helcorr(longg,lat,height,ra,dec,jd)

    savestring = str(plate) + '-' + str(mjd) + '-' + str(fiber) +'.pdf'

    try:
        
        wave0 = openfile[4].data[0]
        wave1 = openfile[4].data[1]
        wave2 = openfile[4].data[2]

        flux0 = openfile[1].data[0]
        flux1 = openfile[1].data[1]
        flux2 = openfile[1].data[2]

                
        c = 299792.458
        lamshift = 1 + (vbc/c)

        wave0 = np.asarray(wave0) * lamshift
        wave1 = np.asarray(wave1) * lamshift
        wave2 = np.asarray(wave2) * lamshift

        leftwindow = functions.find_nearest(wave0,16770)
        rightwindow = functions.find_nearest(wave0,16850)
        #print(leftwindow,rightwindow)
        fluxmax = max(flux0[rightwindow:leftwindow])
        fluxmin = min(flux0[rightwindow:leftwindow])
        #print(fluxmin,fluxmax)

        if fluxmin > 0:
            fluxmin = -0.5*fluxmin
        else:
            fluxmin = 1.25*fluxmin
        print(counter)
        plt.figure(figsize=(20,10))
        plt.plot(wave0,flux0)
        plt.plot(wave1,flux1)
        plt.plot(wave2,flux2)

        plt.axvline(16811,ls='dashed',color='red')
        plt.xlim(16770,16850)
        plt.ylim(fluxmin,1.25*fluxmax)
        plt.savefig('/Users/ballanr/Desktop/Confidence Plots/Confidence2 Plots/'+savestring,bbox_inches='tight',dpi=300)
        plt.clf()
        plt.close()
    except:
        print('Fiber ' + str(fiber) + ' has no file...')

'''
mjd = [57790,57792,57793]

for j in range(3):

    filepath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/9659/' + str(mjd[j]) + '/apVisit-apogee2-9659-' + str(mjd[j])
    #print(filepath)

    for i in range(301):
        if i != 0:
            if i < 10:
                fiber = '00' + str(i)
            elif i > 9 and i < 100:
                fiber = '0' + str(i)
            else:
                fiber = str(i)

            fitsfile = filepath + '-' + str(fiber) + '.fits'
            #print(fitsfile)
            try:
                openfile = fits.open(fitsfile)
                header = openfile[0].header

                wave0 = openfile[4].data[0]
                wave1 = openfile[4].data[1]
                wave2 = openfile[4].data[2]

                flux0 = openfile[1].data[0]
                flux1 = openfile[1].data[1]
                flux2 = openfile[1].data[2]

                try:
                    vbc = header['VHELIO']
                except:
                    ra = header['RA']
                    dec = header['DEC'] 
                    jd = header['JD-MID']
                            
                    height = 2788
                    longg = -105.4913
                    lat = 36.4649
                            
                    vbc,hjd = helcorr(longg,lat,height,ra,dec,jd)
                
                c = 299792.458
                lamshift = 1 + (vbc/c)

                wave0 = np.asarray(wave0) * lamshift
                wave1 = np.asarray(wave1) * lamshift
                wave2 = np.asarray(wave2) * lamshift

                plt.figure(figsize=(20,10))
                plt.plot(wave0,flux0)
                plt.plot(wave1,flux1)
                plt.plot(wave2,flux2)

                plt.axvline(16811,ls='dashed',color='red')
                plt.xlim(16770,16850)
                plt.savefig('/Users/ballanr/Desktop/Confidence Plots/9659/9659-' + str(mjd[j]) + '-' + str(fiber) + '.pdf',bbox_inches='tight',dpi=300)
                plt.clf()
                plt.close()
            except:
                print('Fiber ' + str(fiber) + ' has no file...')'''