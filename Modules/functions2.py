############### UPDATED FUNCTIONS ###############


def Br_EqW_Updated(wave,spec,line,centerline):
    
    import functions
    import functions2
    import numpy as np
    
    #regular windows

    L1 = centerline - 301 # ~ 27.42 Angstroms
    L2 = centerline - 150 # ~ 17.21 Angstroms
    R1 = centerline + 150
    R2 = centerline + 301


    Fluxcontinuum = (np.sum(spec[L1:L2])+np.sum(spec[R1:R2])) / (len(spec[L1:L2])+len(spec[R1:R2]))
    EqW1 = 0

    if Fluxcontinuum == 0:

        EqW1 = 0
        EqW1_rounded = 0
        equivs = 0
        
    if Fluxcontinuum != 0:
        

        for i in range(L2,R1):

            trapezoid = (0.5)*(wave[i+1] - wave[i])*(spec[i+1] + spec[i] - (2*Fluxcontinuum))
            EqW1 += trapezoid
        
    #EqW_rounded = round(EqW1/Fluxcontinuum,5)
        equivs = EqW1/Fluxcontinuum
    
    return equivs,Fluxcontinuum
def Br_Error_Updated(wave,flux,err,line,center):
    
    import functions
    import numpy as np

    try:
        if len(err) > 0:
            
            wave1 = np.asarray(wave[center-301:center+301])
            flux1 = np.asarray(flux[center-301:center+301])
            err1 = np.asarray(err[center-301:center+301])
            
            L1 = 0
            L2 = 151
            R1 = 451
            R2 = 602
            
            lcheck = np.median(err1[L1:L2])
            for k in range(L1,L2):
                if err1[k] > 10*lcheck:
                    err1[k] = 10*lcheck


            rcheck = np.median(err1[R1:R2])
            for k in range(R1,R2):
                if err1[k] > 10*rcheck:
                    err1[k] = 10*rcheck
            
            Fc = (np.sum(flux1[L1:L2]) + np.sum(flux1[R1:R2])) / (len(flux1[L1:L2]) + len(flux1[R1:R2]))
            dFc = (np.sum(err1[L1:L2]) + np.sum(err1[R1:R2])) / (len(err1[L1:L2]) + len(err1[R1:R2]))
            

            fluxavg = 0 
            
            for i in range(L2,R1):

                fluxavg += flux1[i] - Fc
                
            
            N = len(wave1[L2:R1])
            fluxavg = fluxavg / N
            del_lam = wave1[2] - wave1[1]    
            A = N*del_lam*fluxavg

            dA = 0
            for k in range(L2,R1):


                dA += err[k]**2

            dA = del_lam * np.sqrt(dA)

            dEqW = (1/(Fc**2)) * np.sqrt( (Fc*dA)**2 + (A*dFc)**2 )
            
        else:
            dEqW = 0
    except IndexError:
        dEqW = 0

    return dEqW
def Model_Fitter_Updated(eqwarray,error,length):
    import pandas as pd
    import numpy as np

    allfile = '/Users/ballanr/Desktop/File Outputs/Brackett Decrements/Profile Test.csv'
    openallfile = pd.read_csv(allfile)
    cols = openallfile.columns
    eqwarray = np.asarray(eqwarray)
    eqwarray= eqwarray/eqwarray[0]
    errors = []
    test = [(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0)]
    

    for i in range(length):

        dq = eqwarray[i] * np.sqrt( (error[i]/eqwarray[i])**2 + (error[0]/eqwarray[0])**2 )
        errors.append(dq)

    resid = []
    for i in range(len(cols)):
        
        A = openallfile[cols[i]]
        
        b = eqwarray
        
        A = A / A[0]
        
        r = b - A
        
        r2 = 0
        
        for k in range(length):
            w = 1/(errors[k])**2
            r2 += w * (r[k])**2
            test[k]=(k,w)
        
        r2 = np.sqrt(r2)
        resid.append((i,r2))
    
    x = min(b for (a,b) in resid)
    y = np.where(resid == x)
    model = cols[y[0][0]]
    
    
    if model.startswith('1'):
        dens = float(model[:4])
        temp = int(model[5:])

    else:
        dens = float(model[:3])
        temp = int(model[4:])

    return dens,temp
def Brackett_Ratios_Updated_Grid(lines):


    import matplotlib as mpl
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import itertools
    import seaborn as sb
    
    lines = str(lines)
    filepath = '/Users/ballanr/Desktop/File Outputs/Currently Working On/Model Fit Update '+lines+'s.csv'
    filefile = pd.read_csv(filepath)

    counter = 0 
    mpl.rcParams.update(mpl.rcParamsDefault)

    for index,row in itertools.islice(filefile.iterrows(),0,1):
        counter += 1
        equivs = []
        errors = []

        plateid = int(row['Plate ID'])
        MJD = int(row['MJD'])
        modeldens = str(row['Model Density'])
        modelt = str(row['Model Temp'])
        ffiber = str(row['Fiber'])
        
        if len(str(ffiber)) == 3:
            ffiber = str(ffiber)
        elif len(str(ffiber)) == 2:
            ffiber = '0' + str(ffiber) 
        else:
            ffiber = '00' + str(ffiber)

        for j in range(10):
            number = 11+j
            errorstring = 'Br' + str(number) + ' Error'
            equivstring = 'Br' + str(number) + ' EqW'
            hh = row[errorstring]
            gg = row[equivstring]
            equivs.append(gg)
            errors.append(hh)

        savestring = str(plateid) + '-' + str(MJD) + '-' + str(ffiber) + '.pdf'

        equivs = np.asarray(equivs)
        errors = np.asarray(errors)
        equivs1 = equivs/equivs[0]

        summederrors = []

        for i in range(10):

            dq = equivs1[i] * np.sqrt( (errors[i]/equivs[i])**2 + (errors[0]/equivs[0])**2 )
            summederrors.append(dq)

        temps = ['3750 K','5000 K','7500 K','8750 K','10000 K','12500 K','15000 K']

        d8T3750 = np.asarray([0.05071,0.04067,0.03307,0.02709,0.02237,0.01835,0.01527,0.01281,0.01081,0.009207])
        d8T8750 = np.asarray([0.04492,0.03533,0.02837,0.02307,0.019,0.01573,0.01313,0.01105,0.009374,0.00801])
        d8T15000 = np.asarray([0.02642,0.02042,0.01618,0.01303,0.01067,0.008821,0.007353,0.006188,0.005252,0.004492])
        d10T3750 = np.asarray([0.03951,0.03036,0.02369,0.01874,0.01503,0.01222,0.01006,0.008374,0.007042,0.005979])
        d10T8750 = np.asarray([0.02256,0.01705,0.01321,0.01044,0.008401,0.006864,0.005682,0.00476,0.004027,0.003437])
        d10T15000 = np.asarray([0.01032,0.00755,0.00576,0.004524,0.003634,0.00297,0.002463,0.002068,0.001753,0.001499])
        d12T3750 = np.asarray([0.2249,0.1888,0.1567,0.1297,0.1077,0.09002,0.07573,0.06418,0.05477,0.04704])
        d12T8750 = np.asarray([0.9622,1.027,1.079,1.12,1.149,1.162,1.159,1.139,1.104,1.058])
        d12T15000 = np.asarray([1.082,1.167,1.23,1.267,1.274,1.253,1.208,1.147,1.076,0.9992])

        filename = '/Users/ballanr/Desktop/File Outputs/Brackett Decrements/Density Files/Density ' + modeldens + ' Ratios.csv'
        openfile = pd.read_csv(filename)

        y = openfile[modelt+ ' K']
        y = np.asarray(y)
        y = y/y[0]

        plt.figure(figsize=(20,10))
        plt.title(str(plateid)+'-'+str(MJD)+'-'+ffiber,fontsize=24)
        plt.xticks((np.arange(11,21,1)))
        
        nonfits = 0.4

        plt.plot(np.arange(11,21,1),d8T3750/d8T3750[0],color=sb.xkcd_rgb['red'],alpha=nonfits,ls='dashed',label='Density: 8 Temp: 3750')
        plt.plot(np.arange(11,21,1),d8T8750/d8T8750[0],color=sb.xkcd_rgb['red'],alpha=nonfits,ls='dashdot',label='Density: 8 Temp: 8750')
        plt.plot(np.arange(11,21,1),d8T15000/d8T15000[0],color=sb.xkcd_rgb['red'],alpha=nonfits,ls='solid',label='Density: 8 Temp: 15000')
        plt.plot(np.arange(11,21,1),d10T3750/d10T3750[0],color=sb.xkcd_rgb['forest green'],alpha=nonfits,ls='dashed',label='Density: 10 Temp: 3750')
        plt.plot(np.arange(11,21,1),d10T8750/d10T8750[0],color=sb.xkcd_rgb['forest green'],alpha=nonfits,ls='dashdot',label='Density: 10 Temp: 8750')
        plt.plot(np.arange(11,21,1),d10T15000/d10T15000[0],color=sb.xkcd_rgb['forest green'],alpha=nonfits,ls='solid',label='Density: 10 Temp: 15000')
        plt.plot(np.arange(11,21,1),d12T3750/d12T3750[0],color=sb.xkcd_rgb['purple'],alpha=nonfits,ls='dashed',label='Density: 12 Temp: 3750')
        plt.plot(np.arange(11,21,1),d12T8750/d12T8750[0],color=sb.xkcd_rgb['purple'],alpha=nonfits,ls='dashdot',label='Density: 12 Temp: 8750')
        plt.plot(np.arange(11,21,1),d12T15000/d12T15000[0],color=sb.xkcd_rgb['purple'],alpha=nonfits,ls='solid',label='Density: 12 Temp: 15000')

        plt.plot(np.arange(11,21,1),y,color=sb.xkcd_rgb['black'],label='Best Fit, Density: '+str(modeldens)+' Temp: '+str(modelt))
        plt.scatter(np.arange(11,21,1),y,color=sb.xkcd_rgb['black'])
        
        if int(lines)<20:
            plt.errorbar(np.arange(int(lines),21,1),equivs1[int(lines)-11:],yerr=summederrors[int(lines)-11:],ls='dashed',ecolor='red',barsabove=True,capsize=3,fmt='-o',color='red')
        
        plt.errorbar(np.arange(11,int(lines)+1,1),equivs1[:int(lines)-10],yerr=summederrors[:int(lines)-10],ecolor='red',barsabove=True,capsize=3,fmt='-o',color=sb.xkcd_rgb['blue'],label='Br Emission')

        plt.grid(True)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=4,frameon=True,facecolor=None,framealpha=1,fontsize=16)
        plt.savefig('/Users/ballanr/Desktop/File Outputs/Currently Working On/Decrement Plots/'+lines+'/' + savestring,bbox_inches='tight',dpi=300)
        plt.clf()
        plt.close()
        print('Plotted '+ str(counter)+' fits')
def New_All_Model_Fitter():

    import functions
    import functions2
    import pandas as pd
    import itertools
    from astropy.io import fits
    from PyAstronomy.pyasl import helcorr

    cols = ['Location ID','2Mass ID', 'Plate ID','MJD','Fiber','Model Density','Model Temp','Br11 EqW',
            'Br12 EqW','Br13 EqW','Br14 EqW','Br15 EqW','Br16 EqW','Br17 EqW','Br18 EqW','Br19 EqW',
            'Br20 EqW','Br11 Error','Br12 Error','Br13 Error','Br14 Error','Br15 Error','Br16 Error',
            'Br17 Error','Br18 Error','Br19 Error','Br20 Error',]

    df = pd.DataFrame(columns = cols)

    filepath = '/Users/ballanr/Desktop/File Outputs/Currently Working On/Test List.csv'
    openfile = pd.read_csv(filepath)

    for index,row in itertools.islice(openfile.iterrows(),0,None):

        ##### GETTING VBC #####
        loc = str(row['Location ID'])
        twomass = str(row['2Mass ID'])
        plate = str(row['Plate'])
        mjd = str(row['MJD'])
        fiber = str(row['Fiber'])
        if len(fiber) == 1:
            fiber = '00' + fiber
        elif len(fiber) == 2:
            fiber = '0' + fiber
        else:
            fiber = fiber

        if int(plate) < 9700:
            
            try:
                server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                filepath = str(plate) + '/' + str(mjd) + '/apVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
                fitsfile = server + filepath
                openfits = fits.open(fitsfile)
                header = openfits[0].header
            except:
                server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                filepath = str(plate) + '/' + str(mjd) + '/apVisit-r8-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
                fitsfile = server + filepath
                openfits = fits.open(fitsfile)
                header = openfits[0].header

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
            server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/lco25m/'
            filepath = str(plate) + '/' + str(mjd) + '/asVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
            
            fitsfile = server + filepath
            openfits = fits.open(fitsfile)

            header = openfits[0].header

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
            

        ##### IMPORTING WAVELENGTH, FLUXES, AND ERRORS #####
        csvfile = '/Users/ballanr/Desktop/File Outputs/DR15/Wave and Flux/'+plate+'-'+mjd+'-'+fiber+'.csv'
        opencsv = pd.read_csv(csvfile)
        wave = opencsv['Wavelength']
        flux = opencsv['Flux']
        errors = opencsv['Error']

        #### CALCULATIONS #####
        lines = [11,12,13,14,15,16,17,18,19,20]
        equivs = []
        error = []

        for i in range(len(lines)):
            observedw,shift,restw = functions.Barycentric_Correction(lines[i],vbc)
            rest = restw*(10**10)
            centerline = functions.find_nearest(wave,rest)
            equiv,fc = functions2.Br_EqW_Updated(wave,flux,lines[i],centerline)
            errs = functions2.Br_Error_Updated(wave,flux,errors,lines[i],centerline)
            equivs.append(equiv)
            error.append(errs)

        ##### MODEL FITTING #####
        dens,temp = functions2.Model_Fitter_Updated(equivs,error,7)

        ##### OUTPUT #####
        data = [loc,twomass,plate,mjd,fiber,dens,temp,equivs[0],equivs[1],equivs[2],equivs[3],
                equivs[4],equivs[5],equivs[6],equivs[7],equivs[8],equivs[9],error[0],error[1],
                error[2],error[3],error[4],error[5],error[6],error[7],error[8],error[9]]

        df.loc[len(df)+1] = data

    df.to_csv('/Users/ballanr/Desktop/File Outputs/Currently Working On/Model Fit Update 17s.csv',index=False)
    df = df.iloc[0:0]
def False_Spectra():


    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import itertools
    from astropy.io import fits

    ##### Get Wavegrid #####



    ##### Generate Random Spectra #####

    mu = 16811
    sigma = 2
    wave = np.linspace(16790, 16830,150)

    flux1 = gauss(wave,mu+2,sigma)
    flux2 = gauss(wave,mu-2,sigma)
    flux3 = (flux1+flux2)
    flux4 = 0.5*(flux1+flux2)

    plt.plot(wave,flux1)
    plt.plot(wave,flux2)
    plt.plot(wave,flux3)
    plt.plot(wave,flux4)
    plt.show()

def gauss(x,mu,sigma): 
    import numpy as np
    return (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mu)/sigma)**2)

