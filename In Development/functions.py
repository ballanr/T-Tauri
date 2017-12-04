def Directory_Walk():

    import os
    import numpy as np
    import pandas as pd

    directory = '/Volumes/CoveyData/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/apo25m/'
    plates = []
    mjd = []
    fibers = []
    path = []

    df = pd.DataFrame(columns=['Plate','MJD','Fiber','Path'])

    for root,dirs,files in os.walk(directory):
        for name in files:
            #x = os.path.join(root,name)
            if len(name) == 30:
                plates.append(name[11:15])
                mjd.append(name[16:21])
                fibers.append(name[22:25])
                path.append(name)
                print(name)
            
            elif len(name) == 31:
                plates.append(name[11:16])
                mjd.append(name[17:22])
                fibers.append(name[23:26])
                path.append(name)
                print(name)
    
    df['Plate'] = plates
    df['MJD'] = mjd
    df['Fiber'] = fibers
    df['Path'] = path

    df.to_csv('/Users/ballanr/Desktop/File Outputs/DR14/DR14 List.csv',index=False)
def apVisit_Catalog_Output(filename,savefile):
    '''
    Notes:
    '''

    import functions
    import csv

    #Setting up Dictionary----------------------------------------------------------------------------|
    apVisit_Dict = {} #Not used right now. Will shift to it eventually.

    #Creating List------------------------------------------------------------------------------------|
    master = []
    with open(filename) as csvfile:

        reader = csv.reader(csvfile,delimiter=',')
        j = 0
        for row in reader:
            loc_id = int(row[1])
            twomass_id = '%s' % row[0]
            #print(loc_id,twomass_id,type(twomass_id))
            x = functions.apStar_to_apVisit(loc_id,twomass_id)
            for i in range(len(x)):
                plate = x[i][0]
                mjd = x[i][1]
                fiber = x[i][2]
                master.append((loc_id,twomass_id,plate,mjd,fiber))
                j += 1
                print('Appending %s to master' %j)
            
    with open(savefile,'w') as savefile:
        
        writer = csv.writer(savefile,delimiter = '\t')
        writer.writerow(('Location ID','2Mass ID','Plate','MJD','Fiber'))
        for i in range(len(master)):
            writer.writerow((master[i][0],master[i][1],
                            master[i][2],master[i][3],master[i][4]))
            print('Writing row %s' %i)
def find_nearest(array,value):
    
    import numpy as np

    array = np.asarray(array)
    index = (np.abs(array-value)).argmin()
    return index
def apStar_to_apVisit(locid,twomassid):
    import apogee.tools.read as apread

    header = apread.apStar(locid,twomassid,ext=0,header=True)

    visits = header[1]['NVISITS']
    array = []
    for i in range(visits):
        x = i+1
        SFILE = 'SFILE%s' % x
        plate = header[1][SFILE][11:15]
        MJD = header[1][SFILE][16:21]
        fiber = header[1][SFILE][22:25]
        array.append((int(plate),int(MJD),fiber,int(visits)))
    return array,locid,twomassid
def Barycentric_Correction(emission_line,vbc):

    #Constants---------------------------------------------------------------------------------------------------------------------------------------|
    
    n = (float(emission_line))**2 #Beginning electron level
    c = 299792.458 #Speed of light (km/s)
    rydberg_inf =  1.0973731568539*(10**7) #Rydberg constant (m^-1)
    electron = 9.10938356*(10**-31) #Mass of electron (kg)
    nucleus = 1.672621898*(10**-27) #Mass of hydrogen nucleus (kg)
    
    #Equations---------------------------------------------------------------------------------------------------------------------------------------|
    
    rydberg_red = rydberg_inf/(1+(electron/nucleus)) #Reduced mass Rydberg constant (m^-1)
    rest_wavelength1 = rydberg_red*((1./16.)-(1./n)) #(m^-1)
    rest_wavelength = 1/rest_wavelength1 #wavelength (m)
    observed_wavelength1 = rest_wavelength*(1-(vbc/c)) #Finding the location of the peak of the observed spectrum
    observed_wavelength = observed_wavelength1*(10**10)
    shift = (rest_wavelength-observed_wavelength1)*(10**10) #Finding the difference between rest and observed wavelengths

    #Returns-----------------------------------------------------------------------------------------------------------------------------------------|

    return observed_wavelength,shift,rest_wavelength
def Max_Flux_Check(x_axis,y_axis,centerline):
    import functions
    c = 299792.458
    v_window = 100
    right_shift = x_axis[centerline]*(1+(v_window/c))
    left_shift = x_axis[centerline]*(1-(v_window/c))
    leftwindow = functions.find_nearest(x_axis,left_shift)
    rightwindow = functions.find_nearest(x_axis,right_shift)
    y_max = max(y_axis[leftwindow:rightwindow])
    z = y_axis.tolist().index(y_max)
    return y_max,z
def Brackett_Ratios(plateid,mjd,fiber):
    
    import csv
    import os
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import numpy as np

    
    #Header------------------------------------------------------------------------------------------------------------------------------------------|

    filename = File_Path(plateid,mjd,fiber)
    main_header = fits.open(filename)
    
    loc_id = main_header[0].header['LOCID']
    twomass_id = main_header[0].header['OBJID']

    #Reading in the Visits file----------------------------------------------------------------------------------------------------------------------|
    

    fileDir = os.path.dirname(os.path.realpath('__file__'))
    filename = os.path.join(fileDir, '../Data/Average Visits.csv')
    visits = os.path.abspath(os.path.realpath(filename))
    
    br_num = np.asarray([11,12,13,14,15,16,17,18,19,20])
    br_value=np.zeros(10)

    with open(visits) as csvfile:

        reader = csv.DictReader(csvfile,delimiter='\t')

        for row in reader:
            if int(row['Location ID'])==loc_id and row['2Mass ID']==twomass_id:
                for i in range(10):
                    num = 11 + i
                    br = 'Br' + str(num) + ' Avg EqW'
                    br_value[i] = float(row[br])

    #Plotting---------------------------------------------------------------------------------------------------------------------------------------|    
    
    fig,ax=plt.subplots(figsize=(16,8))
    ax.tick_params(axis='both', labelsize=20)
    
    plt.plot(br_num,br_value/br_value[0])
    plt.ylabel('Br n>11 / Br 11',fontsize=24)
    plt.xlabel('n',fontsize=24)
    plt.show()
def apVisit_Updated_Catalog(infile,rowstart):

    import functions
    import pandas as pd
    from astropy.io import fits as fits
    import numpy as np
    import itertools

    x = pd.read_csv(infile,delimiter = '\t')
    problems = []
    cols = ['Location ID','2Mass ID', 'Plate ID','MJD','Fiber','confidence1','confidence2','Mean','arearatios','SNR','Br 11 EqW']#,'Br 12 EqW','Br 13 EqW','Br 14 EqW','Br 15 EqW','Br 16 EqW','Br 17 EqW',
            #'Br 18 EqW','Br 19 EqW','Br 20 EqW']
    
    df = pd.DataFrame(columns = cols)
    g = 0
    #rowstart = 10000
    rowend =  rowstart + 79468
    for index,row in itertools.islice(x.iterrows(),rowstart,rowend):
        try:
            g+=1
            loc = row['Location ID']
            twomass = row['2Mass ID']
            print(str(g))
            plateid = row['Plate']
            mjd = row['MJD']
            if len(str(row['Fiber'])) == 3:
                fiber = str(row['Fiber'])
            elif len(str(row['Fiber'])) == 2:
                fiber = '0' + str(row['Fiber']) 
            else:
                fiber = '00' + str(row['Fiber'])
                    
            fitsfile = functions.File_Path(plateid,mjd,fiber)

            #this passes a filepath that astropy can open with fits, circumventing apogee entirely...
            
            openfile = fits.open(fitsfile)
            snr = openfile[0].header['SNR']
            vbc = openfile[0].header['BC']
            c = 299792.458
            lamshift = 1 + (vbc/c)
            
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
            
            openfile.close()
            
            newflux = functions.OH_Skylines_Remover(wave,flux)
            
            #now we run equiv width calc
            lines = [11,12,13,14,15,16,17,18,19,20]
            EqW = []
            equiv_check = 1000
            
            #just for line 11
            equiv_width,fcontinuum,shift,rest_wavelength = functions.Br_EqW(wave,newflux,11,vbc)
            equiv_check = equiv_width
            #for all the lines
            '''
            for i in range(10):
                
                equiv_width,fcontinuum,shift,rest_wavelength = functions.Br_EqW(wave,newflux,lines[i],vbc)
                EqW.append(equiv_width)
                
                if i == 0:
                    equiv_check = equiv_width
            '''
                
            wave = np.asarray(wave) * lamshift #check which direction shift is going
            
            if equiv_check > 0:
                
                confidence1,confidence2,Mean,arearatios = functions.Confidence_Level(wave,flux,snr)
                #confidence3 = confidence1/confidence2

                data = [int(loc),twomass,int(plateid),int(mjd),fiber,confidence1,confidence2,Mean,arearatios,snr,equiv_width]#EqW[0],EqW[1],EqW[2],EqW[3],EqW[4],EqW[5],EqW[6],EqW[7],EqW[8],EqW[9]]

                df.loc[len(df)+1] = data

                df1 = pd.DataFrame(wave,columns=['Wavelength'])
                df1['Flux'] = newflux
                filename = str(plateid) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
                df1.to_csv('/Users/ballanr/Desktop/File Outputs/Wave and Flux/'+filename,index=False)
                df1 = df1.iloc[0:0]

        except KeyError:
            print('Row '+str(g)+' has no BC value...')
            problems.append((loc,twomass,'KeyError'))
            openfile.close()
        
        except FileNotFoundError:
            print('Row '+str(g)+' doesn\'t exist...')
            problems.append((loc,twomass,'FileNotFound'))

    df.to_csv('/Users/ballanr/Desktop/File Outputs/'+str(rowstart)+'-'+str(rowend)+' End.csv',index=False)
    df = df.iloc[0:0]

    probs = pd.DataFrame(problems,columns = ['Location ID','2Mass ID','Problem Type'])
    probs.to_csv('/Users/ballanr/Desktop/File Outputs/'+str(rowstart)+'-'+str(rowend)+' Problems.csv',index=False)
def Br_EqW(wave,spec,line,vbc):
    
    import functions
    import numpy as np
    
    #wave = np.asarray(wave)
    observed_wavelength,shift,rest_wavelength = functions.Barycentric_Correction(line,vbc)
    rest = rest_wavelength*(10**10)
    centerline = functions.find_nearest(wave,rest)
    
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
    
    return equivs,Fluxcontinuum,shift,rest_wavelength,centerline
def Confidence_Level(wave,flux,restwave):

    import numpy as np
    import functions

    center = functions.find_nearest(wave,restwave)
    L1 = center - 301 # ~ 56 Angstroms
    L2 = center - 150 # ~ 35 Angstroms
    R1 = center + 150
    R2 = center + 301

    #method 1 
    leftwindow = np.mean(flux[L1:L2])
    rightwindow = np.mean(flux[R1:R2])
    cmean = (0.5)*(leftwindow + rightwindow)
    linemean = np.mean(flux[L2:R1])
    lefterr = np.std(flux[L1:L2])
    righterr = np.std(flux[R1:R2])
    cerr = (0.5)*(lefterr + righterr)
    confidence1 = (linemean - cmean)/cerr
    
          
    return confidence1
def skylines_cleaner(wave,flux):
    
    import functions
    import numpy as np

    windows = []
    wave = np.asarray(wave)

    removal_beginning = [15185.0,15240.0,15286.0,15330.0,15393.6,15429.9,15460.8,15473.3,15499.6,15508.0,15516.0,
                        15538.5,15545.0,15568.9,15595.9,15652.9,15701.3,15758.9,15861.5,15867.1,15874.3,15888.8,
                        15971.1,16028.8,16078.2,16126.4,16193.1,16233.4,16269.5,16278.8,16301.1,16314.8,16339.1,
                        16349.3,16356.3,16358.6,16387.1,16413.5,16475.1,16476.9,16478.3,16501.4,16527.4,16530.4,
                        16542.7,16552.9,16558.7,16585.4,16609.4,16611.6,16654.1,16658.2,16666.4,16688.2,16691.2,
                        16701.7,16707.7,16711.8,16717.7,16723.1,16725.1,16730.9,16753.4,16755.9,16761.5,16764.8,
                        16839.5,16852.8,16877.2,16883.0,16888.8,16890.3,16899.6,16902.4,16907.1,16908.6,16909.9,16913.7]

    removal_ending = [15189.0,15242.6,15289.0,15334.7,15396.6,15433.9,15463.8,15475.3,15502.3,15511.0,
                    15518.8,15543.4,15547.6,15571.4,15599.3,15658.9,15704.3,15761.8,15863.6,15871.0,
                    15877.0,15892.7,15974.4,16033.1,16081.6,16131.0,16196.5,16237.9,16271.4,16280.7,
                    16304.7,16318.7,16343.0,16353.7,16357.2,16361.7,16390.3,16416.2,16476.3,16477.9,
                    16479.7,16503.7,16528.9,16531.1,16543.5,16555.1,16560.2,16587.5,16610.7,16612.8,
                    16655.0,16661.0,16667.3,16690.4,16693.9,16703.8,16710.5,16712.5,16719.7,16724.5,
                    16725.9,16734.9,16755.2,16756.8,16762.5,16765.9,16841.6,16853.3,16877.8,16884.5,
                    16889.5,16891.2,16900.2,16905.5,16907.7,16909.3,16910.6,16914.3]

    for i in range(len(removal_beginning)):

        lwindow = functions.find_nearest(wave,removal_beginning[i])
        rwindow = functions.find_nearest(wave,removal_ending[i])

        windows.append((lwindow,rwindow))

    for i in range(len(removal_beginning)):

        lwindowelement = int(windows[i][0])
        rwindowelement = int(windows[i][1])
        leftwindow = wave[windows[i][0]]
        rightwindow = wave[windows[i][1]]
        leftflux = flux[windows[i][0]]
        rightflux = flux[windows[i][1]]

        slope = (rightflux - leftflux) / (rightwindow - leftwindow)

        for k in range(rwindowelement - lwindowelement):

            fluxvalue = slope*(wave[lwindowelement+k] - leftwindow) + leftflux
            flux[lwindowelement+k] = fluxvalue

    return flux
def DR13_Brackett_Catalog():

    import functions
    import pandas as pd
    from astropy.io import fits as fits
    import numpy as np
    import itertools
    from PyAstronomy.pyasl import helcorr as helcorr

    dr13table = pd.read_csv('/Users/ballanr/Desktop/File Outputs/DR13/Visits.csv',delimiter='\t')
    problems = []
    cols = ['Location ID','2Mass ID', 'Plate ID','MJD','Fiber','S/R','Model Density','Model Temp','Overall Confidence',
            'Br11 Error','Br12 Error','Br13 Error','Br14 Error','Br15 Error','Br16 Error','Br17 Error','Br18 Error',
            'Br19 Error','Br20 Error','Br11 EqW','Br12 EqW','Br13 EqW','Br14 EqW','Br15 EqW','Br16 EqW','Br17 EqW',
            'Br18 EqW','Br19 EqW','Br20 EqW']

    rowstart = 500000
    #rowend = rowstart + 100000

    df = pd.DataFrame(columns = cols)
    g = 0

    for index,row in itertools.islice(dr13table.iterrows(),rowstart,None):
        try:
            g+=1
            loc = row['Location ID']
            twomass = row['2Mass ID']
            print(str(g))
            plateid = row['Plate']
            mjd = row['MJD']
            
            if len(str(row['Fiber'])) == 3:
                fiber = str(row['Fiber'])
            elif len(str(row['Fiber'])) == 2:
                fiber = '0' + str(row['Fiber']) 
            else:
                fiber = '00' + str(row['Fiber'])

            fitsfile = '/Volumes/CoveyData/APOGEE_Spectra/python_DR13/dr13/apogee/spectro/redux/r6/apo25m/'+str(plateid)+'/'+str(mjd)+'/apVisit-r6-'+str(plateid)+'-'+str(mjd)+'-'+str(fiber)+'.fits'

            #this passes a filepath that astropy can open with fits, circumventing apogee entirely...
            
            openfile = fits.open(fitsfile)

            header = openfile[0].header

            snr = header['SNR'] 
            ra = header['RA']
            dec = header['DEC'] 
            jd = header['JD-MID']
            
            height = 2788
            longg = -105.4913
            lat = 36.4649

            vbc,hjd = helcorr(longg,lat,height,ra,dec,jd)

            c = 299792.458
            lamshift = 1 + (vbc/c)
            
            fspec = openfile[1]
            ferr = openfile[2]
            fwave = openfile[4]
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
            
            openfile.close()
            
            newflux = functions.skylines_cleaner(wave,flux)
            
            #now we run equiv width calc
            lines = [11,12,13,14,15,16,17,18,19,20]
            EqW = []
            restwaves = []
            econfidences = []
            confidences = []
            equivs_error = []
            equiv_check = 1000
            
            wave = np.asarray(wave) * lamshift #check which direction shift is going

            for i in range(10):
                
                equiv_width,fcontinuum,shift,rest_wavelength,centers = functions.Br_EqW(wave,newflux,lines[i],vbc)
                rest = rest_wavelength*(10**10)
                EqW.append(equiv_width)
                restwaves.append(rest)
                
                dEqW = functions.Br_Error(wave,newflux,error,lines[i],rest)
                equivs_error.append(dEqW)
                
                if i == 0:
                    equiv_check = equiv_width
            
            if equiv_check > 0:
                
                modeldens,modeltemp = functions.Model_Fitter(EqW)

                Confidence = functions.Confidence_Level(wave,flux,restwaves[0])

                data = [int(loc),twomass,int(plateid),int(mjd),fiber,snr,modeldens,modeltemp,Confidence,equivs_error[0],equivs_error[1],
                        equivs_error[2],equivs_error[3],equivs_error[4],equivs_error[5],equivs_error[6],equivs_error[7],
                        equivs_error[8],equivs_error[9],EqW[0],EqW[1],EqW[2],EqW[3],EqW[4],EqW[5],EqW[6],EqW[7],EqW[8],EqW[9]]

                df.loc[len(df)+1] = data

                df1 = pd.DataFrame(wave,columns=['Wavelength'])
                df1['Flux'] = newflux
                df1['Error'] = error
                filename = str(plateid) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
                df1.to_csv('/Users/ballanr/Desktop/File Outputs/DR13/Wave and Flux/'+filename,index=False)
                df1 = df1.iloc[0:0]                        

        except KeyError:
            print('Row '+str(g)+' has no BC value...')
            problems.append((loc,twomass,plateid,mjd,fiber,'KeyError'))
            openfile.close()
        
        except FileNotFoundError:
            print('Row '+str(g)+' doesn\'t exist...')
            problems.append((loc,twomass,plateid,mjd,fiber,'FileNotFound'))


    df.to_csv('/Users/ballanr/Desktop/File Outputs/DR13/Section 6.csv',index=False)
    df = df.iloc[0:0]

    probs = pd.DataFrame(problems,columns = ['Location ID','2Mass ID','Plate ID','MJD','Fiber','Problem Type'])
    probs.to_csv('/Users/ballanr/Desktop/File Outputs/DR13/Section 6 Problems.csv',index=False)
def DR14_Brackett_Catalog():

    #Some of the .fits files in the DR14 folders have NaN values in them and break the loop. Need to
    #find them and fix them before a catalog can be run

    import functions
    import pandas as pd
    from astropy.io import fits as fits
    import numpy as np
    import itertools
    from PyAstronomy.pyasl import helcorr as helcorr

    dr14table = pd.read_csv('/Users/ballanr/Desktop/File Outputs/DR14/DR14 List.csv')
    problems = []
    cols = ['Location ID','2Mass ID', 'Plate ID','MJD','Fiber','S/R','Model Density','Model Temp','Overall Confidence',
            'Br11 Error','Br12 Error','Br13 Error','Br14 Error','Br15 Error','Br16 Error','Br17 Error','Br18 Error',
            'Br19 Error','Br20 Error','Br11 EqW','Br12 EqW','Br13 EqW','Br14 EqW','Br15 EqW','Br16 EqW','Br17 EqW',
            'Br18 EqW','Br19 EqW','Br20 EqW']

    rowstart = 0

    df = pd.DataFrame(columns = cols)
    g = 0

    for index,row in itertools.islice(dr14table.iterrows(),rowstart,rowstart+50000):
        try:
            g+=1
            print(str(g))

            plateid = row['Plate']
            mjd = row['MJD']
            
            if len(str(row['Fiber'])) == 3:
                fiber = str(row['Fiber'])
            elif len(str(row['Fiber'])) == 2:
                fiber = '0' + str(row['Fiber']) 
            else:
                fiber = '00' + str(row['Fiber'])

            endpath = row['Path']
            
            fitsfile = '/Volumes/CoveyData/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/apo25m/'+str(plateid)+'/'+str(mjd)+'/'+endpath

            #this passes a filepath that astropy can open with fits, circumventing apogee entirely...
            
            openfile = fits.open(fitsfile)

            header = openfile[0].header

            ra = header['RA']
            dec = header['DEC'] 
            jd = header['JD-MID']
            loc = header['LOCID']
            twomass = header['OBJID']
            snr = header['SNR']
            
            height = 2788
            longg = -105.4913
            lat = 36.4649
            

            vbc,hjd = helcorr(longg,lat,height,ra,dec,jd)

            c = 299792.458
            lamshift = 1 + (vbc/c)
            
            fspec = openfile[1]
            ferr = openfile[2]
            fwave = openfile[4]
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
            
            openfile.close()
            
            newflux = functions.skylines_cleaner(wave,flux)
            
            #now we run equiv width calc
            lines = [11,12,13,14,15,16,17,18,19,20]
            EqW = []
            restwaves = []
            econfidences = []
            confidences = []
            equivs_error = []
            equiv_check = 1000
            
            wave = np.asarray(wave) * lamshift #check which direction shift is going

            for i in range(10):
                
                equiv_width,fcontinuum,shift,rest_wavelength,centers = functions.Br_EqW(wave,newflux,lines[i],vbc)
                rest = rest_wavelength*(10**10)
                EqW.append(equiv_width)
                restwaves.append(rest)
                
                dEqW = functions.Br_Error(wave,newflux,error,lines[i],rest)
                equivs_error.append(dEqW)
                
                if i == 0:
                    equiv_check = equiv_width
            
            if equiv_check > 0:
                
                modeldens,modeltemp = functions.Model_Fitter(EqW,equivs_error)

                Confidence = functions.Confidence_Level(wave,flux,restwaves[0])
    
                data = [int(loc),twomass,int(plateid),int(mjd),fiber,snr,modeldens,modeltemp,Confidence,equivs_error[0],equivs_error[1],
                        equivs_error[2],equivs_error[3],equivs_error[4],equivs_error[5],equivs_error[6],equivs_error[7],
                        equivs_error[8],equivs_error[9],EqW[0],EqW[1],EqW[2],EqW[3],EqW[4],EqW[5],EqW[6],EqW[7],EqW[8],EqW[9]]

                df.loc[len(df)+1] = data

                df1 = pd.DataFrame(wave,columns=['Wavelength'])
                df1['Flux'] = newflux
                df1['Error'] = error
                filename = str(plateid) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
                df1.to_csv('/Users/ballanr/Desktop/File Outputs/DR14/Wave and Flux/'+filename,index=False)
                df1 = df1.iloc[0:0]                       

        except KeyError:
            print('Row '+str(g)+' has no BC value...')
            problems.append((loc,twomass,plateid,mjd,fiber,'KeyError'))
            openfile.close()
        
        except FileNotFoundError:
            print('Row '+str(g)+' doesn\'t exist...')
            problems.append((loc,twomass,plateid,mjd,fiber,'FileNotFound'))
        except IndexError:
            print('Index '+str(g)+' broke...')
            problems.append((loc,twomass,plateid,mjd,fiber,'IndexError')) 

    df.to_csv('/Users/ballanr/Desktop/File Outputs/DR14/DR14 Piece 1.csv',index=False)
    df = df.iloc[0:0]

    probs = pd.DataFrame(problems,columns = ['Location ID','2Mass ID','Plate ID','MJD','Fiber','Problem Type'])
    probs.to_csv('/Users/ballanr/Desktop/File Outputs/DR14/DR14 Problems 1.csv',index=False)
def DR15_Brackett_Catalog():

    import functions
    import pandas as pd
    from astropy.io import fits as fits
    import numpy as np
    import itertools
    from PyAstronomy.pyasl import helcorr as helcorr

    dr15table = pd.read_csv('/Users/ballanr/Desktop/File Outputs/DR15/appended list.csv')
    problems = []
    cols = ['Location ID','2Mass ID', 'Plate ID','MJD','Fiber','S/R','Model Density','Model Temp','Overall Confidence',
            'Br11 Error','Br12 Error','Br13 Error','Br14 Error','Br15 Error','Br16 Error','Br17 Error','Br18 Error',
            'Br19 Error','Br20 Error','Br11 EqW','Br12 EqW','Br13 EqW','Br14 EqW','Br15 EqW','Br16 EqW','Br17 EqW',
            'Br18 EqW','Br19 EqW','Br20 EqW']

    rowstart = 0

    df = pd.DataFrame(columns = cols)
    g = 0

    for index,row in itertools.islice(dr15table.iterrows(),rowstart,None):
        try:
            g+=1
            #loc = row['LOCID']
            #twomass = row['2MASS_ID']
            print(str(g))
            plateid = row['Plate']
            mjd = row['MJD']
            #snr = row['S/N']
            if len(str(row['Fiber'])) == 3:
                fiber = str(row['Fiber'])
            elif len(str(row['Fiber'])) == 2:
                fiber = '0' + str(row['Fiber']) 
            else:
                fiber = '00' + str(row['Fiber'])

            server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
            filepath = row['filepath']
            fitsfile = server + str(plateid) + '/' + str(mjd) + '/' + str(filepath)

            #this passes a filepath that astropy can open with fits, circumventing apogee entirely...
            
            openfile = fits.open(fitsfile)

            header = openfile[0].header

            loc = header['LOCID']
            twomass = header['OBJID']
            try:
                snr = header['SNR']
            except:
                snr = 0
                print('Row '+str(g)+' has no SNR value...')
                problems.append((loc,twomass,plateid,mjd,fiber,'NaN'))
            
            ra = header['RA']
            dec = header['DEC'] 
            jd = header['JD-MID']
            #telescope = header['TELESCOP']
            
            height = 2788
            longg = -105.4913
            lat = 36.4649

            '''if telescope == 'apo25m':
                height = 2788
                longg = -105.4913
                lat = 36.4649
            else:
                height = 2380
                longg = -70.413336
                lat = -29.05256'''

            vbc,hjd = helcorr(longg,lat,height,ra,dec,jd)

            c = 299792.458
            lamshift = 1 + (vbc/c)
            
            fspec = openfile[1]
            ferr = openfile[2]
            fwave = openfile[4]
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
            
            openfile.close()
            
            newflux = functions.skylines_cleaner(wave,flux)
            
            #now we run equiv width calc
            lines = [11,12,13,14,15,16,17,18,19,20]
            EqW = []
            restwaves = []
            econfidences = []
            confidences = []
            equivs_error = []
            equiv_check = 1000
            
            wave = np.asarray(wave) * lamshift #check which direction shift is going

            for i in range(10):
                
                equiv_width,fcontinuum,shift,rest_wavelength,centers = functions.Br_EqW(wave,newflux,lines[i],vbc)
                rest = rest_wavelength*(10**10)
                EqW.append(equiv_width)
                restwaves.append(rest)
                
                dEqW = functions.Br_Error(wave,newflux,error,lines[i],rest)
                equivs_error.append(dEqW)
                
                if i == 0:
                    equiv_check = equiv_width
            
            if equiv_check > 0:
                
                Confidence = functions.Confidence_Level(wave,newflux,restwaves[0])
            
                if Confidence > 1.25:

                    modeldens,modeltemp = functions.Model_Fitter(EqW,equivs_error)

                    data = [int(loc),twomass,int(plateid),int(mjd),fiber,snr,modeldens,modeltemp,Confidence,equivs_error[0],equivs_error[1],
                            equivs_error[2],equivs_error[3],equivs_error[4],equivs_error[5],equivs_error[6],equivs_error[7],
                            equivs_error[8],equivs_error[9],EqW[0],EqW[1],EqW[2],EqW[3],EqW[4],EqW[5],EqW[6],EqW[7],EqW[8],EqW[9]]

                    df.loc[len(df)+1] = data

                    df1 = pd.DataFrame(wave,columns=['Wavelength'])
                    df1['Flux'] = newflux
                    df1['Error'] = error
                    filename = str(plateid) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
                    df1.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Wave and Flux/'+filename,index=False)
                    df1 = df1.iloc[0:0]                        

        except KeyError:
            print('Row '+str(g)+' has no BC value...')
            problems.append((loc,twomass,plateid,mjd,fiber,'KeyError'))
            openfile.close()
        
        except FileNotFoundError:
            print('Row '+str(g)+' doesn\'t exist...')
            problems.append((loc,twomass,plateid,mjd,fiber,'FileNotFound'))
        except IndexError:
            print('Index '+str(g)+' broke...')
            problems.append((loc,twomass,plateid,mjd,fiber,'IndexError')) 

    df.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/blah.csv',index=False)
    df = df.iloc[0:0]

    probs = pd.DataFrame(problems,columns = ['Location ID','2Mass ID','Plate ID','MJD','Fiber','Problem Type'])
    probs.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/blah Problems.csv',index=False)
def Br_Error(wave,flux,err,line,rest):
    
    import functions
    import numpy as np

    try:
        if len(err) > 0:
            center = functions.find_nearest(wave,rest)
            
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
def Brackett_Ratios_Updated(plate,mjd,fiber):

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import itertools
    
    #finding the equivalent widths and errors for a given visit
    filefile = pd.read_csv('/Users/ballanr/Desktop/File Outputs/DR15/Catalog-test.csv')

    if len(str(fiber)) == 3:
        Fiber = str(fiber)
    elif len(str(fiber)) == 2:
        Fiber = '0' + str(fiber) 
    else:
        Fiber = '00' + str(fiber)
    equivs = []
    errors = []
    for index,row in itertools.islice(filefile.iterrows(),0,None):
        plateid = int(row['Plate ID'])
        MJD = int(row['MJD'])
        ffiber = str(row['Fiber'])
        
        if len(str(ffiber)) == 3:
            ffiber = str(ffiber)
        elif len(str(ffiber)) == 2:
            ffiber = '0' + str(ffiber) 
        else:
            ffiber = '00' + str(ffiber)

        if plateid == plate and MJD == mjd and ffiber == Fiber:

            for j in range(10):
                number = 11+j
                errorstring = 'Br' + str(number) + ' Error'
                equivstring = 'Br' + str(number) + ' EqW'
                hh = row[errorstring]
                gg = row[equivstring]
                equivs.append(gg)
                errors.append(hh)


    equivs = np.asarray(equivs)
    errors = np.asarray(errors)
    equivs1 = equivs/equivs[0]

    summederrors = []

    for i in range(10):

        dq = equivs1[i] * np.sqrt( (errors[i]/equivs[i])**2 + (errors[0]/equivs[0])**2 )
        summederrors.append(dq)
        #print(dq)

    temps = ['3750 K','5000 K','7500 K','8750 K','10000 K','12500 K','15000 K']

    #loop through densities
    densityrange = np.arange(8,12.6,0.2)
    for i in range(len(densityrange)):
        
        density = densityrange[i]
        
        filename = '/Users/ballanr/Desktop/File Outputs/Brackett Decrements/Density Files/Density ' + str(8+(i*0.2)) + ' Ratios.csv'
        openfile = pd.read_csv(filename)

        y = []

        for j in range(len(temps)):

            y = openfile[temps[j]]
                    
            y = np.asarray(y)
            y = y/ y[0]
            plt.figure(figsize=(20,10))
            plt.plot(np.arange(11,21,1),y,color='purple',label = temps[j])
            plt.scatter(np.arange(11,21,1),equivs1,color='purple',label='_nolabel_')
            plt.errorbar(np.arange(11,21,1),equivs1,summederrors,ecolor='red',fmt='--o',elinewidth = 1,capsize = 5,label='Br Emission Lines')
            plt.xticks(np.arange(11,21,1))
            plt.legend(bbox_to_anchor=(1,1),fontsize=18)
            plt.grid(True,color='grey',ls='dashed',alpha=0.7)
            plt.title('Density ' + str(density))
            plt.xlabel('Brackett Lines',fontsize = 18)
            plt.ylabel('$Br_{N\geq 11}/Br_{11}$',fontsize = 18)
            #plt.show()
            savestring = '/Users/ballanr/Desktop/File Outputs/Brackett Decrements/test/D' + str(density) + 'T' + str(temps[j]) + '.pdf'
            plt.savefig(savestring,bbox_inches='tight',dpi=300)
            plt.close()
def Model_Fitter(eqwarray,error):
    import pandas as pd
    import numpy as np

    allfile = '/Users/ballanr/Desktop/File Outputs/Brackett Decrements/Profile Test.csv'
    openallfile = pd.read_csv(allfile)
    cols = openallfile.columns
    eqwarray = np.asarray(eqwarray)
    eqwarray= eqwarray/eqwarray[0]
    errors = []
    test = [(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0)]

    for i in range(10):

        dq = eqwarray[i] * np.sqrt( (error[i]/eqwarray[i])**2 + (error[0]/eqwarray[0])**2 )
        errors.append(dq)

    resid = []
    for i in range(len(cols)):
        
        A = openallfile[cols[i]]
        
        b = eqwarray
        
        A = A / A[0]
        
        r = b - A
        
        r2 = 0
        
        for k in range(len(r)):
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
def Brackett_Ratios_Updated_Grid():

    import matplotlib as mpl
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import itertools
    import seaborn as sb
    
    filepath = '/Users/ballanr/Desktop/File Outputs/DR15/pre-DR15 Cutoff.csv'
    filefile = pd.read_csv(filepath)

    
    mpl.rcParams.update(mpl.rcParamsDefault)

    for index,row in itertools.islice(filefile.iterrows(),0,None):
        
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
        plt.plot(np.arange(11,21,1),y,color=sb.xkcd_rgb['black'],label='Best Fit, Density: '+str(modeldens)+' Temp: '+str(modelt))
        plt.scatter(np.arange(11,21,1),y,color=sb.xkcd_rgb['black'])
        plt.errorbar(np.arange(11,21,1),equivs1,yerr=summederrors,ecolor='red',barsabove=True,
                    capsize=3,fmt='-o',color=sb.xkcd_rgb['blue'],label='Br Emission')

        plt.plot(np.arange(11,21,1),d8T3750/d8T3750[0],color=sb.xkcd_rgb['red'],ls='dashed',label='Density: 8 Temp: 3750')
        plt.plot(np.arange(11,21,1),d8T8750/d8T8750[0],color=sb.xkcd_rgb['red'],ls='dashdot',label='Density: 8 Temp: 8750')
        plt.plot(np.arange(11,21,1),d8T15000/d8T15000[0],color=sb.xkcd_rgb['red'],ls='solid',label='Density: 8 Temp: 15000')
        plt.plot(np.arange(11,21,1),d10T3750/d10T3750[0],color=sb.xkcd_rgb['forest green'],ls='dashed',label='Density: 10 Temp: 3750')
        plt.plot(np.arange(11,21,1),d10T8750/d10T8750[0],color=sb.xkcd_rgb['forest green'],ls='dashdot',label='Density: 10 Temp: 8750')
        plt.plot(np.arange(11,21,1),d10T15000/d10T15000[0],color=sb.xkcd_rgb['forest green'],ls='solid',label='Density: 10 Temp: 15000')
        plt.plot(np.arange(11,21,1),d12T3750/d12T3750[0],color=sb.xkcd_rgb['purple'],ls='dashed',label='Density: 12 Temp: 3750')
        plt.plot(np.arange(11,21,1),d12T8750/d12T8750[0],color=sb.xkcd_rgb['purple'],ls='dashdot',label='Density: 12 Temp: 8750')
        plt.plot(np.arange(11,21,1),d12T15000/d12T15000[0],color=sb.xkcd_rgb['purple'],ls='solid',label='Density: 12 Temp: 15000')
        plt.grid(True)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=4,frameon=True,facecolor=None,framealpha=1,fontsize=16)
        plt.savefig('/Users/ballanr/Desktop/File Outputs/DR15/DR15 Brackett Decrement Plots/' + savestring,bbox_inches='tight',dpi=300)
        plt.clf()
        plt.close()
def Aitoff(dr):
    from astropy.coordinates import SkyCoord  # High-level coordinates
    from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
    from astropy.coordinates import Angle, Latitude, Longitude  # Angles
    import astropy.units as u
    import pandas as pd
    import numpy as np
    import itertools
    from astropy.io import fits
    import matplotlib.pyplot as plt


    if dr == 13:
        glong = []
        glat = []
        error = []

        file1 = '/Users/ballanr/Desktop/File Outputs/DR13/DR13 Protostars.csv'
        openfile = pd.read_csv(file1)

        for index,row in itertools.islice(openfile.iterrows(),0,None):
            try:
                plate = row['Plate ID']
                mjd = row['MJD']
                if len(str(row['Fiber'])) == 3:
                    fiber = str(row['Fiber'])
                elif len(str(row['Fiber'])) == 2:
                    fiber = '0' + str(row['Fiber']) 
                else:
                    fiber = '00' + str(row['Fiber'])

                server = '/Volumes/CoveyData/APOGEE_Spectra/python_DR13/dr13/apogee/spectro/redux/r6/apo25m/'
                filepath = str(plate) + '/' + str(mjd) + '/' + 'apVisit-r6-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'
               
                fitsfile = fits.open(server+filepath)

                Glon = fitsfile[0].header['GLON']
                Glat = fitsfile[0].header['GLAT']

                glong.append(Glon)
                glat.append(Glat)
            
            except KeyError:
                #print('This has no coords...')
                error.append((plate,mjd,fiber))

        glong = np.asarray(glong)
        glat = np.asarray(glat)

        glong1 = []
        glat1 = []
        error1 = []

        file1 = '/Users/ballanr/Desktop/File Outputs/DR13/DR13.csv'
        openfile = pd.read_csv(file1)

        for index,row in itertools.islice(openfile.iterrows(),0,None):
            try:
                plate = row['Plate ID']
                mjd = row['MJD']
                if len(str(row['Fiber'])) == 3:
                    fiber = str(row['Fiber'])
                elif len(str(row['Fiber'])) == 2:
                    fiber = '0' + str(row['Fiber']) 
                else:
                    fiber = '00' + str(row['Fiber'])

                server = '/Volumes/CoveyData/APOGEE_Spectra/python_DR13/dr13/apogee/spectro/redux/r6/apo25m/'
                filepath = str(plate) + '/' + str(mjd) + '/' + 'apVisit-r6-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'
                
                fitsfile = fits.open(server+filepath)

                Glon = fitsfile[0].header['GLON']
                Glat = fitsfile[0].header['GLAT']

                glong1.append(Glon)
                glat1.append(Glat)

            except KeyError:
                #print('This has no coords...')
                error1.append((plate,mjd,fiber))

        glong1 = np.asarray(glong1)
        glat1 = np.asarray(glat1)

    elif dr == 14:
        
        glong = []
        glat = []
        error = []

        file1 = '/Users/ballanr/Desktop/File Outputs/DR14/DR14 Protostars.csv'
        openfile = pd.read_csv(file1)

        for index,row in itertools.islice(openfile.iterrows(),0,None):
            try:
                plate = row['Plate ID']
                mjd = row['MJD']
                if len(str(row['Fiber'])) == 3:
                    fiber = str(row['Fiber'])
                elif len(str(row['Fiber'])) == 2:
                    fiber = '0' + str(row['Fiber']) 
                else:
                    fiber = '00' + str(row['Fiber'])


                if plate > 9700:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/lco25m/'
                    filepath = str(plate) + '/' + str(mjd) + '/' + 'asVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'
                else:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                    filepath = str(plate) + '/' + str(mjd) + '/' + 'apVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'

                fitsfile = fits.open(server+filepath)

                Glon = fitsfile[0].header['GLON']
                Glat = fitsfile[0].header['GLAT']

                glong.append(Glon)
                glat.append(Glat)
            except KeyError:
                #print('This has no coords...')
                error.append((plate,mjd,fiber))

        glong = np.asarray(glong)
        glat = np.asarray(glat)

        glong1 = []
        glat1 = []
        error1 = []

        file1 = '/Users/ballanr/Desktop/File Outputs/DR15/forMdot_analysis.csv'
        openfile = pd.read_csv(file1)

        for index,row in itertools.islice(openfile.iterrows(),0,None):
            try:
                plate = row['PLATE']
                mjd = row['MJD']
                if len(str(row['FIB'])) == 3:
                    fiber = str(row['FIB'])
                elif len(str(row['FIB'])) == 2:
                    fiber = '0' + str(row['FIB']) 
                else:
                    fiber = '00' + str(row['FIB'])

                if plate > 9700:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/lco25m/'
                    filepath = str(plate) + '/' + str(mjd) + '/' + 'asVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'
                else:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                    filepath = str(plate) + '/' + str(mjd) + '/' + 'apVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'

                fitsfile = fits.open(server+filepath)

                Glon = fitsfile[0].header['GLON']
                Glat = fitsfile[0].header['GLAT']

                glong1.append(Glon)
                glat1.append(Glat)

            except KeyError:
                #print('This has no coords...')
                error1.append((plate,mjd,fiber))

        glong1 = np.asarray(glong1)
        glat1 = np.asarray(glat1)
            
        
    elif dr == 15:
        
        RA = []
        DEC = []
        error = []

        file1 = '/Users/ballanr/Desktop/File Outputs/DR15/pre-DR15 Cutoff.csv'
        openfile = pd.read_csv(file1)

        for index,row in itertools.islice(openfile.iterrows(),0,None):
            try:
                plate = row['Plate ID']
                mjd = row['MJD']
                if len(str(row['Fiber'])) == 3:
                    fiber = str(row['Fiber'])
                elif len(str(row['Fiber'])) == 2:
                    fiber = '0' + str(row['Fiber']) 
                else:
                    fiber = '00' + str(row['Fiber'])


                if plate > 9700:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/lco25m/'
                    filepath = str(plate) + '/' + str(mjd) + '/' + 'asVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'
                else:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                    filepath = str(plate) + '/' + str(mjd) + '/' + 'apVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'

                fitsfile = fits.open(server+filepath)

                rra = fitsfile[0].header['RA']
                ddec = fitsfile[0].header['DEC']

                RA.append(rra)
                DEC.append(ddec)

            except KeyError:
                #print('This has no coords...')
                error.append((plate,mjd,fiber))

        RA = np.asarray(RA)
        DEC = np.asarray(DEC)

        RA1 = []
        DEC1 = []
        error1 = []

        file1 = '/Users/ballanr/Desktop/File Outputs/DR15/forMdot_analysis.csv'
        openfile = pd.read_csv(file1)

        for index,row in itertools.islice(openfile.iterrows(),0,None):
            try:
                plate = row['PLATE']
                mjd = row['MJD']
                if len(str(row['FIB'])) == 3:
                    fiber = str(row['FIB'])
                elif len(str(row['FIB'])) == 2:
                    fiber = '0' + str(row['FIB']) 
                else:
                    fiber = '00' + str(row['FIB'])

                if plate > 9700:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/lco25m/'
                    filepath = str(plate) + '/' + str(mjd) + '/' + 'asVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'
                else:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                    filepath = str(plate) + '/' + str(mjd) + '/' + 'apVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + fiber + '.fits'

                fitsfile = fits.open(server+filepath)

                rra1 = fitsfile[0].header['RA']
                ddec1 = fitsfile[0].header['DEC']

                RA1.append(rra1)
                DEC1.append(ddec1)

            except KeyError:
                #print('This has no coords...')
                error1.append((plate,mjd,fiber))

        RA1 = np.asarray(RA1)
        DEC1 = np.asarray(DEC1)
    
    fig=plt.figure(figsize=(20,10))

    c1 = SkyCoord(ra=RA*u.degree , dec = DEC*u.degree, frame='icrs')
    c1 = c1.galactic
    c2 = SkyCoord(ra=RA1*u.degree , dec = DEC1*u.degree, frame='icrs')
    c2 = c2.galactic

    longg=c1.l.wrap_at(180*u.deg).radian
    latg=c1.b.radian
    longg1=c2.l.wrap_at(180*u.deg).radian
    latg1=c2.b.radian

    ax = plt.subplot(111,projection='aitoff')
    ax.tick_params(axis='y', labelsize=20)
    ax.tick_params(axis='x',labelsize=16)
    ax.grid(True)

    ax.scatter(longg1,latg1,s=2,color='green',alpha=0.1)
    ax.scatter(longg,latg,s=100,color='red',marker='*',edgecolor='black',linewidth=0.5)

    #plt.show()
    plt.savefig('/Users/ballanr/Desktop/Aitoff3.png',bbox_inches='tight',dpi=300)
def Protostar_Spectra_Plotter():

    import pandas as pd
    import matplotlib.pyplot as plt
    import itertools

    filepath = '/Users/ballanr/Desktop/File Outputs/DR15/pre-DR15 Cutoff.csv'
    openfile = pd.read_csv(filepath)

    for index,row in itertools.islice(openfile.iterrows(),0,None): 

        plateid = int(row['Plate ID'])
        MJD = int(row['MJD'])
        fiber = str(row['Fiber'])
        
        if len(str(fiber)) == 3:
            fiber = str(fiber)
        elif len(str(fiber)) == 2:
            fiber = '0' + str(fiber) 
        else:
            fiber = '00' + str(fiber)
        
        folder = '/Users/ballanr/Desktop/File Outputs/DR15/Wave and Flux/'
        wavefile = folder + str(plateid) + '-' + str(MJD) + '-' + str(fiber) + '.csv'
        tempfile = pd.read_csv(wavefile)

        wave = tempfile['Wavelength']
        flux = tempfile['Flux']
        error = tempfile['Error']

        maxflux = max(flux) + 100
        minflux = min(flux) - 100

        plt.figure(figsize=(20,10))
        plt.title(str(plateid)+ '-' + str(MJD) + '-' + str(fiber),fontsize=20)
        plt.xlabel('Wavelength',fontsize=16)
        plt.ylabel('Flux',fontsize=16)
        plt.errorbar(wave,flux,yerr=error,ecolor='red')
        plt.ylim(minflux,maxflux)
        
        savestring = str(plateid)+ '-' + str(MJD) + '-' + str(fiber) + '.png'

        plt.savefig('/Users/ballanr/Desktop/File Outputs/DR15/Plots/Spectra/'+savestring,bbox_inches='tight',dpi=300)

        plt.clf()
        plt.close()
def Single_Spectra_Plotter(plateid,MJD,fiber):

        import pandas as pd
        import matplotlib.pyplot as plt


        folder = '/Users/ballanr/Desktop/File Outputs/DR15/Wave and Flux/'
        wavefile = folder + str(plateid) + '-' + str(MJD) + '-' + str(fiber) + '.csv'
        tempfile = pd.read_csv(wavefile)

        wave = tempfile['Wavelength']
        flux = tempfile['Flux']
        error = tempfile['Error']

        maxflux = max(flux) + 100
        minflux = min(flux) - 100

        plt.figure(figsize=(20,10))
        plt.title(str(plateid)+ '-' + str(MJD) + '-' + str(fiber),fontsize=20)
        plt.xlabel('Wavelength',fontsize=16)
        plt.ylabel('Flux',fontsize=16)
        #plt.errorbar(wave,flux,yerr=error,ecolor='red')
        plt.plot(wave,flux)
        plt.ylim(0,5000)
        plt.xlim(16770,16850)

        #plt.show()
        plt.savefig('/Users/ballanr/Desktop/Category 1 Zoomed.pdf',bbox_inches='tight',dpi=300)
def DR15_Uniques():
    import functions
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import pandas as pd
    import numpy as np
    import itertools
    from PyAstronomy.pyasl import helcorr as helcorr

    cols = ['Location ID','2Mass ID','Plate','MJD','Fiber','RA','DEC','SNR','Density','Temp','Confidence','Br11 EqW','Br11 EqW Error']
    df = pd.DataFrame(columns = cols)

    filepath = '/Users/ballanr/Desktop/File Outputs/DR15/Full Uniques.csv'
    openfile = pd.read_csv(filepath)

    g=0

    for index,row in itertools.islice(openfile.iterrows(),0,None):

            g+=1
            print(g)
            loc = int(row['Location ID'])
            twomass = row['2Mass ID']
            plate = int(row['Plate'])
            mjd = int(row['MJD'])
            fiber = int(row['Fiber'])

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

            snr = header['SNR']
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
            EqW = []
            restwaves = []
            confidences = []
            equivs_error = []
            equiv_check = 1000
                
            wave = np.asarray(wave) * lamshift #check which direction shift is going            

            for k in range(10):
                equiv_width,fcontinuum,shift,rest_wavelength,centers = functions.Br_EqW(wave,newflux,lines[k],vbc)
                rest = rest_wavelength*(10**10)
                EqW.append(equiv_width)
                restwaves.append(rest)
                dEqW = functions.Br_Error(wave,newflux,error,lines[k],rest)
                equivs_error.append(dEqW)

            modeldens,modeltemp = functions.Model_Fitter(EqW,equivs_error)

            Confidence = functions.Confidence_Level(wave,newflux,restwaves[0])

            data = [int(loc),twomass,int(plate),int(mjd),fiber,ra,dec,snr,modeldens,modeltemp,Confidence,EqW[0],equivs_error[0]]

            df.loc[len(df)+1] = data

    df.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/somelist1.csv',index=False)
    df = df.iloc[0:0]
def Brackett_Decrement_Modified(equivs,errors,density,temp,folder,title):

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sb

    modeldens = str(density)
    modelt = str(temp)

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
    yerrors = np.asarray(summederrors)

    plt.figure(figsize=(20,10))
    plt.title(str(title) + ' Decrement',fontsize=20)
    plt.xticks((np.arange(11,21,1)))
    plt.plot(np.arange(11,21,1),y,color=sb.xkcd_rgb['black'],label='Best Fit, Density: '+str(modeldens)+' Temp: '+str(modelt))
    plt.scatter(np.arange(11,21,1),y,color=sb.xkcd_rgb['black'])
    plt.errorbar(np.arange(11,21,1),equivs1,yerr=yerrors,ecolor='red',barsabove=True,capsize=3,fmt='-o',color=sb.xkcd_rgb['blue'],label='Br Emission')

    plt.plot(np.arange(11,21,1),d8T3750/d8T3750[0],color=sb.xkcd_rgb['red'],ls='dashed',label='Density: 8 Temp: 3750')
    plt.plot(np.arange(11,21,1),d8T8750/d8T8750[0],color=sb.xkcd_rgb['red'],ls='dashdot',label='Density: 8 Temp: 8750')
    plt.plot(np.arange(11,21,1),d8T15000/d8T15000[0],color=sb.xkcd_rgb['red'],ls='solid',label='Density: 8 Temp: 15000')
    plt.plot(np.arange(11,21,1),d10T3750/d10T3750[0],color=sb.xkcd_rgb['forest green'],ls='dashed',label='Density: 10 Temp: 3750')
    plt.plot(np.arange(11,21,1),d10T8750/d10T8750[0],color=sb.xkcd_rgb['forest green'],ls='dashdot',label='Density: 10 Temp: 8750')
    plt.plot(np.arange(11,21,1),d10T15000/d10T15000[0],color=sb.xkcd_rgb['forest green'],ls='solid',label='Density: 10 Temp: 15000')
    plt.plot(np.arange(11,21,1),d12T3750/d12T3750[0],color=sb.xkcd_rgb['purple'],ls='dashed',label='Density: 12 Temp: 3750')
    plt.plot(np.arange(11,21,1),d12T8750/d12T8750[0],color=sb.xkcd_rgb['purple'],ls='dashdot',label='Density: 12 Temp: 8750')
    plt.plot(np.arange(11,21,1),d12T15000/d12T15000[0],color=sb.xkcd_rgb['purple'],ls='solid',label='Density: 12 Temp: 15000')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=4,frameon=True,facecolor=None,framealpha=1,fontsize=16)
    
    plt.savefig(folder+title+' Decrement.pdf',bbox_inches='tight',dpi=300)

    import matplotlib as mpl
    mpl.rcParams.update(mpl.rcParamsDefault)
def DR15_Plots():

    import functions
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import pandas as pd
    import numpy as np
    import itertools
    from PyAstronomy.pyasl import helcorr as helcorr
    import os
    from matplotlib.backends.backend_pdf import PdfPages


    filepath = '/Users/ballanr/Desktop/File Outputs/DR15/pre-DR15 Catalog1.csv'
    openfile = pd.read_csv(filepath)
    g = 0

    for index,row in itertools.islice(openfile.iterrows(),0,None):

        g += 1
        print(g)

        twomass = row['2Mass ID']
        plate = row['Plate']
        mjd = row['MJD']
        fiber = row['Fiber']
        density = row['Density']
        temp = row['Temp']

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

        title = str(plate) + '-' + str(mjd) + '-' + str(fiber)

        #Calculate Equivs and Errors for decrement plot
        lines = [11,12,13,14,15,16,17,18,19,20]
        EqW = []
        Errors = []
        rest_waves = []
        centerlines = []
        for m in range(len(lines)):
            equivs,Fluxcontinuum,shift,rest_wavelength,centerline = functions.Br_EqW(wave,newflux,lines[m],vbc)
            rest = rest_wavelength*(10**10)
            dEqW = functions.Br_Error(wave,newflux,error,lines[m],rest)
            EqW.append(equivs)
            Errors.append(dEqW)
            rest_waves.append(rest)
            centerlines.append(centerline)

        NGC2264 = [6103,6227]
        IC348 = [6218,6219,6220,6221,6222,6223,7073,7074,7075,7079]
        NGC1333 = [6224,6225,6226,7770,7771,7772]
        orion_A = [7220,7221,7222,7223,7224,7225,7226,7227,7228,7229,7230,7231,7232,7233,7234,9659,9660,9661]
        lam_ori = [8879,8880,8881,8882,8883,8884,8885,8886,9482]
        orion_B = [8890,8891,8892,8893,8894,8895,8896,8897,8898,8899]
        orion_OB1AB = [8900,8901,8902,8903,8904,8905,8906,9468,9469,9470,9471,9472,9473,9474,9475,9476,9477]
        pleiades = [8888,9257]
        taurus = [9258,9259,9287,9288]
        perseus = [9662]
        W3_4 = [9245,9246,9247,9248,9249]
        W5 = [9542]
        J305 = [9753]

        if int(plate) in NGC2264:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/NGC_2264/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in IC348:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/IC_348/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in NGC1333:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/NGC_1333/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in orion_A:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/Orion_A/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in lam_ori:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/Lambda_Ori/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in orion_B:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/Orion_B/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in orion_OB1AB:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/Orion_OB1AB/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)
        
        elif int(plate) in pleiades:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/Pleiades/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)
            
        elif int(plate) in taurus:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/Taurus/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in perseus:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/Perseus/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in W3_4:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/W3_4/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in W5:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/W5/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        elif int(plate) in J305:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/305-00/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        else:
            file_path = '/Users/ballanr/Desktop/File Outputs/DR15/Unknown/Plots/'+str(twomass)+'/'
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)
        
        #Br11 Plot
        with PdfPages(file_path+title+ ' Spectra.pdf') as pdf:
            
            #Full Spectra plot
            minflux = min(newflux) - 100
            maxflux = max(flux) + 100

            plt.figure(figsize=(20,10))
            plt.ylabel('Flux',fontsize=18)
            plt.xlabel('Wavelength',fontsize=18)
            plt.title(title + ' Spectra',fontsize=20)
            plt.plot(wave,newflux)
            for k in range(len(rest_waves)):
                plt.axvline(rest_waves[k],ls='dashed',color='red',linewidth=1)
            plt.ylim(minflux,maxflux)
            pdf.savefig(bbox_inches='tight',dpi=300)
            plt.close()

            #Brackett Lines Plot
            for i in range(10):
                L1 = rest_waves[i] - 27.42 # ~ 27.42 Angstroms
                L2 = rest_waves[i] - 17.21 # ~ 17.21 Angstroms
                R1 = rest_waves[i] + 17.24
                R2 = rest_waves[i] + 27.42
                restR = rest_waves[i] + 40
                restL = rest_waves[i] - 40
                max1 = functions.find_nearest(wave,restL)
                max2 = functions.find_nearest(wave,restR)
            
                maxx = (1.1)*max(newflux[max1:max2])
                minx = (0.9)*min(newflux[max1:max2])

                
                plt.figure(figsize=(20,10))
                
                plt.plot(wave,newflux)
                plt.fill_between((L1,L2),minx-100,maxx+100,color='green',alpha=0.5)
                plt.fill_between((R1,R2),minx-100,maxx+100,color='green',alpha=0.5)
                plt.axvline(rest_waves[i],ls='dashed',color='red')

                plt.ylabel('Flux',fontsize=18)
                plt.xlabel('Wavelength',fontsize=18)
                plt.title(title + ' Br'+str(11+i)+' Line',fontsize=20)
                plt.ylim(minx,maxx)
                plt.xlim(restL,restR)
                plt.xticks(np.linspace(restL,restR,10))

                pdf.savefig(bbox_inches='tight',dpi=300)
                plt.close()

        #Decrement Plot
        functions.Brackett_Decrement_Modified(EqW,Errors,density,temp,file_path,title)

        plt.clf()
        plt.close()
def Catalog_Field_Seperator():

    import pandas as pd
    import itertools

    NGC2264 = [6103,6227]
    IC348 = [6218,6219,6220,6221,6222,6223,7073,7074,7075,7079]
    NGC1333 = [6224,6225,6226,7070,7071,7072]
    orion_A = [7220,7221,7222,7223,7224,7225,7226,7227,7228,7229,7230,7231,7232,7233,7234,9659,9660,9661]
    lam_ori = [8879,8880,8881,8882,8883,8884,8885,8886,9482]
    orion_B = [8890,8891,8892,8893,8894,8895,8896,8897,8898,8899]
    orion_OB1AB = [8900,8901,8902,8903,8904,8905,8906,9468,9469,9470,9471,9472,9473,9474,9475,9476,9477]
    pleiades = [8888,9257]
    taurus = [9258,9259,9287,9288]
    perseus = [9662]
    W3_4 = [9245,9246,9247,9248,9249]
    W5 = [9542]
    J305 = [9753]

    filepath = '/Users/ballanr/Desktop/File Outputs/DR15/pre-DR15 Catalog1.csv'
    openfile = pd.read_csv(filepath)

    cols = ['Location ID','2Mass ID','Plate','MJD','Fiber','RA','DEC','SNR','Density',
                'Temp','Confidence','Br11 EqW','Br11 EqW Error','Standard Deviation']

    dNGC2264 = pd.DataFrame(columns = cols)
    dIC348 = pd.DataFrame(columns = cols)
    dNGC1333 = pd.DataFrame(columns = cols)
    dorion_A = pd.DataFrame(columns = cols)
    dlam_ori = pd.DataFrame(columns = cols)
    dorion_B = pd.DataFrame(columns = cols)
    dorion_OB1AB = pd.DataFrame(columns = cols)
    dpleiades = pd.DataFrame(columns = cols)
    dtaurus = pd.DataFrame(columns = cols)
    dperseus = pd.DataFrame(columns = cols)
    dW3_4 = pd.DataFrame(columns = cols)
    dW5 = pd.DataFrame(columns = cols)
    dJ305 = pd.DataFrame(columns = cols)

    for index,row in itertools.islice(openfile.iterrows(),0,None):

        if (row['Plate']) in NGC2264:
            dNGC2264.loc[len(dNGC2264)+1] = row
        elif (row['Plate']) in IC348:
            dIC348.loc[len(dIC348)+1] = row
        elif (row['Plate']) in NGC1333:
            dNGC1333.loc[len(dNGC1333)+1] = row
        elif (row['Plate']) in orion_A:
            dorion_A.loc[len(dorion_A)+1] = row
        elif (row['Plate']) in lam_ori:
            dlam_ori.loc[len(dlam_ori)+1] = row
        elif (row['Plate']) in orion_B:
            dorion_B.loc[len(dorion_B)+1] = row
        elif (row['Plate']) in orion_OB1AB:
            dorion_OB1AB.loc[len(dorion_OB1AB)+1] = row
        elif (row['Plate']) in pleiades:
            dpleiades.loc[len(dpleiades)+1] = row
        elif (row['Plate']) in taurus:
            dtaurus.loc[len(dtaurus)+1] = row
        elif (row['Plate']) in perseus:
            dperseus.loc[len(dperseus)+1] = row
        elif (row['Plate']) in W3_4:
            dW3_4.loc[len(dW3_4)+1] = row
        elif (row['Plate']) in W5:
            dW5.loc[len(dW5)+1] = row
        elif (row['Plate']) in J305:
            dJ305.loc[len(dJ305)+1] = row

    dNGC2264.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/NGC_2264 Catalog.csv',index=False)
    dIC348.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/IC_348 Catalog.csv',index=False)
    dNGC1333.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/NGC_1333 Catalog.csv',index=False)
    dorion_A.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Orion_A_Catalog.csv',index=False)
    dlam_ori.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Lambda_Ori_Catalog.csv',index=False)
    dorion_B.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Orion_B_Catalog.csv',index=False)
    dorion_OB1AB.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Orion_OB1AB_Catalog.csv',index=False)
    dpleiades.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Pleiades_Catalog.csv',index=False)
    dtaurus.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Taurus_Catalog.csv',index=False)
    dperseus.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Perseus_Catalog.csv',index=False)
    dW3_4.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/W3_4_Catalog.csv',index=False)
    dW5.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/W5_Catalog.csv',index=False)
    dJ305.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/305-00_Catalog.csv',index=False)
def DR14skylines_cleaner():
    
    import functions
    import numpy as np
    import pandas as pd
    from astropy.io import fits
    import itertools

    filepath = '/Users/ballanr/Desktop/mdwarves.csv'
    openfile = pd.read_csv(filepath)

    for index,row in itertools.islice(openfile.iterrows(),0,None): 
        
        cols = ['Wavelength']
        df = pd.DataFrame(columns=cols)

        plate = row['Plate']
        MJD = row['MJD']
        fiber = row['Fiber']

        if len(str(row['Fiber'])) == 3:
            fiber = str(row['Fiber'])
        elif len(str(row['Fiber'])) == 2:
            fiber = '0' + str(row['Fiber']) 
        else:
            fiber = '00' + str(row['Fiber'])

        server = '/Volumes/CoveyData/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/apo25m/'
        filepath = str(plate) + '/' + str(MJD) + '/apVisit-r8-' + str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.fits'
        fitsfile = server + filepath
        filename = str(plate)+'-'+str(MJD)+'-'+str(fiber)

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
        
        c = 299792.458
        lamshift = 1 + (vbc/c)
            
        fspec = openfits[1]
        fwave = openfits[4]
            
        wave = []
        flux = []

            
        for i in range(len(fwave.data[2])):
            wave.append(fwave.data[2][-i-1])
            flux.append(fspec.data[2][-i-1])

        for j in range(len(fwave.data[1])):
            wave.append(fwave.data[1][-j-1])
            flux.append(fspec.data[1][-j-1])

        for k in range(len(fwave.data[0])):
            wave.append(fwave.data[0][-k-1])
            flux.append(fspec.data[0][-k-1])

            
        openfits.close()        
                   
        windows = []
        

        removal_beginning = [15185.0,15240.0,15286.0,15330.0,15393.6,15429.9,15460.8,15473.3,15499.6,15508.0,15516.0,
                            15538.5,15545.0,15568.9,15595.9,15652.9,15701.3,15758.9,15861.5,15867.1,15874.3,15888.8,
                            15971.1,16028.8,16078.2,16126.4,16193.1,16233.4,16269.5,16278.8,16301.1,16314.8,16339.1,
                            16349.3,16356.3,16358.6,16387.1,16413.5,16475.1,16476.9,16478.3,16501.4,16527.4,16530.4,
                            16542.7,16552.9,16558.7,16585.4,16609.4,16611.6,16654.1,16658.2,16666.4,16688.2,16691.2,
                            16701.7,16707.7,16711.8,16717.7,16723.1,16725.1,16730.9,16753.4,16755.9,16761.5,16764.8,
                            16839.5,16852.8,16877.2,16883.0,16888.8,16890.3,16899.6,16902.4,16907.1,16908.6,16909.9,16913.7]

        removal_ending = [15189.0,15242.6,15289.0,15334.7,15396.6,15433.9,15463.8,15475.3,15502.3,15511.0,
                        15518.8,15543.4,15547.6,15571.4,15599.3,15658.9,15704.3,15761.8,15863.6,15871.0,
                        15877.0,15892.7,15974.4,16033.1,16081.6,16131.0,16196.5,16237.9,16271.4,16280.7,
                        16304.7,16318.7,16343.0,16353.7,16357.2,16361.7,16390.3,16416.2,16476.3,16477.9,
                        16479.7,16503.7,16528.9,16531.1,16543.5,16555.1,16560.2,16587.5,16610.7,16612.8,
                        16655.0,16661.0,16667.3,16690.4,16693.9,16703.8,16710.5,16712.5,16719.7,16724.5,
                        16725.9,16734.9,16755.2,16756.8,16762.5,16765.9,16841.6,16853.3,16877.8,16884.5,
                        16889.5,16891.2,16900.2,16905.5,16907.7,16909.3,16910.6,16914.3]

        for i in range(len(removal_beginning)):

            lwindow = functions.find_nearest(wave,removal_beginning[i])
            rwindow = functions.find_nearest(wave,removal_ending[i])

            windows.append((lwindow,rwindow))

        for i in range(len(removal_beginning)):

            lwindowelement = int(windows[i][0])
            rwindowelement = int(windows[i][1])
            leftwindow = wave[windows[i][0]]
            rightwindow = wave[windows[i][1]]
            leftflux = flux[windows[i][0]]
            rightflux = flux[windows[i][1]]

            slope = (rightflux - leftflux) / (rightwindow - leftwindow)

            for k in range(rwindowelement - lwindowelement):

                fluxvalue = slope*(wave[lwindowelement+k] - leftwindow) + leftflux
                flux[lwindowelement+k] = fluxvalue

        
        wave = np.asarray(wave) * lamshift 
        df['Wavelength'] = wave
        df['Flux'] = flux
        df.to_csv('/Users/ballanr/Desktop/JessicasDwarves/'+filename+'.csv',index=False)
        df = df.iloc[0:0]
def Br20_Confidence(file,filecreation,br20eval):

    import functions
    import pandas as pd
    import itertools
    import numpy as np
    from astropy.io import fits
    from PyAstronomy.pyasl import helcorr

    openfile = pd.read_csv(file)
    listnum = 0

    if filecreation == 1:
        for index,row in itertools.islice(openfile.iterrows(),0,None):
            plate = row['Plate']
            MJD = row['MJD']
            if len(str(row['Fiber'])) == 3:
                fiber = str(row['Fiber'])
            elif len(str(row['Fiber'])) == 2:
                fiber = '0' + str(row['Fiber']) 
            else:
                fiber = '00' + str(row['Fiber'])

            serverpath = '/Users/ballanr/Desktop/File Outputs/DR15/Wave and Flux/'
            filestring = str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.csv'

            try:

                fitsfile = pd.read_csv(serverpath+filestring)
                print('File exists!')

            except: 

                print('Creating file!')

                if int(plate) < 9700:
                    try:
                        server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                        filepath = str(plate) + '/' + str(MJD) + '/apVisit-apogee2-' + str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.fits'
                        fitsfile = server + filepath
                        openfits = fits.open(fitsfile)
                        header = openfits[0].header
                    except:
                        server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                        filepath = str(plate) + '/' + str(MJD) + '/apVisit-r8-' + str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.fits'
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
                    filepath = str(plate) + '/' + str(MJD) + '/asVisit-apogee2-' + str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.fits'
                
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
            
                c = 299792.458
                lamshift = 1 + (vbc/c)

                fspec = openfits[1]
                ferr = openfits[2]
                fwave = openfits[4]
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

                openfits.close()
                
                newflux = functions.skylines_cleaner(wave,flux)

                wave = np.asarray(wave) * lamshift

                df1 = pd.DataFrame(wave,columns=['Wavelength'])
                df1['Flux'] = newflux
                df1['Error'] = error
                filename = str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.csv'
                df1.to_csv('/Users/ballanr/Desktop/File Outputs/DR15/Wave and Flux/'+filename,index=False)
                df1 = df1.iloc[0:0]   


    if br20eval == 1:

        ones = []
        twos = []
        threes = []
        fours = []
        others = []
    
        for index,row in itertools.islice(openfile.iterrows(),0,None):

            plate = row['Plate']
            MJD = row['MJD']
            if len(str(row['Fiber'])) == 3:
                fiber = str(row['Fiber'])
            elif len(str(row['Fiber'])) == 2:
                fiber = '0' + str(row['Fiber']) 
            else:
                fiber = '00' + str(row['Fiber'])
            category = row['Photospheric Contamination*']
            br11 = row['Br11 EqW']
            confidence = row['Confidence'] 

            serverpath = '/Users/ballanr/Desktop/File Outputs/DR15/Wave and Flux/'
            filestring = str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.csv'

            fluxfile = pd.read_csv(serverpath+filestring)

            wave = fluxfile['Wavelength']
            flux = fluxfile['Flux']

            if int(plate) < 9700:
                try:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                    filepath = str(plate) + '/' + str(MJD) + '/apVisit-apogee2-' + str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.fits'
                    fitsfile = server + filepath
                    openfits = fits.open(fitsfile)
                    header = openfits[0].header
                except:
                    server = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
                    filepath = str(plate) + '/' + str(MJD) + '/apVisit-r8-' + str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.fits'
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
                filepath = str(plate) + '/' + str(MJD) + '/asVisit-apogee2-' + str(plate) + '-' + str(MJD) + '-' + str(fiber) + '.fits'
            
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
            
            c = 299792.458
            lamshift = 1 + (vbc/c)
           
            #print(listnum)
            #listnum +=1

            observed_wavelength,shift,rest_wavelength = functions.Barycentric_Correction(20,vbc)
            rest = rest_wavelength*(10**10)
            centerline = functions.find_nearest(wave,rest)

            try:
                #regular windows

                L1 = centerline - 301 # ~ 27.42 Angstroms
                L2 = centerline - 150 # ~ 17.21 Angstroms
                R1 = centerline + 150
                R2 = centerline + 301

                a = np.asarray(flux[L2:R1])
                b = np.asarray(flux[L1:L2])
                c = np.asarray(flux[R1:R2])
                #Equivs

                avgc = (np.sum(b) + np.sum(c))/(len(b)+len(c))
                EqW1 = 0

                if avgc == 0:

                    EqW1 = 0
                    EqW1_rounded = 0
                    equivs = 0
            
                if avgc != 0:

                    for i in range(L2,R1):

                        trapezoid = (0.5)*(wave[i+1] - wave[i])*(flux[i+1] + flux[i] - (2*avgc))
                        EqW1 += trapezoid
            
                equivs = EqW1/avgc
        

                #CALCULATIONS
                try:
                    fmax = max(a)
                except:
                    fmax = 0
                try:
                    continuum = np.std(b+c)
                except:
                    continuum = np.std(c)

                brcheck = fmax/continuum

            except:
                brcheck = 0
                equivs = 0
            
            if category == 1:
                ones.append((br11,confidence,brcheck,equivs))
            elif category == 2:
                twos.append((br11,confidence,brcheck,equivs))
            elif category == 3:
                threes.append((br11,confidence,brcheck,equivs))
            elif category == 4:
                fours.append((br11,confidence,brcheck,equivs))
            else:
                others.append((br11,confidence,brcheck,equivs))

        onesarray = np.asarray(ones)
        twosarray = np.asarray(twos)
        threesarray = np.asarray(threes)
        foursarray = np.asarray(fours)
        othersarray = np.asarray(others)
        return onesarray,twosarray,threesarray,foursarray,othersarray
        #print(len(br20))