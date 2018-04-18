def Master_Catalog(full_list):

    import pandas as pd
    import itertools

    ##### Read in file list and create csv file
    openfile = pd.read_csv(full_list)

    cols = ['Location ID', '2Mass ID', 'Plate', 'MJD', 'Fiber', 'RA', 'DEC','Density','Temp','Chi','Br11 EqW',
            'Br12 EqW','Br13 EqW','Br14 EqW','Br15 EqW','Br16 EqW','Br17 EqW','Br18 EqW','Br19 EqW','Br20 EqW','Br11 Error',
            'Br12 Error','Br13 Error','Br14 Error','Br15 Error','Br16 Error','Br17 Error','Br18 Error','Br19 Error','Br20 Error']
    df = pd.DataFrame(columns = cols)

    df_fails = pd.DataFrame(columns = ['Location ID', '2Mass ID','Plate','MJD','Fiber'])

    ##### Open fits to get: 2MID, location, plate, mjd, and fiber
    for index,row in itertools.islice(openfile.iterrows(),0,None):
        
        try:
        
            loc = row['Location ID']
            twomass = row['2Mass ID']

            plate = row['Plate']
            mjd = row['MJD']
            fiber = row['Fiber']

            if len(str(fiber)) == 2:
                fiber = '0' + str(fiber)
            elif len(str(fiber)) == 1:
                fiber = '00' + str(fiber)
            else:
                fiber = str(fiber)

            ra = row['RA']
            dec = row['DEC']
            #byhand = row['byhand']
            #maxline = row['max order']

            csvfilepath = '/Users/ballanr/Desktop/Research/DR15/Spectra Files/Emitters/' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
            csvfile = pd.read_csv(csvfilepath)

            cwave = csvfile['Wavelength']
            cflux = csvfile['Flux']
            cerrors = csvfile['Error']
            csnr = csvfile['SNR']

            print(plate,mjd,fiber)

            ##### Calculate equivalent width (EqW) and error
            equivs, errors = Equivalent_Width(cwave,cflux,cerrors,csnr)

            ##### Calculate confidence score

            ##### Calculate decrement
            index, chi = Decrement_Model(equivs,errors)

            ##### Output: 2MID, location, plate, mjd, fiber, density, temperature, chi, equivs, equiv errors
            if index.startswith('1'):
                density = index[:4]
                temp = index[5:]
            else:
                density = index[:3]
                temp = index[4:]

            ##### Save catalog
            data = [loc,twomass,plate,mjd,fiber,ra,dec,density,temp,chi,equivs[0],equivs[1],
                        equivs[2],equivs[3],equivs[4],equivs[5],equivs[6],equivs[7],equivs[8],equivs[9],errors[0],errors[1],
                        errors[2],errors[3],errors[4],errors[5],errors[6],errors[7],errors[8],errors[9]]
            df.loc[len(df)+1] = data

        except:
            data = [loc,twomass,plate,mjd,fiber]
            df_fails.loc[len(df_fails)+1] = data

    ''' Output '''

    ##### Output to local
    df.to_csv('/Users/ballanr/Desktop/Research/Master_test.csv',index=False)
    df_fails.to_csv('/Users/ballanr/Desktop/Research/Master List Failures.csv',index=False)

    ##### Output to server
    df.to_csv('/Volumes/CoveyData/APOGEE_Spectra/Richard/Master_List.csv',index=False)
    df_fails.to_csv('/Volumes/CoveyData/APOGEE_Spectra/Richard/Master_List_Failures.csv',index=False)

def File_Creator(plate,mjd,fiber,version):

    ''' Importing Packages '''

    import pandas as pd
    from astropy.io import fits
    import numpy as np

    ''' Initialize data frame '''

    cols = ['Wavelength','Flux','Error']
    df = pd.DataFrame(columns = cols)

    ''' Extract file path: plate ID, MJD, and fiber '''
    # Will need to change depending on how the server list is made

    plate = str(plate)
    mjd = str(mjd)

    # The fiber is given in three digits but if it leads with 0s they get dropped,
    # so we need to append them back on

    if len(str(fiber)) == 3:
        fiber = str(fiber)
    elif len(str(fiber)) == 2:
        fiber = '0' + str(fiber) 
    else:
        fiber = '00' + str(fiber)

    ''' Open .fits file and get information from headers '''
    
    wavedata,fluxdata,ferrordata,vbc,ra,dec = Fits_Info(plate,mjd,fiber,version)

    ''' Stitch together one array for wavelength, one for flux, and one for error '''

    wave = []
    flux = []
    error = []

    # The fits files are ordered backwards and in seperate arrays. So we reorder from left to right
    # and put everything in one array for ease of use
    for k in range(3):
        j = 2 - k
        for i in range(len(fluxdata[j])):
            wave.append(wavedata[j][-(i+1)])
            flux.append(fluxdata[j][-(i+1)])
            error.append(ferrordata[j][-(i+1)])

    # Convert to np arrays
    wave = np.asarray(wave)
    flux = np.asarray(flux)
    error = np.asarray(error)

    ''' Barycentric correction '''
    
    # apVisit spectra are given with no shift applied so we need to do that sometimes the shift isn't quite what would make sense visually, 
    # either too small or too large
    c = 299792.458
    wave_shift = 1 + (vbc/c)
    wave = wave * wave_shift

    ''' Clean skylines '''

    cwave,cflux,cerror = Skylines_Cleaner(wave,flux,error,wave_shift)

    ''' Chipgaps '''

    chipwave,chipflux,chiperror = Chipgap(cwave,cflux,cerror,wave_shift)

    ''' SNR '''
    SNR = []
    for g in range(len(chipwave)):
       snr = chipflux[g] / chiperror[g]
       SNR.append(snr) 

    ''' Save csv file '''

    savestring = str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
    df['Wavelength'] = chipwave
    df['Flux'] = chipflux
    df['Error'] = chiperror
    df['SNR'] = SNR
    df.to_csv('/Users/ballanr/Desktop/Research/DR15/Spectra Files/All/' + savestring, index=False)
    #df.to_csv('/Users/ballanr/Desktop/Research/DR15/Spectra Files/' + savestring, index=False)
    df = df.iloc[0:0]

def Fits_Info(plate,mjd,fiber,version):

    from astropy.io import fits
    from PyAstronomy.pyasl import helcorr

    # Now we stitch together the server path to each file. Careful though as there is a folder for the APO telescope and one for LCO that   
    # contain different field. So we will check against LCO as it contains fewer plates

    lco = [9753,9761,9762,9856,9906,9907,9908]
    apopath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/'
    lcopath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/lco25m/'

    # There are also some new t9 files but can ignore for now
    #t9apopath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/apo25m/9069/57822/apVisit-t9-9069-57822-002.fits'
    int_plate = int(plate)

    # Setting path to file
    if int_plate in lco:
        filepath = lcopath + str(plate) + '/' + str(mjd) + '/asVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
    else:
        #filepath = apopath + str(plate) + '/' + str(mjd) + '/apVisit-apogee2-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
        filepath = apopath + str(plate) + '/' + str(mjd) + '/apVisit-' + str(version) + '-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
    
    # Now we can open the fits file on the server
    try:
        fitsfile = fits.open(filepath)
        header = fitsfile[0].header
    except:
        # filepath = apopath + str(plate) + '/' + str(mjd) + '/apVisit-r8-' + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.fits'
        # fitsfile = fits.open(filepath)
        print('Bad file path!')

    ''' Get: wavegrid, flux, error, VHELIO, RA, DEC, TELESCOP '''

    # Grabbing the data for the wavelength and flux
    wavedata = fitsfile[4].data
    fluxdata = fitsfile[1].data
    ferrordata = fitsfile[2].data

    # Getting location information and the barycentric correction
    ra = header['RA']
    dec = header['DEC']
    jd = header['JD-MID']
    telescope = header['TELESCOP']

    # Sometimes there is no reported BC so we need to get one from the helcorr package
    try:
        vbc = header['BC']

    except:
        if telescope == 'apo25m':
            vbc,hjd = helcorr(-105.4913,36.4649,2788,ra,dec,jd)
        else:
            vbc,hjd = helcorr(-70.413336,-29.05256,2380,ra,dec,jd)

    fitsfile.close()
    
    return(wavedata,fluxdata,ferrordata,vbc,ra,dec)

def Skylines_Cleaner(wave,flux,error,wave_shift):

    import numpy as np

    # List of skylines provided by Kevin, mostly OH lines
    lremoval = [15185.0,15240.0,15286.0,15330.0,15393.6,15429.9,15460.8,15473.3,15499.6,15508.0,15516.0,15538.5,15545.0,15568.9,15595.9,
                        15652.9,15701.3,15758.9,15861.5,15867.1,15874.3,15888.8,15971.1,16028.8,16078.2,16126.4,16193.1,16233.4,16269.5,16278.8,
                        16301.1,16314.8,16339.1,16349.3,16356.3,16358.6,16387.1,16413.5,16475.1,16476.9,16478.3,16501.4,16527.4,16530.4,16542.7,
                        16552.9,16558.7,16585.4,16609.4,16611.6,16654.1,16658.2,16666.4,16688.2,16691.2,16701.7,16707.7,16711.8,16717.7,16723.1,
                        16725.1,16730.9,16753.4,16755.9,16761.5,16764.8,16839.5,16852.8,16877.2,16883.0,16888.8,16890.3,16899.6,16902.4,16907.1,
                        16908.6,16909.9,16913.7]

    rremoval = [15189.0,15242.6,15289.0,15334.7,15396.6,15433.9,15463.8,15475.3,15502.3,15511.0,15518.8,15543.4,15547.6,15571.4,15599.3,
                    15658.9,15704.3,15761.8,15863.6,15871.0,15877.0,15892.7,15974.4,16033.1,16081.6,16131.0,16196.5,16237.9,16271.4,16280.7,
                    16304.7,16318.7,16343.0,16353.7,16357.2,16361.7,16390.3,16416.2,16476.3,16477.9,16479.7,16503.7,16528.9,16531.1,16543.5,
                    16555.1,16560.2,16587.5,16610.7,16612.8,16655.0,16661.0,16667.3,16690.4,16693.9,16703.8,16710.5,16712.5,16719.7,16724.5,
                    16725.9,16734.9,16755.2,16756.8,16762.5,16765.9,16841.6,16853.3,16877.8,16884.5,16889.5,16891.2,16900.2,16905.5,16907.7,
                    16909.3,16910.6,16914.3]

    emission_lines = [16811.179, 16411.74, 16113.778, 15884.944, 15705.015, 15560.761, 15443.2, 15346.043, 15264.768, 15196.056]

    # The lines are given in Earth's frame and since we moved to the barycentric frame we need to shift these lines there as well
    lremoval = np.asarray(lremoval) * wave_shift
    rremoval = np.asarray(rremoval) * wave_shift
    emission_lines = np.asarray(emission_lines) * wave_shift

    # Loop over skyline arrays
    for i in range(len(lremoval)):
        # Setting up conditional for whether we are dealing with an emission line or not
        emission = False

        for b in range(10):
            ldifference = np.sqrt((lremoval[i] - emission_lines[b])**2)
            rdifference = np.sqrt((lremoval[i] - emission_lines[b])**2)

            if ldifference < 6 or rdifference < 6:
                emission = True

        ''' Setting up windows '''
        # Setting windows to average on either side of the skyline
        left = Find_Nearest(wave,lremoval[i])
        right = Find_Nearest(wave,rremoval[i])
        width = right - left

        ''' Non-emission '''
        if not emission:

            ''' Medians '''

            window_size = 3 #5,7
            L1 = left - window_size
            L2 = left
            R1 = right
            R2 = right + window_size

            leftwindow = flux[L1:L2]
            rightwindow = flux[R1:R2]

            ''' Median Error Propagation '''

            median = 0.5 * (leftwindow[2] + rightwindow[0])
            med_error = np.sqrt(error[L1 + 2]**2 + error[R1 + 1]**2)

            ''' Changing flux and error arrays '''

            # Setting flux and error values of skyline to median values
            for k in range(width):
                flux[L2 + k] = median
                error[L2 + k] = med_error

        ''' Emission Line Interpolation'''
        if emission:
            # Setting values to the median doesn't work in emission features as it chops a good deal of area out of the emission feature. So, we need 
            # to do interpolation in the emission windows instead

            # Calculating slope and slope error for interpolation
            delx = wave[right] - wave[left]
            dely = flux[right] - flux[left]
            slope = dely / delx

            slope_err = slope * (np.sqrt(error[left]**2 + error[right]**2) / dely)

            # Applying slope and error
            for a in range(width):
                flux[left + a] = slope*(wave[left + a] - wave[left]) + flux[left]
                error[left + a] = slope_err
           
    return(wave,flux,error)

def Find_Nearest(array,value): 
    
    import numpy as np

    ''' Find the element in the array that is closest to the value '''

    array = np.asarray(array)
    index = (np.abs(array-value)).argmin()
    return index

def Chipgap(wave,flux,errors,wave_shift):

    import numpy as np

    # unshifted chip gap = 16435.00768 -> 16472.31973
    del_x = 0.130676344
    
    ''' Brackett 12 Gap '''

    # Unshifted elements of the Br12 chipgap 
    br12lshift = 16432.0050673565 * wave_shift
    br12rshift = 16473.99556 * wave_shift

    br12lchip = Find_Nearest(wave,br12lshift) + 1

    insert_array1 = np.asarray(np.arange(br12lshift,br12rshift,del_x))
    insert_array2 = np.asarray([1 for i in range(322)])
    wave = np.insert(wave,br12lchip,insert_array1)
    flux = np.insert(flux,br12lchip,insert_array2)
    errors = np.insert(errors,br12lchip,insert_array2)

    br12rchip = Find_Nearest(wave,br12rshift) + 1

    br12chiprange = flux[br12lchip:br12rchip]

    #Changing values to garbage number
    for k in range(len(br12chiprange)):
        flux[br12lchip + k] = -9999

    ''' Brackett 15 Gap '''

    # Unshifted elements of the Br12 chipgap 
    br15lshift =  15808.56714 * wave_shift
    br15rshift =  15858.28309 * wave_shift

    br15lchip = Find_Nearest(wave,br15lshift) +1

    insert_array3 = np.asarray(np.arange(br15lshift,br15rshift,del_x))
    insert_array4 = np.asarray([1 for i in range(381)])
    wave = np.insert(wave,br15lchip,insert_array3)
    flux = np.insert(flux,br15lchip,insert_array4)
    errors = np.insert(errors,br15lchip,insert_array4)

    br15rchip = Find_Nearest(wave,br15rshift) + 1

    br15chiprange = flux[br15lchip:br15rchip]

    #Changing values to garbage number
    for k in range(len(br15chiprange)):
        flux[br15lchip + k] = -9999

    return(wave,flux,errors)

def Rest_Emission(emission_line):

    ''' Calcultaing reduced mass Rydberg constant '''

    n = (float(emission_line))**2 #Beginning electron level
    c = 299792.458 #Speed of light (km/s)
    rydberg_inf =  1.0973731568539*(10**7) #Rydberg constant (m^-1)
    electron = 9.10938356*(10**-31) #Mass of electron (kg)
    nucleus = 1.672621898*(10**-27) #Mass of hydrogen nucleus (kg)
    
    ''' Calculating rest emission wavelength '''
    
    rydberg_red = rydberg_inf/(1+(electron/nucleus)) #Reduced mass Rydberg constant (m^-1)
    rest_wavelength1 = rydberg_red*((1./16.)-(1./n)) #(m^-1)
    rest_wavelength = (1/rest_wavelength1)*(10**10) #wavelength (m^-10)

    return rest_wavelength

def Equivalent_Width(wave,flux,errors,snr):

    import pandas as pd
    import numpy as np

    ''' Setting local continuum windows '''
    # Files have already been shifted so no need to take that into account here
    equivs = []
    equiverr = []
    # Setting windows
    for k in range(10):
        line = 11 + k
        rest = Rest_Emission(line)
        elrest = Find_Nearest(wave,rest)
        
        # if line == 12:
        #     L1 = elrest - 281
        #     L2 = elrest - 140
        #     R1 = elrest + 110
        #     R2 = elrest + 140
        # elif line == 14:
        #     L1 = elrest - 171
        #     L2 = elrest - 80
        #     R1 = elrest + 110
        #     R2 = elrest + 140
        # else:
        #     L1 = elrest - 281
        #     L2 = elrest - 140
        #     R1 = elrest + 140
        #     R2 = elrest + 281

        L1 = elrest - 241
        L2 = elrest - 120
        R1 = elrest + 120
        R2 = elrest + 241

        ''' Calculating local continuum via averages '''

        # local_left = np.sum(flux[L1:L2]) / len(flux[L1:L2])
        # local_right = np.sum(flux[R1:R2]) / len(flux[R1:R2])
        # local = 0.5 * (local_left + local_right)

        ''' Calculating local continuum via median '''
        # Setting up windows
        left_win = np.asarray(flux[L1:L2])
        right_win = np.asarray(flux[R1:R2])
        left_win_err = np.asarray(errors[L1:L2])
        right_win_err = np.asarray(errors[R1:R2])
        left_snr = np.asarray(snr[L1:L2])
        right_snr = np.asarray(snr[R1:R2])

        # Only import non-junk data
        left_win_err = left_win_err[left_win != -9999]
        right_win_err = right_win_err[right_win != -9999]
        left_snr = left_snr[left_win != -9999]
        right_snr = right_snr[right_win != -9999]
        left_win = left_win[left_win != -9999]   
        right_win = right_win[right_win != -9999]

        # Only import data above an SNR of 10
        left_win_err = left_win_err[left_snr > 10]
        right_win_err = right_win_err[right_snr > 10]
        left_win = left_win[left_snr > 10]   
        right_win = right_win[right_snr > 10]
        left_snr = left_snr[left_snr > 10]
        right_snr = right_snr[right_snr > 10]

        # Median
        local_range = np.concatenate((left_win,right_win))
        local = np.median(local_range)     

        ''' Calculating equivalent widths '''

        num = 0
        for j in range(len(flux[L2:R1])):
            num += 0.5 * (wave[L2 + j + 1] - wave[L2 + j]) * ((flux[L2 + j + 1] + flux[L2 + j]) - 2*local)

        EqW = num / local
        equivs.append(EqW)
        
        ''' Error propagation for averages''' 

        # Calculating the error of the local continuum average
        # local_left_err = np.asarray(errors[L1:L2])**2
        # local_right_err = np.asarray(errors[R1:R2])**2

        # local_left_err = np.sum(local_left_err) / len(errors[L1:L2])
        # local_right_err = np.sum(local_right_err) / len(errors[R1:R2])
        # local_err = np.sqrt((local_left_err + local_right_err) / 2)

        # Calculating error of continuum window median
        med_win = np.concatenate((left_win_err,right_win_err))
        med_win = np.sort(med_win)

        if len(med_win)%2 == 0:
            left_val = (len(med_win) / 2) - 0.5
            right_val = (len(med_win) / 2) + 0.5
            left_err = med_win[int(left_val)]
            right_err = med_win[int(right_val)]
            local_err = np.sqrt((left_err**2 + right_err**2) / 2)

        else:
            element = (len(med_win) - 1) / 2
            local_err = med_win[int(element)]

        # Calculating the error of the EqW calculation
        eqw_err = np.asarray(errors[L2:R1])**2
        delta_y = np.sqrt(np.sum(eqw_err) + 2 * (local_err**2))

        # for i in range(len(errors[L2:R1])):
        #     y = (errors[L2 + i + 1])**2 + (errors[L2 + i])**2 + 2*(local_err**2)
        #     delta_y += y
        
        frac1 = (delta_y / (2*num))**2
        frac2 = (local_err / local)**2
        delta_eqw = EqW * np.sqrt(frac1 + 2*frac2)

        equiverr.append(delta_eqw)
        #print(local,local_err,2*num,delta_y)

    return (equivs,equiverr)

def Decrement_Model(equivs,equiverr):

    import pandas as pd
    import numpy as np

    # Reading in the density and temperature models
    modelpath = '/Users/ballanr/Desktop/Research/DR15/Density_Temp_Files/Profile Test.csv'
    openmodel = pd.read_csv(modelpath)
    cols = openmodel.columns
    headers = cols.tolist()

    equivs = np.asarray(equivs)
    equiverr = np.asarray(equiverr)

    # Array for storing chi information
    Chi = []

    # Calculations
    for i in range(len(cols)):

        chi_squared = 0

        probs = openmodel[cols[i]]

        for k in range(len(probs)):
            # Calculating numerator
            ratio = equivs[k] / equivs[0]
            numerator = ratio - (probs[k] / probs[0])
            numerator = numerator**2

            # Calculating denominator
            sigma = ratio * np.sqrt((equiverr[k] / equivs[k])**2 + (equiverr[0] / equivs[0])**2)
            denominator = (np.sqrt(0.02) * ratio)**2 + sigma**2 + (0.1**2)
            denominator = np.sqrt(denominator)

            # Adding to Chi
            chi_squared += numerator / denominator

        # Append information to Chi
        if headers[i].startswith('1'):
            Chi.append((headers[i][0:4],headers[i][6:],chi_squared))
        
        else:
            Chi.append((headers[i][0:3],headers[i][4:],chi_squared))

    x = min(c for (a,b,c) in Chi)

    for index, item in enumerate(Chi):
        if item[2] == x:
            index1 = headers[index]

    return(index1,x)

def oswalk():

    import os
    import pandas as pd

    serverpath = '/Volumes/CoveyData/APOGEE_Spectra/preDR15/apogee/spectro/redux/visits/'
    cols = ['Ap2s']
    cols1 = ['R8s']
    cols2 = ['T9s']
    df = pd.DataFrame(columns = cols)
    df1 = pd.DataFrame(columns = cols1)
    df2 = pd.DataFrame(columns = cols2)
    r8s = []
    t9s = []
    ap2s = []
    for root, dirs, files in os.walk(serverpath, topdown=False):
        for element in files:
            if len(element) > 10:
                try:
                    check = str(element[8:10])
                    if check == 't9':
                        if len(element) == 31 or len(element) == 30:
                            t9s.append(element)
                        else:
                            print('Too long!')
                    elif check == 'r8':
                        if len(element) == 30:
                            r8s.append(element)
                        else:
                            print('Too long!')
                    elif check == 'ap':
                        if len(element) == 35:
                            ap2s.append(element)
                        else:
                            print('Too long!')      
                except:
                    print('Didn\'t work!')

    df['Ap2s'] = ap2s
    df1['R8s'] = r8s
    df2['T9s'] = t9s
    
    
    df.to_csv('/Users/ballanr/Desktop/Research/DR15/Ap2.csv',index=False)
    df1.to_csv('/Users/ballanr/Desktop/Research/DR15/R8.csv',index=False)
    df2.to_csv('/Users/ballanr/Desktop/Research/DR15/T9.csv',index=False)

    df.to_csv('/Volumes/CoveyData/APOGEE_Spectra/Richard/DR15/Ap2.csv',index=False)
    df1.to_csv('/Volumes/CoveyData/APOGEE_Spectra/Richard/DR15/R8.csv',index=False)
    df2.to_csv('/Volumes/CoveyData/APOGEE_Spectra/Richard/DR15/T9.csv',index=False)

def KDE_Plot(filepath):

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import itertools

    ''' Read in file '''

    openfile = pd.read_csv(filepath)
    data1 = pd.DataFrame(columns = ['Density','Temp'])
    

    ##### For some of the data
    # dens = []
    # temp = []

    # for index,row in itertools.islice(openfile.iterrows(),0,None):
    #     if row['Br11 EqW'] > 6:

    #         if row['Max Order Line'] >= 17:

    #             dens.append(row['Density'])
    #             temp.append(row['Temp'])

    ##### For all the data
    data1[r'Density (n/cm$^{3}$)'] = openfile['Density']
    data1['Temperature (K)'] = openfile['Temp']

    ''' Plotting '''

    sns.set(style='whitegrid',font_scale=2.5)
    plt.figure(figsize=(13,10))
    p = sns.JointGrid(x = r'Density (n/cm$^{3}$)', y = 'Temperature (K)', data=data1, size=10)
    p = p.plot_joint(plt.scatter,s=5,color='black')
    p = p.plot_joint(sns.kdeplot,shade=True,shade_lowest=False,cmap='RdGy',zorder=0,n_levels=15)
    p = p.plot_marginals(sns.distplot,kde=True)

    ''' Output '''

    ##### Save to local
    plt.savefig('/Users/ballanr/Desktop/Research/DR15/Plots/KDE.pdf',bbox_inches='tight',dpi=300)

    ##### Save to server
    plt.savefig('/Volumes/CoveyData/APOGEE_Spectra/Richard/DR15/Plots/KDE.pdf',bbox_inches='tight',dpi=300)

def Brackett_Decrement_Plot(plate,mjd,fiber):

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    if len(str(fiber)) == 2:

        fiber = '0' + str(fiber)
    
    elif len(str(fiber)) == 1:

        fiber = '00' + str(fiber)

    else:

        fiber = str(fiber)

    serverpath = '/Volumes/CoveyData/APOGEE_Spectra/Richard/DR15/Spectra Files/Emitters/'
    filepath = serverpath + str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
    openfile = pd.read_csv(filepath)

    wave = openfile['Wavelength']
    flux = openfile['Flux']
    error = openfile['Error']
    snr = openfile['SNR']

    ''' Calculating EqW and Errors '''

    equivs, errors = Equivalent_Width(wave,flux,error,snr)

    equivs = np.asarray(equivs)
    errors = np.asarray(errors)

    ##### Normalizing eqws and errors by Br 11
    equivs_scaled = equivs / equivs[0]
    errors_scaled = []

    for i in range(len(equivs_scaled)):

        err = equivs_scaled[i] * np.sqrt((errors[i]/equivs[i])**2 + (errors[0]/equivs[0])**2)
        errors_scaled.append(err)

    ''' Plotting '''

    plt.figure(figsize=(13,10))

    #plt.errorbar(np.arange(11,21,1),equivs,errors,color='green',ecolor='red',capsize=5,label='Original')
    plt.errorbar(np.arange(11,21,1),equivs_scaled,errors_scaled,fmt='.-',ecolor='red',capsize=5,label='Normalized')
    
    plt.ylabel(r'Flux$_{n}$ / Flux$_{11}$',fontsize=20)
    plt.xlabel(r'Brackett $11-20$',fontsize=20)
    plt.xticks(np.arange(11,21,1))
    plt.legend(fontsize=18)

    plt.grid(True,alpha=0.3)

    ''' Output '''

    savestring = str(plate) + '-' + str(mjd) + '-' + str(fiber) + '.pdf'

    ##### Save to local
    #plt.savefig('/Users/ballanr/Desktop/Research/DR15/Plots/Decrements/' + savestring,bbox_inches='tight',dpi=300)

    ##### Save to server
    plt.savefig('/Volumes/CoveyData/APOGEE_Spectra/Richard/DR15/Plots/Decrements/' + savestring,bbox_inches='tight',dpi=300)

    plt.clf()
    plt.close()