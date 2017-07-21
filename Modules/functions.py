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
    from numpy import abs
    index = (abs(array-value)).argmin()
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

def Br_Equiv_Width_Apogee(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    from astropy.io import fits
    import functions

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    #Importing header via astropy--------------------------------------------------------------------------------------------------------------------|
    filename = functions.File_Path(plateid,MJD,fiber)
    main_header = fits.open(filename)

    #Barycentric Correction--------------------------------------------------------------------------------------------------------------------------|
    
    #vbcstring = 'BC' + str(1) Somehow need to figure out which visit this MJD applies to
    vbc = main_header[0].header['BC']
    vhelio = main_header[0].header['VHELIO']
    observed_wavelength,shift,rest_wavelength = functions.Barycentric_Correction(emission_line,vbc)
    
    #Equivalent Width Calculation--------------------------------------------------------------------------------------------------------------------|

    #Finding the centerline and checking that it matches the peak in the window
    centerline = find_nearest(wave,observed_wavelength) #Finds the closest element of wave for our observed peak
    centerline_check = Max_Flux_Check(wave,spec,centerline)
    
    '''if emission_line == 11:
        if spec[centerline] != spec[centerline_check[1]]:
            print('The centerline has changed from ' + str(centerline) + ' to ' + str(centerline_check[1]) + ' with a new flux of ' + str(centerline_check[0]) + 
            ' from ' + str(spec[centerline]) + '.')
            centerline = centerline_check[1]'''
    

    L1 = centerline - 240 # ~ 56 Angstroms
    L2 = centerline - 150 # ~ 35 Angstroms
    R1 = centerline + 150
    R2 = centerline + 240

    Fluxcontinuum = (np.sum(spec[L1:L2])+np.sum(spec[R1:R2])) / (len(spec[L1:L2])+len(spec[R1:R2]))
    EqW1 = 0

    if Fluxcontinuum == 0:

        EqW1 = 0
        EqW1_rounded = 0

    if Fluxcontinuum != 0:

        for i in range(L2,centerline):

            left_area = (wave[i+1]-wave[i])*(spec[i+1]-Fluxcontinuum)-(1./2.)*(wave[i+1]-wave[i])*(spec[i+1]-spec[i])
            EqW1 += left_area

        for i in range(centerline,R1):

            right_area = (wave[i+1]-wave[i])*(spec[i]-Fluxcontinuum)-(1./2.)*(wave[i+1]-wave[i])*(spec[i]-spec[i+1])
            EqW1 += right_area

        EqW_rounded = round(EqW1/Fluxcontinuum,5)
        EqW = EqW1/Fluxcontinuum
    
    return EqW,EqW_rounded,vbc,vhelio,Fluxcontinuum,centerline,shift
 
def Br_Equiv_Width_Plotter(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    import matplotlib.pyplot as plt
    import functions

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    #Values for plotter needed-----------------------------------------------------------------------------------------------------------------------|
    EqW,EqW_rounded,vbc,vhelio,Fluxcontinuum,centerline,shift = functions.Br_Equiv_Width(plateid,MJD,fiber,emission_line)

    #Plot averaged spectrum with EqW-----------------------------------------------------------------------------------------------------------------|
    title = str(plateid)+'-'+str(MJD)+'-'+str(fiber)+'-'+str(emission_line)
    fig,ax = plt.subplots(figsize=(16,8))
    plt.plot(wave+shift,spec,linewidth=2.5,label='Shifted')
    
    plt.axhline(y=Fluxcontinuum,ls='dashed',color='black')
    plt.axvline(x=wave[centerline]+shift,ls='dashed',color='r',label='Rest Emission')
    
    plt.legend(loc=1,prop={'size':18})
    plt.xlabel('Wavelength'+' '+'('+ r'$\AA$'+')', fontsize=24)
    plt.ylabel('Flux (erg s' + r'$^{-1}$'+' cm'+r'$^{-2}$' + r'$\AA^{-1}$'+')', fontsize=24)
    plt.suptitle(title,fontsize = 20)
    plt.xlim(wave[centerline]-40,wave[centerline]+40)
    plt.ylim(Fluxcontinuum-(1/2)*(spec[centerline]-Fluxcontinuum),Fluxcontinuum+2*(spec[centerline]-Fluxcontinuum))
    ax.tick_params(axis='both', labelsize=20)   
    #plt.show()
def Unshifted_Plotter(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    import matplotlib.pyplot as plt
    import functions

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    #Values for plotter needed-----------------------------------------------------------------------------------------------------------------------|
    EqW,EqW_rounded,vbc,vhelio,Fluxcontinuum,centerline,shift = functions.Br_Equiv_Width(plateid,MJD,fiber,emission_line)

    #Plot averaged spectrum with EqW-----------------------------------------------------------------------------------------------------------------|
    title = str(plateid)+'-'+str(MJD)+'-'+str(fiber)+'-'+str(emission_line)
    fig,ax = plt.subplots(figsize=(16,8))
    plt.plot(wave,spec,linewidth=2.5,label='Not Shifted')
    
    plt.axhline(y=Fluxcontinuum,ls='dashed',color='black')
    plt.axvline(x=wave[centerline]+shift,ls='dashed',color='r',label='Rest Emission')
    
    plt.legend(loc=1,prop={'size':18})
    plt.xlabel('Wavelength'+' '+'('+ r'$\AA$'+')', fontsize=24)
    plt.ylabel('Flux (erg s' + r'$^{-1}$'+' cm'+r'$^{-2}$' + r'$\AA^{-1}$'+')', fontsize=24)
    plt.suptitle(title,fontsize = 20)
    plt.xlim(wave[centerline]-40,wave[centerline]+40)
    plt.ylim(Fluxcontinuum-(1/2)*(spec[centerline]-Fluxcontinuum),Fluxcontinuum+2*(spec[centerline]-Fluxcontinuum))
    ax.tick_params(axis='both', labelsize=20)   
    #plt.show()
def Skyline_Plotter(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import functions

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|
    
    spec = apread.apVisit(plateid,MJD,fiber,ext=5,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    #Values for plotter needed-----------------------------------------------------------------------------------------------------------------------|
    EqW,EqW_rounded,vbc,vhelio,Fluxcontinuum,centerline,shift = functions.Br_Equiv_Width(plateid,MJD,fiber,emission_line)

    #Plot averaged spectrum with EqW-----------------------------------------------------------------------------------------------------------------|
    title = str(plateid)+'-'+str(MJD)+'-'+str(fiber)+'-'+str(emission_line)
    fig,ax = plt.subplots(figsize=(16,8))
    plt.plot(wave,spec,linewidth=2.5,label='Shifted')
    
    plt.axhline(y=Fluxcontinuum,ls='dashed',color='black')
    plt.axvline(x=wave[centerline]+shift,ls='dashed',color='r',label='Rest Emission')
    
    
    plt.xlabel('Wavelength'+' '+'('+ r'$\AA$'+')', fontsize=24)
    plt.ylabel('Flux (erg s' + r'$^{-1}$'+' cm'+r'$^{-2}$' + r'$\AA^{-1}$'+')', fontsize=24)
    plt.suptitle(title,fontsize = 20)
    plt.xlim(wave[centerline]-40,wave[centerline]+40)
    plt.ylim(0,3*Fluxcontinuum)
    ax.tick_params(axis='both', labelsize=20)   
    #plt.show()

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
    

def File_Path(plateid,mjd,fiber):

    import os

    #Creating file path------------------------------------------------------------------------------------------------------------------------------|
    '''
    Notes:
        - Need to rework this to include Location ID from master list
    '''
    server = '/Volumes/CoveyData-1/APOGEE_Spectra/python_DR13/dr13/apogee/spectro/redux/r6/apo25m/'
    #server = '/Users/ballanr/Desktop/SummerResearch/dr13/apogee/spectro/redux/r6/apo25m/'
    plate = str(plateid)
    MJD = str(mjd)
    fiber_num = str(fiber)
    therest = 'apVisit-r6-'
    dashes = '-'
    f = '.fits'

    endname = therest + plate + dashes + MJD + dashes + fiber_num + f
    filename = os.path.join(server,plate,MJD,endname)

    return filename
    
def Balmer_Decrement_Plot():

    import matplotlib.pyplot as plt
    import numpy as np

    #Importing Decrement File------------------------------------------------------------------------------------------------------------------------|
    data = np.loadtxt('/Users/ballanr/Desktop/Fwd__Bracket_Decrement/x74.txt',delimiter = '\t',unpack=True)
    print(data[0])
    T = [3750,5000,7500,8750,10000,12500,15000]
    for i in range(len(T)):
        plt.plot(np.linspace(8,12.4,num=23),data[i],label = str(T[i]) + ' K')
        plt.scatter(np.linspace(8,12.4,num=23),data[i])
    plt.xticks(np.arange(8,12.6,0.2))
    for k in range(5):
        x = 8 + k
        plt.axvline(x,ls = 'dashed',linewidth = 0.5,color = 'black')
        g = 0.1 + (0.1*k)
        plt.axhline(g,ls = 'dashed', linewidth = 0.5, color = 'black')
    plt.xlabel('Log $n_e$')
    plt.ylabel('Transition Probabilities')
    plt.legend()
    plt.show()
    
def Probs_Calc(n,T):

    import numpy as np

    #Constants------------------------------------------------------------------------------------------------------------------------|
    k = 8.6173303*(10**(-5)) #eV/K
   
    #Calculations---------------------------------------------------------------------------------------------------------------------|
    n_i = 4**2
    n_k = n**2
    E_1 = -13.6/n_i
    E_2 = -13.6/n_k
    exponent = -(E_2 - E_1)/(k*T)
    boltz = 2*(n_k/n_i)*np.exp(exponent)

    return boltz

def Probs_Plots():
    import numpy as np
    import functions
    import matplotlib.pyplot as plt
    
    T = [3750,5000,7500,8750,10000,12500,15000]
    for j in range(len(T)):
        xx = []
        yy = []
        for i in range(20):
            x = 1 + i
            y = functions.Probs_Calc(x,T[j])
            y = np.log(y)
            #plt.scatter(x,y)
            xx.append(x)
            yy.append(y)
        plt.plot(xx,yy,label = str(T[j]) + ' K')
        plt.scatter(xx,yy)
        plt.xlabel('N')
        plt.ylabel('Probability')
    plt.legend(bbox_to_anchor=(1.25,1))
    #plt.xlim(2,20)
    plt.xticks(np.arange(2,21,1))
    #plt.ylim(-2.75,1)
    plt.show()

def Saha(n,T,density):

    '''
    This isn't really useful right now...this looks at the ionization states,
    not the transition probabilities...
    '''

    import numpy as np

    k_J = 1.380648*(10**(-23)) #Joules
    k_e = 8.617*(10**(-5)) # eVs
    
    n_i = n**2
    n_k = (n+1)**2
    n_e = density
    h = 6.626*(10**(-34))
    m_e = 9.109*(10**(-31))
    E_1 = -13.6/n_i
    E_2 = -13.6/n_k
    E = E_2 - E_1

    part1 = (1/n_e)
    #print(part1)
    part2 = (2*(n_k/n_i))
    part2 = 1
    #print(part2)
    part3 = (((2*np.pi*m_e*k_J*T)/(h**2))**(3/2))
    #print(part3)
    part4 = (np.exp(E_1/(k_e*T)))
    #print(part4)
    
    #A = np.log(part1*part2*part3*part4)
    A = part1*part2*part3*part4
    
    return A
    
    
def Boltzmann(n_upper,n_lower,T):
    
    import numpy as np

    k_e = 8.617*(10**(-5)) # eVs
    

    g_b = 2*(n_upper**2)
    g_a = 2*((n_lower)**2)
    del_E = -3.4+13.6
    E_b = -13.6/(n_upper**2)
    E_a = -13.6/(n_lower**2)

    B = (g_b / g_a)*np.exp(-(E_b - E_a)/(k_e * T))

    return B

def Saha_Boltzmann(n_upper,n_lower,ion,T,density):
    
    import functions
    import numpy as np

    x = functions.Saha(ion,T,density)
    y = functions.Boltzmann(n_upper,n_lower,T)
    
    SB = (y/(1+y))*(1/(1+x))

    return np.log(SB)

def SB_Plotter():

    import functions
    import numpy as np
    import matplotlib.pyplot as plt

    T = [3750,5000,7500,8750,10000,12500,15000]
    gg = [11,12,13,14,15,16,17,18,19,20]
    densities = []
    for i in range(23):
        y = 8 + (0.2*i)
        x = (100**3)*(np.exp(y))
        densities.append(x)
    #Plotter 1
    '''
    yy = []
    for i in range(len(T)):
        y = np.log(functions.Saha_Boltzmann(11,4,1,T[i]))
        yy.append(y)
    plt.plot(T,yy)
    plt.scatter(T,yy)
    '''
    #Plotter 2

    for i in range(len(T)):
        
        for j in range(len(densities)):
            yy = []
            for k in range(len(gg)):
                y = functions.Saha_Boltzmann(gg[k],4,1,T[i],densities[j])
                yy.append(y)
            plt.plot(gg,yy,label='Temp ' + str(T[i])+'density '+str(densities[j]))
            plt.scatter(gg,yy)
    plt.legend(bbox_to_anchor=(1,1))
    plt.show()

def SB_CSV(savefile):

    import csv
    import functions
    import numpy as np
    import matplotlib.pyplot as plt
    
    with open(savefile,'w') as savefile:

        densities = []
        for i in range(23):
            y = 8 + (0.2*i)
            x = (100**3)*(np.exp(y))
            densities.append(x)
        
        T = [3750,5000,7500,8750,10000,12500,15000]
        brackett = [11,12,13,14,15,16,17,18,19,20]

        names = ['Densities','3750','5000','7500','8750','10000','12500','15000']
        writer = csv.DictWriter(savefile,delimiter = '\t',fieldnames = names)
        writer.writeheader()
        b=[]
        for i in range(len(densities)):
            a = []
            for k in range(len(T)):
                y = functions.Saha_Boltzmann(11,4,1,T[k],densities[i])
                a.append(y)
            #writer.writerow({'Densities':densities[i],'3750':a[0],'5000':a[1],'7500':a[2],'8750':a[3],'10000':a[4],'12500':a[5],'15000':a[6]})
            b.append(a)
        for i in range(len(b)):
            for k in range(7):
                b[i][k] = b[i][k]/b[12][3]
        for i in range(len(densities)):
            writer.writerow({'Densities':densities[i],'3750':b[i][0],'5000':b[i][1],'7500':b[i][2],'8750':b[i][3],'10000':b[i][4],'12500':b[i][5],'15000':b[i][6]})

        for i in range(23):
            plt.plot(T,b[i])
            plt.scatter(T,b[i])
        plt.legend(bbox_to_anchor=(1,1))
        plt.show()

def hundred_plotter(filename):

    import csv
    import functions
    import matplotlib.pyplot as plt

    with open(filename) as csvfile:

        reader = csv.reader(csvfile,delimiter = '\t')
        i = 0
        x = []
        for row in reader:
            if i < 102:
                x.append(row)
                i += 1
            else:
                break
    #return x
    print(len(x))
    for k in range(len(x)):
        if k != 0:
            line = 11
            plate = int(x[k][2])
            MJD = int(x[k][3])
            if len(x[k][4]) == 3:
                fiber = x[k][4]
            if len(x[k][4]) == 2:
                fiber = '0'+x[k][4]
            if len(x[k][4]) == 1:
                fiber = '00'+x[k][4]
            y = functions.Br_Equiv_Width_Plotter(plate,MJD,fiber,line)
            print(fiber)
            title = x[k][2]+'-'+x[k][3]+'-'+fiber
            plt.savefig('/Users/ballanr/Desktop/Test/'+title+'-'+str(line)+'.png')
            plt.close()

def twomass_uniqueness(filename,savefile):
    import csv
    import functions
    import math

    uniques = [[str(4614),'2M05412754-0646074']] #Setting up array to check for unique 2M IDs
    nonuniques = [] #The array for storing non-unqiue 2M IDs
    
    g = 1
    gg = 579468

    with open(filename) as csvfile:

        reader = csv.reader(csvfile,delimiter='\t') 
        next(reader, None)  # skip the headers
        for row in reader:
            locid = row[0]
            twomassid = row[1]
            x = functions.Uniqueness_Filter(uniques, locid, twomassid) #Checks the location ID and 2M ID against the growing list of uniques
            
            if x == 1:
                uniques.append((row[0],row[1]))
                #print('Unique visit %s' %g)
                g += 1

            elif x == 0:
                uniques.append((row[0],row[1]))
                nonuniques.append((row[0],row[1]))
                #print('Non-unique visit %s' %g)
                g += 1
            
            percent = (g/gg)*100
            numbertimes = math.floor(percent)

            hashes = '#' * int(numbertimes)
            
            spaces = ' '* (100 - int(numbertimes))

            print('Sorting through files: %.5f done ['%percent + hashes + str(spaces) + ']',end='\r',flush = True)

    with open(savefile,'w') as savefile:
        
        writer = csv.writer(savefile,delimiter = '\t')
        
        writer.writerow(('Location ID','2Mass ID'))
        
        for i in range(len(nonuniques)):
            if i != 0:
                percent = (i/len(nonuniques))*100
                numbertimes = math.floor(percent)
                hashes = '#' * int(numbertimes)
                spaces = ' '* (100 - int(numbertimes))

                writer.writerow((nonuniques[i][1],nonuniques[i][0]))

                print('Writing rows: %.5f done ['%percent + hashes + str(spaces) + ']',end='\r',flush = True)
    
    return uniques,nonuniques

def Uniqueness_Filter(array, locid, twomassid):
    
    state = 0
    
    for a,b in array:

        if a != locid and b != twomassid:
            state = 1
        
        elif a != locid and b == twomassid:
            state = 0
    
    return state

def Equivalent_Width_Error(plateid,MJD,fiber):


    '''
    Notes:

        - Work on this later
    
    '''

    import apogee.tools.read as apread

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)

    filename = functions.File_Path(plateid,MJD,fiber)
    main_header = fits.open(filename)

    vbc = main_header[0].header['BC']
    vhelio = main_header[0].header['VHELIO']
    observed_wavelength,shift,rest_wavelength = Barycentric_Correction(emission_line,vbc)

    '''
    W_lambda = equivalent width, milli-Angstroms
    lambda1, lambda2 = limits of integration for line profile, Angstroms
    l = line signal, flux units
    c = continuum signal, flux units
    del_lambda = linear dispersion/plate scale, milli-Angstroms per data point
    r = residual intensity, l/c
    n_l = number of points between lambda1 and lambda2
    n_c = number of selected continuum points
    SNR = signal-to-noise ratio
    '''
    centerline = find_nearest(wave,observed_wavelength) #Finds the closest element of wave for our observed peak

    L1 = centerline - 240 # ~ 56 Angstroms
    L2 = centerline - 150 # ~ 35 Angstroms
    R1 = centerline + 150
    R2 = centerline + 240

    lambda1 = wave[L2]
    lambda2 = wave[R1]
    n_l = R1 - L2
    del_lambda = (lambda2-lambda1) / n_l

    l = spec[L2:R1]
    leftc = np.sum(spec[L1:L2])
    rightc = np.sum(spec[R1:R2])
    c = (leftc + rightc) / ((L2 - L1) + (R2 - R1))
    
    summed = 0

    for element in l:

        summed += l[element] / c

    W_lambda = (n_l * del_lambda) - del_lambda * sum(l/c)
    sig_squared = n_l * ((del_lambda ** 2)/(SNR ** 2)) * (r/n_c) * (r + n_c)


    '''
    with open(filename) as csvfile:

        reader = csv.reader(csvfile,delimiter='\t')
        
        for row in reader:

            for a,b,c,d in row:

                if a == '''

def Skylines_Handler_Apogee(plateid,MJD,fiber,filename):
    import functions
    import apogee.tools.read as apread
    import numpy as np
    from astropy.io import fits
    '''
    Notes:
    '''

    OH_Lines = [16840.481,16802.368,16414.737,16388.492,16128.608,16079.753,15897.211,
                15890.033,15869.307,15862.489,15702.539,15597.631,15570.159,15546.141,
                15540.945,15539.711,15474.212,15462.125,15432.156,15430.163,15332.402,
                15287.789,15240.954,15187.140]

    averages = []
    windows = []

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)


    #Linear Interpolation-----------------------------------------------------------------------------------------------|

    for i in range(len(OH_Lines)):

        x1 = functions.find_nearest(wave,OH_Lines[i]-1.5)
        x2 = functions.find_nearest(wave,OH_Lines[i]+1.5)
        
        windows.append((x1,x2))
    
    
    f = fits.open(filename,mode='update')

    fspec = f[1]
    
    for i in range(len(windows)):
        
        if i != 14 and i != 18 and i != 15 and i != 19:

            if windows[i][0] > 8192:
                x = int(windows[i][0]) - 8192
                y = int(windows[i][1]) - 8192
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i][0])
                xleft = wave[windows[i][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[0][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            elif windows[i][0] < 8192 and windows[i][0] > 4096:
                x = int(windows[i][0]) - 4096
                y = int(windows[i][1]) - 4096
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i][0])
                xleft = wave[windows[i][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[1][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            else:
                x = int(windows[i][0])
                y = int(windows[i][1])
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i][0])
                xleft = wave[windows[i][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[2][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            print(slope)
            f.flush()
        
        elif i == 14:
            
            #set left bound to be 15s left bound and the right bound to 16s then proceed with calculation
            if windows[i][0] > 8192:
                x = int(windows[i+1][0]) - 8192
                y = int(windows[i][1]) - 8192
                xx = 4096 - x
                yy = 4096 - y
    
                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft) 
    
                for k in range(xx-yy):

                    fspec.data[0][xx-k] = slope*(wave[x3+k]-xleft) + yleft 

            elif windows[i][0] < 8192 and windows[i][0] > 4096:
                x = int(windows[i+1][0]) - 4096
                y = int(windows[i][1]) - 4096
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[1][xx-k] = slope*(wave[x3+k]-xleft) + yleft
            
            else:
                x = int(windows[i+1][0])
                y = int(windows[i][1])
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[2][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            print(slope)
            f.flush()
        elif i == 18:
            
            #set left bound to be 15s left bound and the right bound to 16s then proceed with calculation
            if windows[i][0] > 8192:
                x = int(windows[i+1][0]) - 8192
                y = int(windows[i][1]) - 8192
                xx = 4096 - x
                yy = 4096 - y
    
                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft) 
    
                for k in range(xx-yy):

                    fspec.data[0][xx-k] = slope*(wave[x3+k]-xleft) + yleft 

            elif windows[i][0] < 8192 and windows[i][0] > 4096:
                x = int(windows[i+1][0]) - 4096
                y = int(windows[i][1]) - 4096
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[1][xx-k] = slope*(wave[x3+k]-xleft) + yleft
            
            else:
                x = int(windows[i+1][0])
                y = int(windows[i][1])
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[2][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            print(slope)
            f.flush()
    f.close()

def Skylines_Handler(plateid,MJD,fiber,filename):
    import functions
    import apogee.tools.read as apread
    import numpy as np
    from astropy.io import fits
    '''
    Notes: 
        - Need to change this to work with new output files (the ascending order versions)
        - indexes need to be reversed since the arrays are left to right 
        - can delete this 
    '''

    OH_Lines = [16840.481,16802.368,16414.737,16388.492,16128.608,16079.753,15897.211,
                15890.033,15869.307,15862.489,15702.539,15597.631,15570.159,15546.141,
                15540.945,15539.711,15474.212,15462.125,15432.156,15430.163,15332.402,
                15287.789,15240.954,15187.140]

    averages = []
    windows = []

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)


    #Linear Interpolation-----------------------------------------------------------------------------------------------|

    for i in range(len(OH_Lines)):

        x1 = functions.find_nearest(wave,OH_Lines[i]-1.5)
        x2 = functions.find_nearest(wave,OH_Lines[i]+1.5)
        
        windows.append((x1,x2))
    
    
    f = fits.open(filename,mode='update')

    fspec = f[1]
    
    for i in range(len(windows)):
        
        if i != 14 and i != 18 and i != 15 and i != 19:

            if windows[i][0] > 8192:
                x = int(windows[i][0]) - 8192
                y = int(windows[i][1]) - 8192
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i][0])
                xleft = wave[windows[i][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[0][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            elif windows[i][0] < 8192 and windows[i][0] > 4096:
                x = int(windows[i][0]) - 4096
                y = int(windows[i][1]) - 4096
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i][0])
                xleft = wave[windows[i][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[1][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            else:
                x = int(windows[i][0])
                y = int(windows[i][1])
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i][0])
                xleft = wave[windows[i][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[2][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            print(slope)
            f.flush()
        
        elif i == 14:
            
            #set left bound to be 15s left bound and the right bound to 16s then proceed with calculation
            if windows[i][0] > 8192:
                x = int(windows[i+1][0]) - 8192
                y = int(windows[i][1]) - 8192
                xx = 4096 - x
                yy = 4096 - y
    
                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft) 
    
                for k in range(xx-yy):

                    fspec.data[0][xx-k] = slope*(wave[x3+k]-xleft) + yleft 

            elif windows[i][0] < 8192 and windows[i][0] > 4096:
                x = int(windows[i+1][0]) - 4096
                y = int(windows[i][1]) - 4096
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[1][xx-k] = slope*(wave[x3+k]-xleft) + yleft
            
            else:
                x = int(windows[i+1][0])
                y = int(windows[i][1])
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[2][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            print(slope)
            f.flush()
        elif i == 18:
            
            #set left bound to be 15s left bound and the right bound to 16s then proceed with calculation
            if windows[i][0] > 8192:
                x = int(windows[i+1][0]) - 8192
                y = int(windows[i][1]) - 8192
                xx = 4096 - x
                yy = 4096 - y
    
                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft) 
    
                for k in range(xx-yy):

                    fspec.data[0][xx-k] = slope*(wave[x3+k]-xleft) + yleft 

            elif windows[i][0] < 8192 and windows[i][0] > 4096:
                x = int(windows[i+1][0]) - 4096
                y = int(windows[i][1]) - 4096
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[1][xx-k] = slope*(wave[x3+k]-xleft) + yleft
            
            else:
                x = int(windows[i+1][0])
                y = int(windows[i][1])
                xx = 4096 - x
                yy = 4096 - y

                x3 = int(windows[i+1][0])
                xleft = wave[windows[i+1][0]]
                xright = wave[windows[i][1]]
                yleft = spec[windows[i+1][0]]
                yright = spec[windows[i][1]]
                slope = (yright - yleft) / (xright - xleft)


                for k in range(xx-yy):

                    fspec.data[2][xx-k] = slope*(wave[x3+k]-xleft) + yleft

            print(slope)
            f.flush()
    f.close()



def apVisit_Updated_Catalog(infile):

    import functions
    import pandas as pd
    from astropy.io import fits as fits
    import numpy as np
    import itertools

    x = pd.read_csv(infile,delimiter = '\t')
    problems = []
    cols = ['Location ID','2Mass ID', 'Plate ID','MJD','Fiber','Confidence1','Confidence2','Confidence 3',
            'Fractional','Br 11 EqW','Br 12 EqW','Br 13 EqW','Br 14 EqW','Br 15 EqW','Br 16 EqW','Br 17 EqW',
            'Br 18 EqW','Br 19 EqW','Br 20 EqW']
    
    df = pd.DataFrame(columns = cols)
    g = 0
    rowstart = 400000
    rowend =  rowstart + 10000
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
            
            for i in range(10):
                
                equiv_width,fcontinuum,shift,rest_wavelength = functions.Br_EqW(wave,newflux,lines[i],vbc)
                EqW.append(equiv_width)
                
                if i == 0:
                    equiv_check = equiv_width
                
            wave = np.asarray(wave) * lamshift #check which direction shift is going
            
            if equiv_check > 0:
                
                confidence1,confidence2,Mean = functions.Confidence_Level(wave,flux,snr)
                confidence3 = confidence1/confidence2

                data = [int(loc),twomass,int(plateid),int(mjd),fiber,confidence1,confidence2,confidence3,Mean,
                        EqW[0],EqW[1],EqW[2],EqW[3],EqW[4],EqW[5],EqW[6],EqW[7],EqW[8],EqW[9]]

                df.loc[len(df)+1] = data

                df1 = pd.DataFrame(wave,columns=['Wavelength'])
                df1['Flux'] = newflux
                filename = str(plateid) + '-' + str(mjd) + '-' + str(fiber) + '.csv'
                df1.to_csv('/Users/ballanr/Desktop/File Outputs/Wave and Flux/'+filename,index=False)

        except KeyError:
            print('Row '+str(g)+' has no BC value...')
            problems.append((loc,twomass,'KeyError'))
            openfile.close()
        
        except FileNotFoundError:
            print('Row '+str(g)+' doesn\'t exist...')
            problems.append((loc,twomass,'FileNotFound'))

    probs = pd.DataFrame(problems,columns = ['Location ID','2Mass ID','Problem Type'])
    probs.to_csv('/Users/ballanr/Desktop/File Outputs/'+str(rowstart)+'- End Problems.csv')
    #df.to_csv('/Users/ballanr/Desktop/File Outputs/'+str(rowstart)+'- End Equivs.csv',index=False)
    df.to_csv('/Users/ballanr/Desktop/testtest.csv',index=False)

def Br_EqW(wave,spec,line,vbc):
    
    import functions
    import numpy as np
    
    wave = np.asarray(wave)
    observed_wavelength,shift,rest_wavelength = functions.Barycentric_Correction(line,vbc)
    
    centerline = functions.find_nearest(wave,observed_wavelength)
    
    L1 = centerline - 240 # ~ 56 Angstroms
    L2 = centerline - 151 # ~ 35 Angstroms
    R1 = centerline + 150
    R2 = centerline + 241

    Fluxcontinuum = (np.sum(spec[L1:L2])+np.sum(spec[R1:R2])) / (len(spec[L1:L2])+len(spec[R1:R2]))
    EqW1 = 0

    if Fluxcontinuum == 0:

        EqW1 = 0
        EqW1_rounded = 0

    if Fluxcontinuum != 0:

        for i in range(L2,R1):

            trapezoid = (0.5)*(wave[i+1] - wave[i])*(spec[i+1] + spec[i] - (2*Fluxcontinuum))
            EqW1 += trapezoid

        #EqW_rounded = round(EqW1/Fluxcontinuum,5)
        EqW = EqW1/Fluxcontinuum
    
    return EqW,Fluxcontinuum,shift,rest_wavelength

def OH_Skylines_Remover(wave,flux):
    
    import functions
    import numpy as np
    
    OH_Lines = [15187.14,15240.954,15287.789,15332.402,15430.163,15432.156,15462.125,15474.212,
                15539.711,15540.945,15546.141,15570.159,15597.631,15702.539,15862.489,15869.307,
                15890.033,15897.211,16079.753,16128.608,16388.492,16414.737,16802.368,16840.481]
    
    windows = []
    wave = np.asarray(wave)
    
    for i in range(len(OH_Lines)):
        
        lwindow = functions.find_nearest(wave,OH_Lines[i]-1.5)
        rwindow = functions.find_nearest(wave,OH_Lines[i]+1.5)
        
        windows.append((lwindow,rwindow))
    
    for i in range(len(windows)):
        
        if i != 4 and i != 5 and i != 8 and i != 9:
            
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
        
        elif i == 4:
            
            lwindowelement = int(windows[i][0])
            rwindowelement = int(windows[i+1][1])
            leftwindow = wave[windows[i][0]]
            rightwindow = wave[windows[i+1][1]]
            leftflux = flux[windows[i][0]]
            rightflux = flux[windows[i+1][1]]
        
            slope = (rightflux - leftflux) / (rightwindow - leftwindow)
        
            for k in range(rwindowelement - lwindowelement):
                fluxvalue = slope*(wave[lwindowelement+k] - leftwindow) + leftflux
                flux[lwindowelement+k] = fluxvalue
                
        elif i == 8:
            
            lwindowelement = int(windows[i][0])
            rwindowelement = int(windows[i+1][1])
            leftwindow = wave[windows[i][0]]
            rightwindow = wave[windows[i+1][1]]
            leftflux = flux[windows[i][0]]
            rightflux = flux[windows[i+1][1]]
        
            slope = (rightflux - leftflux) / (rightwindow - leftwindow)

            for k in range(rwindowelement - lwindowelement):
                
                fluxvalue = slope*(wave[lwindowelement+k] - leftwindow) + leftflux
                flux[lwindowelement+k] = fluxvalue
                
    return flux

def Brackett_Ratios_Updated(infile):

    '''
    Notes:
        - Updated for use with pandas
    '''

    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import itertools

    opened = pd.read_csv(infile)

    for index,row in itertools.islice(opened.iterrows(),0,31): #start is two behind what you want, end is one behind
        y = row[9:]
        y = y/y[0]
        plate = int(row[2])
        mjd = int(row[3])
        fiber = int(row[4])        
        print(plate,mjd,fiber)
        x = np.arange(11,21)
        plt.plot(x,y)
        plt.scatter(x,y,s=30)
        plt.grid(True,linestyle='dashed',linewidth=0.5)
        plt.savefig('/Users/ballanr/Desktop/File Outputs/Br11 Plots/Ratios/'+str(plate)+'-'+str(mjd)+'-'+str(fiber)+'.pdf',dpi=300)
        plt.clf()

        if len(str(fiber)) == 1:
            fiber = '00'+str(fiber)
        if len(str(fiber)) == 2:
            fiber = '0'+str(fiber)

        filename = str(plate)+'-'+str(mjd)+'-'+str(fiber)
        g = pd.read_csv('/Users/ballanr/Desktop/File Outputs/Wave and Flux/'+filename+'.csv',index_col=False)
        x1 = g['Wavelength']
        y1 = g['Flux']
        plt.plot(x1,y1,linewidth=1)
        plt.scatter(x1,y1,s=25,color='black')
        plt.xlim(x1[10629],x1[11333])
        plt.ylim(0.5*y1[10629],1.5*y1[10629])
        plt.savefig('/Users/ballanr/Desktop/File Outputs/Br11 Plots/Plots/'+str(plate)+'-'+str(mjd)+'-'+str(fiber)+'.pdf',dpi=300)
        plt.clf()

def Confidence_Level(wave,flux,snr):

    import numpy as np
    import functions

    center = functions.find_nearest(wave,16811.17934)
    L1 = center - 240 # ~ 56 Angstroms
    L2 = center - 151 # ~ 35 Angstroms
    R1 = center + 150
    R2 = center + 241

    #method 1 
    leftwindow = np.mean(flux[L1:L2])
    rightwindow = np.mean(flux[R1:R2])
    cmean = (0.5)*(leftwindow + rightwindow)
    linemean = np.mean(flux[L2:R1])
    lefterr = np.std(flux[L1:L2])
    righterr = np.std(flux[R1:R2])
    cerr = (0.5)*(lefterr + righterr)
    confidence1 = (linemean - cmean)/cerr
        
    #method 2
    n_l = len(wave[L2:R1])
    n_c = len(wave[L1:L2])+len(wave[R1:R2])
    l = linemean
    c = cmean
    dellam = wave[L1+1]-wave[L1]
    r = l/c
    top = n_l*(dellam)**2
    bottom = snr**2
    sig = (top/bottom)*(r/n_c)*(r+n_c)
    confidence2 = np.sqrt(sig)
        
    #method 3
    mean_ratio = linemean / cmean

    return confidence1,confidence2,mean_ratio

