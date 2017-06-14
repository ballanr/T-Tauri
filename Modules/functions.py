def find_nearest(array,value):
    from numpy import abs
    index = (abs(array-value)).argmin()
    return index
def Br_Equiv_Width(loc_id,twomass_id,emission_line):
    import numpy as np
    import apogee.tools.read as apread
    import apogee.spec.plot as splot
    

    Lambda = splot.apStarWavegrid()

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|

    header_all = apread.apStar(loc_id,twomass_id,ext=0,header=True)
    spec_noheader = apread.apStar(loc_id,twomass_id,ext=1,header=False)
    
    #Getting NVISITS and initilization---------------------------------------------------------------------------------------------------------------|
    
    nvisits = header_all[1]['NVISITS']
    Single_StD = [] 
    EqW_array = []
    
    #Doppler Shift Constants-------------------------------------------------------------------------------------------------------------------------|

    n = (float(emission_line))**2 #Beginning electron level
    c = 299792.458 #Speed of light (km/s)
    rydberg_inf =  1.0973731568539*(10**7) #Rydberg constant (m^-1)
    electron = 9.10938356*(10**-31) #Mass of electron (kg)
    nucleus = 1.672621898*(10**-27) #Mass of hydrogen nucleus (kg)
    
    rydberg_red = rydberg_inf/(1+(electron/nucleus)) #Reduced mass Rydberg constant (m^-1)
    rest_wavelength1 = rydberg_red*((1./16.)-(1./n)) #(m^-1)
    rest_wavelength = 1/rest_wavelength1 #wavelength (m)

    #EqW Calculations for NVISITS = 1----------------------------------------------------------------------------------------------------------------|
    
    if nvisits == 1:
            
        spec = spec_noheader #Since there is only one visit there is only one array element
        header = header_all
        vhelio = header[1]['VHELIO']

        observed_wavelength1 = rest_wavelength*(1-(vhelio/c)) #Finding the location of the peak of the observed spectrum
        observed_wavelength = observed_wavelength1*(10**10)
        shift = (rest_wavelength-observed_wavelength1)*(10**10) #Finding the difference between rest and observed wavelengths

        #Generic continuum lines 45 elements wide (~9 Angstroms)
        
        centerline = find_nearest(Lambda,observed_wavelength) #Finds the closest element of Lambda for our observed peak

        L1 = centerline - 135 # ~ 30 Angstroms
        L2 = centerline - 90 # ~ 21 Angstroms
        R1 = centerline + 90
        R2 = centerline + 135
            
        #Equivalent Width Calculation

        Fluxcontinuum = (np.sum(spec[L1:L2])+np.sum(spec[R1:R2])) / (len(spec[L1:L2])+len(spec[R1:R2]))
        EqW = 0
        
        if Fluxcontinuum != 0:

            for i in range(L2,centerline):

                left_area = (Lambda[i+1]-Lambda[i])*(spec[i+1]-Fluxcontinuum)-(1./2.)*(Lambda[i+1]-Lambda[i])*(spec[i+1]-spec[i])
                EqW += left_area

            for i in range(centerline,R1):

                right_area = (Lambda[i+1]-Lambda[i])*(spec[i]-Fluxcontinuum)-(1./2.)*(Lambda[i+1]-Lambda[i])*(spec[i]-spec[i+1])
                EqW += right_area

            EqW_rounded = round(EqW/Fluxcontinuum,5)
            EqW_array.append(EqW_rounded)
            
        if Fluxcontinuum == 0:

            EqW_array.append(0)
    
    #EqW Calculations for NVISITS > 1----------------------------------------------------------------------------------------------------------------|

    if nvisits != 1:

        for i in range(nvisits):

            spec = spec_noheader[2 + i]
    
            header = header_all
            visit_vhelio = 'VHELIO' + str(i+1)
            vhelio = header[1][visit_vhelio]

            observed_wavelength1 = rest_wavelength*(1-(vhelio/c)) #Finding the location of the peak of the observed spectrum
            observed_wavelength = observed_wavelength1*(10**10)
            shift = (rest_wavelength-observed_wavelength1)*(10**10) #Finding the difference between rest and observed wavelengths

            #Generic continuum lines 45 elements wide (~9 Angstroms)

            centerline = find_nearest(Lambda,observed_wavelength) #Finds the closest element of Lambda for our observed peak

            L1 = centerline - 135
            L2 = centerline - 90
            R1 = centerline + 90
            R2 = centerline + 135
            
            #Equivalent Width Calculation

            Fluxcontinuum = (np.sum(spec[L1:L2])+np.sum(spec[R1:R2])) / (len(spec[L1:L2])+len(spec[R1:R2]))
            EqW = 0
        
            if Fluxcontinuum != 0:

                for i in range(L2,centerline):

                    left_area = (Lambda[i+1]-Lambda[i])*(spec[i+1]-Fluxcontinuum)-(1./2.)*(Lambda[i+1]-Lambda[i])*(spec[i+1]-spec[i])
                    EqW += left_area

                for i in range(centerline,R1):

                    right_area = (Lambda[i+1]-Lambda[i])*(spec[i]-Fluxcontinuum)-(1./2.)*(Lambda[i+1]-Lambda[i])*(spec[i]-spec[i+1])
                    EqW += right_area

                EqW_rounded = round(EqW/Fluxcontinuum,5)
                EqW_array.append(EqW_rounded)
            
            if Fluxcontinuum == 0:

                EqW_array.append(0)

    #Error Calculations------------------------------------------------------------------------------------------------------------------------------|
    
    avg = np.mean(EqW_array)
    EqW_avg = round(avg,5)
    StD_avg = round(np.std(EqW_array),5)

    #Calculating the standard dev for each element
    for i in EqW_array:
        squared = (i - avg)**2
        final=round(np.sqrt(squared/len(EqW_array)),5)
        Single_StD.append(final) 
    
    #Returning Values--------------------------------------------------------------------------------------------------------------------------------|

    return EqW_array,Single_StD,EqW_avg,StD_avg