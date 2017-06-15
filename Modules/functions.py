def find_nearest(array,value):
    from numpy import abs
    index = (abs(array-value)).argmin()
    return index
def Br_Equiv_Width(plateid,MJD,fiber,emission_line):
    import numpy as np
    import apogee.tools.read as apread

    #Importing spectrum via apogee-------------------------------------------------------------------------------------------------------------------|

    spec = apread.apVisit(plateid,MJD,fiber,ext=1,header=False)
    wave = apread.apVisit(plateid,MJD,fiber,ext=4,header=False)
    header = apread.apStar(4586,'2M03434449+3143092',ext=0,header=True)
    
    #Barycentric Correction--------------------------------------------------------------------------------------------------------------------------|
    
    vbcstring = 'BC' + str(1) #Somehow need to figure out which visit this MJD applies to
    vbc = header[1][vbcstring]
    observed_wavelength,shift = Barycentric_Correction(emission_line,vbc)
    
    #Equivalent Width Calculation--------------------------------------------------------------------------------------------------------------------|

    #Finding the centerline and checking that it matches the peak in the window
    centerline = find_nearest(wave,observed_wavelength) #Finds the closest element of wave for our observed peak
    centerline_check = Max_Flux_Check(spec,centerline)

    if spec[centerline] != spec[centerline_check[1]]:
        print('The centerline has changed from ' + str(centerline) + ' to ' + str(centerline_check[1]) + ' with a new flux of ' + str(centerline_check[0]) + 
        ' from ' + str(spec[centerline]) + '.')
        centerline = centerline_check[1]
    

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
    
    return EqW,EqW_rounded
 

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

    return observed_wavelength,shift

def Max_Flux_Check(array,centerline):
    
    array = array[centerline-165:centerline+165]
    y = max(array)
    z = array.tolist().index(y)
    z = z + centerline - 165
    return y,z


def Brackett_Ratios():
    #blah