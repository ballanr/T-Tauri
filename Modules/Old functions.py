'''
Gotta look through these and fix them!!!!

'''

def Brackett_Equiv_Width(loc_id,twomass_id,j,EqW_array,emission_dict):

    
    #given a certain spec
    spec_header = apread.apStar(loc_id,twomass_id,ext=0,header=True)
    spec_noheader = apread.apStar(loc_id,twomass_id,ext=1,header=False)
    nvisits = spec_header[1]['NVISITS']
    vhelio = spec_header[1]['VHELIO']
    
    
    #Calculating Equivalent Width
    n=float((11+j)**2)
    c = 299792
    rydberg = 1.0973731568539*(10**7)
    electron = 9.10938356*(10**-31)
    nucleus = 1.672621898*(10**-27)
    fracryd = rydberg/(1+(electron/nucleus))
    
    
    vacuum = fracryd*((1./16.)-(1./n))
    lambda_obs = 1/vacuum
    calculated_point1 = lambda_obs*(1+(-vhelio/c))
    calculated_point2 = calculated_point1*(10**10)
    
    #EqW Calculations for NVISITS = 1
    if nvisits == 1:
            
        spec = spec_noheader
    
        centerline = find_nearest(Lambda,calculated_point2)
        L1 = centerline - 135
        L2 = centerline - 90
        R1 = centerline + 90
        R2 = centerline + 135
            
        #generic continuum lines 35 elements wide 
        lsum= np.sum(spec[L1:L2])/ len(spec[L1:L2])
        rsum = np.sum(spec[R1:R2])/len(spec[R1:R2])
        Fc= (lsum+rsum)/2
        EqW=0
        
        if Fc!=0:
            for i in range(L2,centerline):
                left_area = (Lambda[i+1]-Lambda[i])*(spec[i+1]-Fc)-(1./2.)*(Lambda[i+1]-Lambda[i])*(spec[i+1]-spec[i])
                EqW += left_area
            for i in range(centerline,R1):
                right_area = (Lambda[i+1]-Lambda[i])*(spec[i]-Fc)-(1./2.)*(Lambda[i+1]-Lambda[i])*(spec[i]-spec[i+1])
                EqW += right_area
            EqW_rounded = round(EqW/Fc,5)
            EqW_array.append(EqW_rounded)
            
        if Fc==0:
            EqW_array.append(0)
    
    #EqW Calculations for NVISITS > 1
    if nvisits != 1:
        for i in range(nvisits):
            spec = spec_noheader[2 + i]
    
            centerline = find_nearest(Lambda,calculated_point2)
            L1 = centerline - 135
            L2 = centerline - 90
            R1 = centerline + 90
            R2 = centerline + 135
                
            #generic continuum lines 35 elements wide 
            lsum= np.sum(spec[L1:L2])/ len(spec[L1:L2])
            rsum = np.sum(spec[R1:R2])/len(spec[R1:R2])
            Fc= (lsum+rsum)/2
            EqW=0
        
            if Fc!=0:
                for i in range(L2,centerline):
                    left_area = (Lambda[i+1]-Lambda[i])*(spec[i+1]-Fc)-(1./2.)*(Lambda[i+1]-Lambda[i])*(spec[i+1]-spec[i])
                    EqW += left_area
                for i in range(centerline,R1):
                    right_area = (Lambda[i+1]-Lambda[i])*(spec[i]-Fc)-(1./2.)*(Lambda[i+1]-Lambda[i])*(spec[i]-spec[i+1])
                    EqW += right_area
                EqW_rounded = round(EqW/Fc,5)
                EqW_array.append(EqW_rounded)
            
            if Fc==0:
                EqW_array.append(0)
    
    #Error Calculations
    Single_StD = [] 
    avg = np.mean(EqW_array)

    #Calculating the standard dev for each element
    for i in EqW_array:
        squared = (i - avg)**2
        final=round(np.sqrt(squared/len(EqW_array)),5)
        Single_StD.append(final) 

    
    
    #Average standard dev
    StD_avg = round(np.std(EqW_array),5)
        
    #Average EqW
    EqW_avg = round(np.mean(EqW_array),5)
    
    emission_dict['StD_avg'+str(11+j)]= StD_avg
    emission_dict['EqW_avg'+str(11+j)]=EqW_avg
    emission_dict['EqW_array'+str(11+j)]=EqW_array
    emission_dict['Single_StD'+str(11+j)]=Single_StD
    
    
    #values={'StD_avg'+str(11+j):StD_avg,'EqW_avg'+str(11+j):EqW_avg,'EqW_array'+str(11+j):EqW_array,
    #        'Single_StD'+str(11+j):Single_StD}
    #print(j)
    #print(values)
    #return values




def Brackett_Catalog(loc_id,twomass_id):
    emission_dict={}
    
    #given a certain spec
    spec_header = apread.apStar(loc_id,twomass_id,ext=0,header=True)
    spec_noheader = apread.apStar(loc_id,twomass_id,ext=1,header=False)
    nvisits = spec_header[1]['NVISITS']
    emission_dict['nvisits']=nvisits
    EqW_array = []
    
    
    #Getting HJD for visits
    hjd=[]
    for i in range(nvisits):
        hjd_visit = str('HJD'+str(i+1))
        hjd.append(spec_header[1][hjd_visit])
    emission_dict['HJD']=hjd
    
    #Getting RA
    RA = spec_header[1]['RA']
    emission_dict['RA']=RA
    
    #Getting Dec
    dec = spec_header[1]['DEC']
    emission_dict['DEC']=dec
    
    #Getting J Mag
    J_mag = spec_header[1]['J']
    emission_dict['J']=J_mag
    
    #Getting H Mag
    H_mag = spec_header[1]['H']
    emission_dict['H']=H_mag
    
    #Getting K Mag
    K_mag = spec_header[1]['K']
    emission_dict['K']=K_mag  
    
    #Getting TEffective
    Teff_avg = spec_header[1]['RVTEFF']
    Teff_single=[]
    for i in range(nvisits):
        teff_visit = str('RVTEFF'+str(i+1))
        Teff_single.append(spec_header[1][teff_visit])
    emission_dict['Teff_avg']=Teff_avg
    emission_dict['Teff_single']=Teff_single
    
    #Getting LogG
    Logg_avg = spec_header[1]['RVLOGG']
    Logg_single=[]
    for i in range(nvisits):
        logg_visit = str('RVLOGG'+str(i+1))
        Logg_single.append(spec_header[1][logg_visit])
    emission_dict['LogG_avg']=Logg_avg
    emission_dict['LogG_single']=Logg_single
    
    #Getting SNR
    SNR_avg = spec_header[1]['SNR']
    SNR_single=[]
    for i in range(nvisits):
        SNR_visit = str('SNRVIS'+str(i+1))
        SNR_single.append(spec_header[1][SNR_visit])
    emission_dict['SNR_avg']=SNR_avg
    emission_dict['SNR_single']=SNR_single
    
    #Getting Glong
    Glon = spec_header[1]['GLON']
    emission_dict['Glon']=Glon   
    
    #Getting Glat
    Glat = spec_header[1]['GLAT']
    emission_dict['Glat']=Glat
    
    for j in range(10):
        EqW_array = []
        Brackett_Equiv_Width(loc_id,twomass_id,j,EqW_array,emission_dict)
    
    return emission_dict

def Brackett_EqW_Plot(loc_id, two_massid,number,left_limit,right_limit,bottom_limit,top_limit):
    spec_header = apread.apStar(loc_id, two_massid,ext=0,header=True)
    
    nvisits = spec_header[1]['NVISITS']
    n=number
    #calculate the emission line in a vacuum
    vhelio = spec_header[1]['VHELIO']
    doppler = spec_header[1]['VRAD1']
    bc1 = spec_header[1]['BC1']
    
    
    if nvisits == 1:
        spec1 = apread.apStar(loc_id, two_massid,ext=1,header=False)
    else:
        spec1 = apread.apStar(loc_id, two_massid,ext=1,header=False)[3]

    c = 299792
    rydberg = 1.0973731568539*(10**7)
    electron = 9.10938356*(10**-31)
    nucleus = 1.672621898*(10**-27)
    fracryd = rydberg/(1+(electron/nucleus))
    vacuum = fracryd*((1./16.)-(1./(float(n**2))))
    lambda_obs = 1/vacuum
    #calculated_point1 = lambda_obs*((1-(vhelio/c))/(1-(doppler/c)))
    #calculated_point1 = lambda_obs*(1-(bc1/c))
    #calculated_point1 = lambda_obs*(1-(vhelio/c))
    calculated_point1 = lambda_obs*(1-((vhelio)/c))
    diff2 = (lambda_obs-calculated_point1)*(10**10)
    calculated_point2 = calculated_point1*(10**10)
    
    
    
    centerline = find_nearest((Lambda+diff2),(lambda_obs*(10**10)))
    
    
        
    L1 = centerline - 135
    L2 = centerline - 90
    R1 = centerline + 90
    R2 = centerline + 135
    Lwindow = centerline - 160
    Rwindow = centerline + 160
    
    #Calculate Fc
    lsum= np.sum(spec1[L1:L2])/ len(spec1[L1:L2])
    rsum = np.sum(spec1[R1:R2])/len(spec1[R1:R2])
    Fc= (lsum+rsum)/2
    
    #Calculate Equivalent Width
    EqW=0
    for i in range(L2,R1):
        summ=(Fc*(Lambda[i+1]-Lambda[i]))-((1./2.)*(Lambda[i+1]-Lambda[i])*(spec1[i+1]+spec1[i]))
        EqW = EqW + summ
    EqW = abs(EqW/Fc)
    upper = Lambda[centerline]+(EqW/2)
    lower = Lambda[centerline]-(EqW/2)
    
    #Plot averaged spectrum with EqW
    fig,ax = plt.subplots(figsize=(16,8))
    plt.plot((Lambda+diff2),spec1,linewidth=2.5,label='Shifted')
    plt.plot(Lambda,spec1,linewidth=2.5,label='Unshifted')
    plt.axhline(y=Fc,ls='dashed',color='black')
    plt.axvline(x=Lambda[centerline]+diff2,ls='dashed',color='r',label='Rest Emission')
    plt.axvline(calculated_point2,ls=':',color='r',label='Star Emission')
    plt.legend()
    plt.xlabel('Wavelength'+' '+'('+ r'$\AA$'+')', fontsize=24)
    plt.ylabel('Flux (erg s' + r'$^{-1}$'+' cm'+r'$^{-2}$' + r'$\AA^{-1}$'+')', fontsize=24)
    plt.xlim(left_limit,right_limit)
    plt.ylim(bottom_limit,top_limit)
    
    ax.tick_params(axis='both', labelsize=20)

def Avg_Brackett_Ratios(loc_id,twomass_id):
    datax=[11,12,13,14,15,16,17,18,19,20]
    datay=[]
    datayy=[]

    with open('Average Visits2.csv') as csvfile:
        
        reader = csv.DictReader(csvfile,delimiter='\t')
    
        for row in reader:
            if int(row['Location ID'])==loc_id and row['2Mass ID']==twomass_id:
                datay.append(float(row['Br11 Avg EqW']))
                datayy.append(float(row['Br11 Avg EqW']))
                datay.append(float(row['Br12 Avg EqW']))
                datay.append(float(row['Br13 Avg EqW']))
                datay.append(float(row['Br14 Avg EqW']))
                datay.append(float(row['Br15 Avg EqW']))
                datay.append(float(row['Br16 Avg EqW']))
                datay.append(float(row['Br17 Avg EqW']))
                datay.append(float(row['Br18 Avg EqW']))
                datay.append(float(row['Br19 Avg EqW']))
                datay.append(float(row['Br20 Avg EqW']))

    #print(datay)
    #print(datayy)

    for i in range(len(datay)):
        datay[i]=datay[i]/datayy[0]
    
    #print(datay)A

    plt.plot(datax,datay)
    plt.ylabel('Br n>11 / Br 11')
    plt.xlabel('n')
    plt.title(twomass_id)
    #plt.savefig(twomass_id+' Average Ratio Test.pdf',dpi=300)

def All_Brackett_Ratios(loc_id,twomass_id):
    datax=[11,12,13,14,15,16,17,18,19,20]
    dataxx=[]
    datay=[]
    datayy=[]
    JD=[]
    with open('All Visits2.csv') as csvfile:
        
        reader = csv.DictReader(csvfile,delimiter='\t')
    
        for row in reader:
            if int(row['Location ID'])==loc_id and row['2Mass ID']==twomass_id:
                dataxx.append(row['Visit #'])
                datay.append([float(row['Br11 EqW']),
                            float(row['Br12 EqW']),float(row['Br13 EqW']),
                            float(row['Br14 EqW']),float(row['Br15 EqW']),
                            float(row['Br16 EqW']),float(row['Br17 EqW']),
                            float(row['Br18 EqW']),float(row['Br19 EqW']),
                            float(row['Br20 EqW'])])
                datayy.append(float(row['Br11 EqW']))
                JD.append(row['HJD'])

    #print(datay)
    #print(dataxx)
    #print(datayy)

    for i in range(len(dataxx)):
        #print(i)
        for j in range(len(datay[i])):
            datay[i][j]=datay[i][j]/datayy[i]
    
    #print(datay)
    
    fig,ax=plt.subplots(figsize=(16,8))
    ax.tick_params(axis='both', labelsize=20)
    
    for i in range(len(dataxx)):
        plt.plot(datax,datay[i],label=JD[i])
    #plt.axhline(0,ls='dashed',color='black')
    plt.ylabel('Br n>11 / Br 11',fontsize=24)
    plt.xlabel('n',fontsize=24)
    #plt.title(twomass_id)
    #if len(dataxx)>16:
    #    plt.legend(bbox_to_anchor=(1.05, 1),ncol=2, loc=2, borderaxespad=0.)
    #else:
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #plt.legend(title='Julian Date')
        
        
    #plt.savefig('testgraphics/' + twomass_id+' Ratios Test.pdf',dpi=300,bbox_inches='tight')
    
    #Multiplot stuff

    #length=int(len(dataxx))
    #nrows = int(length/2)
    #print(nrows)
    #if nrows % 2 == 0:
    #    fig,ax = plt.subplots(nrows,nrows)
    #
    #    ax=ax.ravel()
    #
    #    for i in range(len(dataxx)):
    #        ax[i].plot(datax,datay[i])
    
    #else:
    #    
    #    fig,ax = plt.subplots(nrows,nrows+1)
    #
    #    ax=ax.ravel()
    #
    #    for i in range(len(dataxx)):
    #        ax[i].plot(datax,datay[i])
            
    #fig,ax = plt.subplots(nrows+1,nrows+1,sharex='all',sharey='all')
    #number = int((nrows+1)**2)
    #ax=ax.ravel()
    
    #for i in range(len(dataxx)):
    #    ax[i].plot(datax,datay[i])
    
    #for i in range((number-length)):
    #    fig.delaxes(ax.flatten()[i+length])
    
    #plt.savefig('All Ratios Test.pdf',dpi=300)