############### UPDATED FUNCTIONS ###############

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