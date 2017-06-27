#Testing of modules done here

import functions
import matplotlib.pyplot as plt

#x = Functions.Br_Equiv_Width(6218,56168,148,11)
#print(x)

#x = Functions.Brackett_Ratios(6218,56168,148)
#x = functions.Balmer_Decrement_Plot()
#g = functions.Probs_Calc(11,3750)
#g = functions.Probs_Plots()
#g = functions.Saha_Boltzmann(11,4,1,3750,2980957987.04)
#g = functions.Saha(1)
#print(g)

#x = functions.Br_Equiv_Width_Plotter(6218,56168,148,11)
'''for i in range(10):
    x = functions.Br_Equiv_Width_Plotter(6218,56168,148,11+i)
    string = 'Br ' + str(11+i)
    plt.savefig('/Users/ballanr/Desktop/Test/' + string) '''
#g = functions.SB_Plotter()
#g = functions.SB_CSV('test3.csv')

'''x = functions.apStar_to_apVisit(4617,'2M06450343-0034140')
for i in range(x[0][3]):
    x = functions.Br_Equiv_Width(x[i][0],x[i][1],x[i][2],11)
    print(x)'''
#x = functions.apVisit_Catalog_Output(
#    '/Users/ballanr/Desktop/Research/test_visit.csv',
#    'testest.csv')

#x = functions.apVisit_Catalog_Output(
#    '/Users/ballanr/Desktop/Research/test_visits.csv'
#    ,'Visits.csv')

x,locid,twomassid = functions.apStar_to_apVisit(4586,'2M03434449+3143092')
plateid = x[1][0]
MJD = x[1][1]
fiber = x[1][2]
emission_line = 12
x = functions.Br_Equiv_Width(plateid,MJD,fiber,emission_line)
#y = functions.Skyline_Plotter(plateid,MJD,fiber,emission_line)
print(x)