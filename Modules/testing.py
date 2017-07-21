#Testing of modules done here

import functions
import matplotlib.pyplot as plt
import apogee.tools.read as apread
from astropy.io import fits
import pandas as pd

'''
2Mass IDs to use:

4586,'2M03434449+3143092'
4617,'2M06450343-0034140'
4593,'2M05361555+3257145'
4380,'2M18194176-1058093'
4581,'2M06525305-1000270'

'''
loc = '/Users/ballanr/Desktop/SummerResearch/T-Tauri/Modules/Visits.csv'
loc1 = '/Users/ballanr/Desktop/testtest.csv'
#g = functions.apVisit_Updated_Catalog(loc)
g = functions.Brackett_Ratios_Updated(loc1)
#g = functions.Confidence_Level('/Users/ballanr/Desktop/File Outputs/100001-200001 Equivs.csv')

'''
for i in range(6):
    if i != 6:
        rangestart = 100000*i
        rangeend = rangestart + 100000
        print(rangestart,rangeend)
    else:
        rangestart = 500000
        rangeend = None
        print(rangestart,rangeend)
'''
