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
#loc = '/Users/ballanr/Desktop/SummerResearch/T-Tauri/Modules/Visits.csv'
#loc1 = '/Users/ballanr/Desktop/peakstocompare.csv'
#g = functions.apVisit_Updated_Catalog(loc,400000)
#g = functions.Brackett_Ratios_Updated(loc1)
#g = functions.Confidence_Level('/Users/ballanr/Desktop/File Outputs/100001-200001 Equivs.csv')

#functions.DR15_Brackett_Catalog()
#functions.Directory_Walk()
#functions.Brackett_Ratios_Updated_Grid(9659,57790,107)
functions.Aitoff(15)
#functions.DR14_Brackett_Catalog()
#functions.testest()