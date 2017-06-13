#Importing the necessary packages
import apogee.tools.read as apread
import apogee.spec.plot as splot
import numpy as np
import matplotlib.pyplot as plt
Lambda = splot.apStarWavegrid()
import csv
import sys
import functions as f

#Importing the apStar file

'''
Steps:

- Create new directory for files
- Read IDs from csv file
- Import spectra and skylines per element from csv file
- Apply barycentric correction to visit data and skylines
- Calculate EqW
    -In future first mask out skylines 
- Create combined spectrum (mean or median?)
- Measure EqW of combined spectrum
- Save all data in new file

'''

with open(csvname) as csvfile:

    reader = csv.DictReader(csvfile,delimiter='\t')
    
    for row in reader:
        loc_id=int(row['Location ID'])
        twomass_id=row['2Mass ID']

        #functions live here