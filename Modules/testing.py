#Testing of modules done here

import functions
import matplotlib.pyplot as plt

x = functions.twomass_uniqueness(
    '/Users/ballanr/Desktop/Research/test.csv','Non Uniques.csv'
    )
print(len(x[0]),len(x[1]))