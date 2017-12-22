import pandas as pd
from astropy.io import fits
import itertools
from PyAstronomy.pyasl import helcorr as helcorr
import numpy as np
import functions
import functions2

'''for i in range(4):
    k = 17 + i
    functions2.Brackett_Ratios_Updated_Grid(k)
'''
functions2.Brackett_Ratios_Updated_Grid(17)