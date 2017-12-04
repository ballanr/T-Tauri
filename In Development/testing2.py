import functions

import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

#excel = '/Users/ballanr/Desktop/File Outputs/Currently Working On/18/18s.csv'
#excel = '/Users/ballanr/Desktop/File Outputs/Currently Working On/19/19s.csv'
excel = '/Users/ballanr/Desktop/File Outputs/Currently Working On/20/20s.csv'
openfile = pd.read_csv(excel)

sns.set(style='whitegrid',font_scale = 2)


data1 = pd.DataFrame(columns = ['Density','Temp'])
data1['Density'] = openfile['Density']
data1['Temp'] = openfile['Temp']

p = sns.JointGrid(y='Temp',x='Density', data=data1, size=10)

p = p.plot_joint(plt.scatter,s=5,color='black')
#plt.xlim(-5,5)
p = p.plot_joint(sns.kdeplot,shade=True,shade_lowest=False,cmap='RdGy',zorder=0,n_levels=15)

p = p.plot_marginals(sns.distplot,kde=True)

plt.suptitle('W/ Brackett 20')
plt.savefig('/Users/ballanr/Desktop/File Outputs/Currently Working On/20/KDE.pdf',bbox_inches='tight',dpi=300)
plt.show()