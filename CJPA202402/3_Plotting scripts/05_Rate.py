from numpy.lib.stride_tricks import as_strided as stride
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import NullFormatter
from pandas import DataFrame
from pylab import *
np = pd.np
import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter as cc

from matplotlib.backends.backend_pdf import PdfPages

df = pd.read_csv('') # use your path
#sum events every 1Ma
binedges = np.arange(0,df.t.max(),1)
binlabels = np.arange(0,df.t.max()-1,1)
df2=df.groupby(['bootrep',pd.cut(df['t'],bins=binedges ,labels=binlabels,include_lowest=True)]).sum()

df2=df2.drop(['t'],1)
list(df2.columns.values)
ac=df2.groupby(['bootrep','t']).sum()

#accumulation through time
ac= ac.sort_index(level=['bootrep','t'], ascending=[True,False])
accum=ac.groupby(['bootrep']).cumsum()

acrolling=ac.groupby(['bootrep']).apply(lambda x: x.rolling(window=2).sum()/2)
for i in 'BEFDHCG':
    acrolling[f'{i}']=accum[f'{i}']

acrolling[['B_shift_1','E_shift_1','F_shift_1','D_shift_1','H_shift_1','C_shift_1','G_shift_1']] = acrolling.groupby('bootrep')['B','E','F','D','H','C','G'].shift(1)
B_rate= acrolling[['B_insitu','B_colonization','B_colonization_from_H','B_colonization_from_E','B_colonization_from_C','B_colonization_from_D','B_colonization_from_F','B_colonization_from_G']].div(acrolling.B_shift_1, axis=0)
H_rate= acrolling[['H_insitu','H_colonization','H_colonization_from_B','H_colonization_from_E','H_colonization_from_C','H_colonization_from_D','H_colonization_from_F','H_colonization_from_G']].div(acrolling.H_shift_1, axis=0)
E_rate= acrolling[['E_insitu','E_colonization','E_colonization_from_H','E_colonization_from_B','E_colonization_from_C','E_colonization_from_D','E_colonization_from_F','E_colonization_from_G']].div(acrolling.E_shift_1, axis=0)
C_rate= acrolling[['C_insitu','C_colonization','C_colonization_from_H','C_colonization_from_E','C_colonization_from_B','C_colonization_from_D','C_colonization_from_F','C_colonization_from_G']].div(acrolling.C_shift_1, axis=0)
D_rate= acrolling[['D_insitu','D_colonization','D_colonization_from_H','D_colonization_from_E','D_colonization_from_C','D_colonization_from_B','D_colonization_from_F','D_colonization_from_G']].div(acrolling.D_shift_1, axis=0)
F_rate= acrolling[['F_insitu','F_colonization','F_colonization_from_H','F_colonization_from_E','F_colonization_from_C','F_colonization_from_D','F_colonization_from_B','F_colonization_from_G']].div(acrolling.F_shift_1, axis=0)
G_rate= acrolling[['G_insitu','G_colonization','G_colonization_from_H','G_colonization_from_E','G_colonization_from_C','G_colonization_from_D','G_colonization_from_F','G_colonization_from_B']].div(acrolling.G_shift_1, axis=0)

rate=pd.concat([B_rate, E_rate,D_rate,F_rate,H_rate,C_rate,G_rate], axis=1)
acc=rate.groupby(['t']).median().reset_index()

low=rate.groupby(['t']).quantile(0.25)
up=rate.groupby(['t']).quantile(0.75)

start=50
yrate=0.2
alpha=0.1
plt.style.use('seaborn-notebook')
sns.set_style("ticks")
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure(figsize=(4,18), dpi=300)

plt.legend(loc='upper right',labels=['in situ speciation','colonization'])

ax1=plt.subplot(711)  # (nRows, nColumns, axes number to plot)
plt.plot(acc.index,acc['H_insitu'], color='salmon')
plt.fill_between(low.index,low['H_insitu'], up['H_insitu'], color='r', alpha=alpha)
plt.plot(acc.index,acc['H_colonization'], color='limegreen')
plt.fill_between(low.index,low['H_colonization'], up['H_colonization'], color='limegreen', alpha=alpha)
plt.xlim(start, 0); plt.ylim(0,yrate)
plt.xticks(range(start, 0, -5)) 

plt.savefig('', dpi=300, bbox_inches='tight')

