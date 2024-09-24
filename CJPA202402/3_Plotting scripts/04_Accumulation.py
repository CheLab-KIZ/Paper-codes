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
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

df = pd.read_csv('') # use your path
binedges = np.arange(0,df.t.max(),1)
binlabels = np.arange(0,df.t.max()-1,1)
df2=df.groupby(['bootrep',pd.cut(df['t'],bins=binedges ,labels=binlabels,include_lowest=True)]).sum()
df2=df2.drop(['t'],1)
ac= df2.sort_index(level=['bootrep','t'], ascending=[True,False])
ac=ac.groupby(['bootrep']).cumsum()
ac=ac.reset_index()
ac=ac.drop(['bootrep'],1)

acc=ac.groupby(['t']).median().reset_index()
low=ac.groupby(['t']).quantile(0.25)
up=ac.groupby(['t']).quantile(0.75)

start=50
alpha=0.1
plt.style.use('seaborn-notebook')
sns.set_style("ticks")
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure(figsize=(4,18), dpi=300)
plt.legend(loc='upper right',labels=['in situ speciation','Colonization'])

ax3=plt.subplot(711)
plt.legend(loc='upper right',labels=['in situ speciation','Colonization'])
plt.plot(acc.index,acc['H_insitu'], color='salmon')
plt.fill_between( low.index,low['H_insitu'], up['H_insitu'], color='r', alpha=alpha)
plt.plot(acc.index,acc['H_colonization'], color='limegreen')
plt.fill_between( low.index,low['H_colonization'], up['H_colonization'], color='limegreen', alpha=alpha)
plt.xlim(start,0)
plt.xticks(range(start, 0, -5)) 
plt.yscale('log')
plt.ylim(0,pow(10,3))
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position('right')

plt.savefig('', dpi=300, bbox_inches='tight')
