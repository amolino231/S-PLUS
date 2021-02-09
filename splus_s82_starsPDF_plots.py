__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import matplotlib.pyplot as plt

data = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
data += 'S82/Dec2017/data_quality/PDFs/stars/starsPDFs.cat'
z,p = U.get_data(data,(0,1))

plt.figure(1)
plt.clf()
plt.fill_between(z[::5],p[::5]/p[::5].sum(),0.,color='blue',alpha=0.2)
plt.plot(z[::5],p[::5]/p[::5].sum(),'b-',lw=1,alpha=0.9)
plt.xlim(0.004,0.89)
plt.ylabel('Probability Density Function',size=24,labelpad=1)
plt.ylim(0.001,0.1)
plt.xlabel('$z$',size=32,labelpad=-3)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(['Stellar-PDF'],loc='upper right',fontsize=30)


