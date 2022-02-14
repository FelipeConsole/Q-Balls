
import matplotlib.pyplot as plt

import numpy as np
#
font1 = {'family': 'serif',
        'color':  'purple',
        'weight': 'normal',
        'size': 16,
        }
font2 = {'family': 'serif',
        'color':  'purple',
        'weight': 'normal',
        'size': 10,
        }

phi1 = np.genfromtxt('phi_0=0_4.dat',delimiter='  ')
#x = phi1[0:len(phi1):1][:,3]
x = phi1[0:10:1][:,0]
t11 = phi1[0:10:1][:,4]
#
# plt.plot(x,phi's',markersize=1)

#plt.semilogy(x,t11,'s',markersize=1)
#plt.semilogy(omega_a_0_1,Q_N_a_0_1,'s',markersize=1)
#plt.text(0.8, 0.7, r'$V(\phi^2) = \phi^6 -2 \phi^4 + 1.1\phi^2$', fontdict=font2)
#plt.title(r'Boson stars, ($\alpha = 4 \pi G =  0.005$)', fontdict=font1)
plt.ylabel(r'$ Q_{N} $',fontdict=font1)
plt.xlabel(r'$ \omega $',fontdict=font1)
#plt.savefig('BS - Q_N x omega - a = 0_005')

#plt.plot(omega_a_0001,Q_N_a_0001,'s',markersize=1)
#plt.plot(x,phi,'s',markersize=1)
#plt.xlim([0.7,0.78])
#plt.ylim([0,200])
plt.show()
#
