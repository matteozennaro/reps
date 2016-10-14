import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cn
from math import pi

plt.rc('axes', linewidth=2)

plt.rc('xtick', labelsize=24)
plt.rc('xtick.major', size=10, width=1.5)
plt.rc('xtick.minor', size=5, width=1)

plt.rc('ytick', labelsize=24)
plt.rc('ytick.major', size=10, width=1.5)
plt.rc('ytick.minor', size=5, width=1)

font = {'family' : 'serif',
	'serif'  : 'Computer Modern Roman',
	'weight' : 'normal',
	'size'   : 24}

plt.rc('font', **font)

plt.rc('text', usetex=True)

def myticks(x,p):
	if float(x) < 0.01:
		return '$10^{\mbox{-} %.0f}$'%(abs(np.log10(x)))
	elif float(x) > 100.0:
		return '$10^{%.0f}$'%(np.log10(x))
	elif float(x) < 0.1:
		return '$%.2f^{ \, }$'%x
	elif float(x) < 1.0:
		return '$%.1f^{ \, }$'%x
	else:
		return '$%.0f^{ \, }$'%x

file_p_c = 'PK_TABS/power_z18_tk.dat'
file_p_m = 'mnu015_case00_Pc_rescaled_z99.0000.txt'

h = 0.6711
Ob = 0.049
mnu = 0.15
On = mnu/(93.14*h**2)
Oc = 0.3175 - Ob - On
As = 2.13e-9
ns = 0.9624
kpivot = 0.05

k_c,tc,tb = np.loadtxt(file_p_c,usecols=(0,1,2),unpack=True)
k_m,pc_m = np.loadtxt(file_p_m,unpack=True)

pk0 = As * (((k_c*h)/kpivot)**(ns-1.0)) * (2.0*pi**2*k_c*h**4)

t = (Ob/(Ob+Oc))*tb + (Oc/(Ob+Oc))*tc
t = tc

pc_c = pk0*t**2

fig = plt.figure(figsize=(10,10))
ax = plt.axes()

#pc_m_interp = np.interp(k_c,k_m,pc_m)
pc_c_interp = np.interp(k_m,k_c,pc_c)

ax.axhspan(0.0001,0.001,color='0.85')
ax.axhspan(0.0001,0.0002,color='0.7')
# ax.plot(k_c,0.0*k_c,color='k',linestyle=':')
# ax.plot(k_c,abs(pc_m_interp/pc_c-1.0),'r-',linewidth=2)
ax.plot(k_m,abs(pc_c_interp/pc_m-1.0),'r-',linewidth=2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(myticks))
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(myticks))
ax.tick_params(axis='x', which='major', pad=8.)
ax.set_xlabel('$k \, [h^{-1} \, \mathrm{Mpc}]$')
ax.set_ylabel('$|P_c(\mathrm{matteo})/P_c(\mathrm{camb})-1|$')
ax.set_ylim([0.0001,0.01])

#plt.figtext(0.2,0.55,'0.2%')
plt.figtext(0.5,0.915,'$M_\\nu=0.15 \, \mathrm{eV}$ - $w_0 = -1$, $w_a = 0$',linespacing=2,horizontalalignment ='center')
ax.text(1.e-4,1.e-3,'0.1\%')
ax.text(1.e-4,2.e-4,'0.02\%')

# plt.savefig("mnu015_w0-0.9_wa-0.3_Pcb_IC_vs_CAMB.pdf",dpi=2000)
plt.show()

'''
h_file_c = 'mnu000_w-1.00/class_tabs/power_background.dat'
h_file_m = 'mnu000_w-1.00/mnu000_case05_hubble.txt'

z_c,h_c = np.loadtxt(h_file_c,usecols=(0,3),unpack=True)
z_m,h_m = np.loadtxt(h_file_m,usecols=(0,1),unpack=True)

h_c *= (cn.c/1000.0)

h_c=h_c[::-1]
h_m=h_m[::-1]
z_m=z_m[::-1]
z_c=z_c[::-1]

fig = plt.figure(figsize=(10,10))
ax = plt.axes()

h_m_interp = np.interp(z_c,z_m,h_m)

ax.plot((1.0+z_c),abs(h_m_interp/h_c-1.0),'r-',linewidth=2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$1+z$',fontsize=20)
ax.set_ylabel('$|H(\mathrm{matteo})/H(\mathrm{class})-1|$',fontsize=20)
ax.set_ylim([1.0e-8,1.0e-3])
ax.set_xlim([1.,100.])


plt.figtext(0.5,0.915,'DEMNUNII - $M_\\nu=0.00 \, \mathrm{eV}$ - $w_0 = -1$, $w_a = 0$',fontsize=20,linespacing=2,horizontalalignment ='center')

plt.savefig("mnu000_w-1_Hubble_IC_vs_CLASS.pdf",dpi=2000)
plt.show()
'''
