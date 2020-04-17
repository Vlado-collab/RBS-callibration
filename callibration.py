# -*- coding: utf-8 -*-
'''Enter filename with measured data'''

measurements = 'C:/Users/kolesar/Documents/cv/51_Schönaich bei Böblingen/04_depth_profiel_analysis/RBS measurements for depth profile/35_keV/2_Standard'

figure1 = 'C:/Users/kolesar/Documents/cv/51_Schönaich bei Böblingen/04_depth_profiel_analysis/results_FZR/multiplefigures/0073_2_Si_Wafer/35_keV/2_Standard.png'
'''Enter input energy'''
E0 = 1700

'''Enter incident angle'''
theta1 = 170

'''Enter types of atoms Z1, Z2, M1, M2'''
from masses import He
Z1, M1 = He()

from masses import C
Z2, M2 = C()

from masses import Si
Z3, M3 = Si()

from masses import Ni
Z4, M4 = Ni()

from masses import Au
Z5, M5 = Au()

'''Energy calculation'''
from kinematic_factor import kinematic_factor_Mayer
K_C = kinematic_factor_Mayer(theta1,M1, M2)
K_Si = kinematic_factor_Mayer(theta1,M1, M3)
K_Ni = kinematic_factor_Mayer(theta1,M1, M4)
K_Au = kinematic_factor_Mayer(theta1,M1, M5)

E_C = E0*K_C
E_Si = E0*K_Si
E_Ni = E0*K_Ni
E_Au = E0*K_Au
Energy = [E_C, E_Si, E_Ni, E_Au]

'''Set the limits'''
xmin2 = 200
xmax2  = 900
ymin2 = -250
ymax2 = 250

xmin3 = 210
xmax3 = 260
ymin3  = -200
ymax3 = 200

xmin4 = 500
xmax4 = 550

xmin5 = 690
xmax5 = 740

xmin6 = 840
xmax6 = 890

number_of_points = 50

'''Enter the fitting parameters: x0 is position of peak maximum
a is peak amplitude and c is parameter realted with thickness of peak'''
x0_C = 231.50353464
a_C  = 78.37
c_C  = 3.28

x0_Si= 525
a_Si = 70
c_Si = 3.28

x0_Ni= 720
a_Ni = 70
c_Ni = 3.28

x0_Au = 870
a_Au =70
c_Au = 3.28

'''spectrum data'''
from spectrum_reader import get_table2
channel1, counts1 = get_table2(measurements)

'''difference calc.'''
difference1 = []
for i in range(len(channel1)):
    dy1 = counts1[i] - counts1[i-1]
    difference1.append(dy1)
import numpy as np
x3 = np.linspace(channel1[xmin3], channel1[xmax3], number_of_points)
x4 = np.linspace(channel1[xmin4], channel1[xmax4], number_of_points)
x5 = np.linspace(channel1[xmin5], channel1[xmax5], number_of_points)
x6 = np.linspace(channel1[xmin6], channel1[xmax6], number_of_points)

difference_C  = []
difference_Si = []
difference_Ni = []
difference_Au = []
for i in range(xmin3,xmax3,1):
    difference_C.append(difference1[i])
for i in range(xmin4,xmax4,1):
    difference_Si.append(difference1[i])
for i in range(xmin5,xmax5,1):
    difference_Ni.append(difference1[i])
for i in range(xmin6,xmax6,1):
    difference_Au.append(difference1[i])

'''Smoothing function'''
def smooth(difference_C, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(difference_C, box, mode='same')
    return y_smooth
y_C = smooth(difference_C, 3)
y_Si = smooth(difference_Si, 3)
y_Ni = smooth(difference_Ni, 3)
y_Au = smooth(difference_Au, 3)

'''Fitting function'''
def gauss_function_C(x3, a_C, x0_C, c_C):
    return -a_C*np.exp(-(x3-x0_C)**2/(2*c_C**2))
def gauss_function_Si(x4, a_Si, x0_Si, c_Si):
    return -a_Si*np.exp(-(x4-x0_Si)**2/(2*c_Si**2))
def gauss_function_Ni(x5, a_Ni, x0_Ni, c_Ni):
    return -a_Ni*np.exp(-(x5-x0_Ni)**2/(2*c_Ni**2))
def gauss_function_Au(x6, a_Au, x0_Au, c_Au):
    return -a_Au*np.exp(-(x6-x0_Au)**2/(2*c_Au**2))

'''Fitting procedure'''
from scipy.optimize import curve_fit
popt_C, pcov_C = curve_fit(gauss_function_C, x3, y_C, p0 = [1, x0_C, c_C])
popt_Si, pcov_Si = curve_fit(gauss_function_Si, x4, y_Si, p0 = [1, x0_Si, c_Si])
popt_Ni, pcov_Ni = curve_fit(gauss_function_Ni, x5, y_Ni, p0 = [1, x0_Ni, c_Ni])
popt_Au, pcov_Au = curve_fit(gauss_function_Au, x6, y_Au, p0 = [1, x0_Au, c_Au])
yfitted_C = gauss_function_C(x3, *popt_C)
yfitted_Si = gauss_function_Si(x4, *popt_Si)
yfitted_Ni = gauss_function_Ni(x5, *popt_Ni)
yfitted_Au = gauss_function_Au(x6, *popt_Au)

Channel = [popt_C[1], popt_Si[1], popt_Ni[1], popt_Au[1]]
from pylab import *
fit = polyfit(Channel,Energy,1)
fit_fn = poly1d(fit)
q = fit_fn[0]
k = fit_fn[1]

outfile1 = open('C:/Users/kolesar/Documents/cv/51_Schönaich bei Böblingen/04_depth_profiel_analysis/callibration.dat', 'w')
outfile1.write('%12.8e %12.8e\n' % (q,k))
outfile1.close()

'''Figures'''
##Energy = [E_C, E_Si, E_Ni, E_Au]

import matplotlib.pyplot as plt
plt.figure(figsize = (10,6))
plt.subplot(2, 3, 1)
#plt.title(measurements +' 1700 keV')
plt.plot(channel1, counts1, label ='spectrum')
plt,annotate("$E_C =%g keV$ " %(E_C), xy = (300,3000))
plt,annotate("$E_{Si} =%g keV$ " %(E_Si), xy = (300,2500))
plt,annotate("$E_{Ni} =%g keV$ " %(E_Ni), xy = (300,2000))
plt,annotate("$E_{Au} =%g keV$ " %(E_Au), xy = (300,1500))
plt.xticks(np.arange(0, 1201, 300))
plt.yticks(np.arange(1000, 4001, 1000))
plt.legend(prop={'size':12})
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.grid()

plt.subplot(2, 3, 2)
plt.title('Energy per channel')
plt.plot(Channel,Energy,'bo', label = 'calc.')
plt.plot(Channel, fit_fn(Channel),'r-' ,label='fit')
plt.xticks(np.arange(200, 901, 200))
plt.yticks(np.arange(600, 1601, 300))
plt.annotate("$k =%g$" %(fit_fn[1]), xy =(300, 1400))
plt.annotate("$q =%g$" %(fit_fn[0]), xy =(300, 1200))
plt.legend(prop={'size':12}, loc = 4)
plt.xlabel('Channel')
plt.ylabel('Energy [keV]')
plt.grid()

plt.subplot(2, 3, 3)
plt.title('difference_C')
plt.plot(x3,difference_C,'bo-',  label = 'diff')
plt.plot(x3,smooth(difference_C,3),'g-', lw=2, label = 'smooth')
plt.plot(x3, yfitted_C,'r-' ,label='fit')
plt.legend(prop={'size':12})
plt.xlim(xmin3, xmax3)
plt.ylim(ymin3, ymax3)
plt.xlabel('Channel')
plt.ylabel('Yield derivative')
plt.grid()

plt.subplot(2,3,4)
plt.title('difference_Si')
plt.plot(x4, difference_Si, 'bo-', label = 'diff')
plt.plot(x4, smooth(difference_Si,3), 'g-', lw=2, label = 'smooth')
plt.plot(x4, yfitted_Si, 'r-', label='fit')
plt.xlim(xmin4, xmax4)
plt.xlabel('Channel')
plt.ylabel('Yield derivative')
plt.grid()
plt.legend(prop={'size':12})

plt.subplot(2,3,5)
plt.title('difference_Ni')
plt.plot(x5, difference_Ni, 'bo-', label = 'diff')
plt.plot(x5, smooth(difference_Ni,3), 'g-', lw=2, label = 'smooth')
plt.plot(x5, yfitted_Ni, 'r-', label='fit')
plt.xlim(xmin5, xmax5)
plt.xlabel('Channel')
plt.ylabel('Yield derivative')
plt.grid()
plt.legend(prop={'size':12})


plt.subplot(2,3,6)
plt.title('difference_Au')
plt.plot(x6, difference_Au, 'bo-', label = 'diff')
plt.plot(x6, smooth(difference_Au,3), 'g-', lw=2, label = 'smooth')
plt.plot(x6, yfitted_Au, 'r-', label='fit')
plt.xlim(xmin6, xmax6)
plt.xlabel('Channel')
plt.ylabel('Yield derivative')
plt.grid()
plt.legend(prop={'size':12})

plt.tight_layout()
##figure1 = 'Results_FZR/multiplefigures/figure.png'
plt.savefig(figure1)
