##E0 = 1700
##theta1 = 170
##from masses import He
##Z1, M1 = He()
##from masses import As
##Z2, M2 = As()
##print 'M2, Z2=', M2, Z2
import numpy as np

def kinematic_factor_Mayer(theta1, M1, M2):
    first_brack = M1/(M2 + M1)
    theta = np.radians(theta1)
    if M1 < M2:
        second_brack = np.cos(theta) + np.sqrt((M2/M1)**2 - (np.sin(theta))**2)
    else:
        second_brack = np.cos(theta) - np.sqrt((M2/M1)**2 - (np.sin(theta))**2)
    return first_brack**2 * second_brack**2


##print 'kinematc_factor_Mayer=', kinematic_factor_Mayer(theta1, M1, M2)
##K = kinematic_factor_Mayer(theta1, M1, M2)
##E = K*E0
####energy_per_channel = 1.763
##print 'E=', E
####print 'channel=', E/energy_per_channel
