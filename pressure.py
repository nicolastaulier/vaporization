import math
import numpy as np

###############################################################################
def P_out(P,R,g_lw):
    r""" It derives the pressure out of a droplet from the pressure inside a droplet
    
    .. math::
        P_{out} = P_{in} - \frac{2\gamma}{R}
        
    :param P: Pressure :math:`P_{in}` inside the droplet, in Pa.
    :type command: float

    :param R: Pressure :math:`R` of the droplet, in m.
    :type command: float

    :param g: Interfacial tension :math:`\gamma` of the interface between the inside and outside space of the droplet, in N/m.
    :type command: float
        
    :return: :math:`P_{out}` outside the droplet, in Pa.
    :rtype: float        
                
    """
    value = P - 2*g_lw/R
    return value
###############################################################################

###############################################################################
def P_in(P, R, g_lw):
    r""" It derives the pressure inside droplet from the pressure outside a droplet

    .. math::
        P_{in} = P_{out} +  \frac{2\gamma}{R}

    :param P: Pressure :math:`P_{out}` outside the droplet, in Pa.
    :type command: float

    :param R: Pressure :math:`R` of the droplet, in m.
    :type command: float

    :param g: Interfacial tension :math:`\gamma` of the interface between the inside and outside space of the droplet, in N/m.
    :type command: float
        
    :return: :math:`P_{in}` inside the droplet, in Pa.
    :rtype: float        
                
    """
    value = P + 2*g_lw/R
    return value
###############################################################################

###############################################################################
def c_n(n):
    value = 1/np.sqrt(2*np.log(n))
    return value
###############################################################################

###############################################################################
def d_n(n):
    value = np.sqrt(2*np.log(n)) - (math.log(4*math.pi) + np.log(np.log(n)))/(2*np.sqrt(2*np.log(n)))
    return value
###############################################################################

###############################################################################
def p_n(n,P,P1,sigma1):
    mu_n = P1 - sigma1*d_n(n)
    beta_n = sigma1*c_n(n)
    value = 1 - np.exp(-np.exp((P-mu_n)/beta_n))
    return value
###############################################################################

###############################################################################
def P_n(n,P1,sigma):
    value  = P1 - sigma*(d_n(n) - c_n(n)*math.log(math.log(2)))
    return value
###############################################################################
