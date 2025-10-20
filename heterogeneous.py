import math
import numpy
import scipy.constants

import pressure
import surfactant
import tolman


# -*- coding: utf-8 -*-
import math
import numpy
import scipy.constants

import pressure
import surfactant
import tolman

###############################################################################
class surface():

    def __init__(self, PFC, T, f, g_gl, g_lw, g_gw, g_gw_0, k, delta, a, b, c):
        r""" It set the values of the interface tension between PFC liquid/water (:math:`\gamma_{lw}`),
        PFC gas/water (:math:`\gamma_{gw}`) and PFC gas/ PFC liquid (:math:`\gamma_{gl}`).
        The class PFC is asked.

        :param T: The temperature :math:`T`, in K.
        :type command: float

        :param freq: The acoustic frequency :math:`f` of the acoustic wave, in Hz.
        :type command: float

        :param PFC: Type of perfluorocarbone.
        :type command: string

        """
        self.PFC = PFC
        self.T = T
        self.f = f
        self.g_gl = g_gl
        self.g_lw = g_lw
        self.g_gw = g_gw
        self.delta = delta
        self.k = k
        self.g_gw_0 = g_gw_0
        self.ratio = 1
        self.r1_final = 1e-10
        self.r2_final = 1e-10
        self.d_final = 1e-10
        self.deltaP_in = 1e6
        self.m_final = 0.5
        self.g_gw_r2 = 0.03
        self.a = a
        self.b = b
        self.c = c
        if a == 0 and b == 0 and c == 0 :
            self.g_gw_case = 'unknown'
        else:
            self.g_gw_case = 'known'

    def Wc_rigid(self,P_out, R):
        """ Derive the energy to create a nucleus on a solid surface
        """
        deltaP = pressure.P_in(P_out,R,self.g_lw) - self.PFC.Pv
        g_gw_R = tolman.gamma_tolman(R,self.g_gw,self.delta)
        r =  self.r1(deltaP)
        g_gl_r = tolman.gamma_tolman(r,self.g_gl,self.delta)

        m = (self.g_lw - g_gw_R)/g_gl_r

        v = self.v(r,R,m)
        a_1 = self.a_r1(r, R, m)
        a_2 = self.a_r2(r, R, m)
        value2 = (g_gl_r*a_1 + (g_gw_R - self.g_lw)*a_2 - deltaP*v)/(scipy.constants.k*self.T)

        x = R*deltaP/(2*g_gl_r)
        value = 16*math.pi*g_gl_r**3*self.h(x,m)/(3*math.pow(deltaP,2)*scipy.constants.k*self.T)
        #print(value2,value)
        return (value2, m, r)

    def g_gw_r2_c(self, r2, m_phi, R):
        a_gw_R_over_a_gw_r2 = (R/r2)**2*(1 - math.sqrt(1 - (r2/R)**2*(1 - m_phi**2)))/(1 - m_phi)
        G = 1/(scipy.constants.N_A*87e-20)
        G_max = 2e-6
        ratio = a_gw_R_over_a_gw_r2*G/G_max
        value_C = self.g_gw_0 + G_max*scipy.constants.R*self.T*(math.log(1 - ratio) - 0.5*self.k*ratio**2)
        value = value_C/(1 + 2*self.delta/r2)
        return value
    """
    def g_gw_r2(self, g_lw_R, g_gl, m_chi, m_psi, m_phi):
        value = (g_lw_R*m_chi - g_gl*m_psi)/m_phi
        return value
    """

    def r1(self, deltaP, g_gl):
        """ Derive the critical radius of the nucleus inside the droplet

        .. math::
                r^{*}_{1} = \frac{2\gamma_{gl}}{P_{v}-P_{l}} - 2\delta


        :param DeltaP: Pressure :math:`\Delta P`
        :type command: float

        :param g_gl: interfacial tension between the gaseous and liquid perfluorocarbone :math:`\gamma_{gl}`
        :type command: float

        :param delta: Tolman lenght :math:`\delta`
        :type command: float

        :return: The critical radius of nucleation :math:`r_{1}`
        :rtype: float
        """
        value = 2*g_gl/deltaP - 2*self.delta
        return value

    """
    def m_psi(self, r1, r2, m_phi):
        #value = (r2*m + r1)/numpy.sqrt(r2**2 + r1**2 + 2*r1*r2*m)
        #print(round((r2/r1)**2*(1 - m_phi**2) - 1,3), round(1 - m_phi**2,3))
        value = math.sqrt((r2/r1)**2*(1 - m_phi**2) - 1)
        return value
    """
    def m_chi(self, R, r2, m_phi):
        #value = numpy.sqrt(1 - numpy.power(r2/R,2)*(1 - numpy.power(self.m_phi(r1, r2, m),2)))
        value = math.sqrt(1 - (r2/R)**2*(1 - m_phi**2))
        return value

    def m(self,r1, r2, R, g_gl_r1, g_gw_r2):
        a = (g_gl_r1*r2 + g_gw_r2*r1)**2 - self.g_lw**2*(r1**2*r2**2)/R**2
        b = 2*(g_gl_r1*r2 + g_gw_r2*r1)*(g_gl_r1*r1 + g_gw_r2*r2) - 2*self.g_lw**2*r1*r2
        c = (g_gl_r1*r1 + g_gw_r2*r2)**2 + self.g_lw**2*((r1**2*r2**2)/R**2 - r1**2 - r2**2)
        sqrtDelta = math.sqrt(b**2 - 4*a*c)
        m1 = (-b - sqrtDelta)/(2*a)
        m2 = (-b + sqrtDelta)/(2*a)
        #print(m1,m2)
        if abs(m1) > 1:
            if abs(m2) > 1:
                #print('m > 1')
                if m2 < m1:
                    m = m2
                else:
                    m = m1
            else:
                m = m2
        else:
            m = m1
        return m

    def test_m(self,Param, args):
        #print(args)
        #r1, r2, R, g_gl_r1, g_gw_r2, g_lw_R = args
        r1 = args[0][0]
        r2 = args[0][1]
        R = args[0][2]
        g_gl_r1 = args[0][3]
        g_gw_r2 = args[0][4]
        g_lw_R = args[0][5]
        m = Param[0]
        #print(m)
        a = (g_gl_r1*r2 + g_gw_r2*r1)**2 - g_lw_R**2*(r1**2*r2**2)/R**2
        b = 2*(g_gl_r1*r2 + g_gw_r2*r1)*(g_gl_r1*r1 + g_gw_r2*r2) - 2*g_lw_R**2*r1*r2
        c = (g_gl_r1*r1 + g_gw_r2*r2)**2 + g_lw_R**2*((r1**2*r2**2)/R**2 - r1**2 - r2**2)
        value = a*m**2+b*m+c
        #print(value)
        return [value]

    def func_g_gw_P(self):
        value = self.a*self.PFC.Pv**2 + self.b*self.PFC.Pv + self.c
        return value

    def func_g_gw_R(self,R):
        value = self.a*numpy.exp(R*self.b) + self.c #self.a*R**2 + self.b*R + self.c
        return value

    def func_g_gw_R_linear(self,R):
        value = self.a*R + self.b
        return value

    def W_05(self,m,A,g_gl_r):
        """
        .. math::
            W_{0.5} = k_{b}T \frac{\ln(\Pi_{0)A}{2\ln(2)f}


        """
        value = math.log((self.Pi_0(m,g_gl_r)*A)/(2*math.log(2)*self.f))
        return value

    def Pi_0(self,m,g_gl_r):
        r""" It provides the surface nucleation rate :math:`\Pi_{0}` (see https://doi.org/10.1002/aic.690210502):

        .. math::
                \Pi_{0} = (N_A \rho)^{2/3} \left( \frac{1-m}{2} \right) \sqrt{\frac{2\gamma_{lw}}{\pi M}}

        where

        * :math:`N_A` is the Avogadro number,
        * :math:`\rho` is the density of the liquid in which nucleation occurs,
        * :math:`\gamma_{lw}` is the interface tension of the surface where nucleation occurs,
        * :math:`M` is the mass of a molecule of the liquid (that is :math:`M_{w}/N_{A}`, where :math:`M_{w}` is the molecular weigh of the molecule),
        * finally  :math:`m` is the ratio :math:`\frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}} = \cos\theta`.

        :param m:  It is the ratio :math:`m = \frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}} = \cos\theta`.
        :type command: float

        :param g_lw: Interfacial tension :math:`\gamma_{lw}` between the PFC liquid and the water inside the droplet, in N/m.
        :type command: float

        :return: The surface nucleation rate :math:`\Pi_{0}`
        :rtype: float

        """
        value = math.pow(scipy.constants.N_A*self.PFC.rho/self.PFC.Mw,2/3)*(1 - m)*0.5*math.sqrt(2*g_gl_r/(self.PFC.Mw*math.pi/scipy.constants.N_A))
        return value

    def func_g_gw(self,Param,args):
        self.g_gw_r2 = Param[0]
        (R, surfaceType, P_out) = args[0]
        value = self.Wdiff(P_out,R,surfaceType)
        return [value]

    def func_delta_2(self,Param,args):
        self.delta_2 = Param[0]
        (R, surfaceType, P_out) = args[0]
        value = self.Wdiff(P_out,R,surfaceType)
        return [value]

    def d(self, r1, r2, m):
        value = math.sqrt(r1**2 + r2**2 + 2*r1*r2*m)
        self.d_final = value
        return value

    def  m_psi(self, r1, r2, m):
        value = (r1 + r2*m)/self.d(r1, r2, m)
        return value

    def m_phi(self, r1, r2, m):
        value = (r1*m + r2)/self.d(r1, r2, m)
        return value

    def m_chi(self, r1, r2, R, m):
        value = math.sqrt(1 - (r2/R)**2*(1 - self.m_phi(r1, r2, m)**2))
        return value

    def area(self,r,m):
        value = 2*math.pi*r**2*(1-m)
        return value

    def a_r1(self, r1, r2, m):
        value = self.area(r1,self.m_psi(r1,r2,m))
        return value

    def a_r2(self, r1, r2, m):
        value = self.area(r2,self.m_phi(r1,r2,m))
        return value

    def a_R(self, r1, r2, R, m):
        value = self.area(R,self.m_chi(r1, r2, R, m))
        return value

    def v(self, r1, r2, m):
        r""" It derives the volume :math:`v`  of a gas nuclei

        .. math::
            v = \frac{\pi}{3} \left[ R^{3} (2 - 3\cos\phi + \cos^{3}\phi) + r^{3}(2 - 3 \cos\psi +  \cos^{3}\psi) \right]

        Where :math:`R` and :math:`r` are the radius of the droplet and  of the nuclei, respectively, and

        .. math::
            \cos\psi = \frac{Rm+r}{d}

        .. math::
            \cos\phi = \frac{R+rm}{d}

        :param r: Radius :math:`r` of the nuclei, in m.
        :type command: float

        :param R: Radius :math:`R` of the droplet, in m.
        :type command: float

        :return: Volume :math:`v` of a nuclei, in m\ :sup:`3`.
        :rtype: float

        """
        value = math.pi/3*(math.pow(r1,3)*(2 - 3*self.m_psi(r1, r2, m) + math.pow(self.m_psi(r1, r2, m),3)) + math.pow(r2,3)*(2 - 3*self.m_phi(r1, r2, m) + math.pow(self.m_phi(r1, r2, m),3)))
        return value

    def h(self,x,m):
        r""" It derives the function :math:`f`

        .. math::
            f(m,x) = \frac{1}{2}\left\lbrace 1- \left( \frac{1+mx}{g} \right)^{3} - x^{3} \left[2-3\left(\frac{x+m}{g}\right) + \left(\frac{x+m}{g}\right)^{3}\right] -3mx^{2}\left(1-\frac{x+m}{g}\right)\right\rbrace

        with :math:`g = \sqrt{1 + x^{2} + 2mx}`

        .. math::
           h(m,x) = \left(\frac{m x-1}{2}+ x^{2}\right)\sqrt{1 + x^{2}+2mx}    - 1.5 m x^{2}  - x^{3}  + \frac{1}{2}




        :param x: Parameter :math:`x = \frac{R}{r^{*}}`, without unit.
        :type command: x

        :param m: Parameter :math:`m = \frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}}`, without unit.
        :type command: float

        :return: :math:`f`, without unit.
        :rtype: float

        """
        #g = math.sqrt(1 + math.pow(x,2) + 2*m*x)
        #t  =  (x + m)/g
        #value = 0.5*(1 - math.pow((1 + m*x)/g,3) - math.pow(x,3)*(2 - 3*t + math.pow(t,3) ) - 3*m*math.pow(x,2)*(1-t))
        value =   -1.5*m*x**2 + 0.5*m*x*math.sqrt(2*m*x + x**2 + 1) - x**3 + x**2*math.sqrt(2*m*x + x**2 + 1) - 0.5*math.sqrt(2*m*x + x**2 + 1) + 0.5
        #value = (0.5*(m*x-1) + x**2)*math.sqrt(1+x**2+2*m*x) - 1.5*m*x**2 - x**3 + 0.5
        return value

###############################################################################

###############################################################################
class external(surface):

    def __init__(self, PFC, T, f, g_gl, g_lw, g_gw, g_gw_0, k, delta, a, b, c):
        super().__init__(PFC, T, f, g_gl, g_lw, g_gw, g_gw_0, k, delta, a, b, c)

    def r2_func(self, deltaP, R, r1, g_lw, g_gw_r2, g_gl_r1):
        value = 2*self.g_gw_r2/(deltaP + 2*self.g_lw/R)
        #or alternatively value = g_gw_r2/(g_gl_r1/r1 + g_lw/R)
        return value

    def Wc_soft(self, P_out, R, g_lw, g_gl):
        deltaP =  self.PFC.Pv - pressure.P_in(P_out, R, g_lw)
        self.deltaP_in = deltaP
        r1 = self.r1(deltaP, g_gl)
        self.r1_final = r1
        g_gl_r1 = tolman.gamma_tolman(r1,g_gl,self.delta)

        if self.g_gw_case == 'known':
            #self.g_gw_r2 =  self.func_g_gw_P(deltaP)
            self.g_gw_r2 =  self.func_g_gw_R(R)

        r2 = self.r2_func(deltaP, R, r1, g_lw, self.g_gw_r2, g_gl_r1)
        self.r2_final = r2

        m = self.m(r1, r2, R, g_gl_r1, self.g_gw_r2)
        self.m_final = m

        v = self.v(r1, r2, m)
        a_r1 = self.a_r1(r1, r2, m)
        a_r2 = self.a_r2(r1, r2, m)
        a_R = self.a_R(r1, r2, R, m)
        self.ratio  = a_r2/a_R

        value = (g_gl_r1*a_r1 + self.g_gw_r2*a_r2 - self.g_lw*a_R - deltaP*v)/(scipy.constants.k*self.T)
        return (value, m, r1,r2)

    def Wdiff(self,P_out,R,surfaceType):
        r""" It derives the value of the difference :math:`(W^{*}-W_{0.5})` that is equal to:

        .. math::
            \frac{16\pi\gamma^{3}_{gl}}{3P_{in}^{2}} h(m,x) - k_{B}T \ln\left(\frac{\Pi_{0}nA}{2f\ln(2)}\right)

        where

        * :math:`P_{in} = P_{0.5} + \frac{2\gamma_{lw}}{R}` is the acoustic pressure inside the droplet,
        *  :math:`P_{0.5}` is the acoustic pressure outside the droplet,
        * :math:`k_{B}` is the Boltzmann constant,
        * :math:`f` is the acoustic frequency,
        * :math:`T` is the temperature,
        * :math:`A` is the droplet surface (:math:`= 4\pi R^{2}` where :math:`R` is the droplet radius),
        * :math:`\gamma_{gl}`, is the interfacial tension between gaseous PFC and liquid PFC,
        * :math:`\gamma_{lw}`, is the interfacial tension between liquid PFC liquid and water,
        * :math:`\gamma_{gw}`, is the interfacial tension between gaseous PFC liquid and water,
        * :math:`m = \frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}}`
        * :math:`\Pi_{0}` is the surface nucleation rate.

        :param P_out: The acoustic pressure outside the droplet:math:`P_{0.5}`, in Pa.
        :type command: float

        :param R: Radius :math:`R` in meter of a droplet on which surface nucleation occurs, in m
        :type command: float

        :param n: Number :math:`n` of droplets.
        :type command: float

        :param g_lw: Interfacial tension :math:`\gamma_{lw}` between the PFC liquid and the water inside the droplet, in N/m.
        :type command: float

        :param g_gw: Interfacial tension :math:`\gamma_{gw}` between the PFC gas inside the nuclei and the water droplets, in N/m.
        :type command: float

        :param g_gl: Interfacial tension :math:`\gamma_{lw}` between the PFC gas inside the nuclei and the PFC liquid, in N/m.
        :type command: float

        :return:  the difference :math:`(W^{*}-W_{0.5})`
        :rtype: float

        """
        A = 4*math.pi*math.pow(R,2)
        g_lw = self.g_lw
        g_gl = self.g_gl

        if surfaceType == 'rigid':
            (Wc, m, r) = self.Wc_rigid(P_out, R)
            g_gl_r =  tolman.gamma_tolman(r, g_gl, self.delta)
            value = self.W_05(m, A, g_gl_r) - Wc

        elif surfaceType == 'soft':
            (Wc, m, r1, r2) = self.Wc_soft(P_out, R, g_lw, g_gl)
            g_gl_r1 =  tolman.gamma_tolman(r1, g_gl, self.delta)
            value = self.W_05(m, A, g_gl_r1) - Wc
        return value

    def func_P(self,Param,args):
        P_out = Param[0]
        (R, surfaceType) = args[0]
        value = self.Wdiff(P_out, R, surfaceType)
        return [value]
###############################################################################


###############################################################################
class internal(surface):

    def __init__(self, PFC, T, f, g_gl, g_lw, g_gw, g_gw_0, k, delta, a, b, c, gi_gl, gi_lw, gi_gw):
        super().__init__(PFC, T, f, g_gl, g_lw, g_gw, g_gw_0, k, delta, a, b, c)
        self.ge_gl = g_gl
        self.gi_gl = gi_gl
        self.gi_lw = gi_lw
        self.gi_gw = gi_gw
        self.gi_gw_r2 = 0.03

    def r2_func(self,deltaP, R_w, gi_gw_r2):
        value = 2*gi_gw_r2/(deltaP - 2*gi_lw/R_w)
        return value

    def Wc_soft(self, P_out, R, R_w, ge_lw, gi_gl):
        deltaP =  self.PFC.Pv - pressure.P_in(P_out, R, ge_lw)
        self.deltaP_in = deltaP
        r1 = self.r1(deltaP, gi_gl)
        self.r1_final = r1
        gi_gl_r1 = tolman.gamma_tolman(r1,gi_gl,self.delta)

        if self.gi_gw_case == 'known':
            #self.g_gw_r2 =  self.func_g_gw_P(deltaP)
            self.gi_gw_r2 =  self.func_g_gw_R(R_w)

        r2_func(self,deltaP, R_w, gi_gw_r2)
        self.r2_final = r2

        m = self.m(r1, r2, R_w, g_gl_r1, self.g_gw_r2)
        self.m_final = m

        v = self.v(r1, r2, m)
        a_r1 = self.a_r1(r1, r2, m)
        a_r2 = self.a_r2(r1, r2, m)
        a_R = self.a_R(r1, r2, R_w, m)
        self.ratio  = a_r2/a_R
        value = (g_gl_r1*a_r1 + self.g_gw_r2*a_r2 - self.g_lw*a_R - deltaP*v)/(scipy.constants.k*self.T)

        return (value, m, r1,r2)

    def Wdiff_internal(self,P_out,R, R_w,surfaceType):
        r""" It derives the value of the difference :math:`(W^{*}-W_{0.5})` that is equal to:

        .. math::
            \frac{16\pi\gamma^{3}_{gl}}{3P_{in}^{2}} h(m,x) - k_{B}T \ln\left(\frac{\Pi_{0}nA}{2f\ln(2)}\right)

        where

        * :math:`P_{in} = P_{0.5} + \frac{2\gamma_{lw}}{R}` is the acoustic pressure inside the droplet,
        *  :math:`P_{0.5}` is the acoustic pressure outside the droplet,
        * :math:`k_{B}` is the Boltzmann constant,
        * :math:`f` is the acoustic frequency,
        * :math:`T` is the temperature,
        * :math:`A` is the droplet surface (:math:`= 4\pi R^{2}` where :math:`R` is the droplet radius),
        * :math:`\gamma_{gl}`, is the interfacial tension between gaseous PFC and liquid PFC,
        * :math:`\gamma_{lw}`, is the interfacial tension between liquid PFC liquid and water,
        * :math:`\gamma_{gw}`, is the interfacial tension between gaseous PFC liquid and water,
        * :math:`m = \frac{\gamma_{lw}-\gamma_{gw}}{\gamma_{gl}}`
        * :math:`\Pi_{0}` is the surface nucleation rate.

        :param P_out: The acoustic pressure outside the droplet:math:`P_{0.5}`, in Pa.
        :type command: float

        :param R: Radius :math:`R` in meter of a droplet on which surface nucleation occurs, in m
        :type command: float

        :param n: Number :math:`n` of droplets.
        :type command: float

        :param g_lw: Interfacial tension :math:`\gamma_{lw}` between the PFC liquid and the water inside the droplet, in N/m.
        :type command: float

        :param g_gw: Interfacial tension :math:`\gamma_{gw}` between the PFC gas inside the nuclei and the water droplets, in N/m.
        :type command: float

        :param g_gl: Interfacial tension :math:`\gamma_{lw}` between the PFC gas inside the nuclei and the PFC liquid, in N/m.
        :type command: float

        :return:  the difference :math:`(W^{*}-W_{0.5})`
        :rtype: float

        """
        A = 4*math.pi*math.pow(R_w,2)
        gi_lw = self.gi_lw
        gi_gl = self.gi_gl

        (Wc, m, r1, r2) = self.Wc_soft(P_out, R, R_w, gi_lw, gi_gl)
        gi_gl_r1 =  tolman.gamma_tolman(r1, gi_gl, self.delta)
        value = self.W_05(m, A, gi_gl_r1) - Wc
        return value

    def func_P(self,Param,args):
        P_out = Param[0]
        (R, R_w, surfaceType) = args[0]
        value = self.Wdiff_internal(P_out, R, R_w, surfaceType)
        return [value]
###############################################################################









