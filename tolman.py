
###############################################################################
def gamma_tolman(R,g,d):
    r""" It derives the effective interface tension of a curved surface
    
    .. math::
        \gamma(R) = \frac{\gamma(\infty)}{1+\frac{2\delta}{R^2}}
        
    :param R: Radius :math:`R` of the surface curvature, in m.
    :type command: float

    :param g: Interface tension :math:`\gamma` of the interface, in N/m.
    :type command: float

    :param d: Tolman length :math:`\delta`, in m.
    :type command: float
        
    :return: :math:`\gamma(R)`, the effective interface tension of the interface, in N/m.
    :rtype: float        
                
    """
    value = g/(1+2*d/R)
    return value
###############################################################################

###############################################################################
def gamma_infty(R,g,d):
    r""" It derives the effective interface tension of a curved surface

    .. math::
        \gamma(\infty)  = \gamma(R)(1+\frac{2\delta}{R^2})

    :param R: Radius :math:`R` of the surface curvature, in m.
    :type command: float

    :param g: Interface tension :math:`\gamma` of the interface, in N/m.
    :type command: float

    :param d: Tolman length :math:`\delta`, in m.
    :type command: float

    :return: :math:`\gamma(R)`, the effective interface tension of the interface, in N/m.
    :rtype: float

    """
    value = g*(1+2*d/R)
    return value
###############################################################################


"""
###############################################################################
def plotTolman():
    
    marker_size = 12
    markerwidth = 2
    width = 1
    cap = 3
    
    fig = matplotlib.pyplot.figure(num=0,figsize=(7,5),facecolor='w',edgecolor='w')
    ax = fig.add_subplot(111)
    matplotlib.pyplot.subplots_adjust(left=0.150,right=0.97,bottom=0.170,top=0.960,hspace = 0.12,wspace = 0)
    matplotlib.pyplot.tick_params(axis='both',which='both',direction ='in',bottom=True,top=True,left=True,right=True)
    
    R = numpy.arange(0.001,30,0.001)*1e-6

    g = 0.035
    d = 10e-10
    y = gamma_tolman(R,g,d)
    matplotlib.pyplot.plot(R*1e6, y*1e3,'b', linewidth=2)
    d = 1e-10
    y = gamma_tolman(R,g,d)
    matplotlib.pyplot.plot(R*1e6, y*1e3,'b:', linewidth=2)


    g = 0.012
    d = 10e-10
    y = gamma_tolman(R,g,d)
    matplotlib.pyplot.plot(R*1e6, y*1e3,'r', linewidth=2)
    d = 1e-10
    y = gamma_tolman(R,g,d)
    matplotlib.pyplot.plot(R*1e6, y*1e3,'r:', linewidth=2)    
    
    
    matplotlib.pyplot.xlabel(r'$R$ ($\mu$m)',fontsize = 28)
    matplotlib.pyplot.ylabel(r'$\gamma(R)$ (mN/m)',fontsize = 28)
    matplotlib.pyplot.xticks(fontsize = 24)
    matplotlib.pyplot.yticks(fontsize = 24)
    matplotlib.pyplot.xlim(0,0.5)
    matplotlib.pyplot.ylim(0,40)
    #ax.set_xticks([0,0.25,0.5,0.75,1])
    #ax.xaxis.set_minor_locator(MultipleLocator(1))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    matplotlib.pyplot.savefig('../tolman.png',dpi=600)
###############################################################################  
#plotTolman() 
#matplotlib.pyplot.show()
"""
