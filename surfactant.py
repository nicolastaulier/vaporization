
###############################################################################
class surfactant():
    r""" This class provides the surface tensions according to the type of surfactant used.

    :param abbr: surfactant to be considered, could be either FTAC or KRYTOX
    :type command: string

    """
    
    def  __init__(self,abbr):
        r""" 
        """
        if abbr == 'FTAC':
            self.g_gl = 0.012 # in N/m
            self.g_lw = 0.025 # in N/m
            self.g_gw = 0.038 # in N/m
            
        elif abbr == 'KRYTOX':
            self.g_gl = 0.012 # in N/m
            self.g_lw = 0.018 # in N/m
            self.g_gw = 0.066 # in N/m
            
###############################################################################

