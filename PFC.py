
###############################################################################
class PFC():
    r""" This class provides the properties of several perfluorocarbon liquids. 
    Namely, the density :math:`\rho` (in liquid phase), 
    the surface tension :math:`\gamma` (between PCF and air),
    the molecular weight :math:`M_{w}`,
    and the boiling temperature :math:`T_{b}`.

    :param abbr: perfluorocarbon to be considered, could be PFH
    :type command: string

    """
    
    def  __init__(self,abbr):
        r""" The values for PFH are 
        :math:`\rho = 1670` kg/m\ :sup:`3` ,
        :math:`\gamma = 0.056` N/m,  
        :math:`M_{w}=0.338` kg/mol, 
        and :math:`T_{b} = 330.15` K.
        """
        if abbr == 'PFH':
            self.rho = 1670 # in kg/m:sup:`3` at 25°C
            self.gamma_0 = 0.056 # in N/m
            self.Mw = 0.33806 # in kg/mol
            self.Tb = 57 + 271.15 # in K
            self.Pv = 29330 # Pa

        elif abbr == 'PFP':
            self.rho = 1603 # in kg/m:sup:`3` at 25°C
            self.gamma_0 = 0.056 # in N/m
            self.Mw = 0.08805 # in kg/mol
            self.Tb = 29 + 271.15 # in K
            self.Pv = 73990 #Pa

###############################################################################

