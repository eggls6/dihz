import numpy as np
from .seff import *

__all__=['SSHZ']

#######################################
# Single star habitable zone
######################################

def SSHZ(L, teff):
    """ Calculate single star habitable zone limits
    according to Kopparapu et al. (2014).

    Parameters:
    ----------
    L    ... stellar luminosity [solar luminosities]
    teff ... effective stellar temperature [K]

    Returns:
    -------
    dinner ... distance to host star of inner HZ border [au]
    douter ... distance to host star of outer HZ border [au]

    Requires:
    --------
    Functions sinner, souter
    """
    return [np.sqrt(L/seffi(teff)), np.sqrt(L/seffo(teff))]

