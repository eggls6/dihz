#!/bin/python
import numpy as np
from .seff import *

################################
# Circumstellar (S-type systems)
###############################

sqrt=np.sqrt

__all__=['reqb','reqpS','eforcedS','epmaxSe0','epmaxS','eav2S','eav2Se0','PHZ','AHZ']

#######################################
# Insolation equivalent orbit distance
######################################
def reqb(a, e):
    """Insolation equivalent orbit distance

    Parameters:
    ----------
    a ... semimajor axis
    e ... orbital eccentricity
    """
    reqp = a*(1.-e*e)**(0.25)
    return reqp


def eforcedS(ab, eb, ap):
    """Forced eccentricity estimate of planetary orbit in S-type
    binary system.

    Parameters:
    ----------
    ab ... binary orbit semimajor axis [au]
    eb ... binary orbit eccentricity []
    ap ... planetary orbit semimajor axis [au]

    Returns:
    -------
    emax ... maximum eccentricity of circumstellar orbit []
    """
    emax = 5./2.*ap/ab*(eb/(1.-eb*eb))
    return emax

def epmaxS(ab, eb, ap):
    """Maximum eccentricity of planetary orbit in S-type binary system
    assuming the planet stats out on its forced orbit.

    Parameters:
    ----------
    ab ... binary orbit semimajor axis [au]
    eb ... binary orbit eccentricity []
    ap ... planetary orbit semimajor axis [au]

    Returns:
    -------
    emax ... max ecc of planetary orbit in S-type binary system
    """
    emax = eforcedS(ab, eb, ap)
    return emax

def epmaxSe0(ab, eb, ap):
    """Maximum eccentricity of initially circular
    planetary orbit in S-type binary system.

    Parameters:
    ----------
    ab ... binary orbit semimajor axes
    eb ... binary orbit eccentricity
    ap ... planetary orbit semimajor axes

    Returns:
    -------
    emax ... maximum eccentricity of planetary orbit
    """
    emax = 2.*eforcedS(ab, eb, ap)
    return emax

def eav2Se0(ab, eb, ap):
    """Average square eccentricity <e^2> of
    planetary orbit in S-type binary system
    for planets on initially circular orbits.

    Parameters:
    ----------
    ab ... binary orbit semimajor axes
    eb ... binary orbit eccentricity
    ap ... planetary orbit semimajor axes

    Returns:
    -------
    eav2 ... average squared eccentricity of planetary orbit
    """
    eav2 = 2.*eforcedS(ab, eb, ap)**2
    return eav2

def eav2S(ab, eb, ap):
    """Average square eccentricity <e^2> of
    planetary orbit in S-type binary system.

    Parameters:
    ----------
    ab ... binary orbit semimajor axes
    eb ... binary orbit eccentricity
    ap ... planetary orbit semimajor axes


    Returns:
    -------
    eav2 ... average squared eccentricity of planetary orbit
    """
    eav2 = eforcedS(ab, eb, ap)**2
    return eav2

def reqpS(ab, eb, ap):
    """Insolation equivalence radius for a planet orbiting
    a single star in a binary star system (Eggl, 2018)

    Parameters:
    ----------
    ab ... binary orbit semimajor axes [au]
    eb ... binary orbit eccentricity
    ap ... planetary orbit semimajor axes [au]

    Returns:
    -------
    reqp ... equivalent radius [au]
    """
    reqp = ap*(1.-eav2S(ab, eb, ap))**(0.25)
    return reqp

def AHZ(LA, teffA, LB, teffB, ab, eb):
    """Averaged Habitable Zone (AHZ) for S-type
     binary star systems (Eggl, 2018).

    Parameters:
    ----------
    LA     ... luminosity of primary star [Lsun]
    teffA  ... effective temperature of primary star [K]
    LB     ... luminosity of secondary star [Lsun]
    teffB  ... effective temperature of secondary star [K]
    ab     ... binary star orbit semimajor axes [au]
    eb     ... binary star orbit eccentricity

    Returns:
    -------
    ahzi   ... inner edge of the AHZ [au]
    ahzo   ... outer edge of the AHZ [au]

    Requires:
    --------
    import numpy as np
    Functions sinner, souter
    """
    #analytic approximation
    AI = LA/seffi(teffA)
    BI = LB/seffi(teffB)

    AO = LA/seffo(teffA)
    BO = LB/seffo(teffB)

#     apI = sqrt(AI)
#     apO = sqrt(AO)

#     reqpSI = reqpS(ab, eb, apI)
#     reqpSO = reqpS(ab, eb, apO)
#     reqbS = reqb(ab, eb)

#     ahzi = AI/reqpSI**2+BI/(reqbS**2-reqpSI**2)
#     ahzo = AO/reqpSO**2+BO/(reqbS**2-reqpSO**2)

    #apI=np.sqrt(AI)
    #apO=np.sqrt(AO)
    
    #reqpSI=reqpS(ab,eb,apI)
    #reqpSO=reqpS(ab,eb,apO)
    #reqb=reqb(ab,eb)
    
    ahzi=sqrt(AI)*(1.+BI/(ab**2*sqrt(1-eb**2)-AI))
    ahzo=sqrt(AO)*(1.+BO/(ab**2*sqrt(1-eb**2)-AO))

    return [ahzi, ahzo]

def PHZ(LA, teffA, LB, teffB, ab, eb):
    """Permanently Habitable Zone (PHZ) for S-type
    binary star systems (Eggl, 2018).

    Parameters:
    ----------
    LA     ... luminosity of primary star [Lsun]
    teffA  ... effective temperature of primary star [K]
    LB     ... luminosity of secondary star [Lsun]
    teffB  ... effective temperature of secondary star [K]
    ab     ... binary star orbit semimajor axes [au]
    eb     ... binary star orbit eccentricity

    Returns:
    -------
    phzi   ... inner edge of the PHZ [au]
    phzo   ... outer edge of the PHZ [au]

    Requires:
    --------
    import numpy as np
    Functions sinner, souter
    """
    AI = LA/seffi(teffA)
    BI = LB/seffi(teffB)

    AO = LA/seffo(teffA)
    BO = LB/seffo(teffB)

    apI = sqrt(AI)
    apO = sqrt(AO)

    epmaxi = epmaxS(ab, eb, apI)
    epmaxo = epmaxS(ab, eb, apO)

    qpI = apI*(1.-epmaxi)
    # qpO = apO*(1.-epmaxo)

    # apopI = apI*(1.+epmaxi)
    apopO = apO*(1.+epmaxo)

    qb = ab*(1.-eb)
    apob = ab*(1.+eb)

    phzi = AI/qpI+BI*qpI/(qpI-qb)**2
    phzo = AO/apopO+BO*apopO/(apopO-apob)**2

    return [phzi, phzo]