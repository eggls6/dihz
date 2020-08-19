#!/bin/python
import numpy as np
from .seff import *

################################
# Circumbinary (P-type systems)
###############################

sqrt=np.sqrt

__all__=['reqb','reqpP','eforcedP','epmaxPe0','epmaxP','eav2P','eav2Pe0','PHZ','AHZ']

#######################################
# Insolation equivalent orbit distance
######################################


def eforcedP(ab, eb, ap, mA, mB):
    """Forced eccentricity of planetary orbits in P-type binary system.

    Parameters:
    ----------
    ab ... binary orbit semimajor axis [au]
    eb ... binary orbit eccentricity []
    ap ... planetary orbit semimajor axis [au]
    mA ... mass of primary star [Msun]
    mB ... mass of secondary star [Msun]

    Returns:
    -------
    emax ... forced eccentricity of planetary orbit in S-type binary system
    """
    mu = (mB/(mA+mB))
    emax = 5./2.*ab/ap*(1-2*mu)*(4.*eb+3*eb**3)/(4.+6.*eb*eb)
    return emax

def epmaxPe0(ab, eb, ap, mA, mB):
    """Maximum eccentricity of initially circular
    planetary orbit in P-type binary system.

    Parameters:
    ----------
    ab ... binary orbit semimajor axes
    eb ... binary orbit eccentricity
    ap ... planetary orbit semimajor axes
    mA ... mass of primary star [Msun]
    mB ... mass of secondary star [Msun]

    Returns:
    -------
    emax ... maximum eccentricity of planetary orbit
    """
    # max ecc of planetary orbit in S-type binary system
    emax = 2.*eforcedP(ab, eb, ap, mA, mB)
    return emax

def epmaxP(ab, eb, ap, mA, mB):
    """Maximum eccentricity of
    planetary orbit in P-type binary system.

    Parameters:
    ----------
    ab ... binary orbit semimajor axes
    eb ... binary orbit eccentricity
    ap ... planetary orbit semimajor axes
    mA ... mass of primary star [Msun]
    mB ... mass of secondary star [Msun]

    Returns:
    -------
    emax ... maximum eccentricity of planetary orbit
    """

    emax = eforcedP(ab, eb, ap, mA, mB)
    return emax

def eav2P(ab, eb, ap, mA, mB):
    """Average square eccentricity <e^2> of
    planetary orbit in P-type binary system.

    Parameters:
    ----------
    ab ... binary orbit semimajor axes
    eb ... binary orbit eccentricity
    ap ... planetary orbit semimajor axes
    mA ... mass of primary star [Msun]
    mB ... mass of secondary star [Msun]

    Returns:
    -------
    eav2 ... average squared eccentricity of planetary orbit
    """
    efp = eforcedP(ab, eb, ap, mA, mB)
    eav2 = efp**2
    return eav2

def eav2Pe0(ab, eb, ap, mA, mB):
    """Average square eccentricity <e^2> of
    planetary orbit in P-type binary system
    for planets on initially circular orbits.

    Parameters:
    ----------
    ab ... binary orbit semimajor axes
    eb ... binary orbit eccentricity
    ap ... planetary orbit semimajor axes
    mA ... mass of primary star [Msun]
    mB ... mass of secondary star [Msun]

    Returns:
    -------
    eav2 ... average squared eccentricity of planetary orbit
    """
    efp = circumbinary.forcedP(ab, eb, ap, mA, mB)
    eav2 = 2.*efp**2
    return eav2

def reqb(a, e):
    """Insolation equivalent orbit distance

    Parameters:
    ----------
    a ... semimajor axis
    e ... orbital eccentricity
    """
    reqp = a*(1.-e*e)**(0.25)
    return reqp         
         
def reqpP(mA, mB, ab, eb, ap):
    """Insolation equivalence radius for a planet orbiting
     a binary star (Eggl, 2018).

    Parameters:
    ----------
    mA ... mass of primary star [Msun]
    mB ... mass of secondary star [Msun]
    ab ... binary orbit semimajor axes [au]
    eb ... binary orbit eccentricity
    ap ... planetary orbit semimajor axes [au]

    Returns:
    -------
    reqp ... equivalent radius [au]
    """
    reqp = ap*(1-circumbinary.eav2P(ab, eb, ap, mA, mB))**0.25
    return reqp

def PHZ(LA, teffA, mA, LB, teffB, mB, ab, eb):
    """Permanently Habitable Zone (PHZ) for P-type
    binary star systems (Eggl, 2018).

    Parameters:
    ----------
    LA     ... luminosity of primary star [Lsun]
    teffA  ... effective temperature of primary star [K]
    mA     ... mass of primary star [Msun]
    LB     ... luminosity of secondary star [Lsun]
    teffB  ... effective temperature of secondary star [K]
    mB     ... mass of secondary star [Msun]
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
    mu = mB/(mA+mB)

    AI = LA/seffi(teffA)
    BI = LB/seffi(teffB)

    AO = LA/seffo(teffA)
    BO = LB/seffo(teffB)

    qpI = sqrt(AI+BI)
    qpO = sqrt(AO+BO)

    epmaxi = epmaxP(ab, eb, qpI, mA, mB)
    epmaxo = epmaxP(ab, eb, qpO, mA, mB)

    apI=qpI/(1.-epmaxi)
    apO=qpO/(1.+epmaxo)
    
    qb=ab*(1.-eb)
    apob=ab*(1.+eb)
    
    a=qpI-mu*apob
    b=qpI+(1-mu)*apob
    phzi=(AI/a**2+BI/b**2)*apI
        
    a=qpO+mu*apob
    b=qpO-(1-mu)*apob
    phzo=(AO/a**2+BO/b**2)*apO

    return [phzi, phzo]

def AHZ(LA, teffA, mA, LB, teffB,  mB, ab, eb):
    """Averaged Habitable Zone (AHZ) for S-type
     binary star systems (Eggl, 2018).

    Parameters:
    ----------
    LA     ... luminosity of primary star [Lsun]
    teffA  ... effective temperature of primary star [K]
    mA     ... mass of primary star [Msun]
    LB     ... luminosity of secondary star [Lsun]
    teffB  ... effective temperature of secondary star [K]
    mB     ... mass of secondary star [Msun]
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
    mu = mB/(mA+mB)

    AI = LA/seffi(teffA)
    BI = LB/seffi(teffB)

    AO = LA/seffo(teffA)
    BO = LB/seffo(teffB)

    apI = sqrt(AI+BI)
    apO = sqrt(AO+BO)

    reqbP = reqb(ab, eb)

    reqA = reqbP*mu
    reqB = reqbP*(1.-mu)

    ahzi = AI*apI/(apI**2-(mu*reqA)**2)+BI*apI/(apI**2+((1-mu)*reqB)**2)
    ahzo = AO*apO/(apO**2-(mu*reqA)**2)+BO*apO/(apO**2+((1-mu)*reqB)**2)

    return [ahzi, ahzo]