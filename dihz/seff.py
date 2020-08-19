#!/bin/python

### Calculate effective insolation values (S_eff) for Habitable Zones. Kopparapu et al. (2014)

############################
# Seff values
############################


# Solar Effective Temperature [K]
teffsun = 5777.

__all__=['seffi','seffo']


def seffi(teff):
    """Calculate effective insolation (S_eff) following
    Kopparapu et al. (2014): Ruaway Greenhouse limit.

    Parameters:
    -----------
    teff...   [K] effective stellar temperature

    Returns:
    -------
    sinner... [] S_eff for the inner Habitable Zone border
    """
    tstar = teff-teffsun
    tstar2 = tstar*tstar
    tstar3 = tstar2*tstar
    tstar4 = tstar3*tstar

    seff0 = 1.107
    a = 1.332e-4
    b = 1.58e-8
    c = -8.308e-12
    d = -1.931e-15

    sinner = seff0+a*tstar + b*tstar2 + c*tstar3+d*tstar4
    return sinner

def seffo(teff):
    """Calculate effective insolation (S_eff) following
    Kopparapu et al. (2014) Maximum Greenhouse limit.

    Parameters:
    -----------
    teff...   [K] effective stellar temperature
    Returns:
    -------
    souter... S_eff for the inner Habitable Zone border
    """
    tstar = teff-teffsun
    tstar2 = tstar*tstar
    tstar3 = tstar2*tstar
    tstar4 = tstar3*tstar

    # Kopparapu et al 2014 Maximum Greenhouse
    seff0 = 0.356
    a = 6.171e-5
    b = 1.698e-9
    c = -3.198e-12
    d = -5.575e-16
    souter = seff0 + a*tstar + b*tstar2 + c*tstar3+d*tstar4
    return souter
