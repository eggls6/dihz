#!/bin/python
import numpy as np

__all__=['stabilityLimit','hw99S','hw99P']

###########################
# Orbital stability
############################
def stabilityLimit(binary_star_type, mA, mB, ab, eb):
    """Routines for caculating dynamical stability
    for Earth-like planets in binary star systems following
    Holman & Wiegert (1999).

    Parameters:
    -----------
    mA... mass of primary star [Msun]
    mB... mass of secondary star [Msun]
    ab... semimajor axis of binary orbit [au]
    eb... orbital eccentricity of binary

    Calculates:
    -----------
    stability_limit  ... maximum (cicumstellar) or minimum (circumbinary)
                         stable distance of planet on circular orbit from
                         host star(s)
    """

    tp=binary_star_type

    if(ab <= 0):
        raise ValueError('Semimajor axis ab must be >= 0!')

    if(eb < 1 and eb >= 0):
        raise ValueError('Eccentricity eb must be < 1!')

    if(tp in ['S', 's', 'S-type', 'S-Type', 'circumstellar']):
           stability_limit = hw99S(mA, mB, ab, eb)

    elif(tp in ['P', 'p', 'P-type', 'P-Type', 'circumbinary']):
           stability_limit = hw99P(mA, mB, ab, eb)
    else:
        raise ValueError('Binary star type not recognized. \
                              Choose "S" or "P".')

def hw99S(mA, mB, ab, eb):
        """ Circumstellar (S-type) stability limit for binary star systems
        according to Holman & Wiegert (1999)

        Parameters:
        -----------
        mA... mass of primary star [Msun]
        mB... mass of secondary star [Msun]
        ab... semimajor axis of binary orbit [au]
        eb... orbital eccentricity of binary

        Returns:
        -----------
        ap... maximum stable distance of planet on circular orbit from
              host star

        """
        mu = mB/(mA+mB)
        eb2 = eb*eb
        ap = ab*(0.464-0.38*mu-0.631*eb+0.586*mu*eb + 0.15*eb2-0.198*mu*eb2)
        return ap

def hw99P(mA, mB, ab, eb):
        """ Circumbinary (P-type) stability limit for binary star systems
        according to Holman & Wiegert (1999)

        Parameters:
        ----------
        mA... mass of primary star [Msun]
        mB... mass of secondary star [Msun]
        ab... semimajor axis of binary orbit [au]
        eb... orbital eccentricity of binary

        Returns:
        -------
        ap... maximum stable distance of planet on circular orbit from
              host star

        """
        mu = mB/(mA+mB)
        eb2 = eb*eb
        mu2 = mu*mu
        c = [0, 1.6, +5.1, -2.22, 4.12, -4.27, -5.09, 4.61]
        ap = ab*(c[1]+c[2]*eb+c[3]*eb2+c[4]*mu+c[5]*eb*mu +
                 c[6]*mu2+c[7]*eb2*mu2)
        return ap