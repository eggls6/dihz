#!/bin/python
import numpy as np
from scipy.optimize import fsolve

from .seff import *
import circumstellar
import circumbinary 


sqrt = np.sqrt

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

# circumstellar

def ahziS(LA, teffA, LB, teffB, ab, eb, apI):
    AI = LA/seffi(teffA)
    BI = LB/seffi(teffB)

    reqpSI = circumstellar.reqpS(ab, eb, apI)
    reqbI = reqb(ab, eb)

    lhsm1 = AI/reqpSI**2+BI/(reqbI**2-reqpSI**2)-1

    return lhsm1

def ahzoS(LA, teffA, LB, teffB, ab, eb, apO):

    AO = LA/seffo(teffA)
    BO = LB/seffo(teffB)

    reqpSO = circumstellar.reqpS(ab, eb, apO)
    reqbO = reqb(ab, eb)

    lhsm1 = AO/reqpSO**2+BO/(reqbO**2-reqpSO**2)-1

    return lhsm1

def phziS(LA, teffA, LB, teffB, ab, eb, apI):
    AI = LA/seffi(teffA)
    BI = LB/seffo(teffB)

    epmaxi = circumstellar.epmaxS(ab, eb, apI)

    qpI = apI*(1.-epmaxi)
#       apopI = apI*(1.+epmaxi)
    qb = ab*(1.-eb)

    lhsm1 = AI/qpI**2+BI/(qpI-qb)**2-1.

    return lhsm1
###

def phzoS(LA, teffA, LB, teffB, ab, eb, apO):

    AO = LA/seffo(teffA)
    BO = LB/seffo(teffB)

    epmaxo = circumstellar.epmaxS(ab, eb, apO)
#        qpO = apO*(1.-epmaxo)
    apopO = apO*(1.+epmaxo)
    apob = ab*(1.+eb)
    lhsm1 = AO/apopO**2+BO/(apopO-apob)**2-1.

    return lhsm1

def PHZ_A(LA, teffA, LB, teffB, ab, eb):

    [apI_initial_guess, apO_initial_guess] = SSHZ(LA, teffA)

    funci = lambda apI: phziS(LA, teffA, LB, teffB, ab, eb, apI)
    funco = lambda apO: phzoS(LA, teffA, LB, teffB, ab, eb, apO)

    phzi = fsolve(funci, apI_initial_guess)
    phzo = fsolve(funco, apO_initial_guess)

    if(phzi > phzo or phzi < 0 or phzo < 0):
        phzi = 0
        phzo = 0
        print("Error in class: semianalytic: no PHZ for given parameters")

    return [phzi, phzo]


##########################################################
# Circumbinary
##########################################################
def ahziP(LA, teffA, mA, LB, teffB, mB, ab, eb, ap):

    mu = mB/(mA+mB)

    AI = LA/seffi(teffA)
    BI = LB/seffi(teffB)

#       apI=np.sqrt(AI+BI)
    apI = ap

    reqpPI = circumbinary.reqpP(mA, mB, ab, eb, apI)
    reqbP = reqb(ab, eb)

    reqA = reqbP*mu
    reqB = reqbP*(1-mu)

    ahzi = AI/(reqpPI**2-(mu*reqA)**2)+BI/(reqpPI**2+((1-mu)*reqB)**2)-1

    return ahzi

def ahzoP(LA, teffA, mA, LB, teffB, mB, ab, eb, ap):

    mu = mB/(mA+mB)

    AO = LA/seffo(teffA)
    BO = LB/seffo(teffB)

#       apO=np.sqrt(AO+BO)
    apO = ap

    reqpPO = circumbinary.reqpP(mA, mB, ab, eb, apO)
    reqbP = reqb(ab, eb)

    reqA = reqbP*mu
    reqB = reqbP*(1-mu)

    ahzo = AO/(reqpPO**2-(mu*reqA)**2)+BO/(reqpPO**2+((1-mu)*reqB)**2)-1

    return ahzo

def phziP(LA, teffA, mA, LB, teffB, mB, ab, eb, ap):

    mu = mB/(mA+mB)

    AI = LA/seffi(teffA)
    BI = LB/seffi(teffB)

#       apI=np.sqrt(AI+BI)
    apI = ap
    epmaxi = circumbinary.epmaxP(ab, eb, apI)
    qpI = apI*(1.-epmaxi)
#        apopI = apI*(1.+epmaxi)
#        qb = ab*(1.-eb)
    apob = ab*(1.+eb)

    phzi = AI/(qpI-mu*apob)**2+BI/(qpI+(1-mu)*apob)**2-1

    return phzi

def phzoP(LA, teffA, mA, LB, teffB, mB, ab, eb, ap):

    mu = mB/(mA+mB)

    AO = LA/seffo(teffA)
    BO = LB/seffo(teffB)

#       apO=np.sqrt(AO+BO)
    apO = ap
    epmaxo = circumbinary.epmaxP(ab, eb, apO)
#        qpO = apO*(1.-epmaxo)
    apopO = apO*(1.+epmaxo)
#        qb = ab*(1.-eb)
    apob = ab*(1.+eb)

    phzo = AO/(apopO+mu*apob)**2+BO/(apopO-(1-mu)*apob)**2-1

    return phzo
