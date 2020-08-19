import numpy as np
import matplotlib.pyplot as plt
import circumbinary
import circumstellar
import stability

__all__=['circumbinaryhz2D','circumstellarhz2D']

### Plot Dynamically Informed Habitable Zones

sqrt=np.sqrt
linspace=np.linspace
meshgrid=np.meshgrid

def circumbinaryhz2D(LA,LB,teffA,teffB,mA,mB,ab,eb,xmin=-4,xmax=4,ymin=-4,ymax=4,title=''):
    """Plot the circumbinary dynamically informed habitable zones. 
    
    Parameters:
    ----------
    LA     ... luminosity of primary star [Lsun]
    LB     ... luminosity of secondary star [Lsun]
    teffA  ... effective temperature of primary star [K]
    teffB  ... effective temperature of secondary star [K]
    mA     ... mass of primary star [Msun]
    mB     ... mass of secondary star [Msun]
    ab     ... binary star orbit semimajor axes [au]
    eb     ... binary star orbit eccentricity

    Returns:
    -------
    matplotlib pyplot plot 

    """
    xs = linspace(xmin,xmax, 201)
    ys = linspace(ymin,ymax, 201)

    # mesh for contours
    xv,yv = meshgrid(xs,ys)

    # generate the level map
    r = sqrt(xv**2 + yv**2)

    [phzi,phzo]=circumbinary.PHZ(LA,LB,teffA,teffB,mA,mB,ab,eb)
    [ahzi,ahzo]=circumbinary.AHZ(LA,LB,teffA,teffB,mA,mB,ab,eb)

    astab=hw99P(mA,mB,ab,eb)

    plt.figure(figsize=(4,4),dpi=150)

    # plot the contours with two levels only
    # notice the xv, yv parameters
    plt.title(title)
    plt.contourf(xs,ys,r, levels=[0,astab], colors=('#bf00ff'),alpha=0.4,hatch='//')
    plt.contour(xs,ys,r, levels=[astab], colors=('#ac10e0'),linestyles='dashed')

    plt.contourf(xs,ys,r, levels=[ahzi,ahzo], colors=('#EEA700'))
    plt.contourf(xs, ys, r, levels=[phzi,phzo], colors=('blue'))

    # plot the two circles
    plt.contour(xs,ys,r, levels=[ahzi], colors=('k'),linewidths=0.6)
    plt.contour(xs,ys,r, levels=[ahzo], colors=('k'),linewidths=0.6)
    plt.contour(xs,ys,r, levels=[phzi], colors=('k'),linewidths=0.6)
    plt.contour(xs,ys,r, levels=[phzo], colors=('k'),linewidths=0.6)


    plt.text(0,0,'AB',horizontalalignment='center',verticalalignment='center')

    plt.text(xmin+0.5,ymax-0.3,'Averaged Habitable Zone',horizontalalignment='left',verticalalignment='center',color='#EEA700')
    plt.text(xmin+0.5,ymax-0.9,'Permanently Habitable Zone',horizontalalignment='left',verticalalignment='center',color='b')
    plt.text(xmin+0.5,ymin+0.5,'Unstable Orbits',horizontalalignment='left',verticalalignment='center',color='#bf00ff')
    
    plt.xlabel('x [au]')
    plt.ylabel('y [au]')

    plt.show()

    
def circumstellarhz2D(LA,LB,teffA,teffB,mA,mB,ab,eb,xmin=-4,xmax=4,ymin=-4,ymax=4):
    """Plot the circumstellar dynamically informed habitable zones. 
    
    Parameters:
    ----------
    LA     ... luminosity of primary star [Lsun]
    LB     ... luminosity of secondary star [Lsun]
    teffA  ... effective temperature of primary star [K]
    teffB  ... effective temperature of secondary star [K]
    mA     ... mass of primary star [Msun]
    mB     ... mass of secondary star [Msun]
    ab     ... binary star orbit semimajor axes [au]
    eb     ... binary star orbit eccentricity

    Returns:
    -------
    matplotlib pyplot plot 

    """
    
    xs = linspace(xmin,xmax, 201)
    ys = linspace(ymin,ymax, 201)

# mesh for contours
    xv,yv = meshgrid(xs,ys)

# generate the level map
    r = sqrt(xv**2 + yv**2)

    [phzi,phzo]=circumstellar.PHZ(LA,LB,teffA,teffB,ab,eb)
    [ahzi,ahzo]=circumstellar.AHZ(LA,LB,teffA,teffB,ab,eb)


    astab=stability.hw99S(m0,m1,ab,eb)

    plt.figure(figsize=(4,4),dpi=150)

# plot the contours with two levels only
# notice the xv, yv parameters

    plt.contourf(xs,ys,r, levels=[0,astab], colors=('g'),alpha=0.4)
    plt.contour(xs,ys,r, levels=[astab], colors=('g'),linestyles='dashed')

    plt.contourf(xs,ys,r, levels=[ahzi,ahzo], colors=('#EEA700'))
    plt.contourf(xs, ys, r, levels=[phzi,phzo], colors=('blue'))

    plt.contour(xs,ys,r, levels=[ahzi], colors=('k'),linewidths=0.6)
    plt.contour(xs,ys,r, levels=[ahzo], colors=('k'),linewidths=0.6)
    plt.contour(xs,ys,r, levels=[phzi], colors=('k'),linewidths=0.6)
    plt.contour(xs,ys,r, levels=[phzo], colors=('k'),linewidths=0.6)

    plt.text(0,0,'A',horizontalalignment='center',verticalalignment='center')

    plt.text(xmin+0.5,ymax-0.3,'Averaged Habitable Zone',horizontalalignment='left',verticalalignment='center',color='#EEA700')
    plt.text(xmin+0.5,ymax-0.9,'Permanently Habitable Zone',horizontalalignment='left',verticalalignment='center',color='b')
    plt.text(xmin+0.5,ymin+0.5,'Stable Orbits',horizontalalignment='left',verticalalignment='center',color='g')

    plt.xlabel('x [au]')
    plt.ylabel('y [au]')
    plt.show()
    