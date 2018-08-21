

def gauss_area(height, centre, fwhm, dheight=0, dcentre=0, dfwhm=0):
    """

    :param height:
    :param centre:
    :param fwhm:
    :param dheight:
    :param dcentre:
    :param dfwhm:
    :return:
    """
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    dsigma = dfwhm/(2*np.sqrt(2*np.log(2)))
    area = np.abs(height*sigma*np.log(2))
    darea = area*np.sqrt( (dheight/height)**2 + ( dsigma/sigma)**2 )
    return area, darea


def lorentz_area(height, centre, fwhm, dheight=0, dcentre=0, dfwhm=0):
    """

    :param height:
    :param centre:
    :param fwhm:
    :param dheight:
    :param dcentre:
    :param dfwhm:
    :return:
    """
    area = np.abs(np.pi*height*fwhm/2)
    darea = area*np.sqrt( (dheight/height)**2 + ( dfwhm/fwhm)**2 )
    return area, darea


def pvoight_area(height, centre, fwhm, lorfrac, dheight=0, dcentre=0, dfwhm=0, dlorfrac=0):
    """

    :param height:
    :param centre:
    :param fwhm:
    :param lorfrac:
    :param dheight:
    :param dcentre:
    :param dfwhm:
    :param dlorfrac:
    :return:
    """

    # Calculated Voight area = Gaussian + Voight
    sigma = fwhm/(2*np.sqrt(2*np.log(2))) # Gaussian sigma
    dsigma = dfwhm/((2*np.sqrt(2*np.log(2))))
    Garea = np.abs(height*sigma*np.sqrt(2*np.pi))
    Larea = np.pi*height*fwhm/2
    area = lorfrac*Larea + (1-lorfrac)*Garea

    # Error on area
    dGarea = Garea*np.sqrt( (dheight/height)**2 + (dsigma/sigma)**2 )
    dLarea = Larea*np.sqrt( (dheight/height)**2 + (dfwhm/fwhm)**2 )
    dVarea1= (1-lorfrac)*Garea*np.sqrt( (dlorfrac/(1-lorfrac))**2 + (dGarea/Garea)**2 )
    dVarea2= lorfrac*Larea*np.sqrt( (dlorfrac/lorfrac)**2 + (dLarea/Larea)**2 )
    darea = np.sqrt( dVarea1**2 + dVarea2**2 )
    return area, darea