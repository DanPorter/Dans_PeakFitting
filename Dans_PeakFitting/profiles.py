####################################################################################################################
#########################################  FITTING FUNCTIONS  ######################################################
####################################################################################################################


def default_function(x):
    return x


def background(x, bkg=0):
    return bkg * np.ones(len(x))


def straightline(x, gradient=1.0, intercept=0.0):
    """
    Straight line
    """
    return gradient * x + intercept


def slope(x, gradient=1.0, bkg=0):
    """
    sloping background
    """
    centre = np.mean(x)
    return gradient * (x - centre) + bkg


def step(x, height=0, centre=0, bkg=0):
    y = np.zeros(len(x))
    y[x <= centre] = bkg
    y[x > centre] = bkg + height
    return y


def simple(x, height=1, centre=0, fwhm=0.5):
    "Plot an Illustration of simpfit"

    minpos = centre - fwhm
    maxpos = centre + fwhm
    y = np.zeros(len(x))
    y[len(x) // 5:-len(x) // 5] += height / 2
    y[np.logical_and(x > minpos, x < maxpos)] += height / 2
    return y


def gauss(x, height=1, centre=0, fwhm=1):
    return height * np.exp(-np.log(2) * ((x - centre) / (fwhm / 2)) ** 2)


def lorentz(x, height=1, centre=0, fwhm=0.5):
    return height / (1 + ((x - centre) / (fwhm / 2)) ** 2)


def pvoight(x, height=1, centre=0, fwhm=0.5, lorfrac=0.5):
    hwhm = fwhm / 2.0
    ln2 = 0.69314718055994529
    pos = x - centre
    L = lorfrac / (1 + (pos / hwhm) ** 2)
    G = (1 - lorfrac) * np.exp(-ln2 * (pos / hwhm) ** 2)
    return height * (G + L)


def orderparameter(x, Tc=100, beta=0.5, amp=1):
    """
    Generate an order parameter
    """
    # op = amp*np.real(np.power(np.complex(Tc-x),beta))
    op = amp * np.power(Tc - x, beta)
    op[np.isnan(op)] = 0.0
    return op