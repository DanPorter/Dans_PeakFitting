####################################################################################################################
#########################################  ESTIMATION FUNCTIONS  ###################################################
####################################################################################################################

# Parameters
_max_height_factor = 3  # times the max height


def guess_height(y):
    return np.max(y) - np.min(y)


def max_height(y):
    return _max_height_factor * np.max(y)


def min_height(y):
    return np.min(y)


def guess_bkg(y):
    return np.min(y)


def max_bkg(y):
    return np.max(y)


def min_bkg(y):
    return -np.inf


def guess_centre(x, y):
    return x[np.argmax(y)]


def max_centre(x):
    return 2 * np.max(x) - np.mean(x)


def min_centre(x):
    return 2 * np.min(x) - np.mean(x)


def guess_fwhm(x, y, interpolate=False):
    """
    Calculate a simple FWHM from a peak
    """

    if interpolate:
        interx = np.linspace(x[0], x[-1], len(x) * 100)
        intery = np.interp(interx, x, y)
        x, y = interx, intery

    mx = max(y)
    ln = len(y)

    # Peak position
    pkpos = y.argmax()

    # Split into two parts - before and after the peak
    hfxx1 = x[:pkpos + 1]
    hfxx2 = x[pkpos:]

    # Find the half-max positions
    hfmx1 = abs(y[:pkpos + 1] - mx // 2)
    hfmx2 = abs(y[pkpos:] - mx // 2)

    hfpos1 = hfxx1[hfmx1.argmin()]
    hfpos2 = hfxx2[hfmx2.argmin()]

    # Return FWHM
    return abs(hfpos2 - hfpos1)


def max_fwhm(x):
    return 2 * (2 * np.max(x) - np.min(x))


def min_fwhm(x):
    return 0.5 * np.min(np.diff(x))


def guess_slope(x, y):
    return (y[-1] - y[0]) / (x[-1] - x[0])