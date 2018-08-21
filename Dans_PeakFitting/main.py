
import numpy as np
from scipy.optimize import curve_fit # Peak fitting
import inspect
import copy
import DansGeneralProgs as dgp


####################################################################################################################
###############################################  PROFILE CLASS  ####################################################
####################################################################################################################


class Profile:
    """
    Multi-function profile generator, adds together multiple functions into a single profile
    Useful for generating profiles with multiple peaks
    Allows certain parameters to be fixed

    Use:
        fn = profile(function, nfuncs, *args, **kwargs)

        function = a profile function with inputs like gauss(x, a,b,c) where a,b,c are the function parameters
        nfuncs = number of times this function will be repeated
        args = additinal values as default values, in order of input functions
        kwargs = additinal values as default values, in order of input functions
        fn = profile function can be called using all non-fixed input arguments

    Note that the kwargs keywords are based on the input keywords of function, appended
    with a number based on the number of functions stored in profile

    Note that the instantiated profile function, fn, has an __add__ (+) property allowing multiple
    profile functions based on different underlying functions be be joined together. The numeric sequence
    of the inputs will be based on the order of the functions in each profile.

    E.G.:

        x = np.arange(-5,5,0.1)
        profile1 = profile(gauss) # single gaussian
        profile2 = profile(gauss,2) # two gaussian functions
        profile3 = profile(background) # simple background function
        multiprofile = profile1+profile2+profile3

        y = multiprofile(x,
                height1 = 1,
                centre1 = 0,
                fwhm1 = 0.1,
                height2 = 1,
                centre2 =-1,
                fwhm2 = 0.1,
                height3 = 1,
                centre3 = 1,
                fwhm3 = 0.1,
                bkg4 = 0)
        or:
        y = multiprofile(1,0,0.1,1,-1,0.1,1,1,0.1,0)

    Fixed Parameters:
    any arguments given as inputs to profile will be fixed and
    do not need to be called by the resultant function, this allows the function
    to be called by curve_fit:

        profile1 = profile(gauss, centre1=0.3)
        profile1.fix('centre1')
        y = profile1(x, height1=3, fwhm1=0.5)
    or:
        y = profile1(x, 3, 0.5)

    Fitting:
    curve_fit cannot automatically determine the number of input arguments so default
    values must be expressed using the p0 parameter:

        from scipy.optimize import curve_fit
        profile1 = profile(gauss, centre1=0.3)
        popt, pcov = curve_fit(profile1, x, y, p0=[3,0.5])
        yfit = profile1(x, *popt)

    limits:
    Limits can be set so any argument values outside the limit range will instead use a default value

        profile1 = profile(gauss)
        profile1.limit('height1', high_limit=100, low_limit=0, default=1.0)
        y = profile1(x, 100, 0, 0.5)
        y.max() >> 1.0
    """
    def __init__(self, func, nfuncs=1, *args, **kwargs):
        """
        Build lists of argments and functions
        """
        # Generate function
        self.functions = [func ] *nfuncs

        # get the function input arguments
        arg_names, varargs, varkw, defaults = inspect.getargspec(func)

        self.arguments = []
        self.errors = []
        self.input_arguments = []
        self.argument_map = {}
        arg_val = 0
        n_inputargs = len(args)
        inputkwargs = kwargs.keys()
        for n in range(nfuncs):
            profile_arguments = {}
            profile_errors = {}
            for m in range(1 ,len(arg_names)):
                name = arg_names[m]
                newname = nam e +str( n +1)
                self.input_arguments += [newname]
                self.argument_map[newname] = [name, n]

                # Assign argument values
                if arg_val < n_inputargs:
                    profile_arguments[name] = args[arg_val]
                    arg_val += 1
                elif newname in inputkwargs:
                    profile_arguments[name] = kwargs[newname]
                else:
                    profile_arguments[name] = 0
                profile_errors[name] = 0
            self.arguments += [profile_arguments]
            self.errors += [profile_errors]

        # fixed arguments
        self.fixed = []

        # Limits and defaults
        # limits[input_argument] = [high, low, default]
        self.limits = {}
        self.high_limits = {}
        self.low_limits = {}
        self.defaults = {}

    def get(self, arg_name):
        """
        Returns argument value, using numeric argument name
        """
        name, funcno = self.argument_map[arg_name]
        return self.arguments[funcno][name]

    def get_error(self, arg_name):
        """
        Returns argument value, using numeric argument name
        """
        name, funcno = self.argument_map[arg_name]
        return self.errors[funcno][name]

    def set(self, arg_name, value):
        """
        Set the current value of an argument
        """
        name, funcno = self.argument_map[arg_name]
        self.arguments[funcno][name] = value

    def set_error(self, arg_name, value):
        """
        Set the current error value of an argument
        """
        name, funcno = self.argument_map[arg_name]
        self.errors[funcno][name] = value

    def fix(self, arg_names):
        """
        Fix parameter arg_name
        """
        arg_names = list(np.asarray(arg_names).reshape(-1))
        for arg in arg_names:
            if arg in self.input_arguments:
                self.fixed += [arg]
                self.input_arguments.remove(arg)
            else:
                print('%s is not an input argument ' %arg)

    def free(self, arg_names):
        """
        Free parameter arg_name
        """
        arg_names = list(np.asarray(arg_names).reshape(-1))
        for arg in arg_names:
            if arg in self.fixed:
                self.input_arguments += [arg]
                self.fixed.remove(arg)
            else:
                print('%s is not an input argument ' %arg)

    def limit(self, arg_name, high_limit, low_limit, default=None):
        """
        Set limits on input arguments, such that inputting a value outside the limits
        will be reset to a default value
        """
        if default is None:
            # use current value
            default = self.get(arg_name)
        self.limits[arg_name] = [high_limit, low_limit, default]
        self.high_limits[arg_name] = high_limit
        self.low_limits[arg_name] = low_limit
        self.defaults[arg_name] = default

    def applylimits(self):
        """
        Apply any provided limits to the function parameters
        """
        for arg in self.high_limits.keys():
            high_lim = self.high_limits[arg]
            low_lim = self.low_limits[arg]
            def_val = self.defaults[arg]
            name, funcno = self.argument_map[arg]
            val = self.arguments[funcno][name]
            if val > high_lim:
                self.arguments[funcno][name] = def_val
            elif val < low_lim:
                self.arguments[funcno][name] = def_val

    def applydefaults(self):
        """
        Resents arguments to their defaults
        """
        for arg in self.defaults.keys():
            self.set(arg, self.defaults[arg])

    def ninputs(self):
        """
        Returns the number of non-fixed input arguments
        """
        return len(self.input_arguments)

    def inputstring(self):
        """
        Returns string of function inputs
        """
        input_str = ', '.join(['%s=%s ' %(arg ,self.get(arg)) for arg in self.input_arguments])
        return 'fn(x, %s) ' %input_str

    def values(self):
        """
        Returns the current values stored in profile
        """
        return [self.arguments[ self.argument_map[arg][1] ][ self.argument_map[arg][0] ] for arg in self.input_arguments]

    def function_value_errors(self):
        """
        Returns the values and errors of each function as
          [[value1, value2, error1, error2], [value1, value2, error1, error2]]
        """
        output = []
        for n in range(len(self.functions)):
            fun_output = []
            args = self.arguments[n]
            errs = self.errors[n]
            for k in args.keys():
                fun_output += [args[k]]
            for k in errs.keys():
                fun_output += [errs[k]]
            output += [fun_output]
        return output

    def check(self):
        """
        Print profile argument
        """
        print('Profiles:')
        for n in range(len(self.functions)):
            print('%2d %s ' %( n +1 ,self.functions[n].__name__))
            args = self.arguments[n]
            errs = self.errors[n]
            for k in args.keys():
                arg = '%s%d ' %(k , n +1)
                out = '%15s : %s ' %(arg, dgp.stfm(args[k], errs[k]))
                if arg not in self.input_arguments:
                    out += ' (fixed)'
                if arg in self.high_limits.keys():
                    high_lim = self.high_limits[arg]
                    low_lim = self.low_limits[arg]
                    def_val = self.defaults[arg]
                    out += '   high_limit=%s, low_limit=%s, default=%s ' %(high_lim, low_lim, def_val)
                print(out)

    def fit(self, x, y, dy=None):
        """
        Fit data
        """
        # Set dy to 1 if not given
        if dy is None: d y =np.ones(len(y))

        # Remove zeros from x - causes errors in covariance matrix
        xold = 1. 0 *x
        offset = 0.
        if any(np.abs(x ) <0.0001):
            print( 'Zero detected - adding 0.0001 to x values' )
            offset = 0.0001
            x = x + offset
        if any(np.isnan(dy)):
            print( 'Ignoring errors due to NaNs' )
            d y =np.ones(len(y))

        # Handle zero intensities
        y[ y< 0.01] = y[y < 0.01] + 0.01
        dy[dy < 0.01] = dy[dy < 0.01] + 0.01

        # weights
        weights = 1 / dy ** 2

        # Estimate starting parameters
        ini_parameters = self.values()

        # Initial Fit
        try:
            fitvals, covmat = curve_fit(self, x, y, p0=ini_parameters, sigma=dy, absolute_sigma=True)
        except RuntimeError:
            print('fit failed!')
            fitvals = 1 * ini_parameters
            covmat = 0 * np.eye(len(ini_parameters))
        yfit = self.__call__(xold, *fitvals)  # New curve
        self.chi = np.sum((y - yfit) ** 2 / dy)  # Calculate CHI^2

        # Add errors
        errvals = np.sqrt(np.diag(covmat))
        for n, arg in enumerate(self.input_arguments):
            self.set(arg, fitvals[n])
            self.set_error(arg, errvals[n])
        return yfit

    def __add__(self, newprofile):
        """
        Add additional profile
        """
        new = copy.deepcopy(self)
        # copy arguments without reference
        new.functions += copy.deepcopy(newprofile.functions)
        new.arguments += copy.deepcopy(newprofile.arguments)
        new.errors += copy.deepcopy(newprofile.errors)

        # add new profile function to arguments
        # this requires assigning new numbers to each argument
        nprofiles = len(self.functions)
        for name in newprofile.input_arguments:
            name_str = ''.join([i for i in name if not i.isdigit()])
            name_num = int(''.join([i for i in name if i.isdigit()]))
            new_name = name_str + str(name_num + nprofiles)
            new.input_arguments += [new_name]

        new_keys = newprofile.argument_map.keys()
        for name in new_keys:
            name_str = ''.join([i for i in name if not i.isdigit()])
            name_num = int(''.join([i for i in name if i.isdigit()]))
            new_name = name_str + str(name_num + nprofiles)
            new.argument_map[new_name] = copy.deepcopy(newprofile.argument_map[name])
            new.argument_map[new_name][1] += nprofiles

        # copy across fixed parameters
        for name in newprofile.fixed:
            name_str = ''.join([i for i in name if not i.isdigit()])
            name_num = int(''.join([i for i in name if i.isdigit()]))
            new_name = name_str + str(name_num + nprofiles)
            new.fixed += [new_name]

        # copy across limits and defaults
        for arg in newprofile.limits.keys():
            high_lim = copy.deepcopy(newprofile.limits[arg][0])
            low_lim = copy.deepcopy(newprofile.limits[arg][1])
            def_val = copy.deepcopy(newprofile.limits[arg][2])
            name_str = ''.join([i for i in arg if not i.isdigit()])
            name_num = int(''.join([i for i in arg if i.isdigit()]))
            new_name = name_str + str(name_num + nprofiles)
            new.limits[new_name] = [high_lim, low_lim, def_val]
        return new

    def __call__(self, x, *newargs, **newkwargs):
        """
        Call the profile functions
        """
        # loop over input arguments (non-keyword)
        for n in range(len(newargs)):
            arg = self.input_arguments[n]
            name, funcno = self.argument_map[arg]
            self.arguments[funcno][name] = newargs[n]
        # Loop over keyword input arguments
        for arg in newkwargs.keys():
            name, funcno = self.argument_map[arg]
            self.arguments[funcno][name] = newkwargs[arg]
        # Apply limits to arguments
        # Any argument values outside the limits will be reverted to defaults
        self.applylimits()
        # Call each function using the stored argument values
        # Sum all the profile functions together
        return np.sum([fn(x, **ags) for fn, ags in zip(self.functions, self.arguments)], axis=0)
