###############################################################################
#                                                                             #
# Copyright (C) 2003, 2008 Edward d'Auvergne                                  #
#                                                                             #
# This file is part of the minfx optimisation library.                        #
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                             #
###############################################################################

# Python module imports.
from copy import deepcopy
from numpy import dot, sqrt

# Minfx module imports.
from interpolate import cubic_ext, quadratic_fafbga

# Rename the interpolation functions.
quadratic = quadratic_fafbga


def nocedal_wright_interpol(func, args, x, f, g, p, a_init=1.0, mu=0.001, print_flag=0):
    """A line search algorithm based on interpolation.

    Pages 56-57, from 'Numerical Optimization' by Jorge Nocedal and Stephen J. Wright, 1999, 2nd ed.

    Requires the gradient function.


    Function options
    ~~~~~~~~~~~~~~~~

    func       - The function to minimise.
    func_prime - The function which returns the gradient vector.
    args       - The tuple of arguments to supply to the functions func and dfunc.
    x          - The parameter vector at minimisation step k.
    f          - The function value at the point x.
    g          - The function gradient vector at the point x.
    p          - The descent direction.
    a_init     - Initial step length.
    mu         - Constant determining the slope for the sufficient decrease condition (0 < mu < 1).
    """

    # Initialise values.
    i = 1
    f_count = 0
    a0 = {}
    a0['a'] = 0.0
    a0['phi'] = f
    a0['phi_prime'] = dot(g, p)

    # Initialise sequence data.
    a = {}
    a['a'] = a_init
    a['phi'] = func(*(x + a['a']*p,)+args)
    f_count = f_count + 1

    if print_flag:
        print "\n<Line search initial values>"
        print_data("Pre (a0)", i-1, a0)

    # Test for errors.
    if a0['phi_prime'] >= 0.0:
        raise NameError, "The gradient at point 0 of this line search is positive, ie p is not a descent direction and the line search will not work."

    # Check for sufficient decrease.  If so, return a_init.  Otherwise the interval [0, a_init] contains acceptable step lengths.
    if a['phi'] <= a0['phi'] + mu * a['a'] * a0['phi_prime']:
        return a['a'], f_count

    # Backup a_last.
    a_last = deepcopy(a)

    # Quadratic interpolation.
    a_new = - 0.5 * a0['phi_prime'] * a['a']**2 / (a['phi'] - a0['phi'] - a0['phi_prime']*a['a'])
    a['a'] = a_new
    a['phi'] = func(*(x + a['a']*p,)+args)
    f_count = f_count + 1

    # Check for sufficient decrease.  If so, return a['a'].
    if a['phi'] <= a0['phi'] + mu * a['a'] * a0['phi_prime']:
        return a['a'], f_count

    while 1:
        if print_flag:
            print "<Line search iteration i = " + `i` + " >"
            print_data("Initial (a)", i, a)
            print_data("Initial (a_last)", i, a_last)

        beta1 = 1.0 / (a_last['a']**2 * a['a']**2 * (a['a'] - a_last['a']))
        beta2 = a['phi'] - a0['phi'] - a0['phi_prime'] * a['a']
        beta3 = a_last['phi'] - a0['phi'] - a0['phi_prime'] * a_last['a']

        fact_a = beta1 * (a_last['a']**2 * beta2 - a['a']**2 * beta3)
        fact_b = beta1 * (-a_last['a']**3 * beta2 + a['a']**3 * beta3)

        a_new = {}
        a_new['a'] = (-fact_b + sqrt(fact_b**2 - 3.0 * fact_a * a0['phi_prime'])) / (3.0 * fact_a)
        a_new['phi'] = func(*(x + a_new['a']*p,)+args)
        f_count = f_count + 1

        # Check for sufficient decrease.  If so, return a_new['a'].
        if a_new['phi'] <= a0['phi'] + mu * a_new['a'] * a0['phi_prime']:
            return a_new['a'], f_count

        # Safeguarding.
        if a['a'] - a_new['a'] > 0.5 * a['a'] or 1.0 - a_new['a']/a['a'] < 0.9:
            a_new['a'] = 0.5 * a['a']
            a_new['phi'] = func(*(x + a_new['a']*p,)+args)
            f_count = f_count + 1

        # Updating.
        a_last = deepcopy(a)
        a = deepcopy(a_new)


def print_data(text, k, a):
    """Temp func for debugging."""

    print text + " data printout:"
    print "   Iteration:      " + `k`
    print "   a:              " + `a['a']`
    print "   phi:            " + `a['phi']`
    print "   phi_prime:      " + `a['phi_prime']`
