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
# GNU General Public License for more dc2ils.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                             #
###############################################################################

# Python module imports.
from copy import deepcopy
from numpy import dot, sqrt
import pdb

def nocedal_wright_wolfe(func, func_prime, x, f, g, p, a_init=1.0, max_a=1e5, c1=0.001, c2=0.9, tol=1e-10):
    """A line search algorithm implemented using the strong Wolfe conditions.

    Algorithm 3.2, page 59, from 'Numerical Optimization' by Jorge Nocedal and Stephen J. Wright,
    1999, 2nd ed.

    Requires the gradient function.


    Function options
    ~~~~~~~~~~~~~~~~

    func       - The n-dim function to minimise.
    func_prime - The function which returns the gradient vector.
    x          - The parameter vector for starting point of line search.
    f          - The function value f(x) at the point x.
    g          - The function gradient vector f'(x) at the point x.
    p          - The descent direction.
    a_init     - Initial step length.
    max_a      - The maxium value for the step length.

    c1         - Constant determining the slope for the sufficient decrease condition
    c2         - Constant used for the Wolfe curvature condition (0 < c1 < c2 < 1).
    for c1 and c2 0 < c1 < c2 < 1 must hold.
    """

    assert 0 < c1 and c1 < c2 and c2 < 1

    # Initialise values.
    a0 = {} ##ok
    a0['a'] = 0.0                ##alpha_0
    a0['phi'] = f                ##phi(0)
    a0['phi_prime'] = dot(g, p)  ##phi'(0)
    a_last = deepcopy(a0)        ##?
    a_max = {}
    a_max['a'] = max_a           ##alpha_max
    a_max['phi'] = func(x + a_max['a']*p)##phi(alpha_max)
    a_max['phi_prime'] = dot(func_prime(x + a_max['a']*p), p)##phi'(alpha_max)

    # Initialise sequence data.
    a = {}
    a['a'] = a_init ##alpha_1
    a['phi'] = func(x + a['a']*p)##phi(alpha_1)
    a['phi_prime'] = dot(func_prime(x + a['a']*p), p)##phi'(alpha_1)

    i = 0
    while 1:

        ## first if in description of algorithm 3.2
        ## FIXME: second term in disjunction missing
        # Check if the sufficient decrease condition is violated.
        # If so, the interval (a(i-1), a) will contain step lengths satisfying the strong Wolfe conditions.
        ##if not a['phi'] <= a0['phi'] + c1 * a['a'] * a0['phi_prime']:
        if a['phi'] > a0['phi'] + c1 * a['a'] * a0['phi_prime'] or \
               (a['phi'] >= a_last['phi'] and i > 1):
            print 'triggerd first cond'
            return zoom(func, func_prime, x, f, g, p, c1, c2, a0, a_last, a, tol)

        ## second if in description of algorithm 3.2
        # Check if the curvature condition is met and if so, return the step length a which satisfies the strong Wolfe conditions.
        if abs(a['phi_prime']) <= -c2 * a0['phi_prime']:
            print 'triggerd second cond'
            return a['a']

        ## third if in description of algorithm 3.2
        # Check if the gradient is positive.
        # If so, the interval (a(i-1), a) will contain step lengths satisfying the strong Wolfe conditions.
        if a['phi_prime'] >= 0.0:
            print 'triggerd third cond'
            # The arguments to zoom are a followed by a_last, because the function value at a_last will be higher than at a.
            return zoom(func, func_prime, x, f, g, p, c1, c2, a0, a, a_last, tol)

        ## FIXME: review interpolation
        # The strong Wolfe conditions are not met, and therefore interpolation between a and a_max will be used to find a value for a(i+1).
        #a_new = cubic_ext(a['a'], a_max['a'], a['phi'], a_max['phi'], a['phi_prime'], a_max['phi_prime'], full_output=0)
        a_new = a['a'] + 0.25 * (a_max['a'] - a['a'])

        # Update.
        a_last = deepcopy(a)
        a['a'] = a_new
        a['phi'] = func(x + a['a']*p)
        a['phi_prime'] = dot(func_prime(x + a['a']*p,), p)

        ## additional test mentioned in the supplemental text of the description
        # Test if the difference in function values is less than the tolerance.
        if abs(a_last['phi'] - a['phi']) <= tol:
            print 'triggerd fourht cond'
            return a['a']
        i = i + 1

def zoom(func, func_prime, x, f, g, p, c1, c2, a0, a_lo, a_hi, tol):
    """Find the minic1m function value in the open interval (a_lo, a_hi)

    Algorithm 3.3, page 60, from 'Numerical Optimization' by Jorge Nocedal and Stephen J. Wright,
    1999, 2nd ed.
    """

    # Initialize aj.
    aj = {}
    aj_last = deepcopy(a_lo)

    while 1:
        ## FIXME: review interpolation
        # Interpolate to find a trial step length aj between a_lo and a_hi.
        aj_new = quadratic(a_lo['a'], a_hi['a'], a_lo['phi'], a_hi['phi'], a_lo['phi_prime'])
        ##assert a_lo['a'] < aj_new and aj_new < a_hi['a']

        ## FIXME: review saveguarding
        # Safeguarding aj['a']
        aj['a'] = max(aj_last['a'] + 0.66*(a_hi['a'] - aj_last['a']), aj_new)

        # Calculate the function and gradient value at aj['a'].
        aj['phi'] = func(x + aj['a']*p)
        aj['phi_prime'] = dot(func_prime(x + aj['a']*p), p)

        # Check if the sufficient decrease condition is violated.
        if not aj['phi'] <= a0['phi'] + c1 * aj['a'] * a0['phi_prime']:
            a_hi = deepcopy(aj)
        else:
            # Check if the curvature condition is met and if so, return the step length ai which satisfies the strong Wolfe conditions.
            if abs(aj['phi_prime']) <= -c2 * a0['phi_prime']:
                return aj['a']

            # Determine if a_hi needs to be reset.
            if aj['phi_prime'] * (a_hi['a'] - a_lo['a']) >= 0.0:
                a_hi = deepcopy(a_lo)

            a_lo = deepcopy(aj)

        # Test if the difference in function values is less than the tolerance.
        if abs(aj_last['phi'] - aj['phi']) <= tol:
            return aj['a']

        # Update.
        aj_last = deepcopy(aj)

##def quadratic(a, b, fa, fb, ga):
##    """Quadratic interpolation using f(a), f(b), and g(a).

##    The extremum of the quadratic is given by:

##                     1             g(a)
##        aq  =  a  +  - . -------------------------
##                     2   f(a) - f(b) - (a - b)g(a)
##    """

##    denom = fa - fb - (a - b)*ga
##    if denom == 0.0:
##        return 1e99
##    else:
##        return a + 0.5 * ga / denom

def quadratic(a, b, fa, fb, ga):

    num   =(2*fa*a - a*(2*fb + ga*a) + ga*b*b)
    denom = (2*(fa - fb + ga*(-a + b)))
    if denom == 0.0:
        return 1e99
    else:
        return num / denom

