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
from math import sqrt
from numpy import dot
import sys

# Minfx module imports.
from interpolate import cubic_int, cubic_ext, quadratic_fafbga, quadratic_gagb

# Rename the interpolation functions.
cubic = cubic_int
quadratic = quadratic_fafbga
secant = quadratic_gagb


def more_thuente(func, func_prime, x, f, g, p, a_init=1.0, a_min=1e-25, a_max=None, a_tol=1e-10, phi_min=-1e3, mu=0.001, eta=0.1):
    """A line search algorithm from More and Thuente.

    More, J. J., and Thuente, D. J. 1994, Line search algorithms with guaranteed sufficient
    decrease. ACM Trans. Math. Softw. 20, 286-307.

    Function options (FIX)
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

    Internal variables
    ~~~~~~~~~~~~~~~~~~

    a0, the null sequence data structure containing the following keys:
        'a'        - 0
        'phi'        - phi(0)
        'phi_prime'    - phi'(0)

    a, the sequence data structure containing the following keys:
        'a'        - alpha
        'phi'        - phi(alpha)
        'phi_prime'    - phi'(alpha)

    Ik, the interval data structure containing the following keys:
        'a'        - The current interval Ik = [al, au]
        'phi'        - The interval [phi(al), phi(au)]
        'phi_prime'    - The interval [phi'(al), phi'(au)]

    Instead of using the modified function:
        psi(a) = phi(a) - phi(0) - a.phi'(0),
    the function:
        psi(a) = phi(a) - a.phi'(0),
    was used as the phi(0) component has no effect on the results.
    """

    # Initialise values.
    k = 0
    mod_flag = 1
    bracketed = 0

    a0 = {}
    a0['a'] = 0.0
    a0['phi'] = f
    a0['phi_prime'] = dot(g, p)
    if not a_min:
        a_min = 0.0
    if not a_max:
        a_max = 4.0*max(1.0,a_init)
        
    Ik_lim = [0.0, 5.0*a_init]
    width = a_max - a_min
    width2 = 2.0*width

    # Initialise sequence data.
    a = {}
    a['a'] = a_init
    a['phi'] = func(x + a['a']*p)
    a['phi_prime'] = dot(func_prime(x + a['a']*p), p)

    # Initialise interval data.
    Ik = {}
    Ik['a'] = [0.0, 0.0]
    Ik['phi'] = [a0['phi'], a0['phi']]
    Ik['phi_prime'] = [a0['phi_prime'], a0['phi_prime']]


    # Test for errors.
    if a0['phi_prime'] > 0.0:
        raise ArgumentError, "The gradient at point 0 of this line search is positive, ie p is not a descent direction and the line search will not work."
    if a['a'] < a_min:
        raise ArgumentError, "Alpha is less than alpha_min, " + `a['a']` + " > " + `a_min`
    if a['a'] > a_max:
        raise ArgumentError, "Alpha is greater than alpha_max, " + `a['a']` + " > " + `a_max`

    while 1:
        # Test values.
        curv = mu * a0['phi_prime']
        suff_dec = a0['phi'] + a['a'] * curv

        # Modification flag, 0 - phi, 1 - psi.
        if mod_flag:
            if a['phi'] <= suff_dec and a['phi_prime'] >= 0.0:
                mod_flag = 0

        # Test for convergence using the strong Wolfe conditions.
        if a['phi'] <= suff_dec and abs(a['phi_prime']) <= eta * abs(a0['phi_prime']):
            return a['a']

        # Test if limits have been reached.
        if a['a'] == a_min:
            if a['phi'] > suff_dec or a['phi_prime'] >= curv:
                return a['a']
        if a['a'] == a_max:
            if a['phi'] <= suff_dec and a['phi_prime'] <= curv:
                return a['a']

        if bracketed:
            # Test for roundoff error.
            if a['a'] <= Ik_lim[0] or a['a'] >= Ik_lim[1]:
                return a['a'], f_count, g_count
            # Test to see if a_tol has been reached.
            if Ik_lim[1] - Ik_lim[0] <= a_tol * Ik_lim[1]:
                return a['a'], f_count, g_count
            
        # Choose a safeguarded ak in set Ik which is a subset of [a_min, a_max], and update the interval Ik.
        a_new = {}
        if mod_flag and a['phi'] <= Ik['phi'][0] and a['phi'] > suff_dec:
            # Calculate the modified function values and gradients at at, al, and au.
            psi = a['phi'] - curv * a['a']
            psi_l = Ik['phi'][0] - curv * Ik['a'][0]
            psi_u = Ik['phi'][1] - curv * Ik['a'][1]
            psi_prime = a['phi_prime'] - curv
            psi_l_prime = Ik['phi_prime'][0] - curv
            psi_u_prime = Ik['phi_prime'][1] - curv

            a_new['a'], Ik_new, bracketed = update(a, Ik, a['a'], Ik['a'][0], Ik['a'][1], psi, psi_l, psi_u, psi_prime, psi_l_prime, psi_u_prime, bracketed, Ik_lim)
        else:
            a_new['a'], Ik_new, bracketed = update(a, Ik, a['a'], Ik['a'][0], Ik['a'][1], a['phi'], Ik['phi'][0], Ik['phi'][1], a['phi_prime'], Ik['phi_prime'][0], Ik['phi_prime'][1], bracketed, Ik_lim)

        # Bisection step.
        if bracketed:
            size = abs(Ik_new['a'][0] - Ik_new['a'][1])
            if size >= 0.66 * width2:
                a_new['a'] = 0.5 * (Ik_new['a'][0] + Ik_new['a'][1])
            width2 = width
            width = size

        # Limit.
        if bracketed:
            Ik_lim[0] = min(Ik_new['a'][0], Ik_new['a'][1])
            Ik_lim[1] = max(Ik_new['a'][0], Ik_new['a'][1])
        else:
            Ik_lim[0] = a_new['a'] + 1.1 * (a_new['a'] - Ik_new['a'][0])
            Ik_lim[1] = a_new['a'] + 4.0 * (a_new['a'] - Ik_new['a'][0])

        if bracketed:
            if a_new['a'] <= Ik_lim[0] or a_new['a'] >= Ik_lim[1] or Ik_lim[1] - Ik_lim[0] <= a_tol * Ik_lim[1]:
                a_new['a'] = Ik['a'][0]

        # The step must be between a_min and a_max.
        if a_new['a'] < a_min:
            a_new['a'] = a_min
        if a_new['a'] > a_max:
            a_new['a'] = a_max

        # Calculate new values.
        a_new['phi'] = func(x + a_new['a']*p)
        a_new['phi_prime'] = dot(func_prime(x + a_new['a']*p), p)

        # Shift data from k+1 to k.
        k = k + 1
        a = deepcopy(a_new)
        Ik = deepcopy(Ik_new)

def update(a, Ik, at, al, au, ft, fl, fu, gt, gl, gu, bracketed, Ik_lim, d=0.66):
    """Trial value selection and interval updating.

    Trial value selection
    ~~~~~~~~~~~~~~~~~~~~~

    fl, fu, ft, gl, gu, and gt are the function and gradient values at the interval end points al
        and au, and at the trial point at.
    ac is the minimiser of the cubic that interpolates fl, ft, gl, and gt.
    aq is the minimiser of the quadratic that interpolates fl, ft, and gl.
    as is the minimiser of the quadratic that interpolates fl, gl, and gt.

    The trial value selection is divided into four cases.

    Case 1: ft > fl.  In this case compute ac and aq, and set

               / ac,            if |ac - al| < |aq - al|,
        at+ = <
               \ 1/2(aq + ac),  otherwise.


    Case 2: ft <= fl and gt.gl < 0.  In this case compute ac and as, and set

               / ac,            if |ac - at| >= |as - at|,
        at+ = <
               \ as,            otherwise.


    Case 3: ft <= fl and gt.gl >= 0, and |gt| <= |gl|.  In this case at+ is chosen by extrapolating
    the function values at al and at, so the trial value at+ lies outside th interval with at and al
    as endpoints.  Compute ac and as.

        If the cubic tends to infinity in the direction of the step and the minimum of the cubic is
        beyound at, set

                   / ac,            if |ac - at| < |as - at|,
            at+ = <
                   \ as,            otherwise.

        Otherwise set at+ = as.


        Redefine at+ by setting

                   / min{at + d(au - at), at+},        if at > al.
            at+ = <
                   \ max{at + d(au - at), at+},        otherwise,

        for some d < 1.


    Case 4: ft <= fl and gt.gl >= 0, and |gt| > |gl|.  In this case choose at+ as the minimiser of
    the cubic that interpolates fu, ft, gu, and gt.


    Interval updating
    ~~~~~~~~~~~~~~~~~

    Given a trial value at in I, the endpoints al+ and au+ of the updated interval I+ are determined
    as follows:
        Case U1: If f(at) > f(al), then al+ = al and au+ = at.
        Case U2: If f(at) <= f(al) and f'(at)(al - at) > 0, then al+ = at and au+ = au.
        Case U3: If f(at) <= f(al) and f'(at)(al - at) < 0, then al+ = at and au+ = al.
    """

    # Trial value selection.

    # Case 1.
    if ft > fl:
        # The minimum is bracketed.
        bracketed = 1

        # Interpolation.
        ac = cubic(al, at, fl, ft, gl, gt)
        aq = quadratic(al, at, fl, ft, gl)

        # Return at+.
        if abs(ac - al) < abs(aq - al):
            at_new = ac
        else:
            at_new = 0.5*(aq + ac)


    # Case 2.
    elif gt * gl < 0.0:
        # The minimum is bracketed.
        bracketed = 1

        # Interpolation.
        ac = cubic(al, at, fl, ft, gl, gt)
        as = secant(al, at, gl, gt)

        # Return at+.
        if abs(ac - at) >= abs(as - at):
            at_new = ac
        else:
            at_new = as


    # Case 3.
    elif abs(gt) <= abs(gl):
        # Interpolation.
        ac, beta1, beta2 = cubic_ext(al, at, fl, ft, gl, gt, full_output=1)

        if ac > at and beta2 != 0.0:
            # Leave ac as ac.
            pass
        elif at > al:
            # Set ac to the upper limit.
            ac = Ik_lim[1]
        else:
            # Set ac to the lower limit.
            ac = Ik_lim[0]

        as = secant(al, at, gl, gt)

        # Test if bracketed.
        if bracketed:
            if abs(ac - at) < abs(as - at):
                at_new = ac
            else:
                at_new = as

            # Redefine at+.
            if at > al:
                at_new = min(at + d*(au - at), at_new)
            else:
                at_new = max(at + d*(au - at), at_new)
        else:
            if abs(ac - at) > abs(as - at):
                at_new = ac
            else:
                at_new = as

            # Check limits.
            if at_new < Ik_lim[0]:
                at_new = Ik_lim[0]
            if at_new > Ik_lim[1]:
                at_new = Ik_lim[1]

    # Case 4.
    else:
        if bracketed:
            at_new = cubic(au, at, fu, ft, gu, gt)
        elif at > al:
            at_new = Ik_lim[1]
        else:
            at_new = Ik_lim[0]


    # Interval updating algorithm.
    Ik_new = deepcopy(Ik)

    if ft > fl:
        Ik_new['a'][1] = at
        Ik_new['phi'][1] = a['phi']
        Ik_new['phi_prime'][1] = a['phi_prime']
    elif gt*(al - at) > 0.0:
        Ik_new['a'][0] = at
        Ik_new['phi'][0] = a['phi']
        Ik_new['phi_prime'][0] = a['phi_prime']
    else:
        Ik_new['a'][0] = at
        Ik_new['phi'][0] = a['phi']
        Ik_new['phi_prime'][0] = a['phi_prime']
        Ik_new['a'][1] = al
        Ik_new['phi'][1] = Ik['phi'][0]
        Ik_new['phi_prime'][1] = Ik['phi_prime'][0]

    return at_new, Ik_new, bracketed
