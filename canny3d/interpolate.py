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
from math import sqrt


def cubic_int(a, b, fa, fb, ga, gb):
    """Cubic interpolation using f(a), f(b), g(a), and g(b).

    Equations
    ~~~~~~~~~

    f(a) = a'a**3 + b'a**2 + c'a + d'
    f(b) = a'b**3 + b'b**2 + c'b + d'
    g(a) = 3a'a**2 + 2b'a + c'
    g(b) = 3a'b**2 + 2b'b + c'


    Interpolation
    ~~~~~~~~~~~~~

    The extrema are the roots of the quadratic equation:

        3a'*alpha**2 + 2b'*alpha + c' = 0

    The cubic interpolant is given by the formula:

                           g(b) + beta2 - beta1
        ac = b - (b - a) . ---------------------
                           g(b) - g(a) + 2*beta2

    where:
                                  f(a) - f(b)
        beta1 = g(a) + g(b) - 3 . -----------
                                     a - b

        if a < b:
            beta2 = sqrt(beta1**2 - g(a).g(b))
        else:
            beta2 = -sqrt(beta1**2 - g(a).g(b))
    """


    beta1 = ga + gb - 3.0*(fa - fb)/(a - b)
    if a < b:
        beta2 = sqrt(beta1**2 - ga*gb)
    else:
        beta2 = -sqrt(beta1**2 - ga*gb)

    denom = gb - ga + 2.0*beta2
    if denom == 0.0:
        return -1e99
    else:
        return b - (b - a)*(gb + beta2 - beta1)/denom


def cubic_ext(a, b, fa, fb, ga, gb, full_output=0):
    """Cubic Extrapolation using f(a), f(b), g(a), and g(b).

    Extrapolation
    ~~~~~~~~~~~~~

    The extrema are the roots of the quadratic equation:

        3a'*alpha**2 + 2b'*alpha + c' = 0

    The cubic extrapolant is given by the formula:

                           g(b) + beta2 - beta1
        ac = b - (b - a) . ---------------------
                           g(b) - g(a) + 2*beta2

    where:

                                  f(a) - f(b)
        beta1 = g(a) + g(b) - 3 . -----------
                                     a - b

        if a < b:
            beta2 = sqrt(max(0.0, beta1**2 - g(a).g(b)))
        else:
            beta2 = -sqrt(max(0.0, beta1**2 - g(a).g(b)))
    """


    beta1 = ga + gb - 3.0*(fa - fb)/(a - b)
    if a < b:
        beta2 = sqrt(max(0.0, beta1**2 - ga*gb))
    else:
        beta2 = -sqrt(max(0.0, beta1**2 - ga*gb))

    denom = gb - ga + 2.0*beta2
    if denom == 0.0:
        alpha = -1e99
    else:
        alpha = b - (b - a)*(gb + beta2 - beta1)/denom

    if full_output:
        return alpha, beta1, beta2
    else:
        return alpha


def quadratic_fafbga(a, b, fa, fb, ga):
    """Quadratic interpolation using f(a), f(b), and g(a).

    The extremum of the quadratic is given by:

                     1             g(a)
        aq  =  a  +  - . -------------------------
                     2   f(a) - f(b) - (a - b)g(a)
    """

    denom = fa - fb - (a - b)*ga
    if denom == 0.0:
        return 1e99
    else:
        return a + 0.5 * ga / denom


def quadratic_gagb(a, b, ga, gb):
    """Quadratic interpolation using g(a) and g(b).

    The extremum of the quadratic is given by:

               bg(a) - ag(b)
        aq  =  -------------
                g(a) - g(b)
    """

    denom = ga - gb
    if denom == 0.0:
        return 1e99
    else:
        return (b*ga - a*gb) / denom
