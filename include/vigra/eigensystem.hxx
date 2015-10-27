/************************************************************************/
/*                                                                      */
/*                  Copyright 2004 by Ullrich Koethe                    */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/


#ifndef VIGRA_EIGENSYSTEM_HXX
#define VIGRA_EIGENSYSTEM_HXX

#include <algorithm>
#include <complex>
#include "matrix.hxx"
#include "array_vector.hxx"
#include "polynomial.hxx"

namespace vigra
{

namespace linalg
{

namespace detail
{

// code adapted from JAMA
// a and b will be overwritten
template <class T, class C1, class C2>
void
housholderTridiagonalization(MultiArrayView<2, T, C1> &a, MultiArrayView<2, T, C2> &b)
{
    const MultiArrayIndex n = rowCount(a);
    vigra_precondition(n == columnCount(a),
        "housholderTridiagonalization(): matrix must be square.");
    vigra_precondition(n == rowCount(b) && 2 <= columnCount(b),
        "housholderTridiagonalization(): matrix size mismatch.");

    MultiArrayView<1, T, C2> d = b.bindOuter(0);
    MultiArrayView<1, T, C2> e = b.bindOuter(1);

    for(MultiArrayIndex j = 0; j < n; ++j)
    {
        d(j) = a(n-1, j);
    }

    // Householder reduction to tridiagonalMatrix form.

    for(int i = n-1; i > 0; --i)
    {
        // Scale to avoid under/overflow.

        T scale = 0.0;
        T h = 0.0;
        for(int k = 0; k < i; ++k)
        {
            scale = scale + abs(d(k));
        }
        if(scale == 0.0)
        {
            e(i) = d(i-1);
            for(int j = 0; j < i; ++j)
            {
                d(j) = a(i-1, j);
                a(i, j) = 0.0;
                a(j, i) = 0.0;
            }
        }
        else
        {
            // Generate Householder vector.

            for(int k = 0; k < i; ++k)
            {
                d(k) /= scale;
                h += sq(d(k));
            }
            T f = d(i-1);
            T g = VIGRA_CSTD::sqrt(h);
            if(f > 0) {
                g = -g;
            }
            e(i) = scale * g;
            h -= f * g;
            d(i-1) = f - g;
            for(int j = 0; j < i; ++j)
            {
               e(j) = 0.0;
            }

            // Apply similarity transformation to remaining columns.

            for(int j = 0; j < i; ++j)
            {
               f = d(j);
               a(j, i) = f;
               g = e(j) + a(j, j) * f;
               for(int k = j+1; k <= i-1; ++k)
               {
                   g += a(k, j) * d(k);
                   e(k) += a(k, j) * f;
               }
               e(j) = g;
            }
            f = 0.0;
            for(int j = 0; j < i; ++j)
            {
                e(j) /= h;
                f += e(j) * d(j);
            }
            T hh = f / (h + h);
            for(int j = 0; j < i; ++j)
            {
                e(j) -= hh * d(j);
            }
            for(int j = 0; j < i; ++j)
            {
                f = d(j);
                g = e(j);
                for(int k = j; k <= i-1; ++k)
                {
                    a(k, j) -= (f * e(k) + g * d(k));
                }
                d(j) = a(i-1, j);
                a(i, j) = 0.0;
            }
        }
        d(i) = h;
    }

    // Accumulate transformations.

    for(MultiArrayIndex i = 0; i < n-1; ++i)
    {
        a(n-1, i) = a(i, i);
        a(i, i) = 1.0;
        T h = d(i+1);
        if(h != 0.0)
        {
            for(MultiArrayIndex k = 0; k <= i; ++k)
            {
                d(k) = a(k, i+1) / h;
            }
            for(MultiArrayIndex j = 0; j <= i; ++j)
            {
                T g = 0.0;
                for(MultiArrayIndex k = 0; k <= i; ++k)
                {
                    g += a(k, i+1) * a(k, j);
                }
                for(MultiArrayIndex k = 0; k <= i; ++k)
                {
                    a(k, j) -= g * d(k);
                }
            }
        }
        for(MultiArrayIndex k = 0; k <= i; ++k)
        {
            a(k, i+1) = 0.0;
        }
    }
    for(MultiArrayIndex j = 0; j < n; ++j)
    {
        d(j) = a(n-1, j);
        a(n-1, j) = 0.0;
    }
    a(n-1, n-1) = 1.0;
    e(0) = 0.0;
}

// code adapted from JAMA
// de and z will be overwritten
template <class T, class C1, class C2>
bool
tridiagonalMatrixEigensystem(MultiArrayView<2, T, C1> &de, MultiArrayView<2, T, C2> &z)
{
    MultiArrayIndex n = rowCount(z);
    vigra_precondition(n == columnCount(z),
        "tridiagonalMatrixEigensystem(): matrix must be square.");
    vigra_precondition(n == rowCount(de) && 2 <= columnCount(de),
        "tridiagonalMatrixEigensystem(): matrix size mismatch.");

    MultiArrayView<1, T, C2> d = de.bindOuter(0);
    MultiArrayView<1, T, C2> e = de.bindOuter(1);

    for(MultiArrayIndex i = 1; i < n; i++) {
       e(i-1) = e(i);
    }
    e(n-1) = 0.0;

    T f = 0.0;
    T tst1 = 0.0;
    T eps = VIGRA_CSTD::pow(2.0,-52.0);
    for(MultiArrayIndex l = 0; l < n; ++l)
    {
        // Find small subdiagonalMatrix element

        tst1 = std::max(tst1, abs(d(l)) + abs(e(l)));
        MultiArrayIndex m = l;

        // Original while-loop from Java code
        while(m < n)
        {
            if(abs(e(m)) <= eps*tst1)
                break;
            ++m;
        }

        // If m == l, d(l) is an eigenvalue,
        // otherwise, iterate.

        if(m > l)
        {
            int iter = 0;
            do
            {
                if(++iter > 50)
                   return false; // too many iterations

                // Compute implicit shift

                T g = d(l);
                T p = (d(l+1) - g) / (2.0 * e(l));
                T r = hypot(p,1.0);
                if(p < 0)
                {
                    r = -r;
                }
                d(l) = e(l) / (p + r);
                d(l+1) = e(l) * (p + r);
                T dl1 = d(l+1);
                T h = g - d(l);
                for(MultiArrayIndex i = l+2; i < n; ++i)
                {
                   d(i) -= h;
                }
                f = f + h;

                // Implicit QL transformation.

                p = d(m);
                T c = 1.0;
                T c2 = c;
                T c3 = c;
                T el1 = e(l+1);
                T s = 0.0;
                T s2 = 0.0;
                for(int i = m-1; i >= (int)l; --i)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e(i);
                    h = c * p;
                    r = hypot(p,e(i));
                    e(i+1) = s * r;
                    s = e(i) / r;
                    c = p / r;
                    p = c * d(i) - s * g;
                    d(i+1) = h + s * (c * g + s * d(i));

                    // Accumulate transformation.

                    for(MultiArrayIndex k = 0; k < n; ++k)
                    {
                         h = z(k, i+1);
                         z(k, i+1) = s * z(k, i) + c * h;
                        z(k, i) = c * z(k, i) - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e(l) / dl1;
                e(l) = s * p;
                d(l) = c * p;

                // Check for convergence.

            } while(abs(e(l)) > eps*tst1);
        }
        d(l) = d(l) + f;
        e(l) = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.

    for(MultiArrayIndex i = 0; i < n-1; ++i)
    {
        MultiArrayIndex k = i;
        T p = d(i);
        for(MultiArrayIndex j = i+1; j < n; ++j)
        {
            T p1 = d(j);
            if(p < p1)
            {
                k = j;
                p = p1;
            }
        }
        if(k != i)
        {
            std::swap(d(k), d(i));
            for(MultiArrayIndex j = 0; j < n; ++j)
            {
                std::swap(z(j, i), z(j, k));
            }
        }
    }
    return true;
}

// Nonsymmetric reduction to Hessenberg form.

template <class T, class C1, class C2>
void nonsymmetricHessenbergReduction(MultiArrayView<2, T, C1> & H, MultiArrayView<2, T, C2> & V)
{
    //  This is derived from the Algol procedures orthes and ortran,
    //  by Martin and Wilkinson, Handbook for Auto. Comp.,
    //  Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutines in EISPACK.

    int n = rowCount(H);
    int low = 0;
    int high = n-1;
    ArrayVector<T> ort(n);

    for(int m = low+1; m <= high-1; ++m)
    {
        // Scale column.

        T scale = 0.0;
        for(int i = m; i <= high; ++i)
        {
            scale = scale + abs(H(i, m-1));
        }
        if(scale != 0.0)
        {

            // Compute Householder transformation.

            T h = 0.0;
            for(int i = high; i >= m; --i)
            {
                ort[i] = H(i, m-1)/scale;
                h += sq(ort[i]);
            }
            T g = VIGRA_CSTD::sqrt(h);
            if(ort[m] > 0)
            {
                g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;

            // Apply Householder similarity transformation
            // H = (I-u*u'/h)*H*(I-u*u')/h)

            for(int j = m; j < n; ++j)
            {
                T f = 0.0;
                for(int i = high; i >= m; --i)
                {
                    f += ort[i]*H(i, j);
                }
                f = f/h;
                for(int i = m; i <= high; ++i)
                {
                    H(i, j) -= f*ort[i];
                }
            }

            for(int i = 0; i <= high; ++i)
            {
                T f = 0.0;
                for(int j = high; j >= m; --j)
                {
                    f += ort[j]*H(i, j);
                }
                f = f/h;
                for(int j = m; j <= high; ++j)
                {
                    H(i, j) -= f*ort[j];
                }
            }
            ort[m] = scale*ort[m];
            H(m, m-1) = scale*g;
        }
    }

    // Accumulate transformations (Algol's ortran).

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            V(i, j) = (i == j ? 1.0 : 0.0);
        }
    }

    for(int m = high-1; m >= low+1; --m)
    {
        if(H(m, m-1) != 0.0)
        {
            for(int i = m+1; i <= high; ++i)
            {
                ort[i] = H(i, m-1);
            }
            for(int j = m; j <= high; ++j)
            {
                T g = 0.0;
                for(int i = m; i <= high; ++i)
                {
                    g += ort[i] * V(i, j);
                }
                // Double division avoids possible underflow
                g = (g / ort[m]) / H(m, m-1);
                for(int i = m; i <= high; ++i)
                {
                    V(i, j) += g * ort[i];
                }
            }
        }
    }
}


// Complex scalar division.

template <class T>
void cdiv(T xr, T xi, T yr, T yi, T & cdivr, T & cdivi)
{
    T r,d;
    if(abs(yr) > abs(yi))
    {
        r = yi/yr;
        d = yr + r*yi;
        cdivr = (xr + r*xi)/d;
        cdivi = (xi - r*xr)/d;
    }
    else
    {
        r = yr/yi;
        d = yi + r*yr;
        cdivr = (r*xr + xi)/d;
        cdivi = (r*xi - xr)/d;
    }
}

template <class T, class C>
int hessenbergQrDecompositionHelper(MultiArrayView<2, T, C> & H, int n, int l, double eps,
             T & p, T & q, T & r, T & s, T & w, T & x, T & y, T & z)
{
    int m = n-2;
    while(m >= l)
    {
        z = H(m, m);
        r = x - z;
        s = y - z;
        p = (r * s - w) / H(m+1, m) + H(m, m+1);
        q = H(m+1, m+1) - z - r - s;
        r = H(m+2, m+1);
        s = abs(p) + abs(q) + abs(r);
        p = p / s;
        q = q / s;
        r = r / s;
        if(m == l)
        {
            break;
        }
        if(abs(H(m, m-1)) * (abs(q) + abs(r)) <
            eps * (abs(p) * (abs(H(m-1, m-1)) + abs(z) +
            abs(H(m+1, m+1)))))
        {
                break;
        }
        --m;
    }
    return m;
}



// Nonsymmetric reduction from Hessenberg to real Schur form.

template <class T, class C1, class C2, class C3>
bool hessenbergQrDecomposition(MultiArrayView<2, T, C1> & H, MultiArrayView<2, T, C2> & V,
                                     MultiArrayView<2, T, C3> & de)
{

    //  This is derived from the Algol procedure hqr2,
    //  by Martin and Wilkinson, Handbook for Auto. Comp.,
    //  Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    // Initialize
    MultiArrayView<1, T, C3> d = de.bindOuter(0);
    MultiArrayView<1, T, C3> e = de.bindOuter(1);

    int nn = rowCount(H);
    int n = nn-1;
    int low = 0;
    int high = nn-1;
    T eps = VIGRA_CSTD::pow(2.0, sizeof(T) == sizeof(float)
                                     ? -23.0
                                     : -52.0);
    T exshift = 0.0;
    T p=0,q=0,r=0,s=0,z=0,t,w,x,y;
    T norm = vigra::norm(H);

    // Outer loop over eigenvalue index
    int iter = 0;
    while(n >= low)
    {

        // Look for single small sub-diagonal element
        int l = n;
        while (l > low)
        {
            s = abs(H(l-1, l-1)) + abs(H(l, l));
            if(s == 0.0)
            {
                s = norm;
            }
            if(abs(H(l, l-1)) < eps * s)
            {
                break;
            }
            --l;
        }

        // Check for convergence
        // One root found
        if(l == n)
        {
            H(n, n) = H(n, n) + exshift;
            d(n) = H(n, n);
            e(n) = 0.0;
            --n;
            iter = 0;

        // Two roots found

        }
        else if(l == n-1)
        {
            w = H(n, n-1) * H(n-1, n);
            p = (H(n-1, n-1) - H(n, n)) / 2.0;
            q = p * p + w;
            z = VIGRA_CSTD::sqrt(abs(q));
            H(n, n) = H(n, n) + exshift;
            H(n-1, n-1) = H(n-1, n-1) + exshift;
            x = H(n, n);

            // Real pair

            if(q >= 0)
            {
                if(p >= 0)
                {
                    z = p + z;
                }
                else
                {
                    z = p - z;
                }
                d(n-1) = x + z;
                d(n) = d(n-1);
                if(z != 0.0)
                {
                    d(n) = x - w / z;
                }
                e(n-1) = 0.0;
                e(n) = 0.0;
                x = H(n, n-1);
                s = abs(x) + abs(z);
                p = x / s;
                q = z / s;
                r = VIGRA_CSTD::sqrt(p * p+q * q);
                p = p / r;
                q = q / r;

                // Row modification

                for(int j = n-1; j < nn; ++j)
                {
                    z = H(n-1, j);
                    H(n-1, j) = q * z + p * H(n, j);
                    H(n, j) = q * H(n, j) - p * z;
                }

                // Column modification

                for(int i = 0; i <= n; ++i)
                {
                    z = H(i, n-1);
                    H(i, n-1) = q * z + p * H(i, n);
                    H(i, n) = q * H(i, n) - p * z;
                }

                // Accumulate transformations

                for(int i = low; i <= high; ++i)
                {
                    z = V(i, n-1);
                    V(i, n-1) = q * z + p * V(i, n);
                    V(i, n) = q * V(i, n) - p * z;
                }

            // Complex pair

            }
            else
            {
                d(n-1) = x + p;
                d(n) = x + p;
                e(n-1) = z;
                e(n) = -z;
            }
            n = n - 2;
            iter = 0;

        // No convergence yet

        }
        else
        {

            // Form shift

            x = H(n, n);
            y = 0.0;
            w = 0.0;
            if(l < n)
            {
                y = H(n-1, n-1);
                w = H(n, n-1) * H(n-1, n);
            }

            // Wilkinson's original ad hoc shift

            if(iter == 10)
            {
                exshift += x;
                for(int i = low; i <= n; ++i)
                {
                    H(i, i) -= x;
                }
                s = abs(H(n, n-1)) + abs(H(n-1, n-2));
                x = y = 0.75 * s;
                w = -0.4375 * s * s;
            }

            // MATLAB's new ad hoc shift

            if(iter == 30)
            {
                 s = (y - x) / 2.0;
                 s = s * s + w;
                 if(s > 0)
                 {
                      s = VIGRA_CSTD::sqrt(s);
                      if(y < x)
                      {
                          s = -s;
                      }
                      s = x - w / ((y - x) / 2.0 + s);
                      for(int i = low; i <= n; ++i)
                      {
                          H(i, i) -= s;
                      }
                      exshift += s;
                      x = y = w = 0.964;
                 }
            }

            iter = iter + 1;
            if(iter > 60)
                return false;

            // Look for two consecutive small sub-diagonal elements
            int m = hessenbergQrDecompositionHelper(H, n, l, eps, p, q, r, s, w, x, y, z);
            for(int i = m+2; i <= n; ++i)
            {
                H(i, i-2) = 0.0;
                if(i > m+2)
                {
                    H(i, i-3) = 0.0;
                }
            }

            // Double QR step involving rows l:n and columns m:n

            for(int k = m; k <= n-1; ++k)
            {
                int notlast = (k != n-1);
                if(k != m) {
                    p = H(k, k-1);
                    q = H(k+1, k-1);
                    r = (notlast ? H(k+2, k-1) : 0.0);
                    x = abs(p) + abs(q) + abs(r);
                    if(x != 0.0)
                    {
                        p = p / x;
                        q = q / x;
                        r = r / x;
                    }
                }
                if(x == 0.0)
                {
                    break;
                }
                s = VIGRA_CSTD::sqrt(p * p + q * q + r * r);
                if(p < 0)
                {
                    s = -s;
                }
                if(s != 0)
                {
                    if(k != m)
                    {
                        H(k, k-1) = -s * x;
                    }
                    else if(l != m)
                    {
                        H(k, k-1) = -H(k, k-1);
                    }
                    p = p + s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q = q / p;
                    r = r / p;

                    // Row modification

                    for(int j = k; j < nn; ++j)
                    {
                        p = H(k, j) + q * H(k+1, j);
                        if(notlast)
                        {
                            p = p + r * H(k+2, j);
                            H(k+2, j) = H(k+2, j) - p * z;
                        }
                        H(k, j) = H(k, j) - p * x;
                        H(k+1, j) = H(k+1, j) - p * y;
                    }

                    // Column modification

                    for(int i = 0; i <= std::min(n,k+3); ++i)
                    {
                        p = x * H(i, k) + y * H(i, k+1);
                        if(notlast)
                        {
                            p = p + z * H(i, k+2);
                            H(i, k+2) = H(i, k+2) - p * r;
                        }
                        H(i, k) = H(i, k) - p;
                        H(i, k+1) = H(i, k+1) - p * q;
                    }

                    // Accumulate transformations

                    for(int i = low; i <= high; ++i)
                    {
                        p = x * V(i, k) + y * V(i, k+1);
                        if(notlast)
                        {
                            p = p + z * V(i, k+2);
                            V(i, k+2) = V(i, k+2) - p * r;
                        }
                        V(i, k) = V(i, k) - p;
                        V(i, k+1) = V(i, k+1) - p * q;
                    }
                }  // (s != 0)
            }  // k loop
        }  // check convergence
    }  // while (n >= low)

    // Backsubstitute to find vectors of upper triangular form

    if(norm == 0.0)
    {
        return false;
    }

    for(n = nn-1; n >= 0; --n)
    {
        p = d(n);
        q = e(n);

        // Real vector

        if(q == 0)
        {
            int l = n;
            H(n, n) = 1.0;
            for(int i = n-1; i >= 0; --i)
            {
                w = H(i, i) - p;
                r = 0.0;
                for(int j = l; j <= n; ++j)
                {
                    r = r + H(i, j) * H(j, n);
                }
                if(e(i) < 0.0)
                {
                    z = w;
                    s = r;
                }
                else
                {
                    l = i;
                    if(e(i) == 0.0)
                    {
                        if(w != 0.0)
                        {
                            H(i, n) = -r / w;
                        }
                        else
                        {
                            H(i, n) = -r / (eps * norm);
                        }

                    // Solve real equations

                    }
                    else
                    {
                        x = H(i, i+1);
                        y = H(i+1, i);
                        q = (d(i) - p) * (d(i) - p) + e(i) * e(i);
                        t = (x * s - z * r) / q;
                        H(i, n) = t;
                        if(abs(x) > abs(z))
                        {
                            H(i+1, n) = (-r - w * t) / x;
                        }
                        else
                        {
                            H(i+1, n) = (-s - y * t) / z;
                        }
                    }

                    // Overflow control

                    t = abs(H(i, n));
                    if((eps * t) * t > 1)
                    {
                        for(int j = i; j <= n; ++j)
                        {
                            H(j, n) = H(j, n) / t;
                        }
                    }
                }
            }

        // Complex vector

        }
        else if(q < 0)
        {
            int l = n-1;

            // Last vector component imaginary so matrix is triangular

            if(abs(H(n, n-1)) > abs(H(n-1, n)))
            {
                H(n-1, n-1) = q / H(n, n-1);
                H(n-1, n) = -(H(n, n) - p) / H(n, n-1);
            }
            else
            {
                cdiv(0.0,-H(n-1, n),H(n-1, n-1)-p,q, H(n-1, n-1), H(n-1, n));
            }
            H(n, n-1) = 0.0;
            H(n, n) = 1.0;
            for(int i = n-2; i >= 0; --i)
            {
                T ra,sa,vr,vi;
                ra = 0.0;
                sa = 0.0;
                for(int j = l; j <= n; ++j)
                {
                    ra = ra + H(j, n-1) * H(i, j);
                    sa = sa + H(j, n) * H(i, j);
                }
                w = H(i, i) - p;

                if(e(i) < 0.0)
                {
                    z = w;
                    r = ra;
                    s = sa;
                }
                else
                {
                    l = i;
                    if(e(i) == 0)
                    {
                        cdiv(-ra,-sa,w,q, H(i, n-1), H(i, n));
                    }
                    else
                    {
                        // Solve complex equations

                        x = H(i, i+1);
                        y = H(i+1, i);
                        vr = (d(i) - p) * (d(i) - p) + e(i) * e(i) - q * q;
                        vi = (d(i) - p) * 2.0 * q;
                        if((vr == 0.0) && (vi == 0.0))
                        {
                            vr = eps * norm * (abs(w) + abs(q) +
                            abs(x) + abs(y) + abs(z));
                        }
                        cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi, H(i, n-1), H(i, n));
                        if(abs(x) > (abs(z) + abs(q)))
                        {
                            H(i+1, n-1) = (-ra - w * H(i, n-1) + q * H(i, n)) / x;
                            H(i+1, n) = (-sa - w * H(i, n) - q * H(i, n-1)) / x;
                        }
                        else
                        {
                            cdiv(-r-y*H(i, n-1),-s-y*H(i, n),z,q, H(i+1, n-1), H(i+1, n));
                        }
                    }

                    // Overflow control

                    t = std::max(abs(H(i, n-1)),abs(H(i, n)));
                    if((eps * t) * t > 1)
                    {
                        for(int j = i; j <= n; ++j)
                        {
                            H(j, n-1) = H(j, n-1) / t;
                            H(j, n) = H(j, n) / t;
                        }
                    }
                }
            }
        }
    }

    // Back transformation to get eigenvectors of original matrix

    for(int j = nn-1; j >= low; --j)
    {
        for(int i = low; i <= high; ++i)
        {
            z = 0.0;
            for(int k = low; k <= std::min(j,high); ++k)
            {
                z = z + V(i, k) * H(k, j);
            }
            V(i, j) = z;
        }
    }
    return true;
}

} // namespace detail

/** \addtogroup MatrixAlgebra
*/
//@{
    /** Compute the eigensystem of a symmetric matrix.

        \a a is a real symmetric matrix, \a ew is a single-column matrix
        holding the eigenvalues, and \a ev is a matrix of the same size as
        \a a whose columns are the corresponding eigenvectors. Eigenvalues
        will be sorted from largest to smallest.
        The algorithm returns <tt>false</tt> when it doesn't
        converge. It can be applied in-place, i.e. <tt>&a == &ev</tt> is allowed.
        The code of this function was adapted from JAMA.

        <b>\#include</b> \<vigra/eigensystem.hxx\> or<br>
        <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool
symmetricEigensystem(MultiArrayView<2, T, C1> const & a,
            MultiArrayView<2, T, C2> & ew, MultiArrayView<2, T, C3> & ev)
{
    vigra_precondition(isSymmetric(a),
        "symmetricEigensystem(): symmetric input matrix required.");
    const MultiArrayIndex acols = columnCount(a);
    vigra_precondition(1 == columnCount(ew) && acols == rowCount(ew) &&
                       acols == columnCount(ev) && acols == rowCount(ev),
        "symmetricEigensystem(): matrix shape mismatch.");

    ev.copy(a); // does nothing if &ev == &a
    Matrix<T> de(acols, 2);
    detail::housholderTridiagonalization(ev, de);
    if(!detail::tridiagonalMatrixEigensystem(de, ev))
        return false;

    ew.copy(columnVector(de, 0));
    return true;
}

namespace detail{

template <class T, class C2, class C3>
bool
symmetricEigensystem2x2(T a00, T a01, T a11,
            MultiArrayView<2, T, C2> & ew, MultiArrayView<2, T, C3> & ev)
{
    double evec[2]={0,0};

    /* Eigenvectors*/
    if (a01==0){
        if (fabs(a11)>fabs(a00)){
            evec[0]=0.;
            evec[1]=1.;
            ew(0,0)=a11;
            ew(1,0)=a00;
        }
        else if(fabs(a00)>fabs(a11)) {
            evec[0]=1.;
            evec[1]=0.;
            ew(0,0)=a00;
            ew(1,0)=a11;
        }
        else {
            evec[0]=.5* M_SQRT2;
            evec[1]=.5* M_SQRT2;
            ew(0,0)=a00;
            ew(1,0)=a11;
        }
    }
    else{
        double temp=a11-a00;

        double coherence=sqrt(temp*temp+4*a01*a01);
        evec[0]=2*a01;
        evec[1]=temp+coherence;
        temp=std::sqrt(evec[0]*evec[0]+evec[1]*evec[1]);
        if (temp==0){
            evec[0]=.5* M_SQRT2;
            evec[1]=.5* M_SQRT2;
            ew(0,0)=1.;
            ew(1,0)=1.;
        }
        else{
            evec[0]/=temp;
            evec[1]/=temp;

            /* Eigenvalues */
            ew(0,0)=.5*(a00+a11+coherence);
            ew(1,0)=.5*(a00+a11-coherence);
        }
    }
    ev(0,0)= evec[0];
    ev(1,0)= evec[1];
    ev(0,1)=-evec[1];
    ev(1,1)= evec[0];
    return true;
}

template <class T, class C2, class C3>
bool
symmetricEigensystem3x3(T a00, T a01, T a02, T a11, T a12, T a22,
            MultiArrayView<2, T, C2> & ew, MultiArrayView<2, T, C3> & ev)
{
    symmetric3x3Eigenvalues(a00, a01, a02, a11, a12, a22,
                            &ew(0,0), &ew(1,0), &ew(2,0));

    /* Calculate eigen vectors */
    double a1=a01*a12,
           a2=a01*a02,
           a3=sq(a01);

    double b1=a00-ew(0,0),
           b2=a11-ew(0,0);
    ev(0,0)=a1-a02*b2;
    ev(1,0)=a2-a12*b1;
    ev(2,0)=b1*b2-a3;

    b1=a00-ew(1,0);
    b2=a11-ew(1,0);
    ev(0,1)=a1-a02*b2;
    ev(1,1)=a2-a12*b1;
    ev(2,1)=b1*b2-a3;

    b1=a00-ew(2,0);
    b2=a11-ew(2,0);
    ev(0,2)=a1-a02*b2;
    ev(1,2)=a2-a12*b1;
    ev(2,2)=b1*b2-a3;

    /* Eigen vector normalization */
    double l0=norm(columnVector(ev, 0));
    double l1=norm(columnVector(ev, 1));
    double l2=norm(columnVector(ev, 2));

    /* Detect fail : eigenvectors with only zeros */
    double M = std::max(std::max(abs(ew(0,0)), abs(ew(1,0))), abs(ew(2,0)));
    double epsilon = 1e-12*M;
    if(l0<epsilon) { return false; }
    if(l1<epsilon) { return false; }
    if(l2<epsilon) { return false; }

    columnVector(ev, 0) /= l0;
    columnVector(ev, 1) /= l1;
    columnVector(ev, 2) /= l2;

    /* Succes    */
    return true;
}

} // closing namespace detail


    /** Fast computation of the eigensystem of a 2x2 or 3x3 symmetric matrix.

        The function works like \ref symmetricEigensystem(), but uses fast analytic
        formula to avoid iterative computations.

        <b>\#include</b> \<vigra/eigensystem.hxx\> or<br>
        <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool
symmetricEigensystemNoniterative(MultiArrayView<2, T, C1> const & a,
                                 MultiArrayView<2, T, C2> & ew, MultiArrayView<2, T, C3> & ev)
{
    vigra_precondition(isSymmetric(a),
        "symmetricEigensystemNoniterative(): symmetric input matrix required.");
    const MultiArrayIndex acols = columnCount(a);
    if(acols == 2)
    {
        detail::symmetricEigensystem2x2(a(0,0), a(0,1), a(1,1), ew, ev);
        return true;
    }
    if(acols == 3)
    {
        // try the fast algorithm
        if(detail::symmetricEigensystem3x3(a(0,0), a(0,1), a(0,2), a(1,1), a(1,2), a(2,2),
                                           ew, ev))
            return true;
        // fast algorithm failed => fall-back to iterative algorithm
        return symmetricEigensystem(a, ew, ev);
    }
    vigra_precondition(false,
        "symmetricEigensystemNoniterative(): can only handle 2x2 and 3x3 matrices.");
    return false;
}

    /** Compute the eigensystem of a square, but
        not necessarily symmetric matrix.

        \a a is a real square matrix, \a ew is a single-column matrix
        holding the possibly complex eigenvalues, and \a ev is a matrix of
        the same size as \a a whose columns are the corresponding eigenvectors.
        Eigenvalues will be sorted from largest to smallest magnitude.
        The algorithm returns <tt>false</tt> when it doesn't
        converge. It can be applied in-place, i.e. <tt>&a == &ev</tt> is allowed.
        The code of this function was adapted from JAMA.

        <b>\#include</b> \<vigra/eigensystem.hxx\> or<br>
        <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool
nonsymmetricEigensystem(MultiArrayView<2, T, C1> const & a,
         MultiArrayView<2, std::complex<T>, C2> & ew, MultiArrayView<2, T, C3> & ev)
{
    const MultiArrayIndex acols = columnCount(a);
    vigra_precondition(acols == rowCount(a),
        "nonsymmetricEigensystem(): square input matrix required.");
    vigra_precondition(1 == columnCount(ew) && acols == rowCount(ew) &&
                       acols == columnCount(ev) && acols == rowCount(ev),
        "nonsymmetricEigensystem(): matrix shape mismatch.");

    Matrix<T> H(a);
    Matrix<T> de(acols, 2);
    detail::nonsymmetricHessenbergReduction(H, ev);
    if(!detail::hessenbergQrDecomposition(H, ev, de))
        return false;

    for(MultiArrayIndex i = 0; i < acols; ++i)
    {
        ew(i,0) = std::complex<T>(de(i, 0), de(i, 1));
    }
    return true;
}

    /** Compute the roots of a polynomial using the eigenvalue method.

        \a poly is a real polynomial (compatible to \ref vigra::PolynomialView),
        and \a roots a complex valued vector (compatible to <tt>std::vector</tt>
        with a <tt>value_type</tt> compatible to the type <tt>POLYNOMIAL::Complex</tt>) to which
        the roots are appended. The function calls \ref nonsymmetricEigensystem() with the standard
        companion matrix yielding the roots as eigenvalues. It returns <tt>false</tt> if
        it fails to converge.

        <b>\#include</b> \<vigra/eigensystem.hxx\> or<br>
        <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg

        \see polynomialRoots(), vigra::Polynomial
     */
template <class POLYNOMIAL, class VECTOR>
bool polynomialRootsEigenvalueMethod(POLYNOMIAL const & poly, VECTOR & roots, bool polishRoots)
{
    typedef typename POLYNOMIAL::value_type T;
    typedef typename POLYNOMIAL::Real    Real;
    typedef typename POLYNOMIAL::Complex Complex;
    typedef Matrix<T> TMatrix;
    typedef Matrix<Complex> ComplexMatrix;

    int const degree = poly.order();
    double const eps = poly.epsilon();

    TMatrix inMatrix(degree, degree);
    for(int i = 0; i < degree; ++i)
        inMatrix(0, i) = -poly[degree - i - 1] / poly[degree];
    for(int i = 0; i < degree - 1; ++i)
        inMatrix(i + 1, i) = NumericTraits<T>::one();
    ComplexMatrix ew(degree, 1);
    TMatrix ev(degree, degree);
    bool success = nonsymmetricEigensystem(inMatrix, ew, ev);
    if(!success)
        return false;
    for(int i = 0; i < degree; ++i)
    {
        if(polishRoots)
            vigra::detail::laguerre1Root(poly, ew(i,0), 1);
        roots.push_back(vigra::detail::deleteBelowEpsilon(ew(i,0), eps));
    }
    std::sort(roots.begin(), roots.end(), vigra::detail::PolynomialRootCompare<Real>(eps));
    return true;
}

template <class POLYNOMIAL, class VECTOR>
bool polynomialRootsEigenvalueMethod(POLYNOMIAL const & poly, VECTOR & roots)
{
    return polynomialRootsEigenvalueMethod(poly, roots, true);
}

    /** Compute the real roots of a real polynomial using the eigenvalue method.

        \a poly is a real polynomial (compatible to \ref vigra::PolynomialView),
        and \a roots a real valued vector (compatible to <tt>std::vector</tt>
        with a <tt>value_type</tt> compatible to the type <tt>POLYNOMIAL::Real</tt>) to which
        the roots are appended. The function calls \ref polynomialRootsEigenvalueMethod() and
        throws away all complex roots. It returns <tt>false</tt> if it fails to converge.
        The parameter <tt>polishRoots</tt> is ignored (it is only here for syntax compatibility
        with polynomialRealRoots()).


        <b>\#include</b> \<vigra/eigensystem.hxx\> or<br>
        <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg

        \see polynomialRealRoots(), vigra::Polynomial
     */
template <class POLYNOMIAL, class VECTOR>
bool polynomialRealRootsEigenvalueMethod(POLYNOMIAL const & p, VECTOR & roots, bool /* polishRoots */)
{
    typedef typename NumericTraits<typename VECTOR::value_type>::ComplexPromote Complex;
    ArrayVector<Complex> croots;
    if(!polynomialRootsEigenvalueMethod(p, croots))
        return false;
    for(unsigned int i = 0; i < croots.size(); ++i)
        if(croots[i].imag() == 0.0)
            roots.push_back(croots[i].real());
    return true;
}

template <class POLYNOMIAL, class VECTOR>
bool polynomialRealRootsEigenvalueMethod(POLYNOMIAL const & p, VECTOR & roots)
{
    return polynomialRealRootsEigenvalueMethod(p, roots, true);
}


//@}

} // namespace linalg

using linalg::symmetricEigensystem;
using linalg::symmetricEigensystemNoniterative;
using linalg::nonsymmetricEigensystem;
using linalg::polynomialRootsEigenvalueMethod;
using linalg::polynomialRealRootsEigenvalueMethod;

} // namespace vigra

#endif // VIGRA_EIGENSYSTEM_HXX
