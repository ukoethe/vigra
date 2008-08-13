/************************************************************************/
/*                                                                      */
/*                  Copyright 2007 by Ullrich Koethe                    */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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


#ifndef VIGRA_SINGULAR_VALUE_DECOMPOSITION_HXX
#define VIGRA_SINGULAR_VALUE_DECOMPOSITION_HXX

#include "matrix.hxx"
#include "array_vector.hxx"


namespace vigra
{

namespace linalg
{

   /** Singular Value Decomposition.
       \ingroup MatrixAlgebra

   For an m-by-n matrix \a A with m >= n, the singular value decomposition is
   an m-by-n orthogonal matrix \a U, an n-by-n diagonal matrix S, and
   an n-by-n orthogonal matrix \a V so that A = U*S*V'.

   To save memory, this functions stores the matrix \a S in a column vector of
   appropriate length (a diagonal matrix can be obtained by <tt>diagonalMatrix(S)</tt>).
   The singular values, sigma[k] = S(k, 0), are ordered so that
   sigma[0] >= sigma[1] >= ... >= sigma[n-1].

   The singular value decomposition always exists, so this function will
   never fail (except if the shapes of the argument matrices don't match).
   The effective numerical rank of A is returned.

	(Adapted from JAMA, a Java Matrix Library, developed jointly
	by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).

    <b>\#include</b> \<<a href="singular__value__decomposition_8hxx-source.html">vigra/singular_value_decomposition.hxx</a>\> or<br>
    <b>\#include</b> \<<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>\><br>
        Namespaces: vigra and vigra::linalg
   */
template <class T, class C1, class C2, class C3, class C4>
unsigned int
singularValueDecomposition(MultiArrayView<2, T, C1> const & A,
    MultiArrayView<2, T, C2> &U, MultiArrayView<2, T, C3> &S, MultiArrayView<2, T, C4> &V)
{
    typedef T Real;
    typedef MultiArrayShape<2>::type Shape;

    const MultiArrayIndex rows = rowCount(A);
    const MultiArrayIndex cols = columnCount(A);
    vigra_precondition(rows >= cols,
       "singularValueDecomposition(): Input matrix A must be rectangular with rowCount >= columnCount.");
    vigra_precondition(rowCount(S) == cols && columnCount(S) == 1,
       "singularValueDecomposition(): Output S must be column vector with rowCount == columnCount(A).");
    vigra_precondition(rowCount(U) == rows && columnCount(U) == cols,
       "singularValueDecomposition(): Output matrix U must have the same dimensions as input matrix A.");
    vigra_precondition(rowCount(V) == cols && columnCount(V) == cols,
       "singularValueDecomposition(): Output matrix V must be square with n = columnCount(A).");

    MultiArrayIndex m = rows;
    MultiArrayIndex n = cols;
    MultiArrayIndex nu = n;

    U.init(0.0);
    S.init(0.0);
    V.init(0.0);

    ArrayVector<Real> e((unsigned int)n);
    ArrayVector<Real> work((unsigned int)m);
    Matrix<Real> a(A);
    MultiArrayView<1, T, C3> s = S.bindOuter(0);

    MultiArrayIndex i=0, j=0, k=0;

    // Reduce a to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.
    MultiArrayIndex nct = std::min(m-1,n);
    MultiArrayIndex nrt = std::max((MultiArrayIndex)0,n-2);
    for (k = 0; k < std::max(nct,nrt); ++k)
    {
        if (k < nct)
        {
            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s(k).
            // Compute 2-norm of k-th column without under/overflow.
            s(k) = 0.0;
            for (i = k; i < m; ++i)
            {
               s(k) = hypot(s(k), a(i, k));
            }
            if (s(k) != 0.0)
            {
                if (a(k, k) < 0.0)
                {
                    s(k) = -s(k);
                }
                for (i = k; i < m; ++i)
                {
                   a(i, k) /= s(k);
                }
                a(k, k) += 1.0;
            }
            s(k) = -s(k);
        }
        for (j = k+1; j < n; ++j)
        {
            if ((k < nct) && (s(k) != 0.0))
            {
                // Apply the transformation.
                Real t(0.0);
                for (i = k; i < m; ++i)
                {
                    t += a(i, k)*a(i, j);
                }
                t = -t/a(k, k);
                for (i = k; i < m; ++i)
                {
                    a(i, j) += t*a(i, k);
                }
            }

            // Place the k-th row of a into e for the
            // subsequent calculation of the row transformation.

            e[j] = a(k, j);
        }
        if (k < nct)
        {
            // Place the transformation in U for subsequent back
            // multiplication.

            for (i = k; i < m; ++i)
            {
                U(i, k) = a(i, k);
            }
        }
        if (k < nrt)
        {
            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (i = k+1; i < n; ++i)
            {
               e[k] = hypot(e[k],e[i]);
            }
            if (e[k] != 0.0)
            {
                if (e[k+1] < 0.0)
                {
                    e[k] = -e[k];
                }
                for (i = k+1; i < n; ++i)
                {
                    e[i] /= e[k];
                }
                e[k+1] += 1.0;
            }
            e[k] = -e[k];
            if ((k+1 < m) & (e[k] != 0.0))
            {
                // Apply the transformation.
                for (i = k+1; i < m; ++i)
                {
                    work[i] = 0.0;
                }
                for (j = k+1; j < n; ++j)
                {
                    for (i = k+1; i < m; ++i)
                    {
                        work[i] += e[j]*a(i, j);
                    }
                }
                for (j = k+1; j < n; ++j)
                {
                    Real t(-e[j]/e[k+1]);
                    for (i = k+1; i < m; ++i)
                    {
                        a(i, j) += t*work[i];
                    }
                }
            }
            // Place the transformation in V for subsequent
            // back multiplication.
            for (i = k+1; i < n; ++i)
            {
                V(i, k) = e[i];
            }
        }
    }

    // Set up the final bidiagonal matrix of order p.

    MultiArrayIndex p = n;
    if (nct < n)
    {
        s(nct) = a(nct, nct);
    }
    if (m < p)
    {
        s(p-1) = 0.0;
    }
    if (nrt+1 < p)
    {
        e[nrt] = a(nrt, p-1);
    }
    e[p-1] = 0.0;

    // Generate U.
    for (j = nct; j < nu; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            U(i, j) = 0.0;
        }
        U(j, j) = 1.0;
    }
    for (k = nct-1; k >= 0; --k)
    {
        if (s(k) != 0.0)
        {
            for (j = k+1; j < nu; ++j)
            {
                Real t(0.0);
                for (i = k; i < m; ++i)
                {
                    t += U(i, k)*U(i, j);
                }
                t = -t/U(k, k);
                for (i = k; i < m; ++i)
                {
                    U(i, j) += t*U(i, k);
                }
            }
            for (i = k; i < m; ++i )
            {
                U(i, k) = -U(i, k);
            }
            U(k, k) = 1.0 + U(k, k);
            for (i = 0; i < k-1; ++i)
            {
                U(i, k) = 0.0;
            }
        }
        else
        {
            for (i = 0; i < m; ++i)
            {
                U(i, k) = 0.0;
            }
            U(k, k) = 1.0;
        }
    }

    // Generate V.
    for (k = n-1; k >= 0; --k)
    {
        if ((k < nrt) & (e[k] != 0.0))
        {
            for (j = k+1; j < nu; ++j)
            {
                Real t(0.0);
                for (i = k+1; i < n; ++i)
                {
                    t += V(i, k)*V(i, j);
                }
                t = -t/V(k+1, k);
                for (i = k+1; i < n; ++i)
                {
                    V(i, j) += t*V(i, k);
                }
            }
        }
        for (i = 0; i < n; ++i)
        {
            V(i, k) = 0.0;
        }
        V(k, k) = 1.0;
    }

    // Main iteration loop for the singular values.

    MultiArrayIndex pp = p-1;
    int iter = 0;
    Real eps = NumericTraits<Real>::epsilon()*2.0;
    while (p > 0)
    {
        MultiArrayIndex k=0;
        int kase=0;

        // Here is where a test for too many iterations would go.

        // This section of the program inspects for
        // negligible elements in the s and e arrays.  On
        // completion the variables kase and k are set as follows.

        // kase = 1     if s(p) and e[k-1] are negligible and k<p
        // kase = 2     if s(k) is negligible and k<p
        // kase = 3     if e[k-1] is negligible, k<p, and
        //              s(k), ..., s(p) are not negligible (qr step).
        // kase = 4     if e(p-1) is negligible (convergence).

        for (k = p-2; k >= -1; --k)
        {
            if (k == -1)
            {
                break;
            }
            if (abs(e[k]) <= eps*(abs(s(k)) + abs(s(k+1))))
            {
                e[k] = 0.0;
                break;
            }
        }
        if (k == p-2)
        {
            kase = 4;
        }
        else
        {
            MultiArrayIndex ks;
            for (ks = p-1; ks >= k; --ks)
            {
                if (ks == k)
                {
                    break;
                }
                Real t( (ks != p ? abs(e[ks]) : 0.) +
                        (ks != k+1 ? abs(e[ks-1]) : 0.));
                if (abs(s(ks)) <= eps*t)
                {
                    s(ks) = 0.0;
                    break;
                }
            }
            if (ks == k)
            {
               kase = 3;
            }
            else if (ks == p-1)
            {
               kase = 1;
            }
            else
            {
               kase = 2;
               k = ks;
            }
        }
        ++k;

        // Perform the task indicated by kase.

        switch (kase)
        {
          case 1: // Deflate negligible s(p).
          {
              Real f(e[p-2]);
              e[p-2] = 0.0;
              for (j = p-2; j >= k; --j)
              {
                  Real t( hypot(s(j),f));
                  Real cs(s(j)/t);
                  Real sn(f/t);
                  s(j) = t;
                  if (j != k)
                  {
                      f = -sn*e[j-1];
                      e[j-1] = cs*e[j-1];
                  }
                  for (i = 0; i < n; ++i)
                  {
                      t = cs*V(i, j) + sn*V(i, p-1);
                      V(i, p-1) = -sn*V(i, j) + cs*V(i, p-1);
                      V(i, j) = t;
                  }
              }
              break;
          }
          case 2: // Split at negligible s(k).
          {
              Real f(e[k-1]);
              e[k-1] = 0.0;
              for (j = k; j < p; ++j)
              {
                  Real t(hypot(s(j),f));
                  Real cs( s(j)/t);
                  Real sn(f/t);
                  s(j) = t;
                  f = -sn*e[j];
                  e[j] = cs*e[j];
                  for (i = 0; i < m; ++i)
                  {
                      t = cs*U(i, j) + sn*U(i, k-1);
                      U(i, k-1) = -sn*U(i, j) + cs*U(i, k-1);
                      U(i, j) = t;
                  }
              }
              break;
          }
          case 3: // Perform one qr step.
          {
              // Calculate the shift.
              Real scale = std::max(std::max(std::max(std::max(
                      abs(s(p-1)),abs(s(p-2))),abs(e[p-2])),
                      abs(s(k))),abs(e[k]));
              Real sp = s(p-1)/scale;
              Real spm1 = s(p-2)/scale;
              Real epm1 = e[p-2]/scale;
              Real sk = s(k)/scale;
              Real ek = e[k]/scale;
              Real b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
              Real c = (sp*epm1)*(sp*epm1);
              Real shift = 0.0;
              if ((b != 0.0) || (c != 0.0))
              {
                  shift = VIGRA_CSTD::sqrt(b*b + c);
                  if (b < 0.0)
                  {
                      shift = -shift;
                  }
                  shift = c/(b + shift);
              }
              Real f = (sk + sp)*(sk - sp) + shift;
              Real g = sk*ek;

              // Chase zeros.
              for (j = k; j < p-1; ++j)
              {
                  Real t = hypot(f,g);
                  Real cs = f/t;
                  Real sn = g/t;
                  if (j != k)
                  {
                      e[j-1] = t;
                  }
                  f = cs*s(j) + sn*e[j];
                  e[j] = cs*e[j] - sn*s(j);
                  g = sn*s(j+1);
                  s(j+1) = cs*s(j+1);
                  for (i = 0; i < n; ++i)
                  {
                      t = cs*V(i, j) + sn*V(i, j+1);
                      V(i, j+1) = -sn*V(i, j) + cs*V(i, j+1);
                      V(i, j) = t;
                  }
                  t = hypot(f,g);
                  cs = f/t;
                  sn = g/t;
                  s(j) = t;
                  f = cs*e[j] + sn*s(j+1);
                  s(j+1) = -sn*e[j] + cs*s(j+1);
                  g = sn*e[j+1];
                  e[j+1] = cs*e[j+1];
                  if (j < m-1)
                  {
                      for (i = 0; i < m; ++i)
                      {
                          t = cs*U(i, j) + sn*U(i, j+1);
                          U(i, j+1) = -sn*U(i, j) + cs*U(i, j+1);
                          U(i, j) = t;
                      }
                  }
              }
              e[p-2] = f;
              iter = iter + 1;
              break;
          }
          case 4:  // Convergence.
          {
              // Make the singular values positive.
              if (s(k) <= 0.0)
              {
                  s(k) = (s(k) < 0.0 ? -s(k) : 0.0);
                  for (i = 0; i <= pp; ++i)
                  {
                      V(i, k) = -V(i, k);
                  }
              }

              // Order the singular values.

              while (k < pp)
              {
                  if (s(k) >= s(k+1))
                  {
                      break;
                  }
                  Real t = s(k);
                  s(k) = s(k+1);
                  s(k+1) = t;
                  if (k < n-1)
                  {
                      for (i = 0; i < n; ++i)
                      {
                           t = V(i, k+1); V(i, k+1) = V(i, k); V(i, k) = t;
                      }
                  }
                  if (k < m-1)
                  {
                      for (i = 0; i < m; ++i)
                      {
                          t = U(i, k+1); U(i, k+1) = U(i, k); U(i, k) = t;
                      }
                  }
                  ++k;
              }
              iter = 0;
              --p;
              break;
          }
          default:
              vigra_fail("vigra::svd(): Internal error.");
        }
    }
    Real tol = std::max(m,n)*s(0)*eps;
    unsigned int rank = 0;
    for (MultiArrayIndex i = 0; i < n; ++i)
    {
        if (s(i) > tol)
        {
            ++rank;
        }
    }
    return rank; // effective rank
}

} // namespace linalg

using linalg::singularValueDecomposition;

} // namespace vigra

#endif // VIGRA_SINGULAR_VALUE_DECOMPOSITION_HXX
