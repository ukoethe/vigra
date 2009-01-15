/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
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

#ifndef VIGRA_QUADPROG_HXX
#define VIGRA_QUADPROG_HXX

#include <limits>
#include "mathutil.hxx"
#include "matrix.hxx"
#include "linear_solve.hxx"
#include "numerictraits.hxx"
#include "array_vector.hxx"

namespace vigra {

namespace detail {

template <class T, class C1, class C2, class C3>
bool quadprogAddConstraint(MultiArrayView<2, T, C1> & R, MultiArrayView<2, T, C2> & J, MultiArrayView<2, T, C3> & d, int& iq, double& R_norm)
{
    int n=rowCount(J);
#ifdef TRACE_SOLVER
    std::cerr << "Add constraint " << iq << '/';
#endif
    int i, j, k;
    double cc, ss, h, t1, t2, xny;
    
    /*  we have to find the Givens rotation which will reduce the element
        d(j) to zero.
        if it is already zero we don't have to do anything, except of
        decreasing j */  
    for (j = n - 1; j >= iq + 1; j--)
    {
        /*   The Givens rotation is done with the matrix (cc cs, cs -cc).
             If cc is one, then element (j) of d is zero compared with element
             (j - 1). Hence we don't have to do anything. 
             If cc is zero, then we just have to switch column (j) and column (j - 1) 
             of J. Since we only switch columns in J, we have to be careful how we
             update d depending on the sign of gs.
             Otherwise we have to apply the Givens rotation to these columns.
             The i - 1 element of d has to be updated to h. */
        cc = d(j - 1,0);
        ss = d(j,0);
        h = hypot(cc, ss);
        if (h == 0.0)
            continue;
        d(j,0) = 0.0;
        ss = ss / h;
        cc = cc / h;
        if (cc < 0.0)
        {
            cc = -cc;
            ss = -ss;
            d(j - 1,0) = -h;
        }
        else
            d(j - 1,0) = h;
        xny = ss / (1.0 + cc);
        for (k = 0; k < n; k++)
        {
            t1 = J(k,j - 1);
            t2 = J(k,j);
            J(k,j - 1) = t1 * cc + t2 * ss;
            J(k,j) = xny * (t1 + J(k,j - 1)) - t2;
        }
    }
    /* update the number of constraints added*/
    iq++;
    /* To update R we have to put the iq components of the d vector
    into column iq - 1 of R
    */
    for (i = 0; i < iq; i++)
        R(i,iq - 1) = d(i,0);
#ifdef TRACE_SOLVER
  std::cerr << "add iq: " << iq << std::endl;
#endif
  
    if (abs(d(iq - 1,0)) <= NumericTraits<T>::epsilon() * R_norm)
        // problem degenerate
        return false;
    R_norm = std::max<T>(R_norm, abs(d(iq - 1,0)));
    return true;
}

template <class T, class C1, class C2, class C3, class C4>
void quadprogDeleteConstraint(MultiArrayView<2, T, C1> & R, MultiArrayView<2, T, C2> & J, MultiArrayView<2, int, C3> & A, MultiArrayView<2, T, C4> & u,  int p, int& iq, int l)
{
    int n=rowCount(R);
#ifdef TRACE_SOLVER
    std::cerr << "Delete constraint " << l << ' ' << iq;
#endif
    int i, j, k, qq;
    double cc, ss, h, xny, t1, t2;
  
    /* Find the index qq for active constraint l to be removed */
    for (i = p; i < iq; i++)
        if (A(i,0) == l)
        {
            qq = i;
            break;
        }
      
    /* remove the constraint from the active set and the duals */
    for (i = qq; i < iq - 1; i++)
    {
        A(i,0) = A(i + 1,0);
        u(i,0) = u(i + 1,0);
        for (j = 0; j < n; j++)
            R(j,i) = R(j,i + 1);
    }

    A(iq - 1,0) = A(iq,0);
    u(iq - 1,0) = u(iq,0);
    A(iq,0) = 0; 
    u(iq,0) = 0.0;
    for (j = 0; j < iq; j++)
        R(j,iq - 1) = 0.0;
    /* constraint has been fully removed */
    iq--;
#ifdef TRACE_SOLVER
    std::cerr << '/' << iq << std::endl;
#endif 
  
    if (iq == 0)
        return;

    for (j = qq; j < iq; j++)
    {
        cc = R(j,j);
        ss = R(j + 1,j);
        h = hypot(cc, ss);
        if (h == 0.0)
          continue;
        cc = cc / h;
        ss = ss / h;
        R(j + 1,j) = 0.0;
        if (cc < 0.0)
        {
            R(j,j) = -h;
            cc = -cc;
            ss = -ss;
        }
        else
            R(j,j) = h;

        xny = ss / (1.0 + cc);
        for (k = j + 1; k < iq; k++)
        {
            t1 = R(j,k);
            t2 = R(j + 1,k);
            R(j,k) = t1 * cc + t2 * ss;
            R(j + 1,k) = xny * (t1 + R(j,k)) - t2;
        }
        for (k = 0; k < n; k++)
        {
            t1 = J(k,j);
            t2 = J(k,j + 1);
            J(k,j) = t1 * cc + t2 * ss;
            J(k,j + 1) = xny * (J(k,j) + t1) - t2;
        }
    }
}

} // namespace detail

/** \addtogroup Optimization Optimization and Regression
 */
//@{
   /** Solve Quadratic Programming Problem.

     The quadraticProgramming() function implements the algorithm of Goldfarb and Idnani 
     for the solution of a (convex) Quadratic Programming problem
     by means of a dual method.
         
     <b>\#include</b> \<<a href="quadprog_8hxx-source.html">vigra/quadprog.hxx</a>\>
         Namespaces: vigra and vigra::linalg

     <b> Declaration:</b>

     \code
     namespace vigra { 
         template <class T, class C1, class C2, class C3, class C4, class C5, class C6, class C7>
         T 
         quadraticProgramming(MultiArrayView<2, T, C1> const & GG, MultiArrayView<2, T, C2> const & g,  
                              MultiArrayView<2, T, C3> const & CE, MultiArrayView<2, T, C4> const & ce,  
                              MultiArrayView<2, T, C5> const & CI, MultiArrayView<2, T, C6> const & ci, 
                              MultiArrayView<2, T, C7> & x);
     }
     \endcode

     The problem is in the form:

     min 0.5 * x' G x + g' x
     s.t.
         CE x == ce
         CI x >= ci
                  
     The matrix and vectors dimensions are as follows:
         G: n * n
            g: n
                    
            CE: p * n
         ce: p
                    
          CI: m * n
       ci: m

         x: n
     
     The function will return the cost of the solution written in the x vector or
     std::numeric_limits::infinity() if the problem is infeasible. In the latter case
     the value of the x vector is not correct.
     
     References: D. Goldfarb, A. Idnani. A numerically stable dual method for solving
                 strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

     Implementation based on work by Luca Di Gaspero and Angelo Furfaro.
   */
template <class T, class C1, class C2, class C3, class C4, class C5, class C6, class C7>
T 
quadraticProgramming(MultiArrayView<2, T, C1> const & GG, MultiArrayView<2, T, C2> const & g,  
               MultiArrayView<2, T, C3> const & CE, MultiArrayView<2, T, C4> const & ce,  
               MultiArrayView<2, T, C5> const & CI, MultiArrayView<2, T, C6> const & ci, 
               MultiArrayView<2, T, C7> & x)
{
    using namespace linalg;
    typedef typename MultiArrayShape<2>::type Shape;
    
    int i, k, l; /* indices */
    int ip, me, mi;
    int n=rowCount(g);  int p=rowCount(ce);  int m=rowCount(ci);  
    Matrix<T> G(GG.shape()), R(GG.shape()), J(GG.shape());

    Matrix<T> s(m+p,1), z(n,1), r(m + p,1), d(n,1), u(m + p,1);
    Matrix<T> x_old(n,1), u_old(m + p,1);
    T f_value, psi, c1, c2, ss, R_norm;
    const T inf = std::numeric_limits<T>::infinity();
    double t, t1, t2; /* t is the step lenght, which is the minimum of the partial step length t1 
    * and the full step length t2 */
    Matrix<int> A(m + p,1), A_old(m + p,1), iai(m + p,1);
    int q;
    int iq, iter = 0;
    ArrayVector<UInt8> iaexcl(m + p, false);

    me = p; /* number of equality constraints */
    mi = m; /* number of inequality constraints */
    q = 0;  /* size of the active set A (containing the indices of the active constraints) */

    /*
     * Preprocessing phase
     */
    
    /* compute the trace of the original matrix G */
    c1 = trace(GG);
    /* decompose the matrix G in the form L^T L */
    choleskyDecomposition(GG, G);

    R_norm = 1.0; /* this variable will hold the norm of the matrix R */
  
    /* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
    linearSolveLowerTriangular(G, identityMatrix<T>(n), transpose(J));
    c2 = trace(J);
    T epsilon = NumericTraits<T>::epsilon() * sq(J.norm(0));
#ifdef TRACE_SOLVER
    print_matrix("J", J, n);
#endif
  
    /* c1 * c2 is an estimate for cond(G) */
  
    /* 
     * Find the unconstrained minimizer of the quadratic form 0.5 * x G x + g x
     *      x = G^-1 * (-g)
     * this is a feasible point in the dual space
    */
    choleskySolve(G, -g, x);
    /* and compute the current solution value */ 
    f_value = 0.5 * dot(g, x);
#ifdef TRACE_SOLVER
    std::cerr << "Unconstrained solution: " << f_value << std::endl;
    print_vector("x", x, n);
#endif
  
    /* Add equality constraints to the working set A */
    iq = 0;
    for (i = 0; i < me; i++)
    {
        MultiArrayView<2, T, C3> np = rowVector(CE,i);
        d = transpose(np*J);
        z = J.subarray(Shape(0, iq), Shape(n,n))*subVector(d, iq, n);
        linearSolveUpperTriangular(R.subarray(Shape(0, 0), Shape(iq,iq)), 
                                   subVector(d, 0, iq), 
                                   subVector(r, 0, iq));
#ifdef TRACE_SOLVER
        print_matrix("R", R, iq);
        print_vector("z", z, n);
        print_vector("r", r, iq);
        print_vector("d", d, n);
#endif
    
        /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint 
           becomes feasible */
        t2 = 0.0;
        if (abs(dot(z, z)) > epsilon) // i.e. z != 0
            t2 = (-dot(np, x) + ce(i,0)) / dot(z, np);
    
        x += t2 * z;    
        /* set u = u+ */
        u(iq,0) = t2;
        subVector(u, 0, iq) -= t2 * subVector(r, 0, iq);
        
        /* compute the new solution value */
        f_value += 0.5 * (t2 * t2) * dot(z, np);
        A(i,0) = -i - 1;
    
        vigra_precondition(vigra::detail::quadprogAddConstraint(R, J, d, iq, R_norm),
            "quadraticProgramming(): Equality constraints are linearly dependent.");
    }
  
    /* set iai = K \ A */
    for (i = 0; i < mi; i++)
        iai(i,0) = i;
  
l1:    
    iter++;
#ifdef TRACE_SOLVER
    print_vector("x", x, n);
#endif
    /* step 1: choose a violated constraint */
    for (i = me; i < iq; i++)
    {
        ip = A(i,0);
        iai(ip,0) = -1;
    }
    
    /* compute s(x) = CI * x - ci for all elements of K \ A */
    ss = 0.0;
    psi = 0.0; /* this value will contain the sum of all infeasibilities */
    ip = 0; /* ip will be the index of the chosen violated constraint */
    for (i = 0; i < mi; i++)
    {
        iaexcl[i] = true;
        s(i,0) = dot(rowVector(CI, i), x) - ci(i,0);
        psi += std::min(0.0, s(i,0));
    }
#ifdef TRACE_SOLVER
    print_vector("s", s, mi);
#endif
    
    if (abs(psi) <= mi * NumericTraits<T>::epsilon() * c1 * c2* 100.0)
    {
        /* numerically there are not infeasibilities anymore */
        q = iq;
        return f_value;
    }
    
    /* save old values for u and A */
    for (i = 0; i < iq; i++)
    {
        u_old(i,0) = u(i,0);
        A_old(i,0) = A(i,0);
    }
    /* and for x */
    x_old = x;

#if 0
    for (i = 0; i < n; i++)
        x_old(i,0) = x(i,0);
#endif /* #if 0 */

    
l2: /* Step 2: check for feasibility and determine a new S-pair */
    for (i = 0; i < mi; i++)
    {
        if (s(i,0) < ss && iai(i,0) != -1 && iaexcl[i])
        {
            ss = s(i,0);
            ip = i;
        }
    }
    if (ss >= 0.0)
    {
        q = iq;
        return f_value;
    }
    
    /* set np = n(ip) */
    MultiArrayView<2, T, C5> np = rowVector(CI, ip);
    /* set u = (u 0)^T */
    u(iq,0) = 0.0;
    /* add ip to the active set A */
    A(iq,0) = ip;

#ifdef TRACE_SOLVER
    std::cerr << "Trying with constraint " << ip << std::endl;
    print_vector("np", np, n);
#endif
    
l2a:/* Step 2a: determine step direction */
    /* compute z = H np: the step direction in the primal space (through J, see the paper) */
    d = transpose(np*J);
    z = J.subarray(Shape(0, iq), Shape(n,n))*subVector(d, iq, n);
    /* compute N* np (if q > 0): the negative of the step direction in the dual space */
    linearSolveUpperTriangular(R.subarray(Shape(0, 0), Shape(iq,iq)), 
                               subVector(d, 0, iq), 
                               subVector(r, 0, iq));
#ifdef TRACE_SOLVER
    std::cerr << "Step direction z" << std::endl;
    print_vector("z", z, n);
    print_vector("r", r, iq + 1);
    print_vector("u", u, iq + 1);
    print_vector("d", d, n);
    print_ivector("A", A, iq + 1);
#endif
    
    /* Step 2b: compute step length */
    l = 0;
    /* Compute t1: partial step length (maximum step in dual 
       space without violating dual feasibility */
    t1 = inf; /* +inf */
    /* find the index l s.t. it reaches the minimum of u+(x) / r */
    for (k = me; k < iq; k++)
    {
        if (r(k,0) > 0.0)
        {
            if (u(k,0) / r(k,0) < t1)
            {
                t1 = u(k,0) / r(k,0);
                l = A(k,0);
            }
        }
    }
    /* Compute t2: full step length (minimum step in primal space 
       such that the constraint ip becomes feasible */
    if (abs(dot(z, z))  > epsilon) // i.e. z != 0
        t2 = -s(ip,0) / dot(z, np);
    else
        t2 = inf; /* +inf */
  
    /* the step is chosen as the minimum of t1 and t2 */
    t = std::min(t1, t2);
#ifdef TRACE_SOLVER
    std::cerr << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2 << ") ";
#endif
  
    /* Step 2c: determine new S-pair and take step: */
  
    /* case (i): no step in primal or dual space */
    if (t >= inf)
    {
        /* QPP is infeasible */
        // FIXME: unbounded to raise
        q = iq;
        return inf;
    }
    /* case (ii): step in dual space */
    if (t2 >= inf)
    {
        /* set u = u +  t * [-r 1) and drop constraint l from the active set A */
        subVector(u, 0, iq) -= t * subVector(r, 0, iq);
        u(iq,0) += t;
        iai(l,0) = l;
        vigra::detail::quadprogDeleteConstraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
        std::cerr << " in dual space: " << f_value << std::endl;
        print_vector("x", x, n);
        print_vector("z", z, n);
        print_ivector("A", A, iq + 1);
#endif
        goto l2a;
    }
  
    /* case (iii): step in primal and dual space */
  
    x += t * z;
    /* update the solution value */
    f_value += t * dot(z, np) * (0.5 * t + u(iq,0));
    /* u = u + t * (-r 1) */
    subVector(u, 0, iq) -= t * subVector(r, 0, iq);
    u(iq,0) += t;
#ifdef TRACE_SOLVER
    std::cerr << " in both spaces: " << f_value << std::endl;
    print_vector("x", x, n);
    print_vector("u", u, iq + 1);
    print_vector("r", r, iq + 1);
    print_ivector("A", A, iq + 1);
#endif
  
    if (t == t2)
    {
#ifdef TRACE_SOLVER
        std::cerr << "Full step has taken " << t << std::endl;
        print_vector("x", x, n);
#endif
        /* full step has taken */
        /* add constraint ip to the active set*/
        if (!vigra::detail::quadprogAddConstraint(R, J, d, iq, R_norm))
        {
            iaexcl[ip] = false;
            vigra::detail::quadprogDeleteConstraint(R, J, A, u, p, iq, ip);
#ifdef TRACE_SOLVER
            print_matrix("R", R, n);
            print_ivector("A", A, iq);
#endif
            for (i = 0; i < m; i++)
                iai(i,0) = i;
            for (i = 0; i < iq; i++)
            {
                A(i,0) = A_old(i,0);
                iai(A(i,0),0) = -1;
                u(i,0) = u_old(i,0);
            }
            x = x_old;
            goto l2; /* go to step 2 */
        }    
        else
        {
            iai(ip,0) = -1;
        }
#ifdef TRACE_SOLVER
        print_matrix("R", R, n);
        print_ivector("A", A, iq);
#endif
        goto l1;
    }
  
    /* a patial step has taken */
#ifdef TRACE_SOLVER
    std::cerr << "Partial step has taken " << t << std::endl;
    print_vector("x", x, n);
#endif
    /* drop constraint l */
    iai(l,0) = l;
    vigra::detail::quadprogDeleteConstraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
    print_matrix("R", R, n);
    print_ivector("A", A, iq);
#endif
  
    s(ip,0) = dot(rowVector(CI, ip), x) - ci(ip,0);

#ifdef TRACE_SOLVER
    print_vector("s", s, mi);
#endif
    goto l2a;
}

//@}

} // namespace vigra

#endif
