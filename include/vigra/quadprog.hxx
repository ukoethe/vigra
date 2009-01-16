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
/*    OF MERCHANTABILITY, FITNESS FOR activeSet PARTICULAR PURPOSE AND          */
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
bool quadprogAddConstraint(MultiArrayView<2, T, C1> & R, MultiArrayView<2, T, C2> & J, MultiArrayView<2, T, C3> & d, int& activeConstraintCount, double& R_norm)
{
    typedef typename MultiArrayShape<2>::type Shape;
    int n=rowCount(J);
    double cc, ss, h, xny;
    
    /*  we have to find the Givens rotation which will reduce the element
        d(j) to zero.
        if it is already zero we don'step have to do anything, except of
        decreasing j */  
    for (int j = n - 1; j > activeConstraintCount; --j)
    {
        /*   The Givens rotation is done with the matrix (cc cs, cs -cc).
             If cc is one, then element (j) of d is zero compared with element
             (j - 1). Hence we don'step have to do anything. 
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
        for (int k = 0; k < n; k++)
        {
            double t1 = J(k,j - 1);
            double t2 = J(k,j);
            J(k,j - 1) = t1 * cc + t2 * ss;
            J(k,j) = xny * (t1 + J(k,j - 1)) - t2;
        }
    }
    /* update the number of constraints added*/
    activeConstraintCount++;
    /* To update R we have to put the activeConstraintCount components of the d vector
    into column activeConstraintCount - 1 of R
    */
    columnVector(R, Shape(0, activeConstraintCount - 1), activeConstraintCount) = subVector(d, 0, activeConstraintCount);
  
    if (abs(d(activeConstraintCount - 1,0)) <= NumericTraits<T>::epsilon() * R_norm)
        // problem degenerate
        return false;
    R_norm = std::max<T>(R_norm, abs(d(activeConstraintCount - 1,0)));
    return true;
}

template <class T, class C1, class C2, class C3>
void quadprogDeleteConstraint(MultiArrayView<2, T, C1> & R, MultiArrayView<2, T, C2> & J, ArrayVector<int> & activeSet, MultiArrayView<2, T, C3> & u,  int me, int& activeConstraintCount, int l)
{
    typedef typename MultiArrayShape<2>::type Shape;
    int n=rowCount(R);
    int qq;
    double cc, ss, h, xny, dualStep, primalStep;
  
    /* Find the index qq for active constraint l to be removed */
    for (int i = me; i < activeConstraintCount; i++)
        if (activeSet[i] == l)
        {
            qq = i;
            break;
        }
      
    /* remove the constraint from the active set and the duals */
    for (int i = qq; i < activeConstraintCount - 1; i++)
    {
        activeSet[i] = activeSet[i + 1];
        u(i,0) = u(i + 1,0);
        for (int j = 0; j < n; j++)
            R(j,i) = R(j,i + 1);
    }

    activeSet[activeConstraintCount - 1] = activeSet[activeConstraintCount];
    u(activeConstraintCount - 1,0) = u(activeConstraintCount,0);
    activeSet[activeConstraintCount] = 0; 
    u(activeConstraintCount,0) = 0.0;
    columnVector(R, Shape(0, activeConstraintCount - 1), activeConstraintCount).init(NumericTraits<T>::zero());

    /* constraint has been fully removed */
    activeConstraintCount--;
  
    if (activeConstraintCount == 0)
        return;

    for (int j = qq; j < activeConstraintCount; j++)
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
        for (int k = j + 1; k < activeConstraintCount; k++)
        {
            dualStep = R(j,k);
            primalStep = R(j + 1,k);
            R(j,k) = dualStep * cc + primalStep * ss;
            R(j + 1,k) = xny * (dualStep + R(j,k)) - primalStep;
        }
        for (int k = 0; k < n; k++)
        {
            dualStep = J(k,j);
            primalStep = J(k,j + 1);
            J(k,j) = dualStep * cc + primalStep * ss;
            J(k,j + 1) = xny * (J(k,j) + dualStep) - primalStep;
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
     s.step.
         CE x == ce
         CI x >= ci
                  
     The matrix and vectors dimensions are as follows:
         G: n * n
            g: n
                    
            CE: me * n
         ce: me
                    
          CI: mi * n
       ci: mi

         x: n
     
     The function will return the cost of the solution written in the x vector or
     std::numeric_limits::infinity() if the problem is infeasible. In the latter case
     the value of the x vector is not correct.
     
     References: D. Goldfarb, activeSet. Idnani. activeSet numerically stable dual method for solving
                 strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

     Implementation based on work by Luca Di Gaspero and Angelo Furfaro.
   */
template <class T, class C1, class C2, class C3, class C4, class C5, class C6, class C7>
T 
quadraticProgramming(MultiArrayView<2, T, C1> const & G, MultiArrayView<2, T, C2> const & g,  
               MultiArrayView<2, T, C3> const & CE, MultiArrayView<2, T, C4> const & ce,  
               MultiArrayView<2, T, C5> const & CI, MultiArrayView<2, T, C6> const & ci, 
               MultiArrayView<2, T, C7> & x)
{
    using namespace linalg;
    using namespace functor;
    typedef typename MultiArrayShape<2>::type Shape;
    
    int n  = rowCount(g),
        me = rowCount(ce),
        mi = rowCount(ci),
        constraintCount = me + mi;
    vigra_precondition(columnCount(G) == n && rowCount(G) == n,
        "quadraticProgramming(): Matrix shape mismatch between G and g.");
    vigra_precondition(rowCount(x) == n,
        "quadraticProgramming(): Output vector x has illegal shape.");
    vigra_precondition((me > 0 && columnCount(CE) == n && rowCount(CE) == me) || 
                       (me == 0 && columnCount(CE) == 0),
        "quadraticProgramming(): Matrix CE has illegal shape.");
    vigra_precondition((mi > 0 && columnCount(CI) == n && rowCount(CI) == mi) || 
                       (mi == 0 && columnCount(CI) == 0),
        "quadraticProgramming(): Matrix CI has illegal shape.");

    Matrix<T> J = identityMatrix<T>(n);
    {
        Matrix<T> L(G.shape());
        choleskyDecomposition(G, L);
        // find unconstrained minimizer of the quadratic form  0.5 * x G x + g' x
        choleskySolve(L, -g, x);
        // compute the inverse of the factorized matrix G^-1, this is the initial value for J
        linearSolveLowerTriangular(L, transpose(J), transpose(J));
    }
    // current solution value
    T f_value = 0.5 * dot(g, x);
    
    T epsilonZ   = NumericTraits<T>::epsilon() * sq(J.norm(0)),
      epsilonPsi = NumericTraits<T>::epsilon() * trace(G)*trace(J)*100.0,
      inf        = std::numeric_limits<T>::infinity();
    
    Matrix<T> R(n, n), r(constraintCount, 1), u(constraintCount,1);
    T R_norm = NumericTraits<T>::one();
    
    ArrayVector<int> activeSet(constraintCount), A_old(constraintCount), iai(constraintCount);
    for (int i = 0; i < mi; i++)
        iai[i] = i;
    
    int activeConstraintCount = 0;
    // add equality constraints to the active set
    for (MultiArrayIndex i = 0; i < me; i++)
    {
        MultiArrayView<2, T, C3> np = rowVector(CE, i);
        Matrix<T> d = transpose(np*J);
        Matrix<T> z = J.subarray(Shape(0, activeConstraintCount), Shape(n,n))*subVector(d, activeConstraintCount, n);
        linearSolveUpperTriangular(R.subarray(Shape(0, 0), Shape(activeConstraintCount,activeConstraintCount)), 
                                   subVector(d, 0, activeConstraintCount), 
                                   subVector(r, 0, activeConstraintCount));
       /* compute full step length primalStep: i.e., the minimum step in primal space s.step. the contraint 
           becomes feasible */
        T primalStep = 0.0;
        if (squaredNorm(z) > epsilonZ) // z != 0
            primalStep = (-dot(np, x) + ce(i,0)) / dot(z, np);
    
        x += primalStep * z;    
        u(activeConstraintCount,0) = primalStep;
        subVector(u, 0, activeConstraintCount) -= primalStep * subVector(r, 0, activeConstraintCount);
        
        f_value += 0.5 * sq(primalStep) * dot(z, np);
        activeSet[i] = -i - 1;
    
        vigra_precondition(vigra::detail::quadprogAddConstraint(R, J, d, activeConstraintCount, R_norm),
            "quadraticProgramming(): Equality constraints are linearly dependent.");
    }
  
    /* Step 2: check for feasibility and determine a new S-pair */
    // compute values of inequality constraints
    Matrix<T> s = CI * x - ci;
    int constraintToBeAdded = argMin(s);
    T ss = s(constraintToBeAdded, 0);

    int iter = 0, maxIter = 10*mi;
    while(iter++ < maxIter)
    {        
        if (ss >= 0.0)       // all constraints are satisfied
            return f_value;  // => solved!

        MultiArrayView<2, T, C5> np = rowVector(CI, constraintToBeAdded);
        u(activeConstraintCount,0) = 0.0;
        activeSet[activeConstraintCount] = constraintToBeAdded;

        /* Step 2a: determine step direction */
        /* compute step direction in the primal space (through J, see the paper) */
        Matrix<T> d = transpose(np*J);
        Matrix<T> z = J.subarray(Shape(0, activeConstraintCount), Shape(n,n))*subVector(d, activeConstraintCount, n);
        /* compute negative of the step direction in the dual space */
        linearSolveUpperTriangular(R.subarray(Shape(0, 0), Shape(activeConstraintCount,activeConstraintCount)), 
                                   subVector(d, 0, activeConstraintCount), 
                                   subVector(r, 0, activeConstraintCount));

        /* Step 2b: compute step length */
        /* Compute primalStep (minimum step in primal space 
           such that constraintToBeAdded becomes feasible */
        T primalStep = (squaredNorm(z) <= epsilonZ) // i.e. z == 0
                          ? inf
                          : -s(constraintToBeAdded,0) / dot(z, np);
      
        /* Compute dualStep (maximum step in dual space without violating dual feasibility) */
        T dualStep = inf; 
        int constraintToBeRemoved;
        /* find the index that achieves the minimum of u+(x) / r */
        for (int k = me; k < activeConstraintCount; k++)
        {
            if (r(k,0) > 0.0)
            {
                if (u(k,0) / r(k,0) < dualStep)
                {
                    dualStep = u(k,0) / r(k,0);
                    constraintToBeRemoved = activeSet[k];
                }
            }
        }
        
        /* the step is chosen as the minimum of dualStep and primalStep */
        T step = std::min(dualStep, primalStep);
      
        /* Step 2c: determine new S-pair and take step: */
      
        // case (i): no step in primal or dual space possible
        if (step >= inf) // QPP is infeasible 
            return inf;

        /* case (ii): step in dual space */
        if (primalStep >= inf)
        {
            /* set u = u +  step * [-r 1) and drop constraintToBeRemoved from the active set */
            subVector(u, 0, activeConstraintCount) -= step * subVector(r, 0, activeConstraintCount);
            u(activeConstraintCount,0) += step;
            iai[constraintToBeRemoved] = constraintToBeRemoved;
            vigra::detail::quadprogDeleteConstraint(R, J, activeSet, u, me, activeConstraintCount, constraintToBeRemoved);
            continue;
        }
      
        /* case (iii): step in primal and dual space */
      
        x += step * z;
        /* update the solution value */
        f_value += step * dot(z, np) * (0.5 * step + u(activeConstraintCount,0));
        /* u = u + step * (-r 1) */
        subVector(u, 0, activeConstraintCount) -= step * subVector(r, 0, activeConstraintCount);
        u(activeConstraintCount,0) += step;
      
        if (step == primalStep)
        {
            /* add constraintToBeAdded to the active set*/
            vigra::detail::quadprogAddConstraint(R, J, d, activeConstraintCount, R_norm);
            iai[constraintToBeAdded] = -1;
        }
        else
        {
            /* drop constraintToBeRemoved from the active set */
            iai[constraintToBeRemoved] = constraintToBeRemoved;
            vigra::detail::quadprogDeleteConstraint(R, J, activeSet, u, me, activeConstraintCount, constraintToBeRemoved);
        }
        
        for (int i = me; i < activeConstraintCount; i++)
            iai[activeSet[i]] = -1;

        // compute values of inequality constraints
        s = CI * x - ci;
        ss = 0.0;
        constraintToBeAdded;
        /* Step 2: check for feasibility and determine a new S-pair */
        for (int i = 0; i < mi; i++)
        {
            if (s(i,0) < ss && iai[i] != -1)
            {
                ss = s(i,0);
                constraintToBeAdded = i;
            }
        }
    }
    return inf;
}

//@}

} // namespace vigra

#endif
