/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
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
/*    OF MERCHANTABILITY, FITNESS FOR activeSet PARTICULAR PURPOSE AND  */
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
bool quadprogAddConstraint(MultiArrayView<2, T, C1> & R, MultiArrayView<2, T, C2> & J, MultiArrayView<2, T, C3> & d, 
                           int activeConstraintCount, double& R_norm)
{
    typedef typename MultiArrayShape<2>::type Shape;
    int n=columnCount(J);
    linalg::detail::qrGivensStepImpl(0, subVector(d, activeConstraintCount, n),
                                     J.subarray(Shape(activeConstraintCount,0), Shape(n,n)));
    if (abs(d(activeConstraintCount,0)) <= NumericTraits<T>::epsilon() * R_norm) // problem degenerate
        return false;
    R_norm = std::max<T>(R_norm, abs(d(activeConstraintCount,0)));

    ++activeConstraintCount;   
    // add d as a new column to R
    columnVector(R, Shape(0, activeConstraintCount - 1), activeConstraintCount) = subVector(d, 0, activeConstraintCount);  
    return true;
}

template <class T, class C1, class C2, class C3>
void quadprogDeleteConstraint(MultiArrayView<2, T, C1> & R, MultiArrayView<2, T, C2> & J, MultiArrayView<2, T, C3> & u, 
                              int activeConstraintCount,  int constraintToBeRemoved)
{
    typedef typename MultiArrayShape<2>::type Shape;
    
    int newActiveConstraintCount = activeConstraintCount - 1;

    if(constraintToBeRemoved == newActiveConstraintCount)
        return;

    std::swap(u(constraintToBeRemoved,0), u(newActiveConstraintCount,0));
    columnVector(R, constraintToBeRemoved).swapData(columnVector(R, newActiveConstraintCount));
    linalg::detail::qrGivensStepImpl(0, R.subarray(Shape(constraintToBeRemoved, constraintToBeRemoved), 
                                                   Shape(newActiveConstraintCount,newActiveConstraintCount)),
                                        J.subarray(Shape(constraintToBeRemoved, 0), 
                                                   Shape(newActiveConstraintCount,newActiveConstraintCount)));
}

} // namespace detail

/** \addtogroup Optimization Optimization and Regression
 */
//@{
   /** Solve Quadratic Programming Problem.

     The quadraticProgramming() function implements the algorithm described in
     
     D. Goldfarb, A. Idnani: <i>"A numerically stable dual method for solving
                 strictly convex quadratic programs"</i>, Mathematical Programming 27:1-33, 1983. 
     
     for the solution of (convex) quadratic programming problems by means of a primal-dual method.
         
     <b>\#include</b> \<<a href="quadprog_8hxx-source.html">vigra/quadprog.hxx</a>\>
         Namespaces: vigra

     <b>Declaration:</b>

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

     The problem must be specified in the form:

     \f{eqnarray*}
        \mbox{minimize } &\,& \frac{1}{2} \mbox{\bf x}'\,\mbox{\bf G}\, \mbox{\bf x} + \mbox{\bf g}'\,\mbox{\bf x} \\
        \mbox{subject to} &\,& \mbox{\bf C}_E\, \mbox{\bf x} = \mbox{\bf c}_e \\
         &\,& \mbox{\bf C}_I\,\mbox{\bf x} \ge \mbox{\bf c}_i
     \f}            
     Matrix <b>G</b> G must be symmetric positive definite, and matrix <b>C</b><sub>E</sub> must have full row rank. 
     Matrix and vector dimensions must be as follows:
     <ul>
     <li> <b>G</b>: [n * n], <b>g</b>: [n * 1]
     <li> <b>C</b><sub>E</sub>: [me * n], <b>c</b><sub>e</sub>: [me * 1]
     <li> <b>C</b><sub>I</sub>: [mi * n], <b>c</b><sub>i</sub>: [mi * 1]
     <li> <b>x</b>: [n * 1]
     </ul>
     
     The function writes the optimal solution into the vector \a x and returns the cost of this solution. 
     If the problem is infeasible, std::numeric_limits::infinity() is returned. In this case
     the value of vector \a x is undefined.
     
     <b>Usage:</b>
     
     Minimize <tt> f = 0.5 * x'*G*x + g'*x </tt> subject to <tt> -1 &lt;= x &lt;= 1</tt>. 
     The solution is <tt> x' = [1.0, 0.5, -1.0] </tt> with <tt> f = -22.625</tt>.
     \code
      double Gdata[] = {13.0, 12.0, -2.0,
                        12.0, 17.0,  6.0,
                        -2.0,  6.0, 12.0};

      double gdata[] = {-22.0, -14.5, 13.0};

      double CIdata[] = { 1.0,  0.0,  0.0,
                          0.0,  1.0,  0.0,
                          0.0,  0.0,  1.0,
                         -1.0,  0.0,  0.0,
                          0.0, -1.0,  0.0,
                          0.0,  0.0, -1.0};
                        
      double cidata[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};

      Matrix<double> G(3,3, Gdata), 
                     g(3,1, gdata), 
                     CE,             // empty since there are no equality constraints
                     ce,             // likewise
                     CI(7,3, CIdata), 
                     ci(7,1, cidata), 
                     x(3,1);
                   
      double f = quadraticProgramming(G, g, CE, ce, CI, ci, x);
     \endcode
   */
template <class T, class C1, class C2, class C3, class C4, class C5, class C6, class C7>
T 
quadraticProgramming(MultiArrayView<2, T, C1> const & G, MultiArrayView<2, T, C2> const & g,  
               MultiArrayView<2, T, C3> const & CE, MultiArrayView<2, T, C4> const & ce,  
               MultiArrayView<2, T, C5> const & CI, MultiArrayView<2, T, C6> const & ci, 
               MultiArrayView<2, T, C7> & x)
{
    using namespace linalg;
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
        linearSolveLowerTriangular(L, J, J);
    }
    // current solution value
    T f_value = 0.5 * dot(g, x);
    
    T epsilonZ   = NumericTraits<T>::epsilon() * sq(J.norm(0)),
      epsilonPsi = NumericTraits<T>::epsilon() * trace(G)*trace(J)*100.0,
      inf        = std::numeric_limits<T>::infinity();
    
    Matrix<T> R(n, n), r(constraintCount, 1), u(constraintCount,1);
    T R_norm = NumericTraits<T>::one();
    
    // incorporate equality constraints
    for (int i=0; i < me; ++i)
    {
        MultiArrayView<2, T, C3> np = rowVector(CE, i);
        Matrix<T> d = J*transpose(np);
        Matrix<T> z = transpose(J).subarray(Shape(0, i), Shape(n,n))*subVector(d, i, n);
        linearSolveUpperTriangular(R.subarray(Shape(0, 0), Shape(i,i)), 
                                   subVector(d, 0, i), 
                                   subVector(r, 0, i));
        // compute step in primal space so that the constraint becomes satisfied
        T step = (squaredNorm(z) <= epsilonZ) // i.e. z == 0
                     ? 0.0 
                     : (-dot(np, x) + ce(i,0)) / dot(z, np);
    
        x += step * z;    
        u(i,0) = step;
        subVector(u, 0, i) -= step * subVector(r, 0, i);
        
        f_value += 0.5 * sq(step) * dot(z, np);
    
        vigra_precondition(vigra::detail::quadprogAddConstraint(R, J, d, i, R_norm),
            "quadraticProgramming(): Equality constraints are linearly dependent.");
    }
    int activeConstraintCount = me;
  
    // determine optimum solution and corresponding active inequality constraints
    ArrayVector<int> activeSet(mi);
    for (int i = 0; i < mi; ++i)
        activeSet[i] = i;

    int constraintToBeAdded;
    T ss = 0.0;
    for (int i = activeConstraintCount-me; i < mi; ++i)
    {
        T s = dot(rowVector(CI, activeSet[i]), x) - ci(activeSet[i], 0);
        if (s < ss)
        {
            ss = s;
            constraintToBeAdded = i;
        }
    }

    int iter = 0, maxIter = 10*mi;    
    while(iter++ < maxIter)
    {        
        if (ss >= 0.0)       // all constraints are satisfied
            return f_value;  // => solved!

        // determine step direction in the primal space (through J, see the paper)
        MultiArrayView<2, T, C5> np = rowVector(CI, activeSet[constraintToBeAdded]);
        Matrix<T> d = J*transpose(np);
        Matrix<T> z = transpose(J).subarray(Shape(0, activeConstraintCount), Shape(n,n))*subVector(d, activeConstraintCount, n);
        
        // compute negative of the step direction in the dual space
        linearSolveUpperTriangular(R.subarray(Shape(0, 0), Shape(activeConstraintCount,activeConstraintCount)), 
                                   subVector(d, 0, activeConstraintCount), 
                                   subVector(r, 0, activeConstraintCount));

        // determine minimum step length in primal space such that activeSet[constraintToBeAdded] becomes feasible
        T primalStep = (squaredNorm(z) <= epsilonZ) // i.e. z == 0
                          ? inf
                          : -ss / dot(z, np);
      
        // determine maximum step length in dual space that doesn't violate dual feasibility
        // and the corresponding index
        T dualStep = inf; 
        int constraintToBeRemoved;
        for (int k = me; k < activeConstraintCount; ++k)
        {
            if (r(k,0) > 0.0)
            {
                if (u(k,0) / r(k,0) < dualStep)
                {
                    dualStep = u(k,0) / r(k,0);
                    constraintToBeRemoved = k;
                }
            }
        }
        
        // the step is chosen as the minimum of dualStep and primalStep
        T step = std::min(dualStep, primalStep);
      
        // take step and update matrizes
      
        if (step == inf)
        {
            // case (i): no step in primal or dual space possible
            return inf; // QPP is infeasible 
        }
        if (primalStep == inf)
        {
            // case (ii): step in dual space
            subVector(u, 0, activeConstraintCount) -= step * subVector(r, 0, activeConstraintCount);
            vigra::detail::quadprogDeleteConstraint(R, J, u, activeConstraintCount, constraintToBeRemoved);
            --activeConstraintCount;
            std::swap(activeSet[constraintToBeRemoved-me], activeSet[activeConstraintCount-me]);
            continue;
        }
      
        // case (iii): step in primal and dual space      
        x += step * z;
        // update the solution value
        f_value += 0.5 * sq(step) * dot(z, np);
        // u = [u 1]' + step * [-r 1]
        subVector(u, 0, activeConstraintCount) -= step * subVector(r, 0, activeConstraintCount);
        u(activeConstraintCount,0) = step;
      
        if (step == primalStep)
        {
            // add constraintToBeAdded to the active set
            vigra::detail::quadprogAddConstraint(R, J, d, activeConstraintCount, R_norm);
            std::swap(activeSet[constraintToBeAdded], activeSet[activeConstraintCount-me]);
            ++activeConstraintCount;
        }
        else
        {
            // drop constraintToBeRemoved from the active set
            vigra::detail::quadprogDeleteConstraint(R, J, u, activeConstraintCount, constraintToBeRemoved);
            --activeConstraintCount;
            std::swap(activeSet[constraintToBeRemoved-me], activeSet[activeConstraintCount-me]);
        }
        
        // update values of inactive inequality constraints
        ss = 0.0;
        for (int i = activeConstraintCount-me; i < mi; ++i)
        {
            // compute CI*x - ci with appropriate row permutation
            T s = dot(rowVector(CI, activeSet[i]), x) - ci(activeSet[i], 0);
            if (s < ss)
            {
                ss = s;
                constraintToBeAdded = i;
            }
        }
    }
    return inf; // too many iterations
}

//@}

} // namespace vigra

#endif
