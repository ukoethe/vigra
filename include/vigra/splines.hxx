/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_SPLINES_HXX
#define VIGRA_SPLINES_HXX

#include <cmath>
#include "vigra/config.hxx"
#include "vigra/mathutil.hxx"

namespace vigra {

/** \addtogroup MathFunctions Mathematical Functions
*/
//@{
/*! 

    <b>\#include</b> "<a href="splines_8hxx-source.html">vigra/splines.hxx</a>"<br>
    Namespace: vigra
*/
template <int ORDER, class T = double>
class BSpline
{
  public:
  
    typedef T            value_type;  
    typedef T            argument_type;  
    typedef T            first_argument_type;  
    typedef unsigned int second_argument_type;  
    typedef T            result_type; 
    enum               { order = ORDER };

    result_type operator()(argument_type x) const
    {
        T n12 = (ORDER + 1.0) / 2.0;
        return ((n12 + x) * s1_(x + 0.5) + (n21 - x) * s1_(x - 0.5)) / ORDER;
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
        if(derivative_order == 0)
        {
            return operator()(x);
        }
        else
        {
            --derivative_order;
            return s1_(x + 0.5, derivative_order) - s1_(x - 0.5, derivative_order);
        }
    }
    
    value_type operator[](value_type x) const
        { return operator()(x); }
    
    double radius() const
        { return (ORDER + 1) * 0.5; }

  private:
    BSpline<ORDER-1, T> s1_;  
};

/********************************************************/
/*                                                      */
/*                     BSpline<0, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<0, T>
{
  public:
  
    typedef T            value_type;  
    typedef T            argument_type;  
    typedef T            first_argument_type;  
    typedef unsigned int second_argument_type;  
    typedef T            result_type; 
    enum               { order = 0 };

    result_type operator()(argument_type x) const
    {
         return VIGRA_CSTD::fabs(x) <= 0.5 
                       ? 1.0
                       : 0.0;
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
        return 0.0;
    }
    
    value_type operator[](value_type x) const
        { return operator()(x); }
    
    double radius() const
        { return 0.5; }
};

/********************************************************/
/*                                                      */
/*                     BSpline<1, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<1, T>
{
  public:
  
    typedef T            value_type;  
    typedef T            argument_type;  
    typedef T            first_argument_type;  
    typedef unsigned int second_argument_type;  
    typedef T            result_type; 
    enum               { order = 1 };

    result_type operator()(argument_type x) const
    {
         x = VIGRA_CSTD::fabs(x);
         return x < 1.0 
                       ? 1.0 - x
                       : 0.0;
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const;
    
    value_type operator[](value_type x) const
        { return operator()(x); }
    
    double radius() const
        { return 1.0; }
};

template <class T>
T BSpline<1, T>::operator()(T x, unsigned int derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
            return operator()(x);
        case 1:
            return x >= 0.0 ?
                     x <  1.0 ?
                         -1.0 :
                          0.0 :
                     x > -1.0 ?
                          1.0 :
                          0.0;
        default:
            return 0.0;
    }
}

/********************************************************/
/*                                                      */
/*                     BSpline<2, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<2, T>
{
  public:
  
    typedef T            value_type;  
    typedef T            argument_type;  
    typedef T            first_argument_type;  
    typedef unsigned int second_argument_type;  
    typedef T            result_type; 
    enum               { order = 2 };

    result_type operator()(argument_type x) const;
    
    result_type operator()(first_argument_type x, second_argument_type derivative_order) const;
    
    value_type operator[](value_type x) const
        { return operator()(x); }
    
    double radius() const
        { return 1.5; }
};

template <class T>
T BSpline<2, T>::operator()(T x) const
{
     x = VIGRA_CSTD::fabs(x);
     return x < 0.5 ?
                0.75 - x*x :
                x < 1.5 ?
                    0.5 * sq(1.5 - x) :
                    0.0;
}

template <class T>
T BSpline<2, T>::operator()(T x, unsigned int derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
            return operator()(x);
        case 1:
            return x >= -0.5 ?
                     x <= 0.5 ?
                         -2.0 * x :
                          x < 1.5 ?
                             x - 1.5 :
                             0.0
                     x > -1.5 ?
                         x + 1.5 : 
                         0.0;
        case 2:
            return x >= -0.5 ?
                     x <= 0.5 ?
                         -2.0 :
                          x < 1.5 ?
                             1.0 :
                             0.0
                     x > -1.5 ?
                         1.0 : 
                         0.0;
        default:
            return 0.0;
    }
}

/********************************************************/
/*                                                      */
/*                     BSpline<3, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<3, T>
{
  public:
  
    typedef T            value_type;  
    typedef T            argument_type;  
    typedef T            first_argument_type;  
    typedef unsigned int second_argument_type;  
    typedef T            result_type; 
    enum               { order = 3 };

    result_type operator()(argument_type x) const;
    
    result_type operator()(first_argument_type x, second_argument_type derivative_order) const;
    
    result_type dx(argument_type x) const
        { return operator()(x, 1); }
    
    result_type dxx(argument_type x) const
        { return operator()(x, 2); }
    
    value_type operator[](value_type x) const
        { return operator()(x); }
    
    double radius() const
        { return 2.0; }
};

template <class T>
T BSpline<3, T>::operator()(T x) const
{
    x = VIGRA_CSTD::fabs(x);
    if(x < 1.0)
    {
        return 2.0/3.0 + x*x*(-1.0 + 0.5*x);
    }
    else if(x < 2.0)
    {
        x = 2.0 - x;
        return x*x*x/6.0;
    }
    else
        return 0.0;
}

template <class T>
T BSpline<3, T>::operator()(T x, unsigned int derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
            return operator()(x);
        case 1:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                     s*x*(-2.0 + 1.5*x) :
                     x < 2.0 ?
                       -0.5*s*sq(2.0 - x) :
                        0.0;
        }
        case 2:
        {
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                     3.0*x - 2.0 :
                     x < 2.0 ?
                       2.0 - x :
                       0.0:
        }
        case 3:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                       3.0*s :
                       x < 2.0 ?
                         -s :
                         0.0; 
        }
        default:
            return 0.0;
    }
}

typedef BSpline<3, double> CubicBSplineKernel;

/********************************************************/
/*                                                      */
/*                     BSpline<5, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<5, T>
{
  public:
  
    typedef T            value_type;  
    typedef T            argument_type;  
    typedef T            first_argument_type;  
    typedef unsigned int second_argument_type;  
    typedef T            result_type; 
    enum               { order = 5 };

    result_type operator()(argument_type x) const;
    
    result_type operator()(first_argument_type x, second_argument_type derivative_order) const;
    
    result_type dx(argument_type x) const
        { return operator()(x, 1); }
    
    result_type dxx(argument_type x) const
        { return operator()(x, 2); }
    
    result_type dx3(argument_type x) const
        { return operator()(x, 3); }
    
    result_type dx4(argument_type x) const
        { return operator()(x, 4); }
    
    value_type operator[](value_type x) const
        { return operator()(x); }
    
    double radius() const
        { return 3.0; }
};

template <class T>
T BSpline<5, T>::operator()(T x) const
{
    x = VIGRA_CSTD::fabs(x);
    if(x <= 1.0)
    {
        return 0.55 + x*x*(-0.5 + x*x*(0.25 - x/12.0));
    }
    else if(x < 2.0)
    {
        return 17.0/40.0 + x*(0.625 + x*(-1.75 + x*(1.25 + x*(-0.375 + x/24.0))));
    }
    else if(x < 3.0)
    {
        x = 3.0 - x;
        return x*sq(x*x) / 120.0;
    }
    else
        return 0.0;
}

template <class T>
T BSpline<5, T>::operator()(T x, unsigned int derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
            return operator()(x);
        case 1:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return s*x*(-1.0 + x*x*(1.0 - 5.0/12.0*x));
            }
            else if(x < 2.0)
            {
                return s*(0.625 + x*(-3.5 + x*(3.75 + x*(-1.5 + 5.0/24.0*x))));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return s*sq(x*x) / -24.0;
            }
            else
                return 0.0;
        }
        case 2:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return -1.0 + x*x*(3.0 -5.0/3.0*x);
            }
            else if(x < 2.0)
            {
                return -3.5 + x*(7.5 + x*(-4.5 + 5.0/6.0*x));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return x*x*x / 6.0;
            }
            else
                return 0.0;
        }
        case 3:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return s*x*(6.0 - 5.0*x);
            }
            else if(x < 2.0)
            {
                return s*(7.5 + x*(-9.0 + 2.5*x));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return -0.5*s*x*x;
            }
            else
                return 0.0;
        }
        case 4:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return 6.0 - 10.0*x;
            }
            else if(x < 2.0)
            {
                return -9.0 + 5.0*x;
            }
            else if(x < 3.0)
            {
                return 3.0 - x;
            }
            else
                return 0.0;
        }
        case 5:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                       -10.0*s :
                        x < 2.0 ?
                          5.0*s :
                          x < 3.0 ?
                            -s :
                            0.0; 
        }
        default:
            return 0.0;
    }
}

typedef BSpline<5, double> QuinticBSplineKernel;

/********************************************************/
/*                                                      */
/*                      CatmullRomSpline                */
/*                                                      */
/********************************************************/

template <class T = double>
class CatmullRomSpline
{
public:
    typedef T value_type;
    typedef T argument_type;
    typedef T result_type;

    result_type operator()(argument_type x) const;
    
    T operator[] (T x) const
        { return operator()(x); }

    int radius() const
        {return 2;}
};


template <class T>
T CatmullRomSpline<T>::operator()(T x) const
{
    x = VIGRA_CSTD::fabs(x);
    if (x <= 1.0)
    {
        return 1.0 + x * x * (-2.5 + 1.5 * x);
    }
    else if (x >= 2.0)
    {
        return 0.0;
    }
    else
    {
        return 2.0 + x * (-4.0 + x * (2.5 - 0.5 * x));
    }
}


//@}

} // namespace vigra


#endif /* VIGRA_SPLINES_HXX */
