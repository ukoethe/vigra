/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#include <iostream>
#include <algorithm>
#include <cmath>
#include "unittest.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/rgbvalue.hxx"
#include "vigra/sized_int.hxx"

using namespace vigra::functor;

template <class EXPR>
typename ResultTraits0<UnaryFunctor<EXPR> >::Res 
exec(UnaryFunctor<EXPR> f) { return f(); } 

template <class EXPR, class T1>
typename ResultTraits1<UnaryFunctor<EXPR>, T1>::Res 
exec(UnaryFunctor<EXPR> f, T1 const & a1) { return f(a1); } 

template <class EXPR, class T1, class T2>
typename ResultTraits2<UnaryFunctor<EXPR>, T1, T2>::Res 
exec(UnaryFunctor<EXPR> f, T1 const & a1, T2 const & a2) { return f(a1, a2); } 

template <class EXPR, class T1, class T2, class T3>
typename ResultTraits3<UnaryFunctor<EXPR>, T1, T2, T3>::Res 
exec(UnaryFunctor<EXPR> f, T1 const & a1, T2 const & a2, T3 const & a3) 
{ return f(a1, a2, a3); } 

template <class EXPR, class T1>
void 
exec(UnaryAnalyser<EXPR> f, T1 const & a1) { f(a1); } 

template <class EXPR, class T1, class T2>
void 
exec(UnaryAnalyser<EXPR> f, T1 const & a1, T2 const & a2) { f(a1, a2); } 

template <class EXPR, class T1, class T2, class T3>
void 
exec(UnaryAnalyser<EXPR> f, T1 const & a1, T2 const & a2, T3 const & a3) 
{ f(a1, a2, a3); } 

struct FunctorExpressionTest
{    
    void testArg1()
    {
        should(exec(Arg1(), 1.5) == 1.5); 
        should(exec(Arg1(), 1.5, 2) == 1.5); 
        should(exec(Arg1(), 1.5, 2, true) == 1.5); 
    }
    
    void testArg2()
    {
        should(exec(Arg2(), 1, 2.5) == 2.5); 
        should(exec(Arg2(), 1, 2.5, false) == 2.5); 
    }
    
    void testArg3()
    {
        should(exec(Arg3(), 1, true, 3.5) == 3.5); 
    }

    void testParam()
    {
        should(exec(Param(4.5), 1) == 4.5); 
        should(exec(Param(4.5), 1, true) == 4.5); 
        should(exec(Param(4.5), 1, true, 3) == 4.5); 
    }
    
    void testVar()
    {
        double v = 4.5;
        should(exec(Var(v), 1) == 4.5); 
        should(exec(Var(v), 1, true) == 4.5); 
        should(exec(Var(v), 1, true, 3) == 4.5); 
    }
    
    void testAssignment()
    {
        double v = 0.0;
        exec(Var(v) = Arg1(), 1.5);
        should(v == 1.5);
        exec(Var(v) = Arg2(), 1, 2.5);
        should(v == 2.5);
        exec(Var(v) = Arg3(), 1, 2, 3.5);
        should(v == 3.5);
        exec(Var(v) += Arg1(), 2.5);
        should(v == 6.0);
        exec(Var(v) -= Arg1(), 2.5);
        should(v == 3.5);
        exec(Var(v) *= Arg1(), 2.0);
        should(v == 7.0);
        exec(Var(v) /= Arg1(), 2.0);
        should(v == 3.5);
    }
    
    void testIfThen()
    {
        double v = 4.0;
        exec(ifThen(Arg1(), Var(v) = Param(1.5)), false);
        should(v == 4.0);
        exec(ifThen(Arg1(), Var(v) = Param(1.5)), true);
        should(v == 1.5);
        exec(ifThen(Arg2(), Var(v) = Arg1()), 2.5, false);
        should(v == 1.5);
        exec(ifThen(Arg2(), Var(v) = Arg1()), 2.5, true);
        should(v == 2.5);
        exec(ifThen(Arg3(), Var(v) = Arg2()), 1, 3.5, false);
        should(v == 2.5);
        exec(ifThen(Arg3(), Var(v) = Arg2()), 1, 3.5, true);
        should(v == 3.5);
    }
    
    void testIfThenElse()
    {
        should(exec(ifThenElse(Arg1(), Arg2(), Arg3()), true, 2, 3.5) == 2.0); 
        should(exec(ifThenElse(Arg1(), Arg2(), Arg3()), false, 2, 3.5) == 3.5); 
    }
    
    void testUnary()
    {
        shouldEqual(exec(sqrt(Arg1()), 7.0) , VIGRA_CSTD::sqrt(7.0)); 
        shouldEqual(exec(exp(Arg1()), 4.0) , VIGRA_CSTD::exp(4.0)); 
        shouldEqual(exec(log(Arg1()), 4.0) , VIGRA_CSTD::log(4.0)); 
        shouldEqual(exec(log10(Arg1()), 4.0) , VIGRA_CSTD::log10(4.0)); 
        shouldEqual(exec(sin(Arg1()), 4.0) , VIGRA_CSTD::sin(4.0)); 
        shouldEqual(exec(asin(Arg1()), 0.5) , VIGRA_CSTD::asin(0.5)); 
        shouldEqual(exec(cos(Arg1()), 4.0) , VIGRA_CSTD::cos(4.0)); 
        shouldEqual(exec(acos(Arg1()), 0.5) , VIGRA_CSTD::acos(0.5)); 
        shouldEqual(exec(tan(Arg1()), 4.0) , VIGRA_CSTD::tan(4.0)); 
        shouldEqual(exec(atan(Arg1()), 0.5) , VIGRA_CSTD::atan(0.5)); 
        shouldEqual(exec(abs(Arg1()), -0.5) , 0.5); 
//        shouldEqual(exec(norm(Arg1()), -0.5) , 0.5); 
//        shouldEqual(exec(squaredNorm(Arg1()), -0.5) , 0.25); 
        shouldEqual(exec(floor(Arg1()), -0.5) , -1.0); 
        shouldEqual(exec(ceil(Arg1()), 0.5) , 1.0); 
        shouldEqual(exec(-Arg1(), -0.5) , 0.5); 
        shouldEqual(exec(!Arg1(), true) , false); 
        shouldEqual(exec(~Arg1(),(vigra::UInt32)0xff) , (vigra::UInt32)0xffffff00); 
    }
    
    void testBinary()
    {
        should(exec(Arg1() + Arg2(), 1, 2.5) == 3.5); 
        should(exec(Arg1() - Arg2(), 1, 2.5) == -1.5); 
        should(exec(Arg1() * Arg2(), 1, 2.5) == 2.5); 
        should(exec(Arg1() / Arg2(), 1, 2.0) == 0.5); 
        should(exec(Arg1() == Arg2(), 1, 2.0) == false); 
        should(exec(Arg1() != Arg2(), 1, 2.0) == true); 
        should(exec(Arg1() < Arg2(), 1, 2.0) == true); 
        should(exec(Arg1() <= Arg2(), 1, 2.0) == true); 
        should(exec(Arg1() > Arg2(), 1, 2.0) == false); 
        should(exec(Arg1() >= Arg2(), 1, 2.0) == false); 
        should(exec(Arg1() && Arg2(), true, true) == true); 
        should(exec(Arg1() || Arg2(), true, false) == true); 
        should(exec(Arg1() & Arg2(), 0xff, 0xf) == 0xf); 
        should(exec(Arg1() | Arg2(), 0xff, 0xf) == 0xff); 
        should(exec(Arg1() ^ Arg2(), 0xff, 0xf) == 0xf0); 
        should(exec(pow(Arg1(), Arg2()), 2.0, 3) == VIGRA_CSTD::pow(2.0, 3.0)); 
        shouldEqualTolerance(exec(atan2(Arg1(), Arg2()), 2.0, 3.0), VIGRA_CSTD::atan2(2.0, 3.0), 1e-16); 
        should(exec(fmod(Arg1(), Arg2()), 2.0, 3.0) == VIGRA_CSTD::fmod(2.0, 3.0)); 
        should(exec(min(Arg1(), Arg2()), 2, 3.5) == 2.0); 
        should(exec(max(Arg1(), Arg2()), 2, 3.5) == 3.5); 
    }
    
    void testCombine()
    {
        should(exec(sqrt(Arg1() * Arg2() + Arg3()), 1.5, 2, 1) == 2.0);
    }
    
    void testComma()
    {
        double v1 = 0.0;
        int v2 = 0;
        bool v3 = false;
        
        exec((Var(v1) = Arg1(), Var(v2) = Arg2(), Var(v3) = Arg3()),
             1.5, 2, true);
        should(v1 == 1.5);
        should(v2 == 2);
        should(v3 == true);
        
        should(exec((Var(v1) = Arg3(), Arg2()), true, 2.5, 3.5) == 2.5);
        should(v1 == 3.5);
        
        double count = 0.0;
        should(exec((Var(count) += Param(1), Var(count))) == 1.0);
    }
    
    void testApply()
    {
        should(exec(applyFct((double (*)(double))&VIGRA_CSTD::sqrt, Arg1()), 4.0) == 2.0); 
        should(exec(applyFct((double (*)(double, double))&VIGRA_CSTD::pow, 
                             Arg1(), Arg2()), 4.0, 2) == 16.0); 
    }
    
    void testSequence()
    {
        unsigned char data[] = { 1, 2, 3, 4, 5, 250};
        int res[6];
        int sum = 0;
        int count = 0;
        
        std::for_each(data, data+6, 
        (
            Var(count) += Param(1),
            Var(sum) += Arg1()
        ));
        
        should(count == 6); 
        should(sum == 265); 
        
        std::transform(data, data+6, res,
        (
            Var(count) = Var(count) - Param(1),
            Param(2) * Arg1()
        ));
        
        should(count == 0); 
        
        for(int i=0; i<6; ++i) should(res[i] == 2*data[i]);
    }
};

struct FunctorRGBExpressionTest
{    
    vigra::RGBValue<double> v1_5, v2_5, v3_0, v_1_5;
    vigra::RGBValue<int> v0, v1, v2;
    
    FunctorRGBExpressionTest()
    : v1_5(1.5),
      v2_5(2.5),
      v3_0(3.0),
      v_1_5(-1.5),
      v0(0),
      v1(1),
      v2(2)
    {}
    
    void testArg1()
    {
        should(exec(Arg1(), v1_5) == v1_5); 
        should(exec(Arg1(), v1_5, v2) == v1_5); 
        should(exec(Arg1(), v1_5, v2, true) == v1_5); 
    }
    
    void testArg2()
    {
        should(exec(Arg2(), v1, v2_5) == v2_5); 
        should(exec(Arg2(), v1, v2_5, false) == v2_5);
    }
    
    void testArg3()
    {
        should(exec(Arg3(), 1, true, v2_5) == v2_5); 
    }

    void testParam()
    {
        should(exec(Param(v2_5), 1) == v2_5); 
        should(exec(Param(v2_5), 1, true) == v2_5); 
        should(exec(Param(v2_5), 1, true, 3) == v2_5); 
    }
    
    void testVar()
    {
        should(exec(Var(v2_5), 1) == v2_5); 
        should(exec(Var(v2_5), 1, true) == v2_5); 
        should(exec(Var(v2_5), 1, true, 3) == v2_5); 
    }
    
    void testAssignment()
    {
        vigra::RGBValue<double> v(0.0);
        exec(Var(v) = Arg1(), v1_5);
        should(v == v1_5);
        exec(Var(v) = Arg2(), 1, v2_5);
        should(v == v2_5);
        exec(Var(v) = Arg3(), 1, 2, v1_5);
        should(v == v1_5);
        exec(Var(v) += Arg1(), v1);
        should(v == v2_5);
        exec(Var(v) -= Arg1(), v1);
        should(v == v1_5);
        exec(Var(v) *= Arg1(), 2.0);
        should(v == v3_0);
        exec(Var(v) /= Arg1(), 2.0);
        should(v == v1_5);
        exec(Var(v) *= Arg1(), v2);
        should(v == v3_0);
    }
    
    void testIfThen()
    {
        vigra::RGBValue<double> v(v3_0);
        exec(ifThen(Arg1(), Var(v) = Param(v1_5)), false);
        should(v == v3_0);
        exec(ifThen(Arg1(), Var(v) = Param(v1_5)), true);
        should(v == v1_5);
        exec(ifThen(Arg2(), Var(v) = Arg1()), v2_5, false);
        should(v == v1_5);
        exec(ifThen(Arg2(), Var(v) = Arg1()), v2_5, true);
        should(v == v2_5);
        exec(ifThen(Arg3(), Var(v) = Arg2()), 1, v1_5, false);
        should(v == v2_5);
        exec(ifThen(Arg3(), Var(v) = Arg2()), 1, v1_5, true);
        should(v == v1_5);
    }
    
    void testIfThenElse()
    {
        should(exec(ifThenElse(Arg1(), Arg2(), Arg3()), true, v2, v1_5) == v2); 
        should(exec(ifThenElse(Arg1(), Arg2(), Arg3()), false, v2, v1_5) == v1_5); 
    }
    
    void testUnary()
    {
        should(exec(abs(Arg1()), v_1_5) == v1_5); 
        should(exec(floor(Arg1()), vigra::RGBValue<double>(1.4)) == v1); 
        should(exec(ceil(Arg1()), vigra::RGBValue<double>(1.6)) == v2); 
        should(exec(-Arg1(), v_1_5) == v1_5); 
    }
    
    void testBinary()
    {
        should(exec(Arg1() + Arg2(), v1, v1_5) == v2_5); 
        should(exec(Arg1() - Arg2(), v2_5, v1) == v1_5); 
        should(exec(Arg1() * Arg2(), v2, v1_5) == v3_0); 
        should(exec(Arg1() * Arg2(), v1_5, 2.0) == v3_0); 
        should(exec(Arg1() * Arg2(), 2.0, v1_5) == v3_0); 
        should(exec(Arg1() / Arg2(), v3_0, 2.0) == v1_5); 
        should(exec(Arg1() == Arg2(), v1, v1_5) == false); 
        should(exec(Arg1() != Arg2(), v1, v1_5) == true); 
    }
    
    void testCombine()
    {
        should(exec((Arg1() * Arg2() + Arg3()) / Param(2.0), v1_5, v2, v1) == v2);
    }
    
    void testComma()
    {
        exec((Var(v0) = Arg1(), Var(v3_0) = Arg2()),
             v1, v1_5);
        should(v0 == v1);
        should(v3_0 == v1_5);
        
        should(exec((Var(v0) = Arg2(), Arg3()), true, v2, v2_5) == v2_5);
        should(v0 == v2);
    }
    
    void testApply()
    {
        should(exec(applyFct(
          (vigra::RGBValue<double> (*)(vigra::RGBValue<double> const &))&vigra::abs, Arg1()), 
          v_1_5) == v1_5); 
    }
    
    void testSequence()
    {
        vigra::RGBValue<int> data[] = { v0, v1, v2, vigra::RGBValue<int>(3)};
        vigra::RGBValue<int> sum(0);
        int count = 0;
        
        std::for_each(data, data+4, 
        (
            Var(count) += Param(1),
            Var(sum) += Arg1()
        ));
        
        should(count == 4); 
        should(sum == vigra::RGBValue<int>(6)); 
    }  
};

struct FunctorTestSuite
: public vigra::test_suite
{
    FunctorTestSuite()
    : vigra::test_suite("FunctorTestSuite")
    {
        add( testCase(&FunctorExpressionTest::testArg1));
        add( testCase(&FunctorExpressionTest::testArg2));
        add( testCase(&FunctorExpressionTest::testArg3));
        add( testCase(&FunctorExpressionTest::testParam));
        add( testCase(&FunctorExpressionTest::testVar));
        add( testCase(&FunctorExpressionTest::testAssignment));
        add( testCase(&FunctorExpressionTest::testIfThen));
        add( testCase(&FunctorExpressionTest::testIfThenElse));
        add( testCase(&FunctorExpressionTest::testUnary));
        add( testCase(&FunctorExpressionTest::testBinary));
        add( testCase(&FunctorExpressionTest::testCombine));
        add( testCase(&FunctorExpressionTest::testComma));
        add( testCase(&FunctorExpressionTest::testApply));
        add( testCase(&FunctorExpressionTest::testSequence));

        add( testCase(&FunctorRGBExpressionTest::testArg1));
        add( testCase(&FunctorRGBExpressionTest::testArg2));
        add( testCase(&FunctorRGBExpressionTest::testArg3));
        add( testCase(&FunctorRGBExpressionTest::testParam));
        add( testCase(&FunctorRGBExpressionTest::testVar));
        add( testCase(&FunctorRGBExpressionTest::testAssignment));
        add( testCase(&FunctorRGBExpressionTest::testIfThen));
        add( testCase(&FunctorRGBExpressionTest::testIfThenElse));
        add( testCase(&FunctorRGBExpressionTest::testUnary));
        add( testCase(&FunctorRGBExpressionTest::testBinary));
        add( testCase(&FunctorRGBExpressionTest::testCombine));
        add( testCase(&FunctorRGBExpressionTest::testComma));
        add( testCase(&FunctorRGBExpressionTest::testApply));
        add( testCase(&FunctorRGBExpressionTest::testSequence));
    }
};

int main()
{
    FunctorTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
