/************************************************************************/
/*                                                                      */
/*             Copyright 2011-2012 by Ullrich Koethe                    */
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

#include <iostream>
#include <sstream>
#include <map>
#include <set>

#include <unittest.hxx>

//#include <vigra/accessor.hxx>
//#include <vigra/tinyvector.hxx>
//#include <vigra/rgbvalue.hxx>

//#include <vigra/coordinate_iterator.hxx>
//#include <vigra/object_features.hxx>

//#include <vigra/multi_pointoperators.hxx>

//#include <vigra/basicimage.hxx>
//#include <vigra/stdimage.hxx> // BImage
//#include <vigra/inspectimage.hxx> // FindAverageAndVariance

#include <vigra/multi_array.hxx>
//#include <vigra/multi_convolution.hxx>
//#include <vigra/impex.hxx>
//#include <vigra/imageinfo.hxx>
//#include <vigra/functorexpression.hxx>

//#include <vigra/algorithm.hxx>
//#include <vigra/random.hxx>
//#include <vigra/convolution.hxx>
#include <vigra/accumulator.hxx>

namespace std {

template <unsigned int N, class T, class Stride>
ostream & operator<<(ostream & o, vigra::MultiArrayView<N, T, Stride> const & m)
{
    for(vigra::MultiArrayIndex k=0; k<m.size(); ++k)
        o << m[k] << " ";
    return o;
}

} // namespace std

using namespace vigra;

// mask cl.exe shortcomings
#if defined(_MSC_VER)
#pragma warning( disable : 4503 )
#endif

struct AccumulatorTest
{
    AccumulatorTest()
    {}
    
    void testStandardizeTag()
    {
        using namespace vigra::acc;

            // count
        should((IsSameType<StandardizeTag<Count>::type,
                           Count>::value));
        should((IsSameType<StandardizeTag<Central<Count> >::type,
                           Count>::value));
        should((IsSameType<StandardizeTag<Principal<Count> >::type,
                           Count>::value));
        should((IsSameType<StandardizeTag<Coord<Count> >::type,
                           Count>::value));
        should((IsSameType<StandardizeTag<Weighted<Count> >::type,
                           Weighted<Count> >::value));
        should((IsSameType<StandardizeTag<CoordWeighted<Count> >::type,
                           Weighted<Count> >::value));
        should((IsSameType<StandardizeTag<Coord<Central<Count> > >::type,
                           Count>::value));
        shouldEqual(Count::name(), "PowerSum<0>");
        shouldEqual(Weighted<Count>::name(), "Weighted<PowerSum<0> >");

        // sum 
        should((IsSameType<StandardizeTag<Sum>::type,
                           Sum>::value));
        should((IsSameType<StandardizeTag<Central<Sum> >::type,
                           Central<Sum> >::value));
        should((IsSameType<StandardizeTag<Weighted<Central<Sum> > >::type,
                           Weighted<Central<Sum> > >::value));
        should((IsSameType<StandardizeTag<Central<Weighted<Sum> > >::type,
                           Weighted<Central<Sum> > >::value));
        shouldEqual(Weighted<Central<Sum> >::name(), "Weighted<Central<PowerSum<1> > >");

            // mean
        typedef DivideByCount<Principal<Sum> > PrincipalMean;
        should((IsSameType<StandardizeTag<Mean>::type,
                           Mean>::value));
        should((IsSameType<StandardizeTag<Principal<Mean> >::type,
                           PrincipalMean>::value));
        should((IsSameType<StandardizeTag<Weighted<Principal<Mean> > >::type,
                           Weighted<PrincipalMean> >::value));
        should((IsSameType<StandardizeTag<Principal<Weighted<Mean> > >::type,
                           Weighted<PrincipalMean> >::value));
        shouldEqual(Weighted<PrincipalMean>::name(), "Weighted<DivideByCount<Principal<PowerSum<1> > > >");

            // moment
        should((IsSameType<StandardizeTag<Moment<2> >::type,
                           DivideByCount<PowerSum<2> > >::value));
        should((IsSameType<StandardizeTag<Central<Moment<2> > >::type,
                           DivideByCount<Central<PowerSum<2> > > >::value));
        should((IsSameType<StandardizeTag<Coord<Central<Moment<2> > > >::type,
                           Coord<DivideByCount<Central<PowerSum<2> > > > >::value));
        should((IsSameType<StandardizeTag<Central<Coord<Moment<2> > > >::type,
                           Coord<DivideByCount<Central<PowerSum<2> > > > >::value));
        shouldEqual(Coord<Central<Moment<2> > >::name(), "Coord<DivideByCount<Central<PowerSum<2> > > >");

        should((IsSameType<StandardizeTag<Principal<Moment<2> > >::type,
                           DivideByCount<Principal<PowerSum<2> > > >::value));
        should((IsSameType<StandardizeTag<Coord<Principal<Moment<2> > > >::type,
                           Coord<DivideByCount<Principal<PowerSum<2> > > > >::value));
        should((IsSameType<StandardizeTag<Principal<Coord<Moment<2> > > >::type,
                           Coord<DivideByCount<Principal<PowerSum<2> > > > >::value));

        should((IsSameType<StandardizeTag<Moment<3> >::type,
                           DivideByCount<PowerSum<3> > >::value));
        should((IsSameType<StandardizeTag<Principal<Moment<3> > >::type,
                           DivideByCount<Principal<PowerSum<3> > > >::value));
        should((IsSameType<StandardizeTag<Central<Moment<3> > >::type,
                           DivideByCount<Central<PowerSum<3> > > >::value));

            // SumOfSquaredDifferences
        should((IsSameType<StandardizeTag<SumOfSquaredDifferences>::type,
                           Central<PowerSum<2> > >::value));
        should((IsSameType<StandardizeTag<Weighted<SSD> >::type,
                           Weighted<Central<PowerSum<2> > > >::value));

            // variance
        typedef DivideByCount<Central<PowerSum<2> > > CentralVarianceDesired;
        typedef DivideByCount<Principal<PowerSum<2> > > PrincipalVarianceDesired;
        should((IsSameType<StandardizeTag<Variance>::type,
                           CentralVarianceDesired >::value));
        should((IsSameType<StandardizeTag<Variance>::type,
                           StandardizeTag<Central<Moment<2> > >::type >::value));
        should((IsSameType<StandardizeTag<Principal<Variance> >::type,
                           PrincipalVarianceDesired >::value));
        should((IsSameType<StandardizeTag<Principal<Variance> >::type,
                           StandardizeTag<Principal<Moment<2> > >::type >::value));
        should((IsSameType<StandardizeTag<Coord<Central<Variance> > >::type,
                           Coord<CentralVarianceDesired> >::value));
        should((IsSameType<StandardizeTag<Central<Coord<Variance> > >::type,
                           Coord<CentralVarianceDesired> >::value));
        should((IsSameType<StandardizeTag<Coord<Principal<Variance> > >::type,
                           Coord<PrincipalVarianceDesired> >::value));
        should((IsSameType<StandardizeTag<Principal<Coord<Variance> > >::type,
                           Coord<PrincipalVarianceDesired> >::value));

            // IsCoordinateFeature
        should(!IsCoordinateFeature<StandardizeTag<Variance>::type>::value);
        should(IsCoordinateFeature<StandardizeTag<Coord<Variance> >::type>::value);
        should(!IsCoordinateFeature<StandardizeTag<Principal<Variance> >::type>::value);
        should(IsCoordinateFeature<StandardizeTag<Principal<Coord<Variance> > >::type>::value);
        should(!IsCoordinateFeature<StandardizeTag<Global<Variance> >::type>::value);
        should(IsCoordinateFeature<StandardizeTag<Global<Coord<Variance> > >::type>::value);

            // IsPrincipalFeature
        should(!IsPrincipalFeature<StandardizeTag<Variance>::type>::value);
        should(!IsPrincipalFeature<StandardizeTag<Coord<Variance> >::type>::value);
        should(IsPrincipalFeature<StandardizeTag<Principal<Variance> >::type>::value);
        should(IsPrincipalFeature<StandardizeTag<Principal<Coord<Variance> > >::type>::value);
        should(!IsPrincipalFeature<StandardizeTag<Global<Variance> >::type>::value);
        should(IsPrincipalFeature<StandardizeTag<Global<Principal<Coord<Variance> > > >::type>::value);

            // std dev
        typedef RootDivideByCount<Central<PowerSum<2> > > CentralStdDevDesired;
        typedef RootDivideByCount<Principal<PowerSum<2> > > PrincipalStdDevDesired;
        should((IsSameType<StandardizeTag<StdDev>::type,
                           CentralStdDevDesired>::value));
        should((IsSameType<StandardizeTag<Principal<StdDev> >::type,
                           PrincipalStdDevDesired>::value));
        should((IsSameType<StandardizeTag<Coord<Central<StdDev> > >::type,
                           Coord<CentralStdDevDesired> >::value));
        should((IsSameType<StandardizeTag<Central<Coord<StdDev> > >::type,
                           Coord<CentralStdDevDesired> >::value));
        should((IsSameType<StandardizeTag<Coord<Principal<StdDev> > >::type,
                           Coord<PrincipalStdDevDesired> >::value));
        should((IsSameType<StandardizeTag<Principal<Coord<StdDev> > >::type,
                           Coord<PrincipalStdDevDesired> >::value));

            // skewness
        should((IsSameType<StandardizeTag<Skewness>::type,
                           Skewness >::value));
        should((IsSameType<StandardizeTag<Coord<Skewness> >::type,
                           Coord<Skewness> >::value));
        should((IsSameType<StandardizeTag<Principal<Skewness> >::type,
                           Principal<Skewness> >::value));
        should((IsSameType<StandardizeTag<Coord<Principal<Skewness> > >::type,
                           Coord<Principal<Skewness> > >::value));
        should((IsSameType<StandardizeTag<Principal<Coord<Skewness> > >::type,
                           Coord<Principal<Skewness> > >::value));

            // AbsPowerSum
        should((IsSameType<StandardizeTag<AbsSum>::type,
                           AbsPowerSum<1> >::value));
        should((IsSameType<StandardizeTag<AbsPowerSum<0> >::type,
                           Count>::value));
        should((IsSameType<StandardizeTag<AbsPowerSum<2> >::type,
                           PowerSum<2> >::value));
        should((IsSameType<StandardizeTag<AbsPowerSum<3> >::type,
                           AbsPowerSum<3> >::value));

            // CovarianceEigensystem
        should((IsSameType<StandardizeTag<CovarianceEigensystem>::type,
                           DivideByCount<ScatterMatrixEigensystem> >::value));
        should((IsSameType<StandardizeTag<Central<CovarianceEigensystem> >::type,
                           DivideByCount<ScatterMatrixEigensystem> >::value));
        should((IsSameType<StandardizeTag<Coord<CovarianceEigensystem> >::type,
                           Coord<DivideByCount<ScatterMatrixEigensystem> > >::value));

            // CoordinateSystem 
        should((IsSameType<StandardizeTag<Central<CoordinateSystem> >::type,
                           CoordinateSystem>::value));
        should((IsSameType<StandardizeTag<Principal<CoordinateSystem> >::type,
                           Principal<CoordinateSystem> >::value));

            // RegionRadii
        should((IsSameType<StandardizeTag<RegionRadii>::type,
                           Coord<RootDivideByCount<Principal<PowerSum<2> > > > >::value));
        shouldEqual(StandardizeTag<RegionRadii>::type::name(), "Coord<RootDivideByCount<Principal<PowerSum<2> > > >");

            // HasModifierPriority
        using namespace vigra::acc::detail;
        should((HasModifierPriority<StandardizeTag<Count>::type, AccumulatorPriority>::value));
        should(!(HasModifierPriority<StandardizeTag<Count>::type, AccessDataPriority>::value));
        should((HasModifierPriority<StandardizeTag<Weighted<Count> >::type, WeightingPriority>::value));

            // nested Select
        should((IsSameType<Select<Count, Mean, Select<Sum, Minimum>, Variance>::type,
                           Select<Count, Mean, Sum, Minimum, Variance>::type>::value));
    }

    template <class SOURCE, class REFERENCE>
    void testLongFormImpl(const char * message)
    {
        using namespace vigra::acc;
        using namespace vigra::acc::detail;

        typedef typename StandardizeTag<SOURCE >::type StdSource;
        typedef typename TagLongForm<StdSource, MinPriority>::type LongSource;
        typedef typename StandardizeTagLongForm<LongSource>::type Dest;
            
        shouldMsg((IsSameType<LongSource, REFERENCE >::value), message);
//        shouldMsg((IsSameType<LongSource, REFERENCE >::value), typeid(LongSource).name());
        shouldMsg((IsSameType<StdSource, Dest>::value), message);
    }

    void testTagTransfer()
    {
        using namespace vigra::acc;

#define TEST_LONG_FORM(SOURCE, TARGET) testLongFormImpl<SOURCE, TARGET >(#SOURCE)
#define DM DefaultModifier
        {
            using namespace vigra::acc::detail;


            TEST_LONG_FORM(Minimum, DM<DM<DM<DM<DM<Minimum> > > > >);
            TEST_LONG_FORM(Principal<Minimum>, DM<DM<DM<DM<Principal<Minimum> > > > >);
            TEST_LONG_FORM(Weighted<Coord<Principal<Minimum> > >, DM<Weighted<Coord<DM<Principal<Minimum> > > > >);
            TEST_LONG_FORM(Mean, DM<DM<DM<DivideByCount<DM<Sum> > > > >);
            TEST_LONG_FORM(Coord<Mean>, DM<DM<Coord<DivideByCount<DefaultModifier<Sum> > > > >);
            TEST_LONG_FORM(Count, DM<DM<DM<DM<DM<Count> > > > >);
            TEST_LONG_FORM(Weighted<Count>, DM<Weighted<DM<DM<DM<Count> > > > >);
            TEST_LONG_FORM(Coord<Count>, DM<DM<DM<DM<DM<Count> > > > >);
            TEST_LONG_FORM(Principal<Variance>, DM<DM<DM<DivideByCount<Principal<PowerSum<2> > > > > >);
            TEST_LONG_FORM(Global<Count>, Global<DM<DM<DM<DM<PowerSum<0> > > > > >);
        }
#undef TEST_LONG_FORM
#undef DM
          
          typedef Select<Count, Sum, Mean, Variance>::type Target;

          should((IsSameType<TransferModifiers<Minimum, Target>::type,
                             Target>::value));

          typedef Select<Count, Coord<Sum>, Coord<Mean>, Coord<Variance> >::type Desired1;
          should((IsSameType<TransferModifiers<Coord<Minimum>, Target>::type,
                             Desired1>::value));

          typedef Select<Count, Sum, Mean, Principal<Variance> >::type Desired2;
          should((IsSameType<TransferModifiers<Principal<Minimum>, Target>::type,
                             Desired2>::value));

          typedef Select<Count, Coord<Sum>, Coord<Mean>, Coord<Principal<Variance> > >::type Desired3;
          should((IsSameType<TransferModifiers<Coord<Principal<Minimum> >, Target>::type,
                             Desired3>::value));

          typedef Select<Weighted<Count>, CoordWeighted<Sum>, CoordWeighted<Mean>, CoordWeighted<Variance> >::type Desired4;
          should((IsSameType<TransferModifiers<Coord<Weighted<Minimum> >, Target>::type,
                             Desired4>::value));

          typedef Select<Weighted<Count>, CoordWeighted<Sum>, CoordWeighted<Mean>, CoordWeighted<Principal<Variance> > >::type Desired5;
          should((IsSameType<TransferModifiers<Principal<Coord<Weighted<Minimum> > >, Target>::type,
                             Desired5>::value));

          should((IsSameType<TransferModifiers<Principal<Minimum>, Centralize>::type,
                             PrincipalProjection>::value));
          should((IsSameType<TransferModifiers<Principal<Weighted<Coord<Minimum> > >, Centralize>::type,
                             Weighted<Coord<PrincipalProjection> > >::value));
    }

    void testScalar()
    {
        using namespace vigra::acc;
        
        { 
            typedef AccumulatorChain<double, Select<Count> > A;
            A a;

            shouldEqual(1, a.passesRequired());

            a(1.0);
            a(2.0);
            a(3.0);
            
            shouldEqual(get<Count>(a), 3.0);
            // access to an inactive statistic triggers a static assertion
            // get<Mean>(a);
        }

        {
            typedef AccumulatorChain<double, Select<CovarianceEigensystem, Covariance, UnbiasedVariance, UnbiasedStdDev, 
                                               Variance, StdDev, Minimum, Maximum, Skewness, Kurtosis,
                                               AbsSum, SumOfAbsDifferences, MeanAbsoluteDeviation, 
                                               Principal<Variance>, Principal<CoordinateSystem>
                                              > > A;

            A a;


            shouldEqual(2, a.passesRequired());
            shouldEqual(24, A::staticSize);

            double data[] = { 1.0, 2.0, 3.0, 5.0 };

            for(int k=0; k<4; ++k)
                a(data[k]);

            shouldEqual(get<Count>(a), 4.0);
            shouldEqual(get<Minimum>(a), 1.0);
            shouldEqual(get<Maximum>(a), 5.0);
            shouldEqual(get<Sum>(a), 11.0);
            shouldEqual(get<AbsSum>(a), 11.0);
            shouldEqual(get<Mean>(a), 2.75);
            shouldEqualTolerance(get<UnbiasedVariance>(a), 2.9166666666666665, 1e-15);
            shouldEqualTolerance(get<UnbiasedStdDev>(a), sqrt(2.9166666666666665), 1e-15);
            shouldEqualTolerance(get<Variance>(a), 2.1875, 1e-15);
            shouldEqualTolerance(get<StdDev>(a), sqrt(2.1875), 1e-15);
            shouldEqualTolerance(get<Covariance>(a), 2.1875, 1e-15);

            std::pair<double, double> seigen = get<ScatterMatrixEigensystem>(a);
            shouldEqual(seigen.first, 8.75);
            shouldEqual(seigen.second, 1.0);

            std::pair<double, double> eigen = get<CovarianceEigensystem>(a);
            shouldEqual(eigen.first, 2.1875);
            shouldEqual(eigen.second, 1.0);

            shouldEqual(get<Principal<Variance> >(a), 2.1875);
            shouldEqual(get<Principal<CoordinateSystem> >(a), 1.0);

            for(int k=0; k<4; ++k)
                a.updatePass2(data[k]);

            shouldEqual(get<Count>(a), 4.0);
            shouldEqualTolerance(get<CentralMoment<2> >(a),  2.1875, 1e-15);
            shouldEqualTolerance(get<Skewness>(a), 0.43465075957466565, 1e-15);
            shouldEqualTolerance(get<Kurtosis>(a), -1.1542857142857144, 1e-15);
            shouldEqual(get<SumOfAbsDifferences>(a), 5);
            shouldEqual(get<MeanAbsoluteDeviation>(a), 1.25);
        }

        { 
            DynamicAccumulatorChain<double, Select<Mean, Covariance, StdDev, Minimum, CentralMoment<4> > > a;

            shouldEqual(0, a.passesRequired());

            activate<Count>(a);
            should(isActive<Count>(a));
            should(!isActive<Covariance>(a));
            shouldEqual(1, a.passesRequired());

            a(1.0);
            a(2.0);
            a(3.0);

            shouldEqual(get<Count>(a), 3.0);

            try 
            {
                get<Mean>(a);
                failTest("get<Mean>() failed to throw exception");
            }
            catch(ContractViolation & c) 
            {
                std::string expected("\nPrecondition violation!\nget(accumulator): attempt to access inactive statistic");
                std::string message(c.what());
                should(0 == expected.compare(message.substr(0,expected.size())));
            }

            a.reset();
            shouldEqual(0, a.passesRequired());
            should(!isActive<Count>(a));

            activate<Mean>(a);
            activate<StdDev>(a);
            activate<Covariance>(a);
            activate<CentralMoment<4> >(a);
            //activate<Minimum>(a);
            a.activate("Minimum");

            should(isActive<Count>(a));
            should(isActive<Minimum>(a));
            should(isActive<Sum>(a));
            should(isActive<Mean>(a));
            should(isActive<Variance>(a));
            should(isActive<StdDev>(a));
            should(isActive<CentralMoment<4> >(a));
            should(isActive<Covariance>(a));

            shouldEqual(2, a.passesRequired());

            a.reset();
            a.activateAll();

            a(1.0);
            a(2.0);
            a(3.0);

            shouldEqual(get<Count>(a), 3.0);
            shouldEqual(get<Minimum>(a), 1.0);
            shouldEqual(get<Sum>(a), 6.0);
            shouldEqual(get<Mean>(a), 2.0);
            shouldEqual(get<Variance>(a), 2.0/3.0);
            shouldEqual(get<StdDev>(a), sqrt(2.0/3.0));
            shouldEqual(get<Covariance>(a), 2.0/3.0);

            a.updatePass2(1.0);
            a.updatePass2(2.0);
            a.updatePass2(3.0);

            shouldEqual(get<Count>(a), 3.0);
            shouldEqual(get<CentralMoment<4> >(a), 2.0/3.0);
        }
    }

    void testVector()
    {
        using namespace vigra::acc;

        {
            typedef TinyVector<int, 3> V;
            typedef AccumulatorChain<V, Select<StdDev, Mean, CovarianceEigensystem, Covariance, Minimum, Maximum, CentralMoment<2>,
                                          AbsSum, SumOfAbsDifferences, MeanAbsoluteDeviation, 
                                          Principal<Variance>, Principal<CoordinateSystem>, Principal<Sum>,
                                          Principal<Minimum>, Principal<Maximum>, Principal<Skewness>, Principal<Kurtosis>, Principal<SumOfAbsDifferences>
                                          > > A;
            typedef LookupTag<Mean, A>::value_type W;
            typedef LookupTag<Covariance, A>::value_type Var;

            A a;

            static const int SIZE = 4;
            V d[SIZE] = { V(1,2,3), V(2,3,0), V(3,4,2), V(2,1,2) };

            for(int k=0; k<SIZE; ++k)
                a(d[k]);
            for(int k=0; k<SIZE; ++k)
                a.updatePass2(d[k]);

            shouldEqual(get<Count>(a), 4.0);
            shouldEqual(get<Minimum>(a), V(1,1,0));
            shouldEqual(get<Maximum>(a), V(3,4,3));
            shouldEqual(get<Sum>(a), W(8.0, 10.0, 7.0));
            shouldEqual(get<AbsSum>(a), W(8.0, 10.0, 7.0));
            shouldEqual(get<Mean>(a), W(2.0, 2.5, 7.0/4.0));
            shouldEqual(get<CentralMoment<2> >(a), W(0.5, 1.25, 1.1875));
            shouldEqual(get<Variance>(a),  W(0.5, 1.25, 1.1875));
            shouldEqualTolerance(sqrt(W(0.5, 1.25, 1.1875)), get<StdDev>(a), W(1e-15));
            shouldEqualTolerance(get<SumOfAbsDifferences>(a), W(2.0, 4.0, 3.5), W(1e-15));
            shouldEqualTolerance(get<MeanAbsoluteDeviation>(a), W(0.5, 1.0, 7.0/8.0), W(1e-15));

            double covarianceData[] = { 
                0.5,   0.5,  -0.25,
                0.5,   1.25, -0.375,
               -0.25, -0.375, 1.1875 };
            Var covariance(3,3, covarianceData);
            shouldEqual(get<Covariance>(a).shape(), Shape2(3,3));
            shouldEqual(get<Covariance>(a), covariance);
            std::pair<W const &, Var const &> eigen = get<CovarianceEigensystem>(a);
            W ew(1.8181423035878563, 0.87335382939336145, 0.24600386701878226); 
            shouldEqualTolerance(ew, eigen.first, W(1e-15));
            shouldEqualTolerance(ew, get<Principal<Variance> >(a), W(1e-15));

            double eigenvectorData[] = {
                -0.38281255664062192, -0.19398130489891852, -0.90323075668844639,
                -0.71942795069852928, -0.55075408575738086,  0.42319423528123123,
                 0.57954979961344579,- 0.81181351945583191, -0.071280006991795777 };
            Var ev(3,3, eigenvectorData),
                eps(3, 3, 1e-14);
            shouldEqualTolerance(ev, eigen.second, eps);
            shouldEqualTolerance(ev, get<Principal<CoordinateSystem> >(a), eps);

            shouldEqualTolerance(get<Principal<Sum> >(a), W(0.0), W(1e-15));
            shouldEqualTolerance(get<Principal<Minimum> >(a), W(-1.3739261246727945, -1.2230658133989472, -0.6526113546697957), W(1e-15));
            shouldEqualTolerance(get<Principal<Maximum> >(a), W(1.4669637815066938,  1.1452966161690161, 0.60253363030808593), W(1e-15));
            shouldEqualTolerance(get<Principal<Skewness> >(a), W(0.01148108748350361, -0.07581454384153662, -0.09140344434535799), W(1e-14));
            shouldEqualTolerance(get<Principal<Kurtosis> >(a), W(-1.9829394126459396, -1.6241963546875782, -1.6255854346698215), W(1e-14));
            shouldEqualTolerance(get<Principal<SumOfAbsDifferences> >(a), W(5.3819863149157, 3.5369487298822575, 1.8777415203686885), W(1e-14));
        }

        {
            using namespace vigra::multi_math;
            
            typedef MultiArray<1, int> V;
            typedef TinyVector<int, 3> T;
            typedef AccumulatorChain<V::view_type, Select<Covariance, Mean, StdDev, Minimum, Maximum, CentralMoment<2> > > A;
            typedef LookupTag<Mean, A>::value_type W;
            typedef LookupTag<Covariance, A>::value_type Var;

            A a;

            Shape1 s(3);

            V data[] = { V(s, 1), V(s, 2), V(s, 3) };

            a(data[0]);
            a(data[1]);
            a(data[2]);

            a.updatePass2(data[0]);
            a.updatePass2(data[1]);
            a.updatePass2(data[2]);

            shouldEqual(get<Count>(a), 3.0);
            shouldEqual(get<Minimum>(a), V(s, 1));
            shouldEqual(get<Maximum>(a), V(s, 3));
            shouldEqual(get<Sum>(a), W(s, 6.0));
            shouldEqual(get<Mean>(a),  W(s, 2.0));
            shouldEqual(get<CentralMoment<2> >(a),  W(s, 2.0 / 3.0));

            Var covariance(3,3);
            covariance.init(2.0/3.0);
            shouldEqual(get<Covariance>(a), covariance);

            W variance(s, 2.0/3.0);
            shouldEqual(get<Variance>(a), variance);

            W stddev = sqrt(variance);
            shouldEqualTolerance(stddev, get<StdDev>(a), W(1e-15));

            a.reset(1);

            a(V(s, T(0, 2, 4).begin()));
            shouldEqual(get<Minimum>(a), V(s, T(0,1,1).begin()));
            shouldEqual(get<Maximum>(a), V(s, T(3,3,4).begin()));
        }

        {
            typedef TinyVector<double, 2> V;
            static const int SIZE = 20;

            V data[SIZE] = {
                V(1.88417085437108889, 3.10984300178095197),
                V(3.22221249967652135, 4.62610895051767734),
                V(5.02965943418706019, 3.61409282557627254),
                V(1.99001871343947201, 0.44572597781938073),
                V(0.82190895017090260, 1.52581824695525770),
                V(1.79509471114960295, 4.54126165421070915),
                V(0.63954006398369945, 4.03816177019905265),
                V(-1.19055182745611221, 3.05473509195811443),
                V(1.60460514736327031, 2.39320817128161423),
                V(-1.26828508191601941, 3.08007018243650110),
                V(2.67471223051054885, 2.36574957680121889),
                V(2.54777120650106648, 3.00252905176459528),
                V(-0.49276533572213554, -1.13913810037296859),
                V(3.08185249197166877, 2.61911514709572302),
                V(-3.21266705448485190, 0.62656641585875028),
                V(1.25495155317623874, 0.46593677346153228),
                V(0.71219264300245499, 2.68491068466799110),
                V(0.91993307568972671, 1.99693758751466821),
                V(2.11157305527596195, -1.11069145843301786),
                V(0.10242024165277441, 2.44449590189711241)
            };

            typedef AccumulatorChain<V, Select<Mean, Variance, Covariance, Central<Sum>,
                                          Principal<Variance>, Principal<CoordinateSystem>, Principal<Sum>
                                          > > A;
            typedef LookupTag<Covariance, A>::value_type Var;

            A a;

            for(int k=0; k<SIZE; ++k)
                a(data[k]);

            for(int k=0; k<SIZE; ++k)
                a.updatePass2(data[k]);


            shouldEqual(get<Count>(a), (double)SIZE);
            shouldEqualTolerance(get<Mean>(a), V(1.2114173786271469, 2.2192718726495571), V(1e-15));
            shouldEqualTolerance(get<Central<Sum> >(a), V(0.0),  V(1e-14));
            shouldEqualTolerance(get<Principal<Sum> >(a), V(0.0),  V(1e-14));

            Var eps = Var(2,2, 1e-15);
            double covarianceData[] = {
                 3.24260523085696217,  0.85916467806966068,
                 0.85916467806966068,  2.57086707635742417 };
            Var cov(2,2, covarianceData);
            shouldEqualTolerance(cov, get<Covariance>(a), eps);

            V principalVariance(3.82921757948803698, 1.98425472772634914);
            shouldEqualTolerance(get<Principal<Variance> >(a), principalVariance, V(1e-15));

            double eigenvectorData[] = {
                 0.82586108137035807,  0.56387363325286188,
                 0.56387363325286188, -0.82586108137035807 };
            Var ev(2,2, eigenvectorData);
            shouldEqualTolerance(ev, get<Principal<CoordinateSystem> >(a), eps);

                // desired principal projection
            V principal[SIZE] = {
                V(1.05777049122698186, -0.35614008909044692),
                V(3.01778940082927427, -0.85387872117719943),
                V(3.93984027140249449, 1.00107768000535380),
                V(-0.35703922713407876, 1.90373529408056230),
                V(-0.71270006734035740, 0.35306282845985459),
                V(1.79134520751326742, -1.58852073379423842),
                V(0.55333283753554596, -1.82461691686673699),
                V(-1.51259720503357187, -2.04438366093898516),
                V(0.42279656853422592, 0.07806099602999564),
                V(-1.56250828510342998, -2.10913865877495876),
                V(1.29107318397063953, 0.70414314950436341),
                V(1.54530068874035664, 0.10667306729959503),
                V(-3.30113701276602800, 1.81263639313852387),
                V(1.77018064523347562, 0.72447404411252858),
                V(-4.55176376662593363, -1.17927111226656489),
                V(-0.95270625221923966, 1.47255899419368674),
                V(-0.14972883129345838, -0.66605263830794725),
                V(-0.36609398278207039, 0.01925692021370382),
                V(-1.13427498157747753, 3.25766116941665018),
                V(-0.78887968311061840, -0.81133800523773636)
            };

            for(int k=0; k<SIZE; ++k)
                shouldEqualTolerance(principal[k], getAccumulator<PrincipalProjection>(a)(data[k]), V(1e-14));

                // check that statistics of points in principal representation have expected properties
            a.reset();
            for(int k=0; k<SIZE; ++k)
                a(principal[k]);

            shouldEqual(SIZE, get<Count>(a));
            shouldEqualTolerance(V(0.0), get<Mean>(a), V(1e-14));
            shouldEqualTolerance(get<Variance>(a), principalVariance, V(1e-15));
            shouldEqualTolerance(get<Principal<Variance> >(a), principalVariance, V(1e-15));

            Var ref = linalg::identityMatrix<double>(2);
            shouldEqualTolerance(get<Principal<CoordinateSystem> >(a), ref, eps);
            ref(0,0) = principalVariance[0];
            ref(1,1) = principalVariance[1];
            shouldEqualTolerance(get<Covariance>(a), ref, eps);
        }
    }

    void testMerge()
    {
        using namespace vigra::acc;
        
        typedef AccumulatorChain<double, Select<Covariance, StdDev, Minimum, Maximum, Skewness, Kurtosis, 
                                           CentralMoment<3>, CentralMoment<4> > > A;

        A a, b;

        double data[] = { 1.0, 2.0, 3.0, 4.0, 6.0 };

        for(int k=0; k<3; ++k)
            a(data[k]);

        for(int k=0; k<3; ++k)
            a.updatePass2(data[k]);

        for(int k=3; k<5; ++k)
            b(data[k]);

        for(int k=3; k<5; ++k)
            b.updatePass2(data[k]);

        a += b;

        shouldEqual(get<Count>(a), 5.0);
        shouldEqual(get<Minimum>(a), 1.0);
        shouldEqual(get<Maximum>(a), 6.0);
        shouldEqual(get<Sum>(a), 16.0);
        shouldEqualTolerance(get<Mean>(a), 3.2, 1e-15);
        shouldEqualTolerance(get<SSD>(a), 14.8, 1e-15);
        shouldEqualTolerance(get<Variance>(a), 2.96, 1e-15);
        shouldEqualTolerance(get<StdDev>(a), sqrt(2.96), 1e-15);
        shouldEqualTolerance(get<FlatScatterMatrix>(a), 14.8, 1e-15);
        shouldEqualTolerance(get<Covariance>(a), 2.96, 1e-15);
        shouldEqualTolerance(get<CentralMoment<2> >(a), 2.96, 1e-15);
        shouldEqualTolerance(get<CentralMoment<3> >(a), 2.016, 1e-15);
        shouldEqualTolerance(get<CentralMoment<4> >(a), 17.4752, 1e-15);
        shouldEqualTolerance(get<Skewness>(a), 0.395870337343817, 1e-15);
        shouldEqualTolerance(get<Kurtosis>(a), -1.0054784514243973, 1e-15);
    }

    void testCoordAccess()
    {
        using namespace vigra::acc;

        {
            typedef CoupledIteratorType<3>::type Iterator;
            typedef Iterator::value_type Handle;
            typedef Shape3 V;

            typedef AccumulatorChain<Handle, Select<Coord<Maximum>, Coord<Minimum>, Coord<Mean>, Coord<StdDev>, Coord<Covariance>,
                                               Coord<Principal<Variance> >, Coord<Principal<CoordinateSystem> >,
                                               Coord<AbsSum>, Coord<MeanAbsoluteDeviation>, Coord<CovarianceEigensystem>
                                          > > A;

            typedef LookupTag<Coord<Mean>, A>::value_type W;
            typedef LookupTag<Coord<Covariance>, A>::value_type Var;
            
            A a;

            Iterator i = createCoupledIterator(V(4,4,4));

            a(*(i+V(1,2,3)));
            a(*(i+V(2,3,1)));
            a(*(i+V(3,1,2)));

            a.updatePass2(*(i+V(1,2,3)));
            a.updatePass2(*(i+V(2,3,1)));
            a.updatePass2(*(i+V(3,1,2)));

            shouldEqual(get<Count>(a), 3.0);
            shouldEqual(get<Coord<Minimum> >(a), V(1));
            shouldEqual(get<Coord<Maximum> >(a), V(3));
            shouldEqual(get<Coord<Sum> >(a), W(6.0));
            shouldEqual(get<Coord<AbsSum> >(a), W(6.0));
            shouldEqual(get<Coord<Mean> >(a), W(2.0));
            shouldEqual(get<Coord<CentralMoment<2> > >(a), W(2.0/3.0));
            shouldEqual(get<Coord<Variance> >(a), W(2.0/3.0));
            shouldEqual(get<Coord<SumOfAbsDifferences> >(a), W(2.0));
            shouldEqual(get<Coord<MeanAbsoluteDeviation> >(a), W(2.0/3.0));

            W stddev = sqrt( W(2.0/3.0));
            shouldEqualTolerance(stddev, get<Coord<StdDev> >(a), W(1e-15));

            double covarianceData[] = { 
                2.0/3.0, -1.0/3.0, -1.0/3.0,
               -1.0/3.0,  2.0/3.0, -1.0/3.0,
               -1.0/3.0, -1.0/3.0,  2.0/3.0 };
            Var covariance(3,3, covarianceData);
            shouldEqual(get<Coord<Covariance> >(a), covariance);

            W sew(3.0, 3.0, 0.0); 
            std::pair<W const &, Var const &> seigen = get<Coord<ScatterMatrixEigensystem> >(a);
            shouldEqualTolerance(sew, seigen.first, W(1e-15));

            W ew(1.0, 1.0, 0.0); 
            std::pair<W const &, Var const &> eigen = get<Coord<CovarianceEigensystem> >(a);
            shouldEqualTolerance(ew, eigen.first, W(1e-15));
            shouldEqualTolerance(ew, get<Coord<Principal<Variance> > >(a), W(1e-15));

            double eigenvectorData[] = {
                -0.7071067811865476, -0.4082482904638629, -0.5773502691896257,
                 0.7071067811865476, -0.4082482904638629, -0.5773502691896257,
                 0.0               ,  0.816496580927726,  -0.5773502691896257 };
            Var ev(3,3, eigenvectorData),
                eps(3, 3, 1e-15);
            shouldEqualTolerance(ev, seigen.second, W(1e-15));
            shouldEqualTolerance(ev, get<Coord<Principal<CoordinateSystem> > >(a), W(1e-15));
        }

        {
            typedef CoupledIteratorType<3, double, double>::type Iterator;
            typedef Iterator::value_type Handle;
            typedef Shape3 V;

            typedef AccumulatorChain<Handle, Select<Mean, Coord<Mean>, Coord<Maximum>, Coord<Minimum>, Weighted<Count>,
                                               Weighted<Mean>, CoordWeighted<Mean>,
                                               ArgMinWeight, ArgMaxWeight,
                                               Coord<ArgMinWeight>, Coord<ArgMaxWeight>, WeightArg<2>, DataArg<1>
                                          > > A;
            
            A a;

            typedef LookupTag<Coord<Mean>, A>::value_type W;
            
            MultiArray<3, double> data(Shape3(4,4,4), 1.0);
            data(1,2,3) = 0.5;
            data(3,1,2) = 4.0;
            MultiArray<3, double> weights(Shape3(4,4,4), 1.0);
            weights(1,2,3) = 0.5;
            weights(3,1,2) = 0.25;

            Iterator i = createCoupledIterator(data, weights);

            a(*(i+V(1,2,3)));
            a(*(i+V(2,3,1)));
            a(*(i+V(3,1,2)));

            shouldEqual(get<Count>(a), 3.0);
            shouldEqual(get<Coord<Minimum> >(a), V(1));
            shouldEqual(get<Coord<Maximum> >(a), V(3));
            shouldEqual(get<Coord<Mean> >(a), W(2.0));
            shouldEqualTolerance(get<Weighted<Mean> >(a), 1.2857142857142858, 1e-15);
            shouldEqualTolerance(get<Mean>(a), 1.8333333333333333, 1e-15);
            W coordWeightedMean(1.8571428571428572, 2.4285714285714284,  1.7142857142857142);
            shouldEqualTolerance(coordWeightedMean, get<CoordWeighted<Mean> >(a), W(1e-15));
            shouldEqual(4.0, get<ArgMinWeight>(a));
            shouldEqual(1.0, get<ArgMaxWeight>(a));
            shouldEqual(V(3,1,2), get<Coord<ArgMinWeight> >(a));
            shouldEqual(V(2,3,1), get<Coord<ArgMaxWeight> >(a));
        }
    }

    void testIndexSpecifiers()
    {
        using namespace vigra::acc;

        typedef CoupledIteratorType<3, double, double>::type Iterator;
        typedef Iterator::value_type Handle;
        typedef Shape3 V;
        
        typedef AccumulatorChain<Handle, Select<WeightArg<1>, DataArg<2>, Mean, Coord<Mean>, Coord<Maximum>, Coord<Minimum>, Weighted<Count>, Weighted<Mean>, CoordWeighted<Mean>, ArgMinWeight, ArgMaxWeight, Coord<ArgMinWeight>, Coord<ArgMaxWeight> > >  A;
        A a;
        
        typedef LookupTag<Coord<Mean>, A>::value_type W;
        
        MultiArray<3, double> data(Shape3(4,4,4), 1.0);
        data(1,2,3) = 0.5;
        data(3,1,2) = 4.0;
        MultiArray<3, double> weights(Shape3(4,4,4), 1.0);
        weights(1,2,3) = 0.5;
        weights(3,1,2) = 0.25;
        
        Iterator i = createCoupledIterator(weights, data);
        
        a(*(i+V(1,2,3)));
        a(*(i+V(2,3,1)));
        a(*(i+V(3,1,2)));
        
        shouldEqual(get<Count>(a), 3.0);
        shouldEqual(get<Coord<Minimum> >(a), V(1));
        shouldEqual(get<Coord<Maximum> >(a), V(3));
        shouldEqual(get<Coord<Mean> >(a), W(2.0));
        shouldEqualTolerance(get<Weighted<Mean> >(a), 1.2857142857142858, 1e-15);
        shouldEqualTolerance(get<Mean>(a), 1.8333333333333333, 1e-15);
        W coordWeightedMean(1.8571428571428572, 2.4285714285714284,  1.7142857142857142);
        shouldEqualTolerance(coordWeightedMean, get<CoordWeighted<Mean> >(a), W(1e-15));
        shouldEqual(4.0, get<ArgMinWeight>(a));
        shouldEqual(1.0, get<ArgMaxWeight>(a));
        shouldEqual(V(3,1,2), get<Coord<ArgMinWeight> >(a));
        shouldEqual(V(2,3,1), get<Coord<ArgMaxWeight> >(a));
        
    }
  
    void testHistogram()
    {
        static const int SIZE = 30, HSIZE = 10;
        int data[SIZE] = {4, 3, 2, 2, 2, 0, 3, 6, 8, 8, 4, 0, 2, 0, 2, 8, 7, 8, 6, 0, 9, 3, 7, 0, 9, 5, 9, 9, 2, 4};
        // the same sorted:
        // int data[SIZE] = {0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9};

        using namespace vigra::acc;
        {

            typedef AccumulatorChain<int, Select<StandardQuantiles<UserRangeHistogram<HSIZE> >, StandardQuantiles<AutoRangeHistogram<HSIZE> >,
                                                 StandardQuantiles<IntegerHistogram<HSIZE> >, StandardQuantiles<IntegerHistogram<0> >,
                                                 DivideByCount<UserRangeHistogram<HSIZE> >, StandardQuantiles<UserRangeHistogram<HSIZE+2> >,
                                                 StandardQuantiles<UserRangeHistogram<3> >, StandardQuantiles<IntegerHistogram<HSIZE+2> >
                                    > > A;
            A a;

            getAccumulator<UserRangeHistogram<HSIZE> >(a).setMinMax(-0.5, 9.5);
            getAccumulator<UserRangeHistogram<HSIZE+2> >(a).setMinMax(-1.5, 10.5);
            getAccumulator<UserRangeHistogram<3> >(a).setMinMax(-20.0, 30.0);  // all data in one bin
            getAccumulator<IntegerHistogram<0> >(a).setBinCount(HSIZE);

            shouldEqual(HSIZE, get<IntegerHistogram<HSIZE> >(a).size());
            shouldEqual(HSIZE, get<UserRangeHistogram<HSIZE> >(a).size());
            shouldEqual(HSIZE, get<AutoRangeHistogram<HSIZE> >(a).size());
            shouldEqual(HSIZE, get<IntegerHistogram<0> >(a).size());

            for(int k=0; k<SIZE; ++k)
                a(data[k]);

            for(int k=0; k<SIZE; ++k)
                a.updatePass2(data[k]);

            double h[HSIZE] = { 5.0, 0.0, 6.0, 3.0, 3.0, 1.0, 2.0, 2.0, 4.0, 4.0 };

            shouldEqualSequence(h, h+HSIZE, get<IntegerHistogram<HSIZE> >(a).begin());
            shouldEqualSequence(h, h+HSIZE, get<UserRangeHistogram<HSIZE> >(a).begin());
            shouldEqualSequence(h, h+HSIZE, get<AutoRangeHistogram<HSIZE> >(a).begin());
            shouldEqualSequence(h, h+HSIZE, get<IntegerHistogram<0> >(a).begin());

            double density[HSIZE] = { 5.0/30.0, 0.0, 6.0/30.0, 3.0/30.0, 3.0/30.0, 1.0/30.0, 2.0/30.0, 2.0/30.0, 4.0/30.0, 4.0/30.0 };
            shouldEqualSequence(density, density+HSIZE, get<DivideByCount<UserRangeHistogram<HSIZE> > >(a).begin());

            typedef LookupTag<StandardQuantiles<UserRangeHistogram<HSIZE> >, A>::value_type QuantileVector;
            static const int QSIZE = QuantileVector::static_size;
            shouldEqual(QSIZE, 7);

            double quser[QSIZE] = { 0.0, 0.3, 1.9166666666666666, 3.833333333333333, 7.625, 8.625, 9.0 };
            shouldEqualSequenceTolerance(quser, quser+QSIZE, get<StandardQuantiles<UserRangeHistogram<HSIZE> > >(a).begin(), 1e-15);
            shouldEqualSequenceTolerance(quser, quser+QSIZE, get<StandardQuantiles<UserRangeHistogram<HSIZE+2> > >(a).begin(), 1e-15);

            double q_onebin[QSIZE] = { 0.0, 0.9, 2.25, 4.5, 6.75, 8.1, 9.0 };
            shouldEqualSequenceTolerance(q_onebin, q_onebin+QSIZE, get<StandardQuantiles<UserRangeHistogram<3> > >(a).begin(), 1e-14);
            
            double qauto[QSIZE] = { 0.0, 0.54, 2.175, 3.9, 7.3125, 8.325, 9.0 };
            shouldEqualSequenceTolerance(qauto, qauto+QSIZE, get<StandardQuantiles<AutoRangeHistogram<HSIZE> > >(a).begin(), 1e-15);
            
            double qint[QSIZE] = { 0.0, 0.0, 2.0, 4.0, 7.75, 9.0, 9.0 };
            shouldEqualSequence(qint, qint+QSIZE, get<StandardQuantiles<IntegerHistogram<HSIZE> > >(a).begin());
            shouldEqualSequence(qint, qint+QSIZE, get<StandardQuantiles<IntegerHistogram<0> > >(a).begin());

                // repeat test with negated data => quantiles should be negated, but otherwise the same as before
            a.reset();

            getAccumulator<UserRangeHistogram<HSIZE> >(a).setMinMax(-9.5, 0.5);
            getAccumulator<UserRangeHistogram<HSIZE+2> >(a).setMinMax(-10.5, 1.5);
            getAccumulator<UserRangeHistogram<3> >(a).setMinMax(-30.0, 20.0);
            getAccumulator<IntegerHistogram<0> >(a).setBinCount(HSIZE);

            for(int k=0; k<SIZE; ++k)
                a(-data[k]);

            for(int k=0; k<SIZE; ++k)
                a.updatePass2(-data[k]);

            std::reverse(h, h+HSIZE);
            shouldEqualSequence(h, h+HSIZE, get<UserRangeHistogram<HSIZE> >(a).begin());
            shouldEqualSequence(h, h+HSIZE, get<AutoRangeHistogram<HSIZE> >(a).begin());

            double hneg[HSIZE] = { 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            shouldEqualSequence(hneg, hneg+HSIZE, get<IntegerHistogram<HSIZE> >(a).begin());
            shouldEqualSequence(hneg, hneg+HSIZE, get<IntegerHistogram<0> >(a).begin());

            std::reverse(quser, quser+QSIZE);
            shouldEqualSequenceTolerance(quser, quser+QSIZE, (-get<StandardQuantiles<UserRangeHistogram<HSIZE> > >(a)).begin(), 1e-14);
            
            std::reverse(qauto, qauto+QSIZE);
            shouldEqualSequenceTolerance(qauto, qauto+QSIZE, (-get<StandardQuantiles<AutoRangeHistogram<HSIZE> > >(a)).begin(), 1e-14);

                // repeat test with data shifted by one (test behavior of IntegerHistogram with empty bins at the ends)
            a.reset();

            getAccumulator<UserRangeHistogram<HSIZE> >(a).setMinMax(-0.5, 9.5);
            getAccumulator<UserRangeHistogram<HSIZE+2> >(a).setMinMax(-1.5, 10.5);
            getAccumulator<UserRangeHistogram<3> >(a).setMinMax(-20.0, 30.0);  // all data in one bin
            getAccumulator<IntegerHistogram<0> >(a).setBinCount(HSIZE+2);

            for(int k=0; k<SIZE; ++k)
                a(1+data[k]);

            for(int k=0; k<SIZE; ++k)
                a.updatePass2(1+data[k]);

            shouldEqualSequence(qint, qint+QSIZE, (get<StandardQuantiles<IntegerHistogram<HSIZE+2> > >(a)-QuantileVector(1.0)).begin());
            shouldEqualSequence(qint, qint+QSIZE, (get<StandardQuantiles<IntegerHistogram<0> > >(a)-QuantileVector(1.0)).begin());
        }

        {
            typedef AccumulatorChain<int, Select<UserRangeHistogram<0>, AutoRangeHistogram<0>, IntegerHistogram<0>
                                    > > A;
            A a;

            a.setHistogramOptions(HistogramOptions().setMinMax(-0.5, 9.5).setBinCount(HSIZE));

            shouldEqual(HSIZE, get<UserRangeHistogram<0> >(a).size());
            shouldEqual(HSIZE, get<AutoRangeHistogram<0> >(a).size());
            shouldEqual(HSIZE, get<IntegerHistogram<0> >(a).size());

            extractFeatures(data, data+SIZE, a);

            double h[HSIZE] = { 5.0, 0.0, 6.0, 3.0, 3.0, 1.0, 2.0, 2.0, 4.0, 4.0 };

            shouldEqualSequence(h, h+HSIZE, get<UserRangeHistogram<0> >(a).begin());
            shouldEqualSequence(h, h+HSIZE, get<AutoRangeHistogram<0> >(a).begin());
            shouldEqualSequence(h, h+HSIZE, get<IntegerHistogram<0> >(a).begin());

            try 
            {
                A b;
                extractFeatures(data, data+SIZE, b);
                failTest("extractFeatures() failed to throw exception");
            }
            catch(ContractViolation & c) 
            {
                std::string expected("\nPrecondition violation!\nUserRangeHistogram::update(): setMinMax(...) has not been called.");
                std::string message(c.what());
                shouldEqual(expected, message.substr(0,expected.size()));
            }

            try 
            {
                A b;
                getAccumulator<UserRangeHistogram<0> >(b).setMinMax(-2.5, 12.5);
                failTest("extractFeatures() failed to throw exception");
            }
            catch(ContractViolation & c) 
            {
                std::string expected("\nPrecondition violation!\nRangeHistogramBase::setMinMax(...): setBinCount(...) has not been called.");
                std::string message(c.what());
                shouldEqual(expected, message.substr(0,expected.size()));
            }

            try 
            {
                A b;
                getAccumulator<UserRangeHistogram<0> >(b).setBinCount(HSIZE+2);
                getAccumulator<UserRangeHistogram<0> >(b).setMinMax(-2.5, 12.5);
                extractFeatures(data, data+SIZE, b);
                failTest("extractFeatures() failed to throw exception");
            }
            catch(ContractViolation & c) 
            {
                std::string expected("\nPrecondition violation!\nRangeHistogramBase::setMinMax(...): setBinCount(...) has not been called.");
                std::string message(c.what());
                shouldEqual(expected, message.substr(0,expected.size()));
            }

            try 
            {
                A b;
                b.setHistogramOptions(HistogramOptions().setBinCount(HSIZE));
                getAccumulator<UserRangeHistogram<0> >(b).setMinMax(-2.5, 12.5);

                extractFeatures(data, data+SIZE, b);
                a.merge(b);

                failTest("extractFeatures() failed to throw exception");
            }
            catch(ContractViolation & c) 
            {
                std::string expected("\nPrecondition violation!\nRangeHistogramBase::operator+=(): cannot merge histograms with different data mapping.");
                std::string message(c.what());
                shouldEqual(expected, message.substr(0,expected.size()));
            }

            A b;
            b.setHistogramOptions(HistogramOptions().setBinCount(HSIZE).setMinMax(-0.5, 9.5));

            extractFeatures(data, data+SIZE, b);
            a.merge(b);

            TinyVector<double, HSIZE> h2 = TinyVector<double, HSIZE>(h)*2.0;
            shouldEqualSequence(h2.begin(), h2.end(), get<UserRangeHistogram<0> >(a).begin());
            shouldEqualSequence(h2.begin(), h2.end(), get<AutoRangeHistogram<0> >(a).begin());
            shouldEqualSequence(h2.begin(), h2.end(), get<IntegerHistogram<0> >(a).begin());
        }
    }

    template <class TAG, class A>
    static inline typename acc::LookupDependency<TAG, A>::reference
    getAccumulatorIndirectly(A & a)
    {
        typedef typename acc::LookupDependency<TAG, A>::Tag StandardizedTag;
        typedef typename acc::LookupDependency<TAG, A>::reference reference;
        return acc::detail::CastImpl<StandardizedTag, typename A::Tag, reference>::exec(a);
    }

    void testLabelDispatch()
    {
        using namespace vigra::acc;
        {
            typedef CoupledIteratorType<2, int>::type Iterator;
            typedef Iterator::value_type Handle;
            typedef Shape2 V;

            typedef Select<Count, Coord<Sum>, Global<Count>, Global<Coord<Minimum> >, LabelArg<1>, DataArg<1> > Selected;
            typedef AccumulatorChainArray<Handle, Selected> A;

            should((IsSameType<acc::detail::ConfigureAccumulatorChainArray<Handle, Selected>::GlobalTags, 
                               TypeList<Count,TypeList<Coord<Minimum>,TypeList<DataArg<1>, TypeList<LabelArg<1>, void> > > > >::value));
            should((IsSameType<acc::detail::ConfigureAccumulatorChainArray<Handle, Selected>::RegionTags, 
                               TypeList<Count,TypeList<Coord<Sum>,TypeList<DataArg<1>, void> > > >::value));

            typedef LookupTag<Count, A>::type RegionCount;
            typedef LookupDependency<Global<Count>, RegionCount>::type GlobalCountViaRegionCount;

            should(!(IsSameType<RegionCount, LookupTag<Global<Count>, A>::type>::value));
            should((IsSameType<GlobalCountViaRegionCount, LookupTag<Global<Count>, A>::type>::value));

            MultiArray<2, int> labels(Shape2(3,2));
            labels(2,0) = labels(2,1) = 1;
            Iterator i     = createCoupledIterator(labels),
                     start = i,   
                     end   = i.getEndIterator();

            A a;

            shouldEqual(1, a.passesRequired());

            a.setMaxRegionLabel(1);

            shouldEqual(a.maxRegionLabel(), 1);
            shouldEqual(a.regionCount(), 2);
            should((&getAccumulator<Count, A>(a, 0) != &getAccumulator<Count, A>(a, 1)));
   
            LookupTag<Count, A>::reference rc = getAccumulator<Count>(a, 0);
            LookupTag<Global<Count>, A>::reference gc = getAccumulator<Global<Count> >(a);
            should((&gc == &getAccumulatorIndirectly<Global<Count> >(rc)));
            should((&gc == &getAccumulatorIndirectly<Global<Count> >(getAccumulator<Count>(a, 1))));

            for(; i < end; ++i)
                a(*i);
            
            shouldEqual(4, get<Count>(a, 0));
            shouldEqual(2, get<Count>(a, 1));
            shouldEqual(6, get<Global<Count> >(a));

            shouldEqual(V(2,2), get<Coord<Sum> >(a, 0));
            shouldEqual(V(4,1), get<Coord<Sum> >(a, 1));
            shouldEqual(V(0,0), get<Global<Coord<Minimum> > >(a));

            A b;
            b.ignoreLabel(0);

            i = start;

            for(; i < end; ++i)
                b(*i);
            
            shouldEqual(0, get<Count>(b, 0));
            shouldEqual(2, get<Count>(b, 1));
            shouldEqual(2, get<Global<Count> >(b));

            shouldEqual(V(0,0), get<Coord<Sum> >(b, 0));
            shouldEqual(V(4,1), get<Coord<Sum> >(b, 1));
            shouldEqual(V(2,0), get<Global<Coord<Minimum> > >(b));
        }

        {
            typedef CoupledIteratorType<2, double, int>::type Iterator;
            typedef Iterator::value_type Handle;

            typedef AccumulatorChainArray<Handle, Select<Count, AutoRangeHistogram<3>, GlobalRangeHistogram<3>,
                                                         Global<Count>, Global<AutoRangeHistogram<3> >, DataArg<1>, LabelArg<2>
                                          > > A;

            double d[] = { 1.0, 3.0, 3.0,
                           1.0, 2.0, 5.0 };
            MultiArrayView<2, double> data(Shape2(3,2), d);

            MultiArray<2, int> labels(Shape2(3,2));
            labels(2,0) = labels(2,1) = 1;

            Iterator i     = createCoupledIterator(data, labels),
                     start = i,   
                     end   = i.getEndIterator();

            A a;
            shouldEqual(a.regionCount(), 0);
            shouldEqual(2, a.passesRequired());

            for(; i < end; ++i)
                a(*i);
            
            shouldEqual(a.maxRegionLabel(), 1);
            shouldEqual(a.regionCount(), 2);
            shouldEqual(4, get<Count>(a, 0));
            shouldEqual(2, get<Count>(a, 1));
            shouldEqual(6, get<Global<Count> >(a));

            shouldEqual(1, get<Minimum>(a, 0));
            shouldEqual(3, get<Minimum>(a, 1));
            shouldEqual(1, get<Global<Minimum> >(a));

            shouldEqual(3, get<Maximum>(a, 0));
            shouldEqual(5, get<Maximum>(a, 1));
            shouldEqual(5, get<Global<Maximum> >(a));

            for(i = start; i < end; ++i)
                a.updatePass2(*i);
            
            shouldEqual(4, get<Count>(a, 0));
            shouldEqual(2, get<Count>(a, 1));
            shouldEqual(6, get<Global<Count> >(a));

            typedef TinyVector<double, 3> V;

            shouldEqual(V(2,1,1), get<AutoRangeHistogram<3> >(a, 0));
            shouldEqual(V(1,0,1), get<AutoRangeHistogram<3> >(a, 1));
            shouldEqual(V(3,1,0), get<GlobalRangeHistogram<3> >(a, 0));
            shouldEqual(V(0,1,1), get<GlobalRangeHistogram<3> >(a, 1));
            shouldEqual(V(3,2,1), get<Global<AutoRangeHistogram<3> > >(a));
        }

        {
            typedef CoupledIteratorType<2, double, int>::type Iterator;
            typedef Iterator::value_type Handle;

            typedef DynamicAccumulatorChainArray<Handle, Select<Count, Coord<Mean>, GlobalRangeHistogram<3>,
                                                                AutoRangeHistogram<3>, 
                                                                Global<Count>, Global<Coord<Mean> >, 
                                                                StandardQuantiles<GlobalRangeHistogram<3> >, 
                                                                LabelArg<2>, DataArg<1>
                                                 > > A;

            A a;

            shouldEqual(0, a.passesRequired());

            should(!isActive<Count>(a));
            should(!isActive<Coord<Sum> >(a));
            should(!isActive<GlobalRangeHistogram<3> >(a));

            should(!isActive<Global<Count> >(a));
            should(!isActive<Global<Minimum> >(a));
            should(!isActive<Global<Coord<Sum> > >(a));

            activate<Count>(a);
            should(isActive<Count>(a));
            should(!isActive<Global<Count> >(a));

            //activate<Global<Count> >(a);
            a.activate("Global<PowerSum<0> >");

            should(isActive<Count>(a));
            should(isActive<Global<Count> >(a));

            //activate<Coord<Mean> >(a);
            a.activate("Coord<DivideByCount<PowerSum<1> > >");
            should(isActive<Coord<Mean> >(a));
            should(isActive<Coord<Sum> >(a));
            should(!isActive<Global<Coord<Sum> > >(a));
            should(!isActive<Global<Coord<Mean> > >(a));

            activate<Global<Coord<Mean> > >(a);
            should(isActive<Global<Coord<Sum> > >(a));
            should(isActive<Global<Coord<Mean> > >(a));
            should(!isActive<GlobalRangeHistogram<3> >(a));
            should(!isActive<AutoRangeHistogram<3> >(a));
            should(!isActive<Global<Minimum> >(a));

            shouldEqual(1, a.passesRequired());

            activate<GlobalRangeHistogram<3> >(a);
            a.activate("AutoRangeHistogram<3>");

            should(isActive<GlobalRangeHistogram<3> >(a));
            should(isActive<AutoRangeHistogram<3> >(a));
            should(isActive<Global<Minimum> >(a));

            shouldEqual(2, a.passesRequired());

            MultiArray<2, double> data(Shape2(3,2));
            data(0,0) = 0.1;
            data(2,0) = 1.0;
            data(2,1) = 0.9;
            MultiArray<2, int> labels(Shape2(3,2));
            labels(2,0) = labels(2,1) = 1;
            Iterator i     = createCoupledIterator(data, labels),
                     start = i,   
                     end   = i.getEndIterator();

            for(; i < end; ++i)
                a(*i);
            
            for(i = start; i < end; ++i)
                a.updatePass2(*i);
            
            shouldEqual(a.maxRegionLabel(), 1);
            shouldEqual(4, get<Count>(a, 0));
            shouldEqual(2, get<Count>(a, 1));
            shouldEqual(6, get<Global<Count> >(a));

            typedef TinyVector<double, 2> V;

            shouldEqual(V(0.5, 0.5), get<Coord<Mean> >(a, 0));
            shouldEqual(V(2, 0.5), get<Coord<Mean> >(a, 1));
            shouldEqual(V(1, 0.5), get<Global<Coord<Mean> > >(a));

            should(getAccumulator<GlobalRangeHistogram<3> >(a,0).scale_ == getAccumulator<GlobalRangeHistogram<3> >(a,1).scale_);
            should(getAccumulator<GlobalRangeHistogram<3> >(a,0).scale_ != getAccumulator<AutoRangeHistogram<3> >(a,0).scale_);
            should(getAccumulator<GlobalRangeHistogram<3> >(a,1).scale_ != getAccumulator<AutoRangeHistogram<3> >(a,1).scale_);
            
            typedef TinyVector<double, 3> W;
            shouldEqual(W(4, 0, 0), get<GlobalRangeHistogram<3> >(a,0));
            shouldEqual(W(0, 0, 2), get<GlobalRangeHistogram<3> >(a,1));
            shouldEqual(W(3, 0, 1), get<AutoRangeHistogram<3> >(a,0));
            shouldEqual(W(1, 0, 1), get<AutoRangeHistogram<3> >(a,1));

            A b;
            b.activateAll();

            extractFeatures(start, end, b);
            
            shouldEqual(W(4, 0, 0), get<GlobalRangeHistogram<3> >(b,0));
            shouldEqual(W(0, 0, 2), get<GlobalRangeHistogram<3> >(b,1));

            a += b;
            
            shouldEqual(a.maxRegionLabel(), 1);
            shouldEqual(8, get<Count>(a, 0));
            shouldEqual(4, get<Count>(a, 1));
            shouldEqual(12, get<Global<Count> >(a));
            shouldEqual(V(0.5, 0.5), get<Coord<Mean> >(a, 0));
            shouldEqual(V(2, 0.5), get<Coord<Mean> >(a, 1));
            shouldEqual(V(1, 0.5), get<Global<Coord<Mean> > >(a));
            shouldEqual(W(8, 0, 0), get<GlobalRangeHistogram<3> >(a,0));
            shouldEqual(W(0, 0, 4), get<GlobalRangeHistogram<3> >(a,1));
            shouldEqual(W(4, 0, 0), get<GlobalRangeHistogram<3> >(b,0));
            shouldEqual(W(0, 0, 2), get<GlobalRangeHistogram<3> >(b,1));

            TinyVector<int, 2> labelMapping(2, 3);
            a.merge(b, labelMapping);
            shouldEqual(a.maxRegionLabel(), 3);
            shouldEqual(8, get<Count>(a, 0));
            shouldEqual(4, get<Count>(a, 1));
            shouldEqual(4, get<Count>(a, 2));
            shouldEqual(2, get<Count>(a, 3));
            shouldEqual(18, get<Global<Count> >(a));
            shouldEqual(V(0.5, 0.5), get<Coord<Mean> >(a, 0));
            shouldEqual(V(2, 0.5), get<Coord<Mean> >(a, 1));
            shouldEqual(V(0.5, 0.5), get<Coord<Mean> >(a, 2));
            shouldEqual(V(2, 0.5), get<Coord<Mean> >(a, 3));
            shouldEqual(V(1, 0.5), get<Global<Coord<Mean> > >(a));

            A c;
            c.activateAll();
            c.setHistogramOptions(HistogramOptions().regionAutoInit());
            extractFeatures(start, end, c);

            shouldEqual(getAccumulator<GlobalRangeHistogram<3> >(c,0).scale_, getAccumulator<AutoRangeHistogram<3> >(c,0).scale_);
            shouldEqual(getAccumulator<GlobalRangeHistogram<3> >(c,1).scale_, getAccumulator<AutoRangeHistogram<3> >(c,1).scale_);
            
            shouldEqual(W(3, 0, 1), get<GlobalRangeHistogram<3> >(c,0));
            shouldEqual(W(1, 0, 1), get<GlobalRangeHistogram<3> >(c,1));
            shouldEqual(W(3, 0, 1), get<AutoRangeHistogram<3> >(c,0));
            shouldEqual(W(1, 0, 1), get<AutoRangeHistogram<3> >(c,1));

            c.merge(c, TinyVector<int, 2>(3, 2));

            shouldEqual(c.maxRegionLabel(), 3);
            shouldEqual(get<Count>(c, 0), 4);
            shouldEqual(get<Count>(c, 1), 2);
            shouldEqual(get<Count>(c, 2), 2);
            shouldEqual(get<Count>(c, 3), 4);
            shouldEqual(get<Global<Count> >(c), 12);

            shouldEqual(W(3, 0, 1), get<AutoRangeHistogram<3> >(c,0));
            shouldEqual(W(1, 0, 1), get<AutoRangeHistogram<3> >(c,1));
            shouldEqual(W(1, 0, 1), get<AutoRangeHistogram<3> >(c,2));
            shouldEqual(W(3, 0, 1), get<AutoRangeHistogram<3> >(c,3));

            c.merge(1, 2);

            shouldEqual(c.maxRegionLabel(), 3);
            shouldEqual(get<Count>(c, 0), 4);
            shouldEqual(get<Count>(c, 1), 4);
            shouldEqual(get<Count>(c, 2), 0);
            shouldEqual(get<Count>(c, 3), 4);
            shouldEqual(get<Global<Count> >(c), 12);

            shouldEqual(W(3, 0, 1), get<AutoRangeHistogram<3> >(c,0));
            shouldEqual(W(2, 0, 2), get<AutoRangeHistogram<3> >(c,1));
            shouldEqual(W(0, 0, 0), get<AutoRangeHistogram<3> >(c,2));
            shouldEqual(W(3, 0, 1), get<AutoRangeHistogram<3> >(c,3));
        }
    }
};

struct FeaturesTestSuite : public vigra::test_suite
{
    FeaturesTestSuite()
        : vigra::test_suite("FeaturesTestSuite")
    {
        add(testCase(&AccumulatorTest::testStandardizeTag));
        add(testCase(&AccumulatorTest::testTagTransfer));
        add(testCase(&AccumulatorTest::testScalar));
        add(testCase(&AccumulatorTest::testVector));
        add(testCase(&AccumulatorTest::testMerge));
        add(testCase(&AccumulatorTest::testCoordAccess));
        add(testCase(&AccumulatorTest::testHistogram));
        add(testCase(&AccumulatorTest::testLabelDispatch));
        add(testCase(&AccumulatorTest::testIndexSpecifiers));
    }
};

int main(int argc, char** argv)
{
    FeaturesTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
