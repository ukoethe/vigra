/************************************************************************/
/*                                                                      */
/*            Copyright 2011-2012 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyanalysis_PyArray_API
#define NO_IMPORT_ARRAY

#include "pythonaccumulator.hxx"

namespace python = boost::python;

namespace vigra
{

namespace acc1 
{

AliasMap defineAliasMap()
{
    AliasMap res;
    res["Coord<DivideByCount<PowerSum<1> > >"] = "RegionCenter";
    res["Coord<RootDivideByCount<Principal<PowerSum<2> > > >"] = "RegionRadii";
    res["Coord<Principal<CoordinateSystem> >"] = "RegionAxes";
    res["DivideByCount<Central<PowerSum<2> > >"] = "Variance";
    res["DivideUnbiased<Central<PowerSum<2> > >"] = "UnbiasedVariance";
    res["DivideByCount<Principal<PowerSum<2> > >"] = "Principal<Variance>";
    res["DivideByCount<FlatScatterMatrix>"] = "Covariance";
    res["DivideByCount<PowerSum<1> >"] = "Mean";
    res["PowerSum<1>"] = "Sum";
    res["PowerSum<0>"] = "Count";
    res["Principal<CoordinateSystem>"] = "PrincipalAxes";
    res["AutoRangeHistogram<0>"] = "Histogram";
    res["GlobalRangeHistogram<0>"] = "Histogram";
    res["StandardQuantiles<AutoRangeHistogram<0> >"] = "Quantiles";
    res["StandardQuantiles<GlobalRangeHistogram<0> >"] = "Quantiles";
    res["Weighted<Coord<DivideByCount<PowerSum<1> > > >"] = "Weighted<RegionCenter>";
    res["Weighted<Coord<RootDivideByCount<Principal<PowerSum<2> > > > >"] = "Weighted<RegionRadii>";
    res["Weighted<Coord<Principal<CoordinateSystem> > >"] = "Weighted<RegionAxes>";
    return res;
}

AliasMap createTagToAlias(ArrayVector<std::string> const & names)
{
    static const AliasMap aliases = defineAliasMap();
    AliasMap res;
    for(unsigned int k=0; k<names.size(); ++k)
    {
            // lookup alias names
        AliasMap::const_iterator a = aliases.find(names[k]);
        std::string alias = (a == aliases.end())
                               ? names[k]
                               : a->second;
                               
            // treat FlatScatterMatrix and ScatterMatrixEigensystem as internal,
            // i.e. use names only when they don't contain these strings
        if(alias.find("ScatterMatrixEigensystem") == std::string::npos &&
           alias.find("FlatScatterMatrix") == std::string::npos)
             res[names[k]] = alias;
    }
    return res;   
}

AliasMap createAliasToTag(AliasMap const & tagToAlias)
{
    AliasMap res;
    for(AliasMap::const_iterator k = tagToAlias.begin(); k != tagToAlias.end(); ++k)
        res[normalizeString(k->second)] = normalizeString(k->first);
    return res;
}

ArrayVector<std::string> createSortedNames(AliasMap const & tagToAlias)
{
    ArrayVector<std::string> res;
    for(AliasMap::const_iterator k = tagToAlias.begin(); k != tagToAlias.end(); ++k)
        res.push_back(k->second);
    std::sort(res.begin(), res.end());
    return res;
}

} // namespace acc1

void defineGlobalAccumulators()
{
    using namespace python;
    using namespace vigra::acc1;

    docstring_options doc_options(true, true, false);
    
    PythonAccumulatorBase::definePythonClass();
    PythonRegionAccumulatorBase::definePythonClass();
    
    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, Covariance, 
                   Principal<Variance>, Principal<Skewness>, Principal<Kurtosis>,
                   Principal<CoordinateSystem>,
                   Minimum, Maximum, Principal<Minimum>, Principal<Maximum>
                   > VectorAccumulators;

    definePythonAccumulatorMultiband<3, float, VectorAccumulators>();
    definePythonAccumulatorMultiband<4, float, VectorAccumulators>();
    
    definePythonAccumulator<TinyVector<float, 3>, VectorAccumulators>();

    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, 
                   UnbiasedVariance, UnbiasedSkewness, UnbiasedKurtosis,
                   Minimum, Maximum, StandardQuantiles<AutoRangeHistogram<0> > 
                   > ScalarAccumulators;
    definePythonAccumulatorSingleband<float, ScalarAccumulators>();
}

void defineSinglebandRegionAccumulators();
void defineMultibandRegionAccumulators();

void defineAccumulators()
{
    NumpyArrayConverter<NumpyArray<1, npy_uint32> >();
    NumpyArrayConverter<NumpyArray<1, float> >();
    NumpyArrayConverter<NumpyArray<1, double> >();
    NumpyArrayConverter<NumpyArray<2, MultiArrayIndex> >();
    NumpyArrayConverter<NumpyArray<2, float> >();
    NumpyArrayConverter<NumpyArray<2, double> >();
    NumpyArrayConverter<NumpyArray<3, float> >();
    NumpyArrayConverter<NumpyArray<3, double> >();
    
    defineGlobalAccumulators();
    defineMultibandRegionAccumulators();
    defineSinglebandRegionAccumulators();
}

// TODO:
//  * nested Select
//  * Multiband support
//  * implement PythonAccumulatorArray::merge()
//  * check that merge skips inactive accumulators
//  * implement label remapping in merge()
//  * is there a good implementation of merge for histogramms with different mapping?
//  * multiband histograms
//  * ensure that accumulators promote float arguments to double
//  * general refactoring
//  * better names for PrincipalRadii, PrincipalCoordSystem, MomentsOfInertia, CoordSystemOfInertia
//  * tests and docu
//  * speed-up compilation of Python bindings

} // namespace vigra
