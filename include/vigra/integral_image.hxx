/************************************************************************/
/*                                                                      */
/*    Copyright 2012-2013 by Ullrich Koethe and Anna Kreshuk            */
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

#ifndef VIGRA_INTEGRALIMAGE_HXX
#define VIGRA_INTEGRALIMAGE_HXX

#include <vigra/multi_array.hxx>
#include <vigra/multi_pointoperators.hxx>
#include <vigra/utilities.hxx>
#include <vigra/functorexpression.hxx>

namespace vigra {

template <unsigned int N, class T1, class S1, class T2, class S2, class FUNCTOR>
void 
cumulativeSum(MultiArrayView<N, T1, S1> const & image, 
              MultiArrayView<N, T2, S2> out,
              int axis,
              FUNCTOR const & functor)
{
    typedef typename MultiArrayShape<N>::type ShapeN;
    ShapeN offset = ShapeN::unitVector(axis);
    
    MultiCoordinateIterator<N> i(image.shape()),
                               end(i.getEndIterator());
    for(; i != end; ++i)
    {
        if((*i)[axis] == 0)
        {
            out[*i] = functor(image[*i]);
        }
        else
        {
            out[*i] = functor(image[*i]) + out[*i - offset];
        }
    }
}

template <unsigned int N, class T1, class S1, class T2, class S2, class FUNCTOR>
void 
integralMultiArrayImpl(MultiArrayView<N, T1, S1> const & array, 
                       MultiArrayView<N, T2, S2> intarray,
                       FUNCTOR const & f)
{
    vigra_precondition(array.shape() == intarray.shape(),
        "integralMultiArray(): shape mismatch between input and output.");
        
    cumulativeSum(array, intarray, 0, f);
    
    for(int axis=1; axis < N; ++axis)
        cumulativeSum(intarray, intarray, axis, functor::Identity());
}

template <class T1, class S1, class T2, class S2, class FUNCTOR>
void 
integralMultiArrayImpl(MultiArrayView<2, T1, S1> const & image, 
                       MultiArrayView<2, T2, S2> intimage,
                       FUNCTOR const & functor)
{
    vigra_precondition(image.shape() == intimage.shape(),
        "integralMultiArray(): shape mismatch between input and output.");
        
    int width = image.shape(0);
    int height = image.shape(1);

    T2 s = T2();
    for (int x=0; x<width; ++x)
    {
        s += functor(image(x, 0));
        intimage(x, 0) = s;
    }
    for (int y=1; y<height; ++y)
    {
        s = T2();
        for (int x=0; x<width; ++x)
        {
            s += functor(image(x, y));
            intimage(x, y) = s + intimage(x, y-1);
        }
    }
}

template <class T1, class S1, class T2, class S2, class FUNCTOR>
void 
integralMultiArrayImpl(MultiArrayView<3, T1, S1> const & volume, 
                       MultiArrayView<3, T2, S2> intvolume,
                       FUNCTOR const & functor)
{
    vigra_precondition(volume.shape() == intvolume.shape(),
        "integralMultiArray(): shape mismatch between input and output.");
        
    int nx = volume.shape(0);
    int ny = volume.shape(1);
    int nz = volume.shape(2);

    //this vector will store s2(x, y-1) for all values of x
    MultiArray<1, T2> s2_temp(nx);
    
    T2 s1 = T2();
    T2 s2 = T2();

    for (int iy=0; iy<ny; ++iy)
    {
        s1 = T2();
        for (int ix=0; ix<nx; ++ix)
        {
            s1 += functor(volume(ix, iy, 0));
            s2 = s2_temp(ix) + s1;
            s2_temp(ix) = s2;
            intvolume(ix, iy, 0) = s2;
        }
    }
    
    for (int iz=1; iz<nz; ++iz)
    {    
        s2_temp = T2();
        
        for (int iy=0; iy<ny; ++iy)
        {
            s1 = T2();
            for (int ix=0; ix<nx; ++ix)
            {
                s1 += functor(volume(ix, iy, iz));
                s2 = s2_temp(ix) + s1;
                s2_temp(ix) = s2;
                intvolume(ix, iy, iz) = s2 + intvolume(ix, iy, iz-1);
            }
        }
    }
}

template <unsigned int N, class T1, class S1, class T2, class S2, class FUNCTOR>
inline void 
integralMultiArray(MultiArrayView<N, T1, S1> const & array, 
                   MultiArrayView<N, T2, S2> intarray,
                   FUNCTOR const & f)
{
    integralMultiArrayImpl(array, intarray, f);
}

template <unsigned int N, class T1, class S1, class T2, class S2, class FUNCTOR>
inline void 
integralMultiArray(MultiArrayView<N, Multiband<T1>, S1> const & array, 
                   MultiArrayView<N, Multiband<T2>, S2> intarray,
                   FUNCTOR const & f)
{
    for(int channel=0; channel < array.shape(N-1); ++channel)
        integralMultiArrayImpl(array.bindOuter(channel), intarray.bindOuter(channel), f);
}

template <unsigned int N, class T1, class S1, class T2, class S2>
inline void 
integralMultiArray(MultiArrayView<N, T1, S1> const & array, 
                   MultiArrayView<N, T2, S2> intarray)
{
    integralMultiArray(array, intarray, functor::Identity());
}

template <unsigned int N, class T1, class S1, class T2, class S2>
inline void 
integralMultiArray(MultiArrayView<N, Multiband<T1>, S1> const & array, 
                   MultiArrayView<N, Multiband<T2>, S2> intarray)
{
    integralMultiArray(array, intarray, functor::Identity());
}

template <unsigned int N, class T1, class S1, class T2, class S2>
inline void 
integralMultiArraySquared(MultiArrayView<N, T1, S1> const & array, 
                          MultiArrayView<N, T2, S2> intarray)
{
    using namespace functor;
    integralMultiArray(array, intarray, sq(Arg1()));
}

template <unsigned int N, class T1, class S1, class T2, class S2>
inline void 
integralMultiArraySquared(MultiArrayView<N, Multiband<T1>, S1> const & array, 
                          MultiArrayView<N, Multiband<T2>, S2> intarray)
{
    using namespace functor;
    integralMultiArray(array, intarray, sq(Arg1()));
}

} // namespace vigra
    
#endif // VIGRA_INTEGRALIMAGE_HXX
