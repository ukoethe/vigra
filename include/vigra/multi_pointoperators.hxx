//-- -*- c++ -*-
/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Christian-Dennis Rahn                */
/*                        and Ullrich Koethe                            */
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

#ifndef VIGRA_MULTI_POINTOPERATORS_H
#define VIGRA_MULTI_POINTOPERATORS_H

 //#define MC_SHOW_CONVOLUTION       //-- should not be uncommented


#include <vigra/initimage.hxx>
#include <vigra/copyimage.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/metaprogramming.hxx>



namespace vigra
{

/********************************************************/
/*                                                      */
/*                    initMultiarray                    */
/*                                                      */
/********************************************************/

template <class Iterator, class Accessor, class VALUETYPE>
inline void
initMultiarrayImpl(Iterator s, Iterator send, Accessor a,  VALUETYPE v, MetaInt<0>)
{
    initLine(s, send, a, v);
}
    
template <class Iterator, class Accessor, class VALUETYPE, int N>
void
initMultiarrayImpl(Iterator s, Iterator send, Accessor a,  VALUETYPE v, MetaInt<N>)
{
    for(; s != send; ++s)
    {
        initMultiarrayImpl(s.begin(), s.end(), a, v, MetaInt<N-1>());
    }
}
    
template <class Iterator, class Accessor, class VALUETYPE>
inline void
initMultiarray(Iterator s, Iterator send, Accessor a,  VALUETYPE v)
{
    initMultiarrayImpl(s, send, a, v, MetaInt<Iterator::level>());
}
    
template <class Iterator, class Shape, class Accessor, class VALUETYPE>
inline 
void
initMultiarray(triple<Iterator, Shape, Accessor> const & s, VALUETYPE v)
{
    initMultiarray(s.first, s.first + s.second[s.second.size()-1], s.third, v);
}
    
/********************************************************/
/*                                                      */
/*                    copyMultiarray                    */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiarrayImpl(SrcIterator s, SrcIterator send, SrcAccessor src,
               DestIterator d, DestAccessor dest, MetaInt<0>)
{
    copyLine(s, send, src, d, dest);
}
    
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, int N>
void
copyMultiarrayImpl(SrcIterator s, SrcIterator send, SrcAccessor src,
                   DestIterator d, DestAccessor dest, MetaInt<N>)
{
    for(; s != send; ++s, ++d)
    {
        copyMultiarrayImpl(s.begin(), s.end(), src, d.begin(), dest, MetaInt<N-1>());
    }
}
    
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiarray(SrcIterator s, 
               SrcIterator send, SrcAccessor src,
               DestIterator d, DestAccessor dest)
{    
    copyMultiarrayImpl(s, send, src, d, dest, MetaInt<SrcIterator::level>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiarray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               pair<DestIterator, DestAccessor> const & dest)
{
    
    copyMultiarray(src.first, src.first + src.second[src.second.size()-1], src.third, dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                 transformMultiarray                  */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Functor>
inline void
transformMultiarrayImpl(SrcIterator s, SrcIterator send, SrcAccessor src,
               DestIterator d, DestAccessor dest, Functor const & f, MetaInt<0>)
{
    transformLine(s, send, src, d, dest, f);
}
    
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class Functor, int N>
void
transformMultiarrayImpl(SrcIterator s, SrcIterator send, SrcAccessor src,
                   DestIterator d, DestAccessor dest, 
                   Functor const & f, MetaInt<N>)
{
    for(; s != send; ++s, ++d)
    {
        transformMultiarrayImpl(s.begin(), s.end(), src, d.begin(), dest, 
                                f, MetaInt<N-1>());
    }
}
    
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Functor>
inline void
transformMultiarray(SrcIterator s, SrcIterator send, SrcAccessor src,
                    DestIterator d, DestAccessor dest, Functor const & f)
{    
    transformMultiarrayImpl(s, send, src, d, dest, f, MetaInt<SrcIterator::level>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Functor>
inline void
transformMultiarray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               pair<DestIterator, DestAccessor> const & dest, Functor const & f)
{
    
    transformMultiarray(src.first, src.first + src.second[src.second.size()-1], src.third, 
                        dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*                combineTwoMultiarrays                 */
/*                                                      */
/********************************************************/

template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor,
          class Functor>
inline void
combineTwoMultiarraysImpl(SrcIterator1 s1, SrcIterator1 s1end, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               DestIterator d, DestAccessor dest, Functor const & f, MetaInt<0>)
{
    combineTwoLines(s1, s1end, src1, s2, src2, d, dest, f);
}
    
template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor, 
          class Functor, int N>
void
combineTwoMultiarraysImpl(SrcIterator1 s1, SrcIterator1 s1end, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               DestIterator d, DestAccessor dest, 
                   Functor const & f, MetaInt<N>)
{
    for(; s1 != s1end; ++s1, ++s2, ++d)
    {
        combineTwoMultiarraysImpl(s1.begin(), s1.end(), src1, 
                                  s2.begin(), src2, d.begin(), dest, 
                                  f, MetaInt<N-1>());
    }
}
    
template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor,
          class Functor>
inline void
combineTwoMultiarrays(SrcIterator1 s1, SrcIterator1 s1end, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               DestIterator d, DestAccessor dest, Functor const & f)
{    
    combineTwoMultiarraysImpl(s1, s1end, src1, s2, src2, d, dest, f, 
                              MetaInt<SrcIterator1::level>());
}

template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor, class Functor>
inline void
combineTwoMultiarrays(triple<SrcIterator1, SrcShape, SrcAccessor1> const & src1,
               pair<SrcIterator2, SrcAccessor2> const & src2,
               pair<DestIterator, DestAccessor> const & dest, Functor const & f)
{
    
    combineTwoMultiarrays(
           src1.first, src1.first + src1.second[src1.second.size()-1], src1.third, 
           src2.first, src2.second, dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*               combineThreeMultiarrays                */
/*                                                      */
/********************************************************/

template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class SrcIterator3, class SrcAccessor3,
          class DestIterator, class DestAccessor,
          class Functor>
inline void
combineThreeMultiarraysImpl(SrcIterator1 s1, SrcIterator1 s1end, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               SrcIterator3 s3, SrcAccessor3 src3,
               DestIterator d, DestAccessor dest, Functor const & f, MetaInt<0>)
{
    combineThreeLines(s1, s1end, src1, s2, src2, s3, src3, d, dest, f);
}
    
template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class SrcIterator3, class SrcAccessor3,
          class DestIterator, class DestAccessor, 
          class Functor, int N>
void
combineThreeMultiarraysImpl(SrcIterator1 s1, SrcIterator1 s1end, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               SrcIterator3 s3, SrcAccessor3 src3,
               DestIterator d, DestAccessor dest, 
                   Functor const & f, MetaInt<N>)
{
    for(; s1 != s1end; ++s1, ++s2, ++s3, ++d)
    {
        combineThreeMultiarraysImpl(s1.begin(), s1.end(), src1, 
                                  s2.begin(), src2, s3.begin(), src3, d.begin(), dest, 
                                  f, MetaInt<N-1>());
    }
}
    
template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class SrcIterator3, class SrcAccessor3,
          class DestIterator, class DestAccessor,
          class Functor>
inline void
combineThreeMultiarrays(SrcIterator1 s1, SrcIterator1 s1end, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               SrcIterator3 s3, SrcAccessor3 src3,
               DestIterator d, DestAccessor dest, Functor const & f)
{    
    combineThreeMultiarraysImpl(s1, s1end, src1, s2, src2, s3, src3, d, dest, f, 
                              MetaInt<SrcIterator1::level>());
}

template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class SrcIterator3, class SrcAccessor3,
          class DestIterator, class DestAccessor, class Functor>
inline void
combineThreeMultiarrays(triple<SrcIterator1, SrcShape, SrcAccessor1> const & src1,
               pair<SrcIterator2, SrcAccessor2> const & src2,
               pair<SrcIterator3, SrcAccessor3> const & src3,
               pair<DestIterator, DestAccessor> const & dest, Functor const & f)
{
    
    combineThreeMultiarrays(
           src1.first, src1.first + src1.second[src1.second.size()-1], src1.third, 
           src2.first, src2.second, src3.first, src3.second, dest.first, dest.second, f);
}

};	//-- namespace vigra


#endif	//-- VIGRA_MULTI_POINTOPERATORS_H
