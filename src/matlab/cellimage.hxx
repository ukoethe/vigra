/************************************************************************/
/*                                                                      */
/*        Copyright 2003-2009 by Hans Meine and Ullrich Koethe          */
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

#ifndef CELLIMAGE_HXX
#define CELLIMAGE_HXX

#include <vigra/basicimage.hxx>
#include <vigra/pixelneighborhood.hxx>
#include <functional>

namespace vigra {

namespace cellimage {

enum CellType { CellTypeRegion = 0,
                CellTypeLine = 1,
                CellTypeVertex = 2,
                CellTypeError = 3,
                CellTypeVertexOrLine = 4,
                CellTypeErrorOrLine = 5 };

// -------------------------------------------------------------------
//                          CellPixel, CellImage
// -------------------------------------------------------------------
typedef unsigned int CellLabel;

struct CellPixel
{
private:
    CellLabel typeLabel_;

    friend struct CellPixelSerializer;

public:
    CellPixel() {}
    CellPixel(CellType type, CellLabel label = 0)
    : typeLabel_((label << 2) | type)
    {}

    inline CellType type() const
        { return (CellType)(typeLabel_ & 3); }
    inline void setType(CellType type)
        { typeLabel_ = (label() << 2) | type; }
    inline void setType(CellType type, CellLabel label)
        { typeLabel_ = label << 2 | type; }

    inline CellLabel label() const
        { return typeLabel_ >> 2; }
    inline void setLabel(CellLabel label)
        { typeLabel_ = label << 2 | type(); }
    inline void setLabel(CellLabel label, CellType type)
        { typeLabel_ = label << 2 | type; }

    bool operator==(CellPixel const & rhs) const
        { return typeLabel_ == rhs.typeLabel_; }

    bool operator!=(CellPixel const & rhs) const
        { return typeLabel_ != rhs.typeLabel_; }
};

typedef BasicImage<CellPixel> CellImage;

typedef vigra::NeighborhoodCirculator<CellImage::Iterator, EightNeighborCode>
    CellImageEightCirculator;

// -------------------------------------------------------------------
//                     CellPixel Serialization
// -------------------------------------------------------------------
struct CellPixelSerializer
{
    int operator()(CellPixel const &p) const
    {
        return p.typeLabel_;
    }

    CellPixel operator()(int i) const
    {
        CellPixel result;
        result.typeLabel_ = i;
        return result;
    }
};

// -------------------------------------------------------------------
//                   CellPixel/CellImage Accessors
// -------------------------------------------------------------------
template<class VALUE_TYPE = CellType>
struct TypeAccessor
{
    typedef VALUE_TYPE value_type;
    typedef VALUE_TYPE result_type;

    template<class Iterator>
    value_type operator()(const Iterator &it) const
    {
        return it->type();
    }

    template<class Iterator>
    void set(value_type type, const Iterator &it) const
    {
        it->setType(type);
    }
};

typedef TypeAccessor<unsigned char> TypeAsByteAccessor;

typedef TypeAccessor<> CellTypeAccessor;

struct LabelAccessor
{
    typedef CellLabel value_type;

    template<class Iterator>
    CellLabel operator()(const Iterator &it) const
    {
        return it->label();
    }

    template<class Iterator>
    void set(CellLabel label, const Iterator &it) const
    {
        it->setLabel(label);
    }
};

template<CellType type>
struct LabelWriter
{
    typedef CellLabel value_type;

    template<class Iterator>
    void set(CellLabel label, const Iterator &it) const
    {
        it->setLabel(label, type);
    }
};

template<CellType type>
struct CellTypeEquals : public std::unary_function<CellType, bool>
{
    bool operator()(CellType t) const
    {
        return t == type;
    }

    template<class Iterator>
    bool operator()(const Iterator &it) const
    {
        return it->type() == type;
    }
};

struct CellMask : public std::unary_function<vigra::cellimage::CellPixel, bool>
{
    vigra::cellimage::CellPixel maskPixel_;

    CellMask(vigra::cellimage::CellPixel maskPixel)
    : maskPixel_(maskPixel)
    {}

    template<class Iterator>
    bool operator()(const Iterator &it) const
    {
        return *it == maskPixel_;
    }
};

// -------------------------------------------------------------------
//                        RelabelFunctor (unused!)
// -------------------------------------------------------------------
template<class VALUETYPE>
struct RelabelFunctor
{
    typedef VALUETYPE value_type;
    typedef VALUETYPE argument_type;
    typedef VALUETYPE result_type;

    RelabelFunctor(VALUETYPE oldValue, VALUETYPE newValue)
        : oldValue_(oldValue),
          newValue_(newValue)
    {}

    VALUETYPE operator()(VALUETYPE value) const
    {
        return (value == oldValue_) ? newValue_ : value;
    }

    VALUETYPE oldValue_, newValue_;
};

// -------------------------------------------------------------------
//                              inspectCell
// -------------------------------------------------------------------
// thinking about "typename IteratorTraits<EndIterator>::DefaultAccessor":
// is not needed since we're not implementing srcCellRange here, but
// algorithms.
// srcCellRange can not be implemented that easy, because most VIGRA
// functions expect an ImageIterator, not a std::iterator
template <class EndIterator, class Accessor, class Functor>
void inspectCell(EndIterator endIterator, Accessor a, Functor & f)
{
    for(; endIterator.inRange(); ++endIterator)
        f(a(endIterator));
}

template <class EndIterator, class Functor>
void inspectCell(EndIterator endIterator, Functor & f)
{
    for(; endIterator.inRange(); ++endIterator)
        f(*endIterator);
}

// -------------------------------------------------------------------
//                             transformCell
// -------------------------------------------------------------------
template <class SrcEndIterator, class SrcAccessor,
          class DestEndIterator, class DestAccessor, class Functor>
void transformCell(SrcEndIterator srcEndIterator, SrcAccessor sa,
                   DestEndIterator destEndIterator, DestAccessor da,
                   Functor const & f)
{
    for(; srcEndIterator.inRange(); ++srcEndIterator, ++destEndIterator)
        da.set(f(sa(srcEndIterator)), destEndIterator);
}

template <class SrcEndIterator, class DestEndIterator, class Functor>
void transformCell(SrcEndIterator srcEndIterator,
                   DestEndIterator destEndIterator,
                   Functor const & f)
{
    for(; srcEndIterator.inRange(); ++srcEndIterator, ++destEndIterator)
        *destEndIterator = f(*srcEndIterator);
}

} // namespace cellimage

} // namespace vigra

#endif // CELLIMAGE_HXX
