/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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


#ifndef VIGRA_IMAGEITERATORADAPTER_HXX
#define VIGRA_IMAGEITERATORADAPTER_HXX

#include <iterator>   // iterator tags

namespace vigra {

/** \addtogroup ImageIteratorAdapters Image Iterator Adapters

     Iterate over rows, columns, neighborhoods, contours, and other image subsets
*/
//@{

/********************************************************/
/*                                                      */
/*                      ColumnIterator                  */
/*                                                      */
/********************************************************/

/** \brief Iterator adapter to linearly access colums.

    This iterator may be initialized from any standard ImageIterator,
    a MultibandImageIterator and so on.
    It gives you STL-compatible (random access iterator) access to
    one column of the image. If the underlying iterator is a const iterator,
    the column iterator will also be const (i.e. doesn't allow to change
    the values it points to).
    The iterator gets associated with the accessor of the base iterator.

    Note that image iterators usually have a member <TT>columnIterator()</TT>
    which returns a column iterator optimized for that particular image class.
    ColumnIterator is only necessary if this 'native' column iterator
    is not usable in a particular situation or is not provided.

    <b>\#include</b> \<vigra/imageiteratoradapter.hxx\>

    Namespace: vigra

*/
template <class IMAGE_ITERATOR>
class ColumnIterator : private IMAGE_ITERATOR
{
  public:
        /** the iterator's value type
        */
    typedef typename IMAGE_ITERATOR::value_type value_type;

        /** the iterator's value type
        */
    typedef typename IMAGE_ITERATOR::value_type PixelType;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef typename IMAGE_ITERATOR::reference              reference;

        /** the iterator's index reference type (return type of <TT>iter[n]</TT>)
        */
    typedef typename IMAGE_ITERATOR::index_reference        index_reference;

        /** the iterator's pointer type (return type of <TT>iter.operator->()</TT>)
        */
    typedef typename IMAGE_ITERATOR::pointer                pointer;

        /** the iterator's difference type (argument type of <TT>iter[diff]</TT>)
        */
    typedef typename IMAGE_ITERATOR::difference_type::MoveY difference_type;

        /** the iterator tag (random access iterator)
        */
    typedef std::random_access_iterator_tag                 iterator_category;

        /** the type of the adapted iterator
        */
    typedef IMAGE_ITERATOR Adaptee;

        /** Construct from an the image iterator to be adapted.
       */
    ColumnIterator(IMAGE_ITERATOR  const & i)
    : IMAGE_ITERATOR(i)
    {}

        /** Assignment.
        */
    ColumnIterator & operator=(ColumnIterator  const & i)
    {
        IMAGE_ITERATOR::operator=(i);

        return *this;
    }

        /** Assign a new base iterator.
        */
    ColumnIterator & operator=(IMAGE_ITERATOR  const & i)
    {
        IMAGE_ITERATOR::operator=(i);

        return *this;
    }

    /** @name Navigation */
    //@{
        ///
    ColumnIterator &  operator++()
    {
        ++(this->y);
        return *this;
    }
        ///
    ColumnIterator  operator++(int)
    {
        ColumnIterator ret(*this);
        (this->y)++;
        return ret;
    }

        ///
    ColumnIterator &  operator--()
    {
        --(this->y);
        return *this;
    }

        ///
    ColumnIterator  operator--(int)
    {
        ColumnIterator ret(*this);
        (this->y)--;
        return ret;
    }

        ///
    ColumnIterator &  operator+=(int d)
    {
        this->y += d;
        return *this;
    }

        ///
    ColumnIterator &  operator-=(int d)
    {
        this->y -= d;
        return *this;
    }
    //@}

    /** @name Methods */
    //@{
        /** Construct iterator at a distance.
        */
    ColumnIterator operator+(int d) const
    {
        IMAGE_ITERATOR ret(*this);
        ret.y += d;
        return ColumnIterator(ret);
    }
        /** Construct iterator at a distance.
        */
    ColumnIterator operator-(int d) const
    {
        IMAGE_ITERATOR ret(*this);
        ret.y -= d;
        return ColumnIterator(ret);
    }
        /** Calculate distance.
        */
    int operator-(ColumnIterator const & c) const
    {
        return this->y - c.y;
    }

        /** Equality.
        */
    bool operator==(ColumnIterator const & c) const
    {
        return IMAGE_ITERATOR::operator==(c);
    }

        /** Inequality.
        */
    bool operator!=(ColumnIterator const & c) const
    {
        return IMAGE_ITERATOR::operator!=(c);
    }

        /** Smaller than.
        */
    bool operator<(ColumnIterator const & c) const
    {
        return this->y < c.y;
    }

        /** Access current pixel.
        */
    reference operator*() const
    {
        return IMAGE_ITERATOR::operator*();
    }

        /** Access pixel at distance d.
        */
    index_reference operator[](int d) const
    {
        return IMAGE_ITERATOR::operator()(0, d);
    }

        /** Call member function of current pixel.
        */
    pointer operator->() const
    {
        return IMAGE_ITERATOR::operator->();
    }

        /** Get a reference to the adapted iterator
        */
    Adaptee & adaptee() const { return (Adaptee &)*this; }

    //@}
};

/********************************************************/
/*                                                      */
/*                      RowIterator                     */
/*                                                      */
/********************************************************/

/** \brief Iterator adapter to linearly access row.

    This iterator may be initialized from a standard ImageIterator,
     a MultibandImageIterator and so on.
    It gives you STL-compatible (random access iterator) access to
    one row of the image. If the underlying iterator is a const iterator,
    the row iterator will also be const (i.e. doesn't allow to change
    the values it points to).
    The iterator gets associated with the accessor of the base iterator.

    Note that image iterators usually have a member <TT>rowIterator()</TT>
    which returns a row iterator optimized for that particular image class.
    RowIterator is only necessary if this 'native' row iterator
    is not usable in a particular situation or is not provided.

    <b>\#include</b> \<vigra/imageiteratoradapter.hxx\>

    Namespace: vigra

*/
template <class IMAGE_ITERATOR>
class RowIterator : private IMAGE_ITERATOR
{
  public:
        /** the iterator's value type
        */
    typedef typename IMAGE_ITERATOR::value_type value_type;

        /** the iterator's value type
        */
    typedef typename IMAGE_ITERATOR::value_type PixelType;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef typename IMAGE_ITERATOR::reference              reference;

        /** the iterator's index reference type (return type of <TT>iter[n]</TT>)
        */
    typedef typename IMAGE_ITERATOR::index_reference        index_reference;

        /** the iterator's pointer type (return type of <TT>iter.operator->()</TT>)
        */
    typedef typename IMAGE_ITERATOR::pointer                pointer;

        /** the iterator's difference type (argument type of <TT>iter[diff]</TT>)
        */
    typedef typename IMAGE_ITERATOR::difference_type::MoveY difference_type;

        /** the iterator tag (random access iterator)
        */
    typedef std::random_access_iterator_tag                 iterator_category;

        /** the type of the adapted iterator
        */
    typedef IMAGE_ITERATOR Adaptee;

        /** Construct from an the image iterator to be adapted.
        */
    RowIterator(IMAGE_ITERATOR  const & i)
    : IMAGE_ITERATOR(i)
    {}

        /** Assignment.
        */
    RowIterator & operator=(RowIterator  const & i)
    {
        IMAGE_ITERATOR::operator=(i);

        return *this;
    }

        /** Assign a new base iterator.
        */
    RowIterator & operator=(IMAGE_ITERATOR  const & i)
    {
        IMAGE_ITERATOR::operator=(i);

        return *this;
    }

    /** @name Navigation */
    //@{
        ///
    RowIterator &  operator++()
    {
        ++(this->x);
        return *this;
    }
        ///
    RowIterator  operator++(int)
    {
        RowIterator ret(*this);
        (this->x)++;
        return ret;
    }

        ///
    RowIterator &  operator--()
    {
        --(this->x);
        return *this;
    }

        ///
    RowIterator  operator--(int)
    {
        RowIterator ret(*this);
        (this->x)--;
        return ret;
    }

        ///
    RowIterator &  operator+=(int d)
    {
        this->x += d;
        return *this;
    }

        ///
    RowIterator &  operator-=(int d)
    {
        this->x -= d;
        return *this;
    }
    //@}

    /** @name Methods */
    //@{
        /** Construct iterator at a distance.
        */
    RowIterator operator+(int d) const
    {
        IMAGE_ITERATOR ret(*this);
        ret.x += d;
        return RowIterator(ret);
    }
        /** Construct iterator at a distance.
        */
    RowIterator operator-(int d) const
    {
        IMAGE_ITERATOR ret(*this);
        ret.x -= d;
        return RowIterator(ret);
    }
        /** Calculate distance.
        */
    int operator-(RowIterator const & c) const
    {
        return this->x - c.x;
    }

        /** Equality.
        */
    bool operator==(RowIterator const & c) const
    {
        return IMAGE_ITERATOR::operator==(c);
    }

        /** Inequality.
        */
    bool operator!=(RowIterator const & c) const
    {
        return IMAGE_ITERATOR::operator!=(c);
    }

        /** Smaller than.
        */
    bool operator<(RowIterator const & c) const
    {
        return this->x < c.x;
    }

        /** Access current pixel.
        */
    reference operator*() const
    {
        return IMAGE_ITERATOR::operator*();
    }

        /** Access pixel at distance d.
        */
    index_reference operator[](int d) const
    {
        return IMAGE_ITERATOR::operator()(d, 0);
    }

        /** Call member function of current pixel.
        */
    pointer operator->() const
    {
        return IMAGE_ITERATOR::operator->();
    }

        /** Get a reference to the adapted iterator
        */
    Adaptee & adaptee() const { return (Adaptee &)*this; }

    //@}
};

/********************************************************/
/*                                                      */
/*                     LineIterator                     */
/*                                                      */
/********************************************************/

/** \brief Iterator adapter to iterate along an arbitrary line on the image.

    This iterator may be initialized from a standard ImageIterator,
     a MultibandImageIterator and so on.
    It gives you STL-compatible (forward iterator) access to
    an arbitraty line on the image.
    The iterator gets associated with the accessor of the base iterator.

    <b>\#include</b> \<vigra/imageiteratoradapter.hxx\>

    Namespace: vigra

*/
template <class IMAGE_ITERATOR>
class LineIterator : private IMAGE_ITERATOR
{
  public:
        /** the iterator's value type
        */
    typedef typename IMAGE_ITERATOR::value_type value_type;

        /** the iterator's value type
        */
    typedef typename IMAGE_ITERATOR::value_type PixelType;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef typename IMAGE_ITERATOR::reference              reference;

        /** the iterator's pointer type (return type of <TT>iter.operator->()</TT>)
        */
    typedef typename IMAGE_ITERATOR::pointer                pointer;

        /** the iterator tag (forward iterator)
        */
    typedef std::forward_iterator_tag                       iterator_category;

        /** the type of the adapted iterator
        */
    typedef IMAGE_ITERATOR Adaptee;

        /** Construct from an the image iterator to be adapted.
        */
    LineIterator(IMAGE_ITERATOR  const & start,
                 IMAGE_ITERATOR  const & end)
    : IMAGE_ITERATOR(start), x_(0.0), y_(0.0)
    {
        int dx = end.x - start.x;
        int dy = end.y - start.y;
        int adx = (dx < 0) ? -dx : dx;
        int ady = (dy < 0) ? -dy : dy;
        int dd = (adx > ady) ? adx : ady;
        if(dd == 0) dd = 1;

        dx_ = (double)dx / dd;
        dy_ = (double)dy / dd;
        if(adx > ady) y_ += dy_ / 2.0;
        else          x_ += dx_ / 2.0;
    }

    /** @name Navigation */
    //@{
        ///
    LineIterator &  operator++()
    {
        x_ += dx_;
        if(x_ >= 1.0) {
            x_ -= 1.0;
            ++(this->x);
        }
        else if(x_ <= -1.0) {
            x_ += 1.0;
            --(this->x);
        }
        y_ += dy_;
        if(y_ >= 1.0) {
            y_ -= 1.0;
            ++(this->y);
        }
        else if(y_ <= -1.0) {
            y_ += 1.0;
            --(this->y);
        }
        return *this;
    }
        ///
    LineIterator  operator++(int)
    {
        LineIterator ret(*this);
        operator++();
        return ret;
    }

    //@}

    /** @name Methods */
    //@{
        /** Equality.
       */
    bool operator==(LineIterator const & c) const
    {
        return IMAGE_ITERATOR::operator==(c);
    }

        /** Inequality.
       */
    bool operator!=(LineIterator const & c) const
    {
        return IMAGE_ITERATOR::operator!=(c);
    }

        /** Access current pixel.
       */
    reference operator*() const
    {
        return IMAGE_ITERATOR::operator*();
    }

        /** Call member function for current pixel.
       */
    pointer operator->() const
    {
        return IMAGE_ITERATOR::operator->();
    }

        /** Get a reference to the adapted iterator
       */
    Adaptee & adaptee() const { return (Adaptee &)*this; }

    //@}

  private:

    double x_, y_, dx_, dy_;
};

//@}

} // namespace vigra

#endif // VIGRA_IMAGEITERATORADAPTER_HXX
