/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
 
 
#ifndef VIGRA_IMAGEITERATORADAPTER_HXX
#define VIGRA_IMAGEITERATORADAPTER_HXX

namespace vigra {

/** \addtogroup ImageIteratorAdapters Image Iterator Adapters
    
     Iterate over rows, columns, and other image subsets
*/
//@{

/********************************************************/
/*                                                      */
/*                      ColumnIterator                  */
/*                                                      */
/********************************************************/

/** \brief Iterator adapter to linearly access colums.

    This iterator may be initialized from a standard ImageIterator,
     a MultibandImageIterator and so on. 
    It gives you STL-combatibel (random access iterator) access to 
    one column of the image.
    The iterator gets associated with the accessor of the base iterator.
    
    <b>\#include</b> "<a href="imageiteratoradapter_8hxx-source.html">vigra/imageiteratoradapter.hxx</a>"
    
    Namespace: vigra
    
*/
template <class IMAGE_ITERATOR>
class ColumnIterator : private IMAGE_ITERATOR
{
  public:
        /** the iterator's PixelType
        */
    typedef typename IMAGE_ITERATOR::PixelType PixelType;
    
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
    PixelType & operator*()
    {
        return IMAGE_ITERATOR::operator*(); 
    }
    
        /** Read current pixel.
        */
    PixelType const & operator*() const
    {
        return IMAGE_ITERATOR::operator*(); 
    }
    
        /** Access pixel at distance d.
        */
    PixelType & operator[](int d)
    {
        return IMAGE_ITERATOR::operator()(0, d);
    }

        /** Read pixel at distance d.
        */
    PixelType const & operator[](int d) const
    {
        return IMAGE_ITERATOR::operator()(0, d);
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
    It gives you STL-combatibel (random access iterator) access to 
    one row of the image.
    The iterator gets associated with the accessor of the base iterator.
    
    <b>\#include</b> "<a href="imageiteratoradapter_8hxx-source.html">vigra/imageiteratoradapter.hxx</a>"
    
    Namespace: vigra
    
*/
template <class IMAGE_ITERATOR>
class RowIterator : private IMAGE_ITERATOR
{
  public:
        /** the iterator's PixelType
       */
    typedef typename IMAGE_ITERATOR::PixelType PixelType;
    
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
    PixelType & operator*()
    {
        return IMAGE_ITERATOR::operator*(); 
    }
    
        /** Read current pixel.
        */
    PixelType const & operator*() const
    {
        return IMAGE_ITERATOR::operator*(); 
    }
    
        /** Access pixel at distance d.
        */
    PixelType & operator[](int d)
    {
        return IMAGE_ITERATOR::operator()(d, 0);
    }

        /** Read pixel at distance d.
        */
    PixelType const & operator[](int d) const
    {
        return IMAGE_ITERATOR::operator()(d, 0);
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
    It gives you STL-combatibel (forward iterator) access to 
    an arbitraty line on the image.
    The iterator gets associated with the accessor of the base iterator.
    
    <b>\#include</b> "<a href="imageiteratoradapter_8hxx-source.html">vigra/imageiteratoradapter.hxx</a>"
    
    Namespace: vigra
    
*/
template <class IMAGE_ITERATOR>
class LineIterator : private IMAGE_ITERATOR
{
  public:
        /** the iterator's PixelType
        */
    typedef typename IMAGE_ITERATOR::PixelType PixelType;
    
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
    PixelType & operator*()
    {
        return IMAGE_ITERATOR::operator*(); 
    }
    
        /** Read current pixel.
        */
    PixelType const & operator*() const
    {
        return IMAGE_ITERATOR::operator*(); 
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
