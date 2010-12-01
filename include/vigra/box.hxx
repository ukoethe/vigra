#ifndef VIGRA_BOX_HXX
#define VIGRA_BOX_HXX

#include "metaprogramming.hxx"
#include "numerictraits.hxx"
#include "tinyvector.hxx"

namespace vigra {

namespace detail {

// RangePolicy used for floating point coordinate types
template<class VALUETYPE>
struct EndInsidePolicy
{
    static inline bool isEmptyRange(VALUETYPE b, VALUETYPE e)
    {
        return e < b; // <=
    }

    static inline VALUETYPE pointEnd(VALUETYPE p)
    {
        return p; // +1
    }
};

// RangePolicy used for integer coordinate types
template<class VALUETYPE>
struct EndOutsidePolicy
{
    static inline bool isEmptyRange(VALUETYPE b, VALUETYPE e)
    {
        return e <= b;
    }

    static inline VALUETYPE pointEnd(VALUETYPE p)
    {
        return p+1;
    }
};

} // namespace vigra::detail

    /**
     * Represent an n-dimensional box as a (begin, end) pair.
     * Depending on the value type, end() is considered to be
     * outside the box (as in the STL, for integer types), or
     * inside (for floating point types).  size() will always be
     * end() - begin().
     */
template<class VALUETYPE, unsigned int DIMENSION>
class Box
{
  public:
        /** STL-compatible definition of coordinate valuetype
         */
    typedef VALUETYPE value_type;

        /** Promoted coordinate valuetype, used for volume()
         */
    typedef typename NumericTraits<VALUETYPE>::Promote VolumeType;

        /** Vector type used for begin() and end()
         */
    typedef TinyVector<VALUETYPE, DIMENSION> Vector;

    enum { Dimension = DIMENSION };

  protected:
    Vector begin_, end_;

        /** Range policy (EndInsidePolicy/EndOutsidePolicy, depending on valuetype)
         */
    typedef typename If<typename NumericTraits<VALUETYPE>::isIntegral,
                        detail::EndOutsidePolicy<VALUETYPE>,
                        detail::EndInsidePolicy<VALUETYPE> >::type RangePolicy;

  public:
        /** Construct an empty box (isEmpty() will return true).
         * (Internally, this will initialize all dimensions with the
         * empty range [1..0].)
         */
    Box()
    : begin_(NumericTraits<Vector>::one())
    {}

        /** Construct a box representing the given range.  Depending
         * on the value type, end() is considered to be outside the
         * box (as in the STL, for integer types), or inside (for
         * floating point types).
         */
    Box(Vector const &begin, Vector const &end)
    : begin_(begin), end_(end)
    {}

        /** Construct a box of given size at the origin (i.e. end() ==
         * size()).
         */
    explicit Box(Vector const &size)
    : end_(size)
    {}

        /** Get begin vector (i.e. smallest coordinates for each
         * dimension).  This is the first point (scan-order wise)
         * which is considered to be "in" the box.
         */
    Vector const & begin() const
    {
        return begin_;
    }

        /** Access begin vector (i.e. smallest coordinates for each
         * dimension).  This is the first point (scan-order wise)
         * which is considered to be "in" the box.
         */
    Vector & begin()
    {
        return begin_;
    }

        /** Get end vector (i.e. coordinates higher than begin() in
         * each dimension for non-empty boxes).  This is begin() +
         * size(), and depending on the valuetype (float/int), this is
         * the last point within or the first point outside the box,
         * respectively.
         */
    Vector const & end() const
    {
        return end_;
    }

        /** Access end vector (i.e. coordinates higher than begin() in
         * each dimension for non-empty boxes).  This is begin() +
         * size(), and depending on the valuetype (float/int), this is
         * the last point within or the first point outside the box,
         * respectively.
         */
    Vector & end()
    {
        return end_;
    }

        /** Change begin() without changing end(), changing size()
         * accordingly.
         */
    void setBegin(Vector const &begin)
    {
        begin_ = begin;
    }

        /** Change end() without changing begin(), which will change
         * the size() most probably.
         */
    void setEnd(Vector const &end)
    {
        end_ = end;
    }

        /** Move the whole box so that the given point will be
         * begin() afterwards.
         */
    void moveTo(Vector const &newBegin)
    {
        end_ += newBegin - begin_;
        begin_ = newBegin;
    }

        /** Move the whole box by the given offset.
         * (Equivalent to operator+=)
         */
    void moveBy(Vector const &offset)
    {
        begin_ += offset;
        end_ += offset;
    }

        /** Determine and return the area of this box. That is,
         * if this rect isEmpty(), returns zero, otherwise returns the
         * product of the extents in each dimension.
         */
    VolumeType volume() const
    {
        if(isEmpty())
            return 0;

        VolumeType result(end_[0] - begin_[0]);
        for(unsigned int i = 1; i < DIMENSION; ++i)
            result *= end_[i] - begin_[i];
        return result;
    }

        /** Determine and return the size of this box. The size
         * might be zero or even negative in one or more dimensions,
         * and if so, isEmpty() will return true.
         */
    Vector size() const
    {
        return end_ - begin_;
    }

        /** Resize this box to the given extents. This will
         * change end() only.
         */
    void setSize(Vector const &size)
    {
        end_ = begin_ + size;
    }

        /** Increase the size of the box by the given
         * offset. This will move end() only. (If any of offset's
         * components is negative, the box will get smaller
         * accordingly.)
         */
    void addSize(Vector const &offset)
    {
        end_ += offset;
    }

        /** Adds a border of the given width around the box. That
         * means, begin()'s components are moved by -borderWidth
         * and end()'s by borderWidth. (If borderWidth is
         * negative, the box will get smaller accordingly.)
         */
    void addBorder(VALUETYPE borderWidth)
    {
        for(unsigned int i = 0; i < DIMENSION; ++i)
        {
            begin_[i] -= borderWidth;
            end_[i]   += borderWidth;
        }
    }

        /// equality check
    bool operator==(Box const &r) const
    {
        return (begin_ == r.begin_) && (end_ == r.end_);
    }

        /// inequality check
    bool operator!=(Box const &r) const
    {
        return (begin_ != r.begin_) || (end_ != r.end_);
    }

        /** Return whether this box is considered empty. It is
         * non-empty if all end() coordinates are greater than (or
         * equal, for floating point valuetypes) the corresponding
         * begin() coordinates. Uniting an empty box with something
         * will return the bounding box of the 'something', and
         * intersecting any box with an empty box will again yield an
         * empty box.
         */
    bool isEmpty() const
    {
        for(unsigned int i = 0; i < DIMENSION; ++i)
            if(RangePolicy::isEmptyRange(begin_[i], end_[i]))
                return true;
        return false;
    }

        /** Return whether this box contains the given point.
         * That is, if the point lies within the range [begin, end] in
         * each dimension (excluding end() itself for integer valuetypes).
         */
    bool contains(Vector const &p) const
    {
        for(unsigned int i = 0; i < DIMENSION; ++i)
            if((p[i] < begin_[i]) ||
               RangePolicy::isEmptyRange(p[i], end_[i]))
                return false;
        return true;
    }

        /** Return whether this box contains the given
         * one. <tt>r1.contains(r2)</tt> returns the same as
         * <tt>r1 == (r1|r2)</tt> (but is of course more
         * efficient). That also means, a box (even an empty one!)
         * contains() any empty box.
         */
    bool contains(Box const &r) const
    {
        if(r.isEmpty())
            return true;
        if(!contains(r.begin_))
            return false;
        for(unsigned int i = 0; i < DIMENSION; ++i)
            if(r.end_[i] > end_[i])
                return false;
        return true;
    }

        /** Return whether this box overlaps with the given
         * one. <tt>r1.intersects(r2)</tt> returns the same as
         * <tt>!(r1&r2).isEmpty()</tt> (but is of course much more
         * efficient).
         */
    bool intersects(Box const &r) const
    {
        if(r.isEmpty() || isEmpty())
            return false;
        for(unsigned int i = 0; i < DIMENSION; ++i)
            if(RangePolicy::isEmptyRange(r.begin_[i], end_[i]) ||
               RangePolicy::isEmptyRange(begin_[i], r.end_[i]))
                return false;
        return true;
    }

        /** Modifies this box by including the given point.
         * The result will be the bounding box of the box and the
         * point.  If isEmpty() returns true on the original box, the
         * union will be a box containing only the given point.
         */
    Box &operator|=(Vector const &p)
    {
        if(isEmpty())
        {
            begin_ = p;
            for(unsigned int i = 0; i < DIMENSION; ++i)
                end_[i] = RangePolicy::pointEnd(p[i]);
        }
        else
        {
            for(unsigned int i = 0; i < DIMENSION; ++i)
            {
                if(p[i] < begin_[i])
                    begin_[i] = p[i];
                if(RangePolicy::isEmptyRange(p[i], end_[i]))
                    end_[i] = RangePolicy::pointEnd(p[i]);
            }
        }
        return *this;
    }

        /** Returns the union of this box and the given point.
         * The result will be the bounding box of the box and the
         * point.  If isEmpty() returns true on the original box, the
         * union will be a box containing only the given point.
         */
    Box operator|(Vector const &p) const
    {
        Box result(*this);
        result |= p;
        return result;
    }

        /** Modifies this box by uniting it with the given one.
         * The result will be the bounding box of both boxs. If one of
         * the boxes isEmpty(), the union will be the other one.
         */
    Box &operator|=(Box const &r)
    {
        if(r.isEmpty())
            return *this;
        if(isEmpty())
            return operator=(r);

        for(unsigned int i = 0; i < DIMENSION; ++i)
        {
            if(r.begin_[i] < begin_[i])
                begin_[i] = r.begin_[i];
            if(end_[i] < r.end_[i])
                end_[i] = r.end_[i];
        }
        return *this;
    }

        /** Returns the union of this box and the given one.
         * The result will be the bounding box of both boxs. If one of
         * the boxes isEmpty(), the union will be the other one.
         */
    Box operator|(Box const &r) const
    {
        Box result(*this);
        result |= r;
        return result;
    }

        /** Modifies this box by intersecting it with the given one.
         * The result will be the maximal box contained in both
         * original ones. Intersecting with an empty box will yield
         * again an empty box.
         */
    Box &operator&=(Box const &r)
    {
        if(isEmpty())
            return *this;
        if(r.isEmpty())
            return operator=(r);

        for(unsigned int i = 0; i < DIMENSION; ++i)
        {
            if(begin_[i] < r.begin_[i])
                begin_[i] = r.begin_[i];
            if(r.end_[i] < end_[i])
                end_[i] = r.end_[i];
        }
        return *this;
    }

        /** Intersects this box with the given one.
         * The result will be the maximal box contained in both
         * original ones.  Intersecting with an empty box will yield
         * again an empty box.
         */
    Box operator&(Box const &r) const
    {
        Box result(*this);
        result &= r;
        return result;
    }

        /**
         * Scale box by scalar multiply-assignment.  The same scalar
         * multiply-assignment operation will be performed on both
         * begin() and end().
         */
    Box &operator*=(double scale)
    {
        begin_ *= scale;
        end_   *= scale;
        return *this;
    }

        /**
         * Return box scaled by given factor.  The same scalar
         * multiplication will be performed on both begin() and end().
         */
    Box operator*(double scale)
    {
        Box result(*this);
        result *= scale;
        return result;
    }

        /**
         * Scale box by scalar divide-assignment.  The same scalar
         * divide-assignment operation will be performed on both
         * begin() and end().
         */
    Box &operator/=(double scale)
    {
        begin_ /= scale;
        end_   /= scale;
        return *this;
    }

        /**
         * Return box scaled by inverse of given factor.  The same scalar
         * division will be performed on both begin() and end().
         */
    Box operator/(double scale)
    {
        Box result(*this);
        result /= scale;
        return result;
    }

        /**
         * Translate box by vector addition-assignment.  The same vector
         * addition-assignment operation will be performed on both
         * begin() and end().
         */
    Box &operator+=(const Vector &offset)
    {
        begin_ += offset;
        end_   += offset;
        return *this;
    }

        /**
         * Translate box by vector addition.  The same vector addition
         * operation will be performed on both begin() and end().
         */
    Box operator+(const Vector &offset)
    {
        Box result(*this);
        result += offset;
        return result;
    }

        /**
         * Translate box by vector subtract-assignment.  The same vector
         * subtract-assignment operation will be performed on both
         * begin() and end().
         */
    Box &operator-=(const Vector &offset)
    {
        begin_ -= offset;
        end_   -= offset;
        return *this;
    }

        /**
         * Translate box by vector subtract.  The same vector subtract
         * operation will be performed on both begin() and end().
         */
    Box operator-(const Vector &offset)
    {
        Box result(*this);
        result -= offset;
        return result;
    }
};

} // namespace vigra

#endif // VIGRA_BOX_HXX
