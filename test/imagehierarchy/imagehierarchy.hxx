#ifndef VIGRA_IMAGEHIERARCHY_HXX
#define VIGRA_IMAGEHIERARCHY_HXX

#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/imageiterator.hxx"
#include "vigra/accessor.hxx"
#include "boost/smart_ptr.hpp"
#include "boost/static_assert.hpp"

namespace vigra {

typedef float GrayValue;

class ConstVectorProxy
{
  public:
    typedef GrayValue value_type;
    typedef GrayValue const * iterator;
    typedef GrayValue const * const_iterator;
    
    ConstVectorProxy()
    : data_(0), size_(0)
    {}
    
    ConstVectorProxy(GrayValue const * data, int size)
    : data_(const_cast<value_type *>(data)), size_(size)
    {}
    
    ConstVectorProxy(value_type const & f)
    : data_(const_cast<value_type *>(&f)), size_(1)
    {}
    
    template <int N>
    ConstVectorProxy(TinyVector<value_type, N> const & v)
    : data_(const_cast<TinyVector<value_type, N> &>(v).begin()), size_(N)
    {}
    
    void reset(ConstVectorProxy const & v)
    {
        data_ = v.data_;
        size_ = v.size_;
    }
    
    void reset(value_type const & f)
    {
        data_ = const_cast<value_type *>(&f);
        size_ = 1;
    }
    
    template <int N>
    void reset(TinyVector<value_type, N> const & v)
    {
        data_ = const_cast<TinyVector<value_type, N> &>(v).begin();
        size_ = N;
    }
    
    void resize(int new_size)
        { size_ = new_size; }
    
    iterator begin()
        { return data_; }
    
    iterator end()
        { return data_ + size_; }
    
    const_iterator begin() const
        { return data_; }
    
    const_iterator end() const
        { return data_ + size_; }
    
    value_type const & operator[](int i) const
        { return data_[i]; }
        
    int size() const
        { return size_; }
        
    operator value_type() const
    {
        vigra_precondition(size_ == 1, 
            "ConstVectorProxy::operator value_type(): vector must have size 1.");
        return *data_;
    }
    
    template <int N>
    operator TinyVector<value_type, N>() const
    {
        vigra_precondition(size_ == N, 
            "ConstVectorProxy::operator TinyVector(): size mismatch.");
        TinyVector<value_type, N> res;
        res.init(begin(), end());
        return res;
    }
    
    operator RGBValue<value_type>() const
    {
        vigra_precondition(size_ == 3, 
            "ConstVectorProxy::operator RGBValue(): size mismatch.");
        return RGBValue<value_type>(begin(), end());
    }
    
    bool operator==(ConstVectorProxy const & o) const
    {
        if(size() != o.size())
            return false;
        for(int i=0; i<size(); ++i)
            if ((*this)[i] != o[i])
                return false;
        return true;
    }
     
    bool operator!=(ConstVectorProxy const & o) const
    {
        return !(*this == o);
    }
    
  protected:
  
    ConstVectorProxy & operator=(value_type const & f);
    
    template <int N>
    ConstVectorProxy & operator=(TinyVector<value_type, N> const & v);
    
    value_type * data_;
    int size_;
};

class VectorProxy
: public ConstVectorProxy
{
  public:
    
    typedef GrayValue * iterator;
    typedef GrayValue const * const_iterator;
    
    VectorProxy()
    {}
    
    VectorProxy(GrayValue const * data, int size)
    : ConstVectorProxy(data, size)
    {}
    
    VectorProxy(value_type const & f)
    : ConstVectorProxy(f)
    {}
    
    template <int N>
    VectorProxy(TinyVector<value_type, N> const & v)
    : ConstVectorProxy(v)
    {}
    
    VectorProxy & operator=(ConstVectorProxy const & v)
    {
        vigra_precondition(size_ == v.size(),
           "VectorProxy::operator=(): size mismatch.");
        for(int i=0; i<size_; ++i)
            data_[i] = v[i];
        return *this;
    }
    
    VectorProxy & operator=(value_type const & f)
    {
        vigra_precondition(size_ == 1,
           "VectorProxy::operator=(): size mismatch.");
        *data_ = f;
        return *this;
    }
    
    template <int N>
    VectorProxy & operator=(TinyVector<value_type, N> const & v)
    {
        vigra_precondition(size_ == N,
           "VectorProxy::operator=(): size mismatch.");
        for(int i=0; i<N; ++i)
            data_[i] = v[i];
        return *this;
    }
    
    iterator begin()
        { return data_; }
    
    iterator end()
        { return data_ + size_; }
    
    const_iterator begin() const
        { return data_; }
    
    const_iterator end() const
        { return data_ + size_; }
    
    value_type & operator[](int i)
        { return data_[i]; }
        
    value_type const & operator[](int i) const
        { return data_[i]; }
        
};

template <class IMAGEITERATOR>
class MultibandRowColumnIteratorPolicy
{
  public:
    typedef IMAGEITERATOR                            ImageIterator;
    typedef typename ImageIterator::value_type       value_type;
    typedef int                                      difference_type;
    typedef typename ImageIterator::reference        reference;
    typedef typename ImageIterator::index_reference  index_reference;
    typedef typename ImageIterator::pointer          pointer;
    typedef std::random_access_iterator_tag iterator_category;
    
    
    struct BaseType
    {
        explicit BaseType(pointer c = 0, difference_type o = 0, difference_type s = 0)
        : current_(c), offset_(o), size_(s)
        {}
        
        pointer current_;
        difference_type offset_, size_;
    };
    
    static void initialize(BaseType & d) {}
    
    static reference dereference(BaseType const & d)
        { return reference(d.current_, d.size_); }
    
    static index_reference dereference(BaseType const & d, difference_type n)
    { 
        return index_reference(d.current_+n*d.offset_, d.size_);
    }
    
    static bool equal(BaseType const & d1, BaseType const & d2)
        { return d1.current_ == d2.current_; }
    
    static bool less(BaseType const & d1, BaseType const & d2)
        { return d1.current_ < d2.current_; }
    
    static difference_type difference(BaseType const & d1, BaseType const & d2)
        { return (d1.current_ - d2.current_) / d1.offset_; }
    
    static void increment(BaseType & d)
        { d.current_ += d.offset_; }
    
    static void decrement(BaseType & d)
        { d.current_ -= d.offset_; }
    
    static void advance(BaseType & d, difference_type n)
        { d.current_ += d.offset_*n; }
};

template <class PIXELTYPE, class ITERATOR>
class VariableBandsIteratorBase
{
  protected:
    
    PIXELTYPE * data_;
    int width_, bands_, size_;
    
    VariableBandsIteratorBase(PIXELTYPE * data, int width, int bands, int size)
    : data_(data), width_(width), bands_(bands), size_(size), x(0), y(0)
    {}
    
    VariableBandsIteratorBase()
    : data_(0), width_(0), bands_(0), size_(0), x(0), y(0)
    {}
    
    PIXELTYPE * get()
        { return data_ + ((y * width_) + x)*bands_; }
    
    PIXELTYPE const * get() const
        { return data_ + ((y * width_) + x)*bands_; }
    
    PIXELTYPE * get(int const & dx, int const & dy)
        { return data_ + ((dy + y) * width_ + dx + x)*bands_; }
    
    PIXELTYPE const * get(int const & dx, int const & dy) const
        { return data_ + ((dy + y) * width_ + dx + x)*bands_; }
    
  public:
  
    typedef int MoveX;
    typedef int MoveY;
    
    int x, y;
    
    bool operator==(VariableBandsIteratorBase const & rhs) const
    {
        return (x == rhs.x) && (y == rhs.y);
    }

    bool operator!=(VariableBandsIteratorBase const & rhs) const
    {
        return (x != rhs.x) || (y != rhs.y);
    }
    
    Diff2D operator-(VariableBandsIteratorBase const & rhs) const
    {
        return Diff2D(x - rhs.x, y - rhs.y);
    }
    
    ITERATOR & operator+=(Diff2D const & s)
    {
        x += s.x;
        y += s.y;
        return (ITERATOR &)*this;
    }

    ITERATOR & operator-=(Diff2D const & s)
    {
        x -= s.x;
        y -= s.y;
        return (ITERATOR &)*this;
    }
};

struct VariableBandsIterator
: public VariableBandsIteratorBase<GrayValue, VariableBandsIterator>
{
  public:
  
    typedef VectorProxy          value_type;
    typedef VectorProxy          PixelType;
    typedef VectorProxy          reference;
    typedef VectorProxy          index_reference;
    typedef GrayValue *          pointer;
    typedef Diff2D               difference_type;
    typedef image_traverser_tag  iterator_category;
    typedef IteratorAdaptor<MultibandRowColumnIteratorPolicy<VariableBandsIterator> > 
        row_iterator;
    typedef IteratorAdaptor<MultibandRowColumnIteratorPolicy<VariableBandsIterator> > 
        column_iterator;
    
  
    VariableBandsIterator(GrayValue * data, int width, int bands, int size)
    : VariableBandsIteratorBase<GrayValue, VariableBandsIterator>(data, width, bands, size)
    {}
    
    VariableBandsIterator()
    {}
    
    VariableBandsIterator operator+(Diff2D const & s)
    {
        VariableBandsIterator res(*this);
        res += s;
        return res;
    }

    VariableBandsIterator operator-(Diff2D const & s)
    {
        VariableBandsIterator res(*this);
        res -= s;
        return res;
    }
    
    pointer operator->() const
    {
        return const_cast<pointer>(get());
    }
    
    reference operator*() const
    {
        return reference(get(), size_);
    }
    
    index_reference operator[](Diff2D const & dist) const
    {
        return index_reference(get(dist.x, dist.y), size_);
    }
    
    index_reference operator()(int const & dx, int const & dy) const
    {
        return index_reference(get(dx, dy), size_);
    }
    
    row_iterator rowIterator() const
    {
        return row_iterator(
            row_iterator::BaseType(const_cast<pointer>(get()), bands_, size_));
    }
    
    column_iterator columnIterator() const
    {
        return column_iterator(
            column_iterator::BaseType(const_cast<pointer>(get()), width_*bands_, size_));
    }
};

struct ConstVariableBandsIterator
: public VariableBandsIteratorBase<GrayValue, ConstVariableBandsIterator>
{
  public:
  
    typedef ConstVectorProxy          value_type;
    typedef ConstVectorProxy          PixelType;
    typedef ConstVectorProxy          reference;
    typedef ConstVectorProxy          index_reference;
    typedef GrayValue const *         pointer;
    typedef Diff2D               difference_type;
    typedef image_traverser_tag  iterator_category;
    typedef IteratorAdaptor<MultibandRowColumnIteratorPolicy<ConstVariableBandsIterator> > 
        row_iterator;
    typedef IteratorAdaptor<MultibandRowColumnIteratorPolicy<ConstVariableBandsIterator> > 
        column_iterator;
  
    ConstVariableBandsIterator(GrayValue const * data, int width, int bands, int size)
    : VariableBandsIteratorBase<GrayValue, ConstVariableBandsIterator>(
        const_cast<GrayValue *>(data), width, bands, size)
    {}
    
    ConstVariableBandsIterator()
    {}
    
    ConstVariableBandsIterator operator+(Diff2D const & s)
    {
        ConstVariableBandsIterator res(*this);
        res += s;
        return res;
    }

    ConstVariableBandsIterator operator-(Diff2D const & s)
    {
        ConstVariableBandsIterator res(*this);
        res -= s;
        return res;
    }
    
    pointer operator->() const
    {
        return get();
    }
    
    reference operator*() const
    {
        return reference(get(), size_);
    }
    
    index_reference operator[](Diff2D const & dist) const
    {
        return index_reference(get(dist.x, dist.y), size_);
    }
    
    index_reference operator()(int const & dx, int const & dy) const
    {
        return index_reference(get(dx, dy), size_);
    }
    
    row_iterator rowIterator() const
    {
        return row_iterator(row_iterator::BaseType(get(), bands_, size_));
    }
    
    column_iterator columnIterator() const
    {
        return column_iterator(column_iterator::BaseType(get(), width_*bands_, size_));
    }
};


template <class T>
struct BandsForPixelType;

template <>
struct BandsForPixelType<GrayValue>
    { static const int bands = 1; };

template <>
struct BandsForPixelType<RGBValue<GrayValue> >
    { static const int bands = 3; };

template <>
struct BandsForPixelType<TinyVector<GrayValue, 2> >
    { static const int bands = 2; };

template <>
struct BandsForPixelType<TinyVector<GrayValue, 3> >
    { static const int bands = 3; };

template <>
struct BandsForPixelType<TinyVector<GrayValue, 4> >
    { static const int bands = 4; };
    
class VariableBandsImage
{
    
  public:
    typedef VectorProxy       value_type;
    typedef VectorProxy       PixelType;
    typedef ConstVectorProxy  ConstPixelType;
    typedef VectorProxy       reference;
    typedef ConstVectorProxy  const_reference;
    typedef GrayValue *       pointer;
    typedef GrayValue const * const_pointer;
//    typedef ... iterator;
//    typedef ... ScanOrderIterator;
//    typedef ... const_iterator;
//    typedef ... ConstScanOrderIterator;
    typedef VariableBandsIterator traverser;
    typedef VariableBandsIterator Iterator;
    typedef ConstVariableBandsIterator const_traverser;
    typedef ConstVariableBandsIterator ConstIterator; 

    typedef Diff2D difference_type;
    typedef Diff2D size_type;
    typedef StandardValueAccessor<PixelType> Accessor;
    typedef StandardConstValueAccessor<ConstPixelType> ConstAccessor;
    typedef VariableBandsImage CloneType;
    
    virtual ~VariableBandsImage() {}
    
    virtual CloneType * clone() const = 0;
    virtual VariableBandsImage * shallowCopy() const = 0;
    
    int width() const
        { return roiLR_.x - roiUL_.x; }
    int height() const
        { return roiLR_.y - roiUL_.y; }
    Diff2D size() const 
        { return roiLR_ - roiUL_; }
    Diff2D actualSize() const 
        { return actualLR_; }
    int actualWidth() const
        { return actualSize().x; }
    int actualHeight() const 
        { return actualSize().y; }
    int actualBands() const
        { return actualBands_; }
    int bands() const
        { return bands_; }
    Diff2D roiUpperLeft() const 
        { return roiUL_; }
    
    
    virtual void setROI(Diff2D const & ul, Diff2D const & new_size)
    {
        Diff2D new_ul = roiUL_ + ul;
        Diff2D new_lr;
        new_lr.x = new_size.x > 0 ?
                        new_ul.x + new_size.x :
                        roiLR_.x + new_size.x;
        new_lr.y = new_size.y > 0 ?
                        new_ul.y + new_size.y :
                        roiLR_.y + new_size.y;
                                
        vigra_precondition(new_ul.x >= 0 && new_ul.x < actualWidth() && 
                           new_ul.y >= 0 && new_ul.y < actualHeight(),
            "Image::setROI(): new upper left outside of image.");
        vigra_precondition(new_lr.x > 0 && new_lr.x <= actualWidth() && 
                           new_lr.y > 0 && new_lr.y <= actualHeight(),
            "Image::setROI(): new lower right outside of image.");
        vigra_precondition(new_ul.x < new_lr.x,
            "Image::setROI(): upper left isn't to the left of lower right.");
        vigra_precondition(new_ul.y < new_lr.y,
            "Image::setROI(): upper left isn't above lower right.");
            
        roiUL_ = new_ul;
        roiLR_ = new_lr;
    }
    
    PixelType operator[](Diff2D const & c) 
        { return PixelType(get(c.x, c.y), bands_); }

    PixelType operator()(int x, int y) 
        { return PixelType(get(x,y), bands_); }

    ConstPixelType operator[](Diff2D const & c) const
        { return ConstPixelType(get(c.x, c.y), bands_); }

    ConstPixelType operator()(int x, int y) const
        { return ConstPixelType(get(x,y), bands_); }

    Iterator upperLeft()
        { return Iterator(get(0,0), actualWidth(), actualBands_, bands_); }
        
    Iterator lowerRight()
        { return upperLeft() + size(); }
        
    ConstIterator upperLeft() const
        { return ConstIterator(get(0,0), actualWidth(), actualBands_, bands_); }

    ConstIterator lowerRight() const
        { return upperLeft() + size(); }
    
    Accessor accessor()
        { return Accessor(); }

    ConstAccessor accessor() const
        { return ConstAccessor(); }

  protected:
    VariableBandsImage()
    : data_(0)
    {}
    
    VariableBandsImage(VariableBandsImage const & s)
    : roiUL_(s.roiUL_), roiLR_(s.roiLR_), actualLR_(s.actualLR_),
      bands_(s.bands_), actualBands_(s.actualBands_),
      data_(s.data_)
    {}
    
    VariableBandsImage & operator=(VariableBandsImage const & s)
    {
        roiUL_ = s.roiUL_;
        roiLR_ = s.roiLR_;
        actualLR_ = s.actualLR_;
        bands_ = s.bands_;
        actualBands_ = s.actualBands_;
        data_ = s.data_;
        return *this;
    }
    
    void init(GrayValue * data, Diff2D const & actualLR, int actualBands)
    {
        init(data, actualLR, actualBands, actualBands);
    }
    
    void init(GrayValue * data, Diff2D const & actualLR, int actualBands, int bands)
    {
        roiUL_ = Diff2D(0,0);
        roiLR_ = actualLR;
        actualLR_ = actualLR;
        bands_ = bands;
        actualBands_ = actualBands;
        data_ = data;
    }
    
    static GrayValue * getDataPtr(GrayValue * p)
        { return p; }

    template <class T>
    static GrayValue * getDataPtr(T * p)
    {
        BOOST_STATIC_ASSERT( sizeof(T) ==
               BandsForPixelType<T>::bands*sizeof(GrayValue) );
        return &(*p)[0]; 
    }
        
    GrayValue * get(int const & dx, int const & dy)
        { return data_ + ((roiUL_.y + dy) * actualWidth() + roiUL_.x + dx)*actualBands_; }
    
    GrayValue const * get(int const & dx, int const & dy) const
        { return data_ + ((roiUL_.y + dy) * actualWidth() + roiUL_.x + dx)*actualBands_; }
    
    Diff2D roiUL_, roiLR_, actualLR_;
    int bands_, actualBands_;
    GrayValue * data_;
};

class SelectBandImage;

template <class IMAGE, class ACCESSOR = typename IMAGE::Accessor>
class FixedBandsImage
: public VariableBandsImage
{
  public:
    
    typedef IMAGE InnerImage;
    typedef typename ACCESSOR::value_type value_type;
    typedef value_type          PixelType;
    typedef value_type          ConstPixelType;
    typedef value_type &        reference;
    typedef value_type const &  const_reference;
    typedef value_type *        pointer;
    typedef value_type const *  const_pointer;
    typedef typename IMAGE::Iterator traverser;
    typedef typename IMAGE::Iterator Iterator;
    typedef typename IMAGE::ConstIterator const_traverser;
    typedef typename IMAGE::ConstIterator ConstIterator; 
    typedef typename IMAGE::ScanOrderIterator ScanOrderIterator;
    typedef typename IMAGE::ScanOrderIterator iterator;
    typedef typename IMAGE::ConstScanOrderIterator ConstScanOrderIterator;
    typedef typename IMAGE::ConstScanOrderIterator const_iterator; 

    typedef Diff2D difference_type;
    typedef Diff2D size_type;
    typedef ACCESSOR Accessor;
    typedef ACCESSOR ConstAccessor;
    typedef FixedBandsImage CloneType;
    
    friend class SelectBandImage;
    
    FixedBandsImage(int w, int h)
    : VariableBandsImage(),
      image_(new InnerImage(w,h))
    {
        init(getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
    }
    
    FixedBandsImage(Diff2D const & s)
    : image_(new InnerImage(s))
    {
        init(getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
    }
    
    FixedBandsImage(FixedBandsImage const & s)
    : VariableBandsImage(s),
      image_(s.image_)
    {}
    
    FixedBandsImage(InnerImage * i)
    : image_(i)
    {
        init(getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
    }
    
    FixedBandsImage & operator=(FixedBandsImage const & s)
    {
        if(this != &s)
        {
            image_ = s.image_;
            VariableBandsImage::operator=(s);
        }
        return *this;
    }
    
    virtual ~FixedBandsImage()
        {}
    
    virtual CloneType * clone() const
    {
        InnerImage * newimg = new InnerImage(size());
        copyImage(srcIterRange(upperLeft(), lowerRight()), destImage(*newimg));
        return new FixedBandsImage(newimg); 
    }
    
    virtual FixedBandsImage * shallowCopy() const
    {
        return new FixedBandsImage(*this); 
    }
    
    PixelType & operator[](Diff2D const & c)
        { return image_->operator[](c+roiUL_); }

    PixelType const & operator[](Diff2D const & c) const
        { return image_->operator[](c+roiUL_); }

    PixelType & operator()(int x, int y)
        { return image_->operator()(x+roiUL_.x,y+roiUL_.y); }

    PixelType const & operator()(int x, int y) const
        { return image_->operator()(x+roiUL_.x,y+roiUL_.y); }

    Iterator upperLeft()
        { return image_->upperLeft()+roiUL_; }
        
    Iterator lowerRight()
        { return upperLeft()+size(); }
        
    ConstIterator upperLeft() const
        { return image_->upperLeft()+roiUL_; }

    ConstIterator lowerRight() const
        { return upperLeft()+size(); }

    Accessor accessor() 
        { return Accessor(); }
    
    ConstAccessor accessor() const
        { return ConstAccessor(); }
    
  private:
    boost::shared_ptr<InnerImage> image_;
};

template <class IMAGE>
class FixedRGBImage
: public FixedBandsImage<IMAGE, VectorAccessor<TinyVector<GrayValue, 3> > >
{
  public:
  
    typedef 
        FixedBandsImage<IMAGE, VectorAccessor<TinyVector<GrayValue, 3> > >
        BaseType;
        
    typedef IMAGE InnerImage;
    typedef typename IMAGE::value_type value_type;
    typedef value_type          PixelType;
    typedef value_type          ConstPixelType;
    typedef value_type &        reference;
    typedef value_type const &  const_reference;
    typedef value_type *        pointer;
    typedef value_type const *  const_pointer;
    typedef typename IMAGE::Iterator traverser;
    typedef typename IMAGE::Iterator Iterator;
    typedef typename IMAGE::ConstIterator const_traverser;
    typedef typename IMAGE::ConstIterator ConstIterator; 
    typedef typename IMAGE::ScanOrderIterator ScanOrderIterator;
    typedef typename IMAGE::ScanOrderIterator iterator;
    typedef typename IMAGE::ConstScanOrderIterator ConstScanOrderIterator;
    typedef typename IMAGE::ConstScanOrderIterator const_iterator; 

    typedef Diff2D difference_type;
    typedef Diff2D size_type;
    typedef typename IMAGE::Accessor Accessor;
    typedef typename IMAGE::ConstAccessor ConstAccessor;
    typedef FixedRGBImage CloneType;
    
    friend class SelectBandImage;
    
    FixedRGBImage(int w, int h)
    : BaseType(w,h)
    {}
    
    FixedRGBImage(Diff2D const & s)
    : BaseType(s)
    {}
    
    FixedRGBImage(FixedRGBImage const & s)
    : BaseType(s)
    {}
    
    FixedRGBImage(InnerImage * i)
    : BaseType(i)
    {}
    
    FixedRGBImage & operator=(FixedRGBImage const & s)
    {
        BaseType::operator=(s);
        return *this;
    }
    
    virtual CloneType * clone() const
    {
        InnerImage * newimg = new InnerImage(size());
        copyImage(srcIterRange(upperLeft(), lowerRight()), destImage(*newimg));
        return new FixedRGBImage(newimg); 
    }
    
    virtual FixedRGBImage * shallowCopy() const
    {
        return new FixedRGBImage(*this); 
    }
    
    Accessor accessor() 
        { return Accessor(); }
    
    ConstAccessor accessor() const
        { return ConstAccessor(); }
};

typedef FixedBandsImage<BasicImage<TinyVector<GrayValue, 2> > > Vector2Image;
typedef FixedRGBImage<BasicImage<RGBValue<GrayValue> > > RGBImage;
typedef RGBImage::BaseType Vector3Image;
typedef FixedBandsImage<BasicImage<TinyVector<GrayValue, 4> > > Vector4Image;

template <class PIXELTYPE>
class ConstSelectBandIterator
: public VariableBandsIteratorBase<PIXELTYPE, ConstSelectBandIterator<PIXELTYPE> >
{
  public:
  
    typedef PIXELTYPE          value_type;
    typedef PIXELTYPE          PixelType;
    typedef PIXELTYPE const &  reference;
    typedef PIXELTYPE const &  index_reference;
    typedef PIXELTYPE const *  pointer;
    typedef Diff2D             difference_type;
    typedef image_traverser_tag  iterator_category;
    typedef IteratorAdaptor<ContigousMemoryColumnIteratorPolicy<ConstSelectBandIterator> > 
        row_iterator;
    typedef IteratorAdaptor<ContigousMemoryColumnIteratorPolicy<ConstSelectBandIterator> > 
        column_iterator;
  
    ConstSelectBandIterator(PIXELTYPE const * data, int width, int bands)
    : VariableBandsIteratorBase<PIXELTYPE, ConstSelectBandIterator>(
        const_cast<PIXELTYPE *>(data), width, bands, 1)
    {}
    
    ConstSelectBandIterator()
    {}
    
    ConstSelectBandIterator operator+(Diff2D const & s)
    {
        ConstSelectBandIterator res(*this);
        res += s;
        return res;
    }

    ConstSelectBandIterator operator-(Diff2D const & s)
    {
        ConstSelectBandIterator res(*this);
        res -= s;
        return res;
    }
    
    pointer operator->() const
    {
        return get();
    }
    
    reference operator*() const
    {
        return *get();
    }
    
    index_reference operator[](Diff2D const & dist) const
    {
        return *get(dist.x, dist.y);
    }
    
    index_reference operator()(int const & dx, int const & dy) const
    {
        return *get(dx, dy);
    }
    
    row_iterator rowIterator() const
    {
        return row_iterator(typename row_iterator::BaseType(get(), bands_));
    }
    
    column_iterator columnIterator() const
    {
        return column_iterator(typename column_iterator::BaseType(get(), width_*bands_));
    }
};

template <class PIXELTYPE>
class SelectBandIterator
: public VariableBandsIteratorBase<PIXELTYPE, SelectBandIterator<PIXELTYPE> >
{
  public:
  
    typedef PIXELTYPE          value_type;
    typedef PIXELTYPE          PixelType;
    typedef PIXELTYPE &        reference;
    typedef PIXELTYPE &        index_reference;
    typedef PIXELTYPE *        pointer;
    typedef Diff2D             difference_type;
    typedef image_traverser_tag  iterator_category;
    typedef IteratorAdaptor<ContigousMemoryColumnIteratorPolicy<SelectBandIterator> > 
        row_iterator;
    typedef IteratorAdaptor<ContigousMemoryColumnIteratorPolicy<SelectBandIterator> > 
        column_iterator;
  
    SelectBandIterator(PIXELTYPE * data, int width, int bands)
    : VariableBandsIteratorBase<PIXELTYPE, SelectBandIterator>(data, width, bands, 1)
    {}
    
    SelectBandIterator()
    {}
    
    SelectBandIterator operator+(Diff2D const & s)
    {
        SelectBandIterator res(*this);
        res += s;
        return res;
    }

    SelectBandIterator operator-(Diff2D const & s)
    {
        SelectBandIterator res(*this);
        res -= s;
        return res;
    }
    
    pointer operator->() const
    {
        return const_cast<pointer>(get());
    }
    
    reference operator*() const
    {
        return const_cast<reference>(*get());
    }
    
    index_reference operator[](Diff2D const & dist) const
    {
        return const_cast<index_reference>(*get(dist.x, dist.y));
    }
    
    index_reference operator()(int const & dx, int const & dy) const
    {
        return const_cast<index_reference>(*get(dx, dy));
    }
    
    row_iterator rowIterator() const
    {
        return row_iterator(typename 
            row_iterator::BaseType(const_cast<pointer>(get()), bands_));
    }
    
    column_iterator columnIterator() const
    {
        return column_iterator(typename 
            column_iterator::BaseType(const_cast<pointer>(get()), width_*bands_));
    }
};

class SingleBandImage
: public VariableBandsImage
{
  public:
  
    typedef GrayValue           value_type;
    typedef value_type          PixelType;
    typedef value_type          ConstPixelType;
    typedef value_type &        reference;
    typedef value_type const &  const_reference;
    typedef value_type *        pointer;
    typedef value_type const *  const_pointer;
    typedef SelectBandIterator<PixelType> traverser;
    typedef SelectBandIterator<PixelType> Iterator;
    typedef ConstSelectBandIterator<PixelType> const_traverser;
    typedef ConstSelectBandIterator<PixelType> ConstIterator; 
//    typedef ... ScanOrderIterator;
//    typedef ... iterator;
//    typedef ... ConstScanOrderIterator;
//    typedef ... const_iterator; 

    typedef Diff2D difference_type;
    typedef Diff2D size_type;
    typedef StandardAccessor<PixelType> Accessor;
    typedef StandardConstAccessor<ConstPixelType> ConstAccessor;    
    typedef SingleBandImage CloneType;
    
    virtual ~SingleBandImage()
        {}
    
    virtual CloneType * clone() const = 0;
    virtual SingleBandImage * shallowCopy() const = 0;
    
    PixelType & operator[](Diff2D const & c)
        { return *get(c.x, c.y); }

    PixelType const & operator[](Diff2D const & c) const
        { return *get(c.x, c.y); }

    PixelType & operator()(int x, int y)
        { return *get(x, y); }

    PixelType const & operator()(int x, int y) const
        { return *get(x, y); }

    Iterator upperLeft()
        { return Iterator(get(0,0), actualWidth(), actualBands_); }
        
    Iterator lowerRight()
        { return upperLeft()+size(); }
        
    ConstIterator upperLeft() const
        { return ConstIterator(get(0,0), actualWidth(), actualBands_); }

    ConstIterator lowerRight() const
        { return upperLeft()+size(); }

    Accessor accessor() 
        { return Accessor(); }
    
    ConstAccessor accessor() const
        { return ConstAccessor(); }
        
  protected:
    SingleBandImage()
    {}
    
    SingleBandImage(SingleBandImage const & s)
    : VariableBandsImage(s)
    {}
    
    SingleBandImage & operator=(SingleBandImage const & s)
    {
        VariableBandsImage::operator=(s);
        return *this;
    }
};

class GrayImage
: public SingleBandImage
{
  public:
  
    typedef BasicImage<GrayValue>  InnerImage;
    typedef InnerImage::value_type value_type;
    typedef value_type          PixelType;
    typedef value_type          ConstPixelType;
    typedef value_type &        reference;
    typedef value_type const &  const_reference;
    typedef value_type *        pointer;
    typedef value_type const *  const_pointer;
    typedef InnerImage::Iterator traverser;
    typedef InnerImage::Iterator Iterator;
    typedef InnerImage::ConstIterator const_traverser;
    typedef InnerImage::ConstIterator ConstIterator; 
    typedef InnerImage::ScanOrderIterator ScanOrderIterator;
    typedef InnerImage::ScanOrderIterator iterator;
    typedef InnerImage::ConstScanOrderIterator ConstScanOrderIterator;
    typedef InnerImage::ConstScanOrderIterator const_iterator; 

    typedef Diff2D difference_type;
    typedef Diff2D size_type;
    typedef InnerImage::Accessor Accessor;
    typedef InnerImage::ConstAccessor ConstAccessor;    
    typedef GrayImage CloneType;
    
    GrayImage(int w, int h)
    : image_(new InnerImage(w,h))
    {
        init(getDataPtr(image_->begin()), image_->size(), 1);
    }
    
    GrayImage(Diff2D const & s)
    : image_(new InnerImage(s))
    {
        init(getDataPtr(image_->begin()), image_->size(), 1);
    }
    
    GrayImage(GrayImage const & s)
    : SingleBandImage(s),
      image_(s.image_)
    {}
    
    GrayImage(InnerImage * i)
    : image_(i)
    {
        init(getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
    }
    
    GrayImage & operator=(GrayImage const & s)
    {
        if(this != &s)
        {
            image_ = s.image_;
            SingleBandImage::operator=(s);
        }
        return *this;
    }
    
    virtual ~GrayImage()
        {}
    
    virtual CloneType * clone() const
    {
        InnerImage * newimg = new InnerImage(size());
        copyImage(srcIterRange(upperLeft(), lowerRight()), destImage(*newimg));
        return new GrayImage(newimg); 
    }
    
    virtual GrayImage * shallowCopy() const
    {
        return new GrayImage(*this); 
    }
    
    PixelType & operator[](Diff2D const & c)
        { return image_->operator[](c+roiUL_); }

    PixelType const & operator[](Diff2D const & c) const
        { return image_->operator[](c+roiUL_); }

    PixelType & operator()(int x, int y)
        { return image_->operator()(x+roiUL_.x,y+roiUL_.y); }

    PixelType const & operator()(int x, int y) const
        { return image_->operator()(x+roiUL_.x,y+roiUL_.y); }

    Iterator upperLeft()
        { return image_->upperLeft()+roiUL_; }
        
    Iterator lowerRight()
        { return upperLeft()+size(); }
        
    ConstIterator upperLeft() const
        { return image_->upperLeft()+roiUL_; }

    ConstIterator lowerRight() const
        { return upperLeft()+size(); }

    Accessor accessor() 
        { return Accessor(); }
    
    ConstAccessor accessor() const
        { return ConstAccessor(); }
    
  private:
    boost::shared_ptr<InnerImage> image_;
};

class GrayImage;

class SelectBandImage
: public SingleBandImage
{
    class ImageHandle
    {
      public:
        virtual ~ImageHandle() {}
    };

    template <class IMAGE>
    class ImageHandleImpl
    : public ImageHandle
    {
        boost::shared_ptr<IMAGE> ptr_;

        typedef typename IMAGE::PixelType PixelType;

      public:
        ImageHandleImpl(boost::shared_ptr<IMAGE> const & p)
        : ptr_(p)
        {}
    };
    
  public:
  
    typedef GrayImage CloneType;
  
    template <class IMAGE>
    SelectBandImage(IMAGE const & s, int band)
    : imageHandle_(new ImageHandleImpl<typename IMAGE::InnerImage>(s.image_)),
      band_(band)
    {
        vigra_precondition(band >= 0 && band < s.actualBands(),
            "SelectBandImage(): band out of range.");
        init(getDataPtr(s.image_->begin())+band, s.actualSize(), s.actualBands(), 1);
        setROI(s.roiUpperLeft(), s.size());
    }
    
    SelectBandImage(SelectBandImage const & s)
    : SingleBandImage(s),
      imageHandle_(s.imageHandle_),
      band_(s.band_)
    {}

    SelectBandImage & operator=(SelectBandImage const & s)
    {
        if(this != &s)
        {
            imageHandle_ = s.imageHandle_;
            band_ = s.band_;
            SingleBandImage::operator=(s);
        }
        return *this;
    }
    
    virtual ~SelectBandImage()
        {}
    
    virtual CloneType * clone() const
    {
        CloneType::InnerImage * newimg = new CloneType::InnerImage(size());
        copyImage(srcIterRange(upperLeft(), lowerRight(), accessor()), destImage(*newimg));
        return new CloneType(newimg); 
    }
    
    virtual SelectBandImage * shallowCopy() const
    {
        return new SelectBandImage(*this); 
    }
    
    virtual void setROI(Diff2D const & ul, Diff2D const & new_size)
    {
        VariableBandsImage::setROI(ul, new_size);
    }
    
    virtual void setROI(int band)
    {
        vigra_precondition(band >= 0 && band < actualBands_,
            "SelectbandImage::setROI(): band out of range.");
        data_ = (data_ - band_ + band);
        band_ = band;
    }
        
 private:
    boost::shared_ptr<ImageHandle> imageHandle_;
    int band_;
};
    
/****************************************************************/

#define defineArgumentFactories(Image) \
template <class Accessor> \
inline triple<Image::ConstIterator, Image::ConstIterator, Accessor> \
srcImageRange(Image const & img, Accessor a) \
{ \
    return triple<Image::ConstIterator, Image::ConstIterator, Accessor>( \
            img.upperLeft(), img.lowerRight(), a); \
} \
 \
template <class Accessor> \
inline pair<Image::ConstIterator, Accessor> \
srcImage(Image const & img, Accessor a) \
{ \
    return pair<Image::ConstIterator, Accessor>( \
            img.upperLeft(), a); \
} \
 \
template <class Accessor> \
inline triple<Image::Iterator, Image::Iterator, Accessor> \
destImageRange(Image & img, Accessor a) \
{ \
    return triple<Image::Iterator, Image::Iterator, Accessor>( \
            img.upperLeft(), img.lowerRight(), a); \
} \
 \
template <class Accessor> \
inline pair<Image::Iterator, Accessor> \
destImage(Image & img, Accessor a) \
{ \
    return pair<Image::Iterator, Accessor>( \
            img.upperLeft(), a); \
} \
 \
template <class Accessor> \
inline pair<Image::ConstIterator, Accessor> \
maskImage(Image const & img, Accessor a) \
{ \
    return pair<Image::ConstIterator, Accessor>( \
            img.upperLeft(), a); \
} \
 \
inline triple<Image::ConstIterator, Image::ConstIterator, Image::ConstAccessor> \
srcImageRange(Image const & img) \
{ \
    return triple<Image::ConstIterator, Image::ConstIterator, Image::ConstAccessor>( \
            img.upperLeft(), img.lowerRight(), img.accessor()); \
} \
 \
inline pair<Image::ConstIterator, Image::ConstAccessor> \
srcImage(Image const & img) \
{ \
    return pair<Image::ConstIterator, Image::ConstAccessor>( \
            img.upperLeft(), img.accessor()); \
} \
 \
inline triple<Image::Iterator, Image::Iterator, Image::Accessor> \
destImageRange(Image & img) \
{ \
    return triple<Image::Iterator, Image::Iterator, Image::Accessor>( \
            img.upperLeft(), img.lowerRight(), img.accessor()); \
} \
 \
inline pair<Image::Iterator, Image::Accessor> \
destImage(Image & img) \
{ \
    return pair<Image::Iterator, Image::Accessor>( \
            img.upperLeft(), img.accessor()); \
} \
 \
inline pair<Image::ConstIterator, Image::ConstAccessor> \
maskImage(Image const & img) \
{ \
    return pair<Image::ConstIterator, Image::ConstAccessor>( \
            img.upperLeft(), img.accessor()); \
}

defineArgumentFactories(VariableBandsImage)
defineArgumentFactories(Vector2Image)
defineArgumentFactories(Vector3Image)
defineArgumentFactories(Vector4Image)
defineArgumentFactories(GrayImage)
defineArgumentFactories(RGBImage)
defineArgumentFactories(SelectBandImage)
defineArgumentFactories(SingleBandImage)

template<>
struct IteratorTraits<VariableBandsIterator >
{
    typedef StandardValueAccessor<VectorProxy> DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstVariableBandsIterator >
{
    typedef StandardConstValueAccessor<ConstVectorProxy> DefaultAccessor;
};  

template<>
struct IteratorTraits<SelectBandIterator<GrayValue> >
{
    typedef StandardAccessor<GrayValue> DefaultAccessor;
};  

template<>
struct IteratorTraits<ConstSelectBandIterator<GrayValue> >
{
    typedef StandardConstAccessor<GrayValue> DefaultAccessor;
};  


} // namespace vigra

#endif /* VIGRA_IMAGEHIERARCHY_HXX */
