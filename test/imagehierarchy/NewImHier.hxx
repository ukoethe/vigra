#ifndef VIGRA_IMAGEHIERARCHY_HXX
#define VIGRA_IMAGEHIERARCHY_HXX

#include "vigra/stdimage.hxx" 
#include "vigra/stdimagefunctions.hxx"
#include "vigra/imageiterator.hxx"
#include "vigra/accessor.hxx"
#include "vigra/utilities.hxx"
#include "boost/smart_ptr.hpp" 
#include "boost/static_assert.hpp"

namespace vigra {

typedef float GrayValue;                     

/** Stellt ein Pixel dar, also ein Einband, Zweiband, Dreiband oder Vierband Pixel, also je nachdem um was es
* fuer ein Pixel handelt, soviele Eintraege enthaelt auch der Vectorproxy, also z.B. ein Vierband-Pixel als 
* VectorProxy hat vier Eintraege
*/
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

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION                                      ////AEND_1////////
    template <int N>
    ConstVectorProxy(TinyVector<value_type, N> const & v)
    : data_(const_cast<TinyVector<value_type, N> &>(v).begin()), size_(N)
    {}
#else   
    ConstVectorProxy(TinyVector<value_type, 2> const & v)
    : data_(const_cast<TinyVector<value_type, 2> &>(v).begin()), size_(2)
    {}

    ConstVectorProxy(TinyVector<value_type, 3> const & v)
    : data_(const_cast<TinyVector<value_type, 3> &>(v).begin()), size_(3)
    {}

    ConstVectorProxy(TinyVector<value_type, 4> const & v)
    : data_(const_cast<TinyVector<value_type, 4> &>(v).begin()), size_(4)
    {}
#endif                                                            //////AEND_1//////////////

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
    
    
#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION                                  /////////////AEND_2//////////
    template <int N>
    void reset(TinyVector<value_type, N> const & v)
    {
        data_ = const_cast<TinyVector<value_type, N> &>(v).begin();
        size_ = N;
    }
#else
    void reset(TinyVector<value_type, 2> const & v)
    {
        data_ = const_cast<TinyVector<value_type, 2> &>(v).begin();
        size_ = 2;
    }
    
    void reset(TinyVector<value_type, 3> const & v)
    {
        data_ = const_cast<TinyVector<value_type, 3> &>(v).begin();
        size_ = 3;
    }
    
    void reset(TinyVector<value_type, 4> const & v)
    {
        data_ = const_cast<TinyVector<value_type, 4> &>(v).begin();
        size_ = 4;
    }
#endif                                                                      /////////////AEND_2//////////
    
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
    
#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION                          ///////////AEND_3///////////
    template <int N>
    operator TinyVector<value_type, N>() const
    {
        vigra_precondition(size_ == N, 
            "ConstVectorProxy::operator TinyVector(): size mismatch.");
        TinyVector<value_type, N> res;
        res.init(begin(), end());
        return res;
    }
#else
    operator TinyVector<value_type, 2>() const
    {
        vigra_precondition(size_ == 2, 
            "ConstVectorProxy::operator TinyVector(): size mismatch.");
        TinyVector<value_type, 2> res;
        res.init(begin(), end());
        return res;
    }
    operator TinyVector<value_type, 3>() const
    {
        vigra_precondition(size_ == 3, 
            "ConstVectorProxy::operator TinyVector(): size mismatch.");
        TinyVector<value_type, 3> res;
        res.init(begin(), end());
        return res;
    }
    operator TinyVector<value_type, 4>() const
    {
        vigra_precondition(size_ == 4, 
            "ConstVectorProxy::operator TinyVector(): size mismatch.");
        TinyVector<value_type, 4> res;
        res.init(begin(), end());
        return res;
    }
#endif                                                             ///////////AEND_3///////////
   
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
    
    ConstVectorProxy & operator=(ConstVectorProxy const & v);
  
    ConstVectorProxy & operator=(value_type const & f);
    
    template <int N>
    ConstVectorProxy & operator=(TinyVector<value_type, N> const & v);
    
    value_type * data_;
    int size_;
};

/** Stellt ein Pixel dar, also ein Einband, Zweiband, Dreiband oder Vierband Pixel, also je nachdem um was es
* fuer ein Pixel handelt, soviele Eintraege enthaelt auch der Vectorproxy, also z.B. ein Vierband-Pixel als 
* VectorProxy hat vier Eintraege
*/
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
    
#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION                                      ////AEND_1////////
    template <int N>
    VectorProxy(TinyVector<value_type, N> const & v)
    : ConstVectorProxy(v)
    {}
#else   
    VectorProxy(TinyVector<value_type, 2> const & v)
    : ConstVectorProxy(v)
    {}

    VectorProxy(TinyVector<value_type, 3> const & v)
    : ConstVectorProxy(v)
    {}

    VectorProxy(TinyVector<value_type, 4> const & v)
    : ConstVectorProxy(v)
    {}
#endif                                                                          //////AEND_1//////////////

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION                            /////////AEND_4////////////
    template <int N>
    VectorProxy & operator=(TinyVector<value_type, N> const & v)
    {
        vigra_precondition(size_ == N,
           "VectorProxy::operator=(): size mismatch.");
        for(int i=0; i<N; ++i)
            data_[i] = v[i];
        return *this;
    }
#else
    VectorProxy & operator=(TinyVector<value_type, 2> const & v)
    {
        vigra_precondition(size_ == 2,
           "VectorProxy::operator=(): size mismatch.");
        for(int i=0; i<2; ++i)
            data_[i] = v[i];
        return *this;
    }
    
    VectorProxy & operator=(TinyVector<value_type, 3> const & v)
    {
        vigra_precondition(size_ == 3,
           "VectorProxy::operator=(): size mismatch.");
        for(int i=0; i<3; ++i)
            data_[i] = v[i];
        return *this;
    }
    
    VectorProxy & operator=(TinyVector<value_type, 4> const & v)
    {
        vigra_precondition(size_ == 4,
           "VectorProxy::operator=(): size mismatch.");
        for(int i=0; i<4; ++i)
            data_[i] = v[i];
        return *this;
    }
    
#endif                                                              /////////AEND_4////////////
    
    VectorProxy & operator=(value_type const & v)
    {
        vigra_precondition(size_ == 1,
           "VectorProxy::operator=(): size mismatch.");
        data_[0] = v;
        return *this;
    }
    
    VectorProxy & operator=(VectorProxy const & v)
    {
        vigra_precondition(size_ == v.size(),
           "VectorProxy::operator=(): size mismatch.");
        for(int i=0; i<size_; ++i)
            data_[i] = v[i];
        return *this;
    }
    
    VectorProxy & operator=(ConstVectorProxy const & v)
    {
        vigra_precondition(size_ == v.size(),
           "VectorProxy::operator=(): size mismatch.");
        for(int i=0; i<size_; ++i)
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

template <int N>
inline bool operator==(TinyVector<GrayValue, N> const & l, ConstVectorProxy const & r)
{
    return r == l;
}

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
    typedef std::random_access_iterator_tag          iterator_category;
    
    
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
    
    int x, y;                                    //////////////////////Hier sind die x und y des VariableBandsIterator
    
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

/** 2D Iterator ueber ein VariableBandsImage
*/
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

    VariableBandsIterator operator-(Diff2D const & s)          //////////////////////////////hier die minus-Funktion
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
    
    row_iterator operator[](int dy) const
    {
        return row_iterator(
            row_iterator::BaseType(const_cast<pointer>(get(0, dy)), bands_, size_));
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

/** 2D Iterator ueber ein VariableBandsImage
*/
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
    
    row_iterator operator[](int dy) const
    {
        return row_iterator(row_iterator::BaseType(get(0, dy), bands_, size_));
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

/** Policyklasse des ScanOrderIterators, sie implementiert die Besonderheiten des Durchlaufs ueber eine 2D - Sequenz ,
* wobei die Besonderheiten letzendlich auf Durchlauf einer 1D - Sequenz zurueckgefuehrt werden. ScanOrderIteratorPolicy
* ist dazu da um spaeter (z.B. in VariableBandsImage) als Templateparameter an den "Interface" IteratorAdaptor
* uebergeben zu werden und dann durch  typedef IteratorAdaptor<ScanOrderIteratorPolicy<VariableBandsIterator>> ScanOrderIterator
* den gewoehnte ScanOrderIterator zu sehen bekommen!!! (IteratorAdapter ist die stellt die Schnittstelle zu Verfuegung
* um Mehrdimensionale Sequenzen als eine 1D-Sequenz zu behandeln). 
*/
template <class ImageIterator>
class ScanOrderIteratorPolicy
{
public:

/** Hilfsstruktur, stellt den zugrundeliegenden BaseType zur Verfuegung, speichert in sich den Anfang und das 
* Ende einer 2D-Sequenz (Matrix), die width Membervariable speichert die Anzahl der Elemente, die zum Anfang der
* der Initialisation(Konstruieren eines ROI-Objektes) da waren, also die Spaltenzahl(RegionOfInteress) wird zu
* Erzeugungszeit festgelegt  
*/          
struct ROI
{
    ROI(ImageIterator ul, ImageIterator lr)
    : current(ul), lowerRight(lr), width(lr.x-ul.x)  // initialisiert den Zeiger current zunaechst mit der linken oberen Ecke, den lowerRight entsprechend mit lowerRight Ecke   
    {}

    ROI()
    : width(0) // die Zeiger werden auf null gesetzt => VORSICHT NULLPOINTER !!! die Matrix ist leer => keine Spalten => width = 0
    {}

    ImageIterator current, lowerRight;                                           //der "Anfangsiterator" wird current genannt weil ueber ihn wird entschieden die aktuelle Position des Iterators daher nach ist der Name upperLeft nicht zutreffend, lowerRight bezeichnet wie gewoent das Ende der Sequenz, die in diesem Fall, der Matrix entspricht daher auch der Name   
    int width;                                                                   // speichert Anzahl der Spalten in der Matrix
};

     typedef ROI                                         BaseType;
     typedef typename ImageIterator::value_type          value_type;        // the adpator's value type
     typedef typename ImageIterator::reference           reference;         // the adpator's reference type (result of '*iter')
     typedef typename ImageIterator::index_reference     index_reference;   // the adpator's index_reference type (result of 'iter[n]')
     typedef typename ImageIterator::pointer             pointer;           // the adpator's pointer type (result of 'iter.operator->()')
     typedef std::random_access_iterator_tag             iterator_category; // the adpator's iterator category
     typedef int                                         difference_type;   // the adpator's difference type (result of 'iter1 - iter2', argument of 'iter[

     static void initialize(BaseType & d) {}

     /** Liefert den aktuell referenzierten Objekt.
     */
     static reference dereference(ROI const & r)
     { 
         return *r.current; 
     }

     /** Verschiebt seine aktuelle Position um n Schritte und liefert dann den aktuell referenzierten Objekt.
     */        
     static index_reference dereference(ROI r, difference_type n)
     { 
         advance(r,n);
         return *r.current; 
     }

     /** Haben zwei BaseType, die gleiche aktuelle Zeigerposition so sind sie aequivalent
     */        
     static bool equal(ROI const & l, ROI const & r)
     {
         return l.current == r.current;
     }

     /** Ist der d1 < d2, d.h. in der Matrix d1 steht hoeher als d2 oder stehen sie in der gleichen Zeile
     * und d1 befindet sich links von d2, so wird true zurueckgeliefert
     */
     static bool less(BaseType const & d1, BaseType const & d2)
     { 
         return d1.current.y == d2.current.y ?
                     d1.current.x < d2.current.x :
                     d1.current.y < d2.current.y; 
     }

     /** Liefert die Anzahl der Schritte der einer BaseType (d1) braucht um seine aktuelle Position
     * auf die des anderen (d2) zu veraendern. Dabei ist die Anordnungsposition von d1 und d2 egal, also
     * es kann auch ein negativer Wert zurueckgeliefert werden.
     */
     static difference_type difference(BaseType const & d1, BaseType const & d2)
     { 
         return d1.current.x - d2.current.x + d1.width*(d1.current.y - d2.current.y);         
     }

     /** Inkrementiert die aktuelle Position des ScanOrderIterators, unter Beruecksichtigung
     * des Falles, dass aktuelle Position am rechten Rande der Matrix ist
     */
     static void increment(ROI & r)
     {
         ++r.current.x;
         if(r.current.x == r.lowerRight.x)
         {
             r.current.x -= r.width;
             ++r.current.y;
         }
     }

     /** Dekrementiert die aktuelle Position des ScanOrderIterators, unter Beruecksichtigung
     * des Falles, dass aktuelle Position am linken Rande der Matrix ist
     */
     static void decrement(ROI & r)
     {
         --r.current.x;
         if(r.current.x + r.width == r.lowerRight.x)
         {
             r.current.x += r.width;
             --r.current.y;
         }
     }

     /** verschiebt den Iterator um n Stellen wobei n kann durchaus neagtiv
     * sein, das entspricht dan der Verschiebung nach links. 
     * Precondition : r und n muessen gueltige Werte besitzen, insbesondere soll geachtet werden, dass 
     * n <= width*"height" ist !!! 
     */
     static void advance(ROI & r, difference_type n) //n kann aber Diff2D sein; sollte es nicht int sein
     {
         if(n>0)                                        //ist n positiv so wird so wird die aktuelle Zeigerposition nach "rechts" verschoben
         {
             difference_type dy = (n) / r.width;        // dabei ist dy der Wert um den der Zeiger "in einer Matrix" nach unten verschoben werden soll 
             difference_type dx = (n) % r.width;        // entsprechend ist dx um den der Zeiger nach rechts verschoben werden soll, wobei dx kann ueber die MAtrix hinausragen
             if(dx >= r.lowerRight.x - r.current.x)     // laeuft dx ueber die "Matrix" hinaus so werden dx und dy korrigiert
             {
                 ++dy;                                  // der y-Wert wird einfach erhoeht
                 dx -= r.width;                       // und es wird dann um Anzahl der Zeilenelemente nach LINKS (zurueck) gelaufen
             }
             r.current.x += dx;                         // jetzt enthalten dx und dy korrekte Werte bezueglich denen kann die
             r.current.y += dy;                         // aktuelle Zeigerposition verschoben werden
         }
         else if(n < 0)                                 // ist n negativ so wird die aktuelle Zeigerposition nach "links" verschoben
         {
             n = -n;                                    // laesst sich mit positiver Zahl leichter rechnen
             difference_type dy = n / r.width;          // dabei ist dy der Wert um den der Zeiger "in einer Matrix" nach "oben" verschoben werden soll 
             difference_type dx = n % r.width;          // entsprechend ist dx um den der Zeiger nach "links" verschoben werden soll, wobei dx kann ueber die Matrix hinausragen

             if((r.current.x - r.lowerRight.x) + r.width < dx)// laeuft dx ueber die Kante (r.lowerRight.x - width)(entspricht der upperLeft.x in der Matrix)
             {
                                                        // dx und dy beduerfen die Korrektur weil es durchaus sein kann dass dx < width ist, aber es kann nicht so viele Schritte nach links gegangen werden
                    ++dy;                               // so werden die Differenzen um die die aktuelle Position veraendert werden soll korrigiert
                    dx -= r.width;
             }
             r.current.x -= dx;                         // jetzt enthalten dx und dy korrekte Werte bezueglich denen kann die
             r.current.y -= dy;                         // aktuelle Zeigerposition verschoben werden

         }

     }
};// end of ScanOrderIteratorPolicy

             
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

/** VaribleBandsImage stellt die Basis fuer die ganze Imagehierarchie, alle anderen ImageKlassen erben von ihr.
* Also, SingleBandImage, GrayImage, FixedBansImage usw. sind von VariableBandsImage abgeleitet.
*/    
class VariableBandsImage                          
{
    
  public:
    typedef VectorProxy                                                             value_type;
    typedef VectorProxy                                                             PixelType;
    typedef ConstVectorProxy                                                        ConstPixelType;
    typedef VectorProxy                                                             reference;
    typedef ConstVectorProxy                                                        const_reference;
    typedef GrayValue *                                                             pointer;
    typedef GrayValue const *                                                       const_pointer;
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<VariableBandsIterator> >        iterator;
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<VariableBandsIterator> >        ScanOrderIterator;
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<ConstVariableBandsIterator> >   const_iterator;
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<ConstVariableBandsIterator> >   ConstScanOrderIterator;
    typedef VariableBandsIterator                                                   traverser;
    typedef VariableBandsIterator                                                   Iterator;
    typedef ConstVariableBandsIterator                                              const_traverser;
    typedef ConstVariableBandsIterator                                              ConstIterator; 
    typedef Diff2D                                                                  difference_type;
    typedef Diff2D                                                                  size_type;
    typedef StandardValueAccessor<PixelType>                                        Accessor;
    typedef StandardConstValueAccessor<ConstPixelType>                              ConstAccessor;
    typedef VariableBandsImage                                                      CloneType;
    
    virtual ~VariableBandsImage() {}              
    
    /** Entspricht einer deepCopy(), d.h. es wird eine Kopie vom aufrufendem Objekt in Speicher abgelegt
    *   Erzeugt man einen Clone eines Bildes an dem ROI gesetz ist, dann kriegt man ein Bild
    *   mit der Groesse und Pixelinitialisierung der ROI
    */
    virtual CloneType * clone() const = 0;                                                      
    
    /** Erzeugt eine flache Kopie vom aufrufendem Objekt, d.h. es wird pointer 
    *   auf das aufrufende Objekt zurueckgeliefert.
    *   Erzeugt man eine shallowCopy eines Bildes an dem ROI gesetz ist, dann ist
    *   die shallowCopy eine Copy des gesamten Bildes und nicht nur der ROI
    *   dabei ist an der shallowCopy genau so die ROI gesetzt.
    */ 
    virtual VariableBandsImage * shallowCopy() const = 0;                                       
    
    /** Liefert die Breite der ROI !!! Nicht des Bildes
    */
    int width() const                                   
        { return roiLR_.x - roiUL_.x; } 
    
    /** Liefert die Hoehe der ROI !!! Nicht des Bildes
    */    
    int height() const                                  
        { return roiLR_.y - roiUL_.y; }
    
    /** Liefert die Groesse der ROI !!! Nicht des Bildes
    */    
    Diff2D size() const                                 
        { return roiLR_ - roiUL_; }
        
    /** Liefert die Groesse des Bildes !!! Nicht der ROI
    */
    Diff2D actualSize() const                           
        { return actualLR_; }
    
    /** Liefert die Breite des Bildes !!! Nicht der ROI
    */    
    int actualWidth() const                             
        { return actualSize().x; }
    
    /** Liefert die Hoehe des Bildes !!! Nicht der ROI
    */    
    int actualHeight() const                            
        { return actualSize().y; }
    
    /** Liefert Anzahl der Baender eines Pixels des Bildes !!! Es unterscheidet
    * sich von bands nur fuer SelectBandImage, sonst ist immer actualBands_ == bands_ 
    * erfuelt.
    */    
    int actualBands() const                         
        { return actualBands_; }
    
    /** Liefert Anzahl der Baender eines Pixels der ROI !!! Es unterscheidet
    * sich von actualBands nur fuer SelectBandImage, sonst ist immer actualBands_ == bands_ 
    * erfuelt. Fuer SelectBandImage ist bands_ gleich 1, da die ROI auf ein
    * Band gesetzt ist.
    */    
    int bands() const                                  
        { return bands_; }
    
    /** Liefert die linke obere Ecke der ROI.
    */    
    Diff2D roiUpperLeft() const                        
        { return roiUL_; }
    
    /** Liefert den 1D-iterator auf linke obere Ecke der ROI
    */
    ScanOrderIterator begin()
    {
        return ScanOrderIterator(ScanOrderIteratorPolicy<Iterator>::ROI(upperLeft(), lowerRight()));
    }    
    
    /** Liefert den 1D-iterator hinter der rechten unteren Ecke der ROI
    */
    ScanOrderIterator end()
    {
        return begin() + width()*height();
    }
    
    /** Liefert den ConstIterator auf linke obere Ecke der ROI
    */
    ConstScanOrderIterator begin() const
    {
        return ConstScanOrderIterator(ScanOrderIteratorPolicy<ConstIterator>::ROI(upperLeft(), lowerRight()));
    }    
    
    /** Liefert den ConstIterator hinter der rechten unteren Ecke der ROI
    */
    ConstScanOrderIterator end() const
    {
        return begin() + width()*height();
    }
    
    /** fuellt das komplette Image mit einem Pixel aus, also hinterher haben alle Pixel die gleiche "Farbe",
    * die hier uebergeben wurde
    */
    VariableBandsImage & init(ConstPixelType const & pixel)
    {
        initImage(upperLeft(), lowerRight(), accessor(), pixel);
        
        return *this;
    }    
    
    /** Setzt die ROI (Region of interess) auf dem aufrufendem Bild
    */
    virtual void setROI(Diff2D const & ul, Diff2D const & new_size)             // setROI -methode mit Argumenten Diff2D upperLeft und Diff2D size
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
    
    /** hebt die ROI wieder auf.
    */
    virtual void resetROI()
    {
        roiUL_ = Diff2D(0,0);
        roiLR_ = actualLR_;
    }
    
    /** Stellt den gewohnten Operatrot
    */
    PixelType operator[](Diff2D const & c)                          
        { return PixelType(get(c.x, c.y), bands_); }

    PixelType operator()(int x, int y) 
        { return PixelType(get(x,y), bands_); }

    ConstPixelType operator[](Diff2D const & c) const
        { return ConstPixelType(get(c.x, c.y), bands_); }

    ConstPixelType operator()(int x, int y) const
        { return ConstPixelType(get(x,y), bands_); }

    ScanOrderIterator operator[](int dy)                          
        { return begin()+dy*width(); }

    ConstScanOrderIterator operator[](int dy) const                         
        { return begin()+dy*width(); }

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
        
    /** Gibt an ob es ein Punkt mit Koordinaten (Diff2D) innerhalb des ganzen Bildes existiert
    */    
    bool isInside(Diff2D const & d) const
    {
        return d.x >= 0 && d.y >= 0 &&
               d.x < actualWidth() && d.y < actualHeight();    
    }
    /** Gibt an ob es ein Punkt mit Koordinaten (Diff2D) innerhalb der ROI existiert
    */
    bool isInsideROI(Diff2D const & d) const
    {
        return d.x >= roiUL_.x  && d.y >= roiUL_.y &&
               d.x < roiLR_.x && d.y < roiLR_.y;
    }
    
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
    
    /** Initialisiert die Membervariablen des Images
    */
    void initAdministrationData(GrayValue * data, Diff2D const & actualLR, int actualBands)        
    {
        initAdministrationData(data, actualLR, actualBands, actualBands);
    }
    
    /** Initialisiert die Membervariablen des Images
    */
    void initAdministrationData(GrayValue * data, Diff2D const & actualLR, int actualBands, int bands)
    {
        roiUL_ = Diff2D(0,0);
        roiLR_ = actualLR;
        actualLR_ = actualLR;
        bands_ = bands;
        actualBands_ = actualBands;
        data_ = data;
    }
    
    /** Liefert den Zeiger auf GrayValue
    */
    static GrayValue * getDataPtr(GrayValue * p)          
        { return p; }

    /** Liefert den Zeiger auf GrayValue,                                                            
    */
    template <class T>
    static GrayValue * getDataPtr(T * p)   
    {
        BOOST_STATIC_ASSERT( sizeof(T) ==
               BandsForPixelType<T>::bands*sizeof(GrayValue) );
        return &(*p)[0]; 
    }

    /** Liefert Zeiger auf GrayValue (float) zu einem bestimmten Punkt im Bilde
    */
    GrayValue * get(int dx, int dy)                  
        { return data_ + ((roiUL_.y + dy) * actualWidth() + roiUL_.x + dx)*actualBands_; }
    
    /** Liefert ConstZeiger auf GrayValue (float) zu einem bestimmten Punkt im Bilde
    */
    GrayValue const * get(int dx, int dy) const      
        { return data_ + ((roiUL_.y + dy) * actualWidth() + roiUL_.x + dx)*actualBands_; }
    
    Diff2D roiUL_, roiLR_, actualLR_;// roiUL_ - ist der Abstand von der linken oberen Ecke des Images bis zur linken oberen Ecke von ROI
    int bands_, actualBands_;
    GrayValue * data_;// GrayValue ist float 
};

class SelectBandImage;

template <class IMAGE, class ACCESSOR = typename IMAGE::Accessor>
class FixedBandsImage                                    
: public VariableBandsImage
{
  public:
    
    typedef IMAGE                                                   InnerImage;
    typedef typename ACCESSOR::value_type                           value_type;
    typedef value_type                                              PixelType;
    typedef value_type                                              ConstPixelType;
    typedef value_type &                                            reference;
    typedef value_type const &                                      const_reference;
    typedef value_type *                                            pointer;
    typedef value_type const *                                      const_pointer;
    typedef typename IMAGE::Iterator                                traverser;
    typedef typename IMAGE::Iterator                                Iterator;           //2D Iterator !!! Soll deprecated werden, wurde in traverser umbenannt
    typedef typename InnerImage::ConstIterator                           const_traverser;    
    typedef typename InnerImage::ConstIterator                           ConstIterator;      //2D ConstIterator !!! Soll deprecated werden, wurde in const_traverser umbenannt
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<Iterator> >     ScanOrderIterator;
    typedef ScanOrderIterator                                       iterator;           //1D Iterator entspricht den ScanOrderIterator
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<ConstIterator> > ConstScanOrderIterator;
    typedef ConstScanOrderIterator                                  const_iterator;     //1D ConstIterator entspricht den ConstScanOrderIterator
    typedef Diff2D                                                  difference_type;
    typedef Diff2D                                                  size_type;
    typedef ACCESSOR                                                Accessor;
    typedef ACCESSOR                                                ConstAccessor;
    
#ifndef NO_COVARIANT_RETURN_TYPES                                              ///////////////AEND_6///////////
    typedef FixedBandsImage CloneType;
#else
    typedef VariableBandsImage CloneType;
#endif                                                                      ///////////////AEND_6///////////   

    friend class SelectBandImage;
    /**  DefaultKonstruktor
    */
    FixedBandsImage()
    : VariableBandsImage(),
      image_(new InnerImage(0,0))
    {
        initAdministrationData(getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
    }
    
    FixedBandsImage(int w, int h)                                   
    : VariableBandsImage(),
      image_(new InnerImage(w,h))
    {   
        initAdministrationData(getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
    }
        
    /** Konstruktor, erzeugt ein FixedBandsImage mit width, height, und Pixel_value
    */
    FixedBandsImage(int w, int h, PixelType pixel_value)                                   
    : VariableBandsImage(),
      image_(new InnerImage(w,h, pixel_value))
    {
        initAdministrationData(getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
    }
    
    FixedBandsImage(Diff2D const & s)                                                             
    : VariableBandsImage(),
      image_(new InnerImage(s))
    {
        initAdministrationData(getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
    }
    
    FixedBandsImage(FixedBandsImage const & s)                       
    : VariableBandsImage(s),
      image_(s.image_)
    {}
    
    FixedBandsImage(InnerImage * i)                                                             
    : VariableBandsImage(),
      image_(i)
    {
        initAdministrationData(getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
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
    
    /** Destructor ist leer, weil es wird nur shared_ptr image_, als Klassenvariable angelegt und shared_ptr verwaltet
    * selber Referenzen und gibt den Speicherplatz frei
    */
    virtual ~FixedBandsImage()                                                                          
        {}
    
    /** Entspricht einer deepCopy(), d.h. es wird eine Kopie vom aufrufendem Objekt in Speicher abgelegt
    *   Erzeugt man einen Clone eines Bildes an dem ROI gesetz ist, dann kriegt man ein Bild
    *   mit der Groesse und Pixelinitialisierung der ROI
    */ 
    virtual CloneType * clone() const             
    {
        InnerImage * newimg = new InnerImage(size());
        copyImage(srcIterRange(upperLeft(), lowerRight()), destImage(*newimg));
        return new FixedBandsImage(newimg); 
    }
    
    /** Erzeugt eine flache Kopie vom aufrufendem Objekt, d.h. es wird pointer 
    *   auf das aufrufende Objekt zurueckgeliefert.
    *   Erzeugt man eine shallowCopy eines Bildes an dem ROI gesetz ist, dann ist
    *   die shallowCopy eine Copy des gesamten Bildes und nicht nur der ROI
    *   dabei ist an der shallowCopy genau so die ROI gesetzt.
    */ 
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
      
    PixelType * operator[](int c)                  
        { return image_->operator[](c+roiUL_.y) + roiUL_.x; }              

    PixelType  const * operator[](int c) const         
        { return image_->operator[](c+roiUL_.y) + roiUL_.x; }

    Iterator upperLeft()
        { return image_->upperLeft()+roiUL_; }
        
    Iterator lowerRight()
        { return upperLeft()+size(); }
        
    ConstIterator upperLeft() const
        { return image_->upperLeft()+roiUL_; }

    ConstIterator lowerRight() const
        { return upperLeft()+size(); }

    iterator begin()
        { return ScanOrderIterator(ScanOrderIteratorPolicy<Iterator>::ROI(upperLeft(), lowerRight())); }
        
    iterator end()
        { return begin() + width()*height(); }
        
    const_iterator begin() const
        { return ConstScanOrderIterator(ScanOrderIteratorPolicy<ConstIterator>::ROI(upperLeft(), lowerRight())); }

    const_iterator end() const
        { return begin() + width()*height(); }

    Accessor accessor() 
        { return Accessor(); }
    
    ConstAccessor accessor() const                          
        { return ConstAccessor(); }
    
  protected:
    boost::shared_ptr<InnerImage> image_;                                   // shared_ptr verwaltet selber Referenzen und gibt den Speicherplatz frei
    
};

template <class IMAGE>   
class FixedRGBImage : public FixedBandsImage<IMAGE, VectorAccessor<TinyVector<GrayValue, 3> > >
{
  public:
  
    typedef 
        FixedBandsImage<IMAGE, VectorAccessor<TinyVector<GrayValue, 3> > >  BaseType;
        
    typedef IMAGE                                                           InnerImage;
    typedef typename IMAGE::value_type                                      value_type;
    typedef value_type                                                      PixelType;
    typedef value_type                                                      ConstPixelType;
    typedef value_type &                                                    reference;
    typedef value_type const &                                              const_reference;
    typedef value_type *                                                    pointer;
    typedef value_type const *                                              const_pointer;
    typedef typename IMAGE::Iterator                                        traverser;
    typedef typename IMAGE::Iterator                                        Iterator;           //2D Iterator !!! Soll deprecated werden, wurde in traverser umbenannt
    typedef typename IMAGE::ConstIterator                                   const_traverser;
    typedef typename IMAGE::ConstIterator                                   ConstIterator;      //2D ConstIterator !!! Soll deprecated werden, wurde in const_traverser umbenannt
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<Iterator> >             ScanOrderIterator;
    typedef ScanOrderIterator                                               iterator;           //1D Iterator entspricht den ScanOrderIterator
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<ConstIterator> >        ConstScanOrderIterator;
    typedef ConstScanOrderIterator                                          const_iterator;     //1D ConstIterator entspricht den ConstScanOrderIterator

    typedef Diff2D                                                          difference_type;
    typedef Diff2D                                                          size_type;
    typedef typename IMAGE::Accessor                                        Accessor;
    typedef typename IMAGE::ConstAccessor                                   ConstAccessor;
#ifndef NO_COVARIANT_RETURN_TYPES                                               ///////////////AEND_7////////////
    typedef FixedRGBImage CloneType;
#else
    typedef typename BaseType::CloneType CloneType;
#endif                                                                       ///////////////AEND_7//////////// 
    friend class SelectBandImage;
    
    FixedRGBImage()
    {}
    
    FixedRGBImage(int w, int h)                
    : BaseType(w, h)
    {}
    
    FixedRGBImage(int w, int h, PixelType pix_type)                
    : BaseType(w, h, pix_type)
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
    
    /** Entspricht einer deepCopy(), d.h. es wird eine Kopie vom aufrufendem Objekt in Speicher abgelegt
    *   Erzeugt man einen Clone eines Bildes an dem ROI gesetz ist, dann kriegt man ein Bild
    *   mit der Groesse und Pixelinitialisierung der ROI
    */ 
    virtual CloneType * clone() const                        
    {
        InnerImage * newimg = new InnerImage(size());
        copyImage(srcIterRange(upperLeft(), lowerRight()), destImage(*newimg));
        return new FixedRGBImage(newimg); 
    }
    
    /** Erzeugt eine flache Kopie vom aufrufendem Objekt, d.h. es wird pointer 
    *   auf das aufrufende Objekt zurueckgeliefert.
    *   Erzeugt man eine shallowCopy eines Bildes an dem ROI gesetz ist, dann ist
    *   die shallowCopy eine Copy des gesamten Bildes und nicht nur der ROI
    *   dabei ist an der shallowCopy genau so die ROI gesetzt.
    */ 
    virtual FixedRGBImage * shallowCopy() const               
    {
        return new FixedRGBImage(*this); 
    }
    
    PixelType & operator[](Diff2D const & c)                 
        { return image_->operator[](c+roiUL_); }

    PixelType const & operator[](Diff2D const & c) const
        { return image_->operator[](c+roiUL_); }

    PixelType & operator()(int x, int y)
        { return image_->operator()(x+roiUL_.x,y+roiUL_.y); }

    PixelType const & operator()(int x, int y) const
        { return image_->operator()(x+roiUL_.x,y+roiUL_.y); }
      
    PixelType * operator[](int c)                  
        { return image_->operator[](c+roiUL_.y) + roiUL_.x; }              

    PixelType  const * operator[](int c) const         
        { return image_->operator[](c+roiUL_.y) + roiUL_.x; }

    Iterator upperLeft()
        { return image_->upperLeft()+roiUL_; }
        
    Iterator lowerRight()
        { return upperLeft()+size(); }
        
    ConstIterator upperLeft() const
        { return image_->upperLeft()+roiUL_; }

    ConstIterator lowerRight() const
        { return upperLeft()+size(); }

    iterator begin()
        { return ScanOrderIterator(ScanOrderIteratorPolicy<Iterator>::ROI(upperLeft(), lowerRight())); }
        
    iterator end()
        { return begin() + width()*height(); }
        
    const_iterator begin() const
        { return ConstScanOrderIterator(ScanOrderIteratorPolicy<ConstIterator>::ROI(upperLeft(), lowerRight())); }

    const_iterator end() const
        { return begin() + width()*height(); }

    Accessor accessor()                                       
        { return Accessor(); }
    
    ConstAccessor accessor() const                            
        {return ConstAccessor(); }
        
    
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
    
    row_iterator operator[](int dy) const
    {
        typedef typename row_iterator::BaseType BaseType;
        return row_iterator(BaseType(const_cast<pointer>(get(0, dy)), bands_));   //  braeuchte man hier nicht die Auswahl von bands?
    }
    
    index_reference operator()(int const & dx, int const & dy) const
    {
        return *get(dx, dy);
    }
    
    row_iterator rowIterator() const
    {
        typedef typename row_iterator::BaseType BaseType;
        return row_iterator(BaseType(get(), bands_));
    }
    
    column_iterator columnIterator() const
    {
        typedef typename column_iterator::BaseType BaseType;
        return column_iterator(BaseType(get(), width_*bands_));
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
    
    row_iterator operator[](int dy) const
    {
        typedef typename row_iterator::BaseType BaseType;
        return row_iterator(BaseType(const_cast<pointer>(get(0, dy)), bands_));
    }
    
    index_reference operator()(int const & dx, int const & dy) const
    {
        return const_cast<index_reference>(*get(dx, dy));
    }
    
    row_iterator rowIterator() const
    {
        typedef typename row_iterator::BaseType BaseType;
        return (row_iterator(BaseType(const_cast<pointer>(get()), bands_)));
    }
    
    column_iterator columnIterator() const
    {
        typedef typename column_iterator::BaseType BaseType;
        return (column_iterator(BaseType(const_cast<pointer>(get()), width_*bands_)));
    }
};

class SingleBandImage                          
: public VariableBandsImage                    
{
  public:
  
    typedef GrayValue                                           value_type;         //entsricht float
    typedef value_type                                          PixelType;          // -||- float
    typedef value_type                                          ConstPixelType;     // -||- float
    typedef value_type &                                        reference;          // -||- float const &
    typedef value_type const &                                  const_reference;    // -||- float &
    typedef value_type *                                        pointer;            // -||- float *
    typedef value_type const *                                  const_pointer;      // -||- float const *
    typedef SelectBandIterator<PixelType>                       traverser;          //"normaler" 2D Iterator
    typedef SelectBandIterator<PixelType>                       Iterator;           //2D Iterator !!! Soll deprecated werden, wurde in traverser umbenannt
    typedef ConstSelectBandIterator<PixelType>                  const_traverser;    //"normaler" 2D ConstIterator
    typedef ConstSelectBandIterator<PixelType>                  ConstIterator;      //2D ConstIterator !!! Soll deprecated werden, wurde in const_traverser umbenannt
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<Iterator> > ScanOrderIterator;  // entspricht dem iterator, also gewoehnlicher 1D Iterator 
    typedef ScanOrderIterator                                   iterator;           //1D Iterator entspricht den ScanOrderIterator
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<ConstIterator> > ConstScanOrderIterator;
    typedef ConstScanOrderIterator                              const_iterator;     //1D ConstIterator entspricht den ConstScanOrderIterator
    typedef Diff2D                                              difference_type;
    typedef Diff2D                                              size_type;
    typedef StandardAccessor<PixelType>                         Accessor;
    typedef StandardConstAccessor<ConstPixelType>               ConstAccessor;    
#ifndef NO_COVARIANT_RETURN_TYPES                                                                         ////////AEND_8///////
    typedef SingleBandImage CloneType;
#endif                                                                                                  ////////AEND_8///////   
    virtual ~SingleBandImage()                            
        {}
    
    /** Entspricht einer deepCopy(), d.h. es wird eine Kopie vom aufrufendem Objekt in Speicher abgelegt
    *   Erzeugt man einen Clone eines Bildes an dem ROI gesetz ist, dann kriegt man ein Bild
    *   mit der Groesse und Pixelinitialisierung der ROI
    */ 
    virtual CloneType * clone() const = 0;                                          // dadurch wird gekennzeichent, dass die Funktion "virtual CloneType * clone() const" in der abgeleiteten Klasse implementiert werden soll, das gleiche wie eine abstrakte Funktion   
    
    /** Erzeugt eine flache Kopie vom aufrufendem Objekt, d.h. es wird pointer 
    *   auf das aufrufende Objekt zurueckgeliefert.
    *   Erzeugt man eine shallowCopy eines Bildes an dem ROI gesetz ist, dann ist
    *   die shallowCopy eine Copy des gesamten Bildes und nicht nur der ROI
    *   dabei ist an der shallowCopy genau so die ROI gesetzt.
    */ 
    virtual SingleBandImage * shallowCopy() const = 0;                              // --//-- "virtual SingleBandImage * shallowCopy() const" --//--   
    
    PixelType & operator[](Diff2D const & c)              
        { return *get(c.x, c.y); }

    PixelType const & operator[](Diff2D const & c) const
        { return *get(c.x, c.y); }
 ////////////////////////////////////////////////


    ScanOrderIterator operator[](int dy)                          
        { return begin()+dy*width(); }

    ConstScanOrderIterator operator[](int dy) const                         
        { return begin()+dy*width(); }

////////////////////////////////////////////////

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
    
    iterator begin()
        { return ScanOrderIterator(ScanOrderIteratorPolicy<Iterator>::ROI(upperLeft(), lowerRight())); }
        
    iterator end()
        { return begin() + width()*height(); }
        
    const_iterator begin() const
        { return ConstScanOrderIterator(ScanOrderIteratorPolicy<ConstIterator>::ROI(upperLeft(), lowerRight())); }

    const_iterator end() const
        { return begin() + width()*height(); }
    
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
};// end of SingleBandsImage

/** entspricht BasicImage<float> und stellt die gleiche Operationen zur Verfuegung wie BasicImage,
* also ist ein EinbandImage, der in die Imagehierarchie eingebetet ist.
*/
class GrayImage                                           
: public SingleBandImage
{
  public:
  
    typedef BasicImage<GrayValue>                                      InnerImage;              // entspricht BasicImage<float>
    typedef InnerImage::value_type                                     value_type;              // -||- float
    typedef value_type                                                 PixelType;               // -||- float
    typedef value_type                                                 ConstPixelType;          // -||- float
    typedef value_type &                                               reference;               // -||- float &
    typedef value_type const &                                         const_reference;         // -||- float const &
    typedef value_type *                                               pointer;                 // -||- float *
    typedef value_type const *                                         const_pointer;           // -||- float const *
    typedef InnerImage::Iterator                                       traverser;               //"normaler"(gewoehnlicher) 2D-Iterator
    typedef InnerImage::Iterator                                       Iterator;                // soll deprecated werden !!! Entspricht dem traverser
    typedef InnerImage::ConstIterator                                  const_traverser;         //"normaler"(gewoehnlicher) 2D-ConstIterator
    typedef InnerImage::ConstIterator                                  ConstIterator;           // soll deprecated werden !!! Entspricht dem const_traverser
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<Iterator> >        ScanOrderIterator;       // 1D Iterator ueber das Bild (bzw. ROI), entspricht dem iterator
    typedef ScanOrderIterator                                          iterator;                // "normaler" 1D Iterator ueber das Bild (bzw.ROI)
    typedef IteratorAdaptor<ScanOrderIteratorPolicy<ConstIterator> >   ConstScanOrderIterator;  // 1D ConstIterator ueber das Bild (bzw. ROI), entspricht dem const_iterator
    typedef ConstScanOrderIterator                                     const_iterator;          // "normaler" 1D ConstIterator ueber das Bild (bzw.ROI)

    typedef Diff2D                                                      difference_type;
    typedef Diff2D                                                      size_type;
    typedef InnerImage::Accessor                                        Accessor;               // Accessor aus dem BasicImage<float>
    typedef InnerImage::ConstAccessor                                   ConstAccessor;          // ConstAccessor aus dem BasicImage<float>    
#ifndef NO_COVARIANT_RETURN_TYPES                                                                   //////AEND_9////////
    typedef GrayImage CloneType;
#endif                                                                                            //////AEND_9//////// 
    
    GrayImage()
    : image_(new InnerImage(0,0))                              
    {
        initAdministrationData(0, image_->size(), 1);
    }
    
    GrayImage(int w, int h)                                        
    : image_(new InnerImage(w,h))
    {   
        initAdministrationData((w*h == 0) ? 0 : getDataPtr(image_->begin()), image_->size(), 1);
    }
    
    GrayImage(int w, int h, GrayValue gray_value)               
    : image_(new InnerImage(w,h, gray_value))
    {   
        initAdministrationData((w*h == 0) ? 0 : getDataPtr(image_->begin()), image_->size(), 1);
    }
    
    GrayImage(Diff2D const & s)                                   
    : image_(new InnerImage(s))
    {   
        initAdministrationData((s.x*s.y == 0) ? 0 : getDataPtr(image_->begin()), image_->size(), 1);
    }
    
    GrayImage(GrayImage const & s)                                
    : SingleBandImage(s),
      image_(s.image_)
    {}
    
    GrayImage(InnerImage * i)                              
    : image_(i)
    {   
        initAdministrationData((i->width()*i->height() == 0) ? 0 : getDataPtr(image_->begin()), image_->size(), BandsForPixelType<PixelType>::bands);
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
    
    GrayImage & init(PixelType const & pixel)
    {
        initImage(upperLeft(), lowerRight(), accessor(), pixel);
        return *this;
    }
    
    virtual ~GrayImage()                                        
        {}
    
    /** Entspricht einer deepCopy(), d.h. es wird eine Kopie vom aufrufendem Objekt in Speicher abgelegt
    *   Erzeugt man einen Clone eines Bildes an dem ROI gesetz ist, dann kriegt man ein Bild
    *   mit der Groesse und Pixelinitialisierung der ROI
    */
    virtual CloneType * clone() const
    {   
        InnerImage * newimg = new InnerImage(size());
        copyImage(srcIterRange(upperLeft(), lowerRight()), destImage(*newimg));   // copyImage  kommt ueber imagefunctions.hxx aus copyimage.hxx; srcIterRange kommt aus stdimage.hxx 762; destImage kommt aus stdimage.hxx 733
        return new GrayImage(newimg); 
    }
    
    /** Erzeugt eine flache Kopie vom aufrufendem Objekt, d.h. es wird pointer 
    *   auf das aufrufende Objekt zurueckgeliefert.
    *   Erzeugt man eine shallowCopy eines Bildes an dem ROI gesetz ist, dann ist
    *   die shallowCopy eine Copy des gesamten Bildes und nicht nur der ROI
    *   dabei ist an der shallowCopy genau so die ROI gesetzt.
    */ 
    virtual GrayImage * shallowCopy() const     
    {   
        return new GrayImage(*this); 
    }                                                      
    
    PixelType & operator[](Diff2D const & c)               
        { return image_->operator[](c+roiUL_); }
    
    PixelType const & operator[](Diff2D const & c) const 
        { return image_->operator[](c+roiUL_); }
      
    PixelType * operator[](int c)                  
        { return image_->operator[](c+roiUL_.y) + roiUL_.x; }              

    PixelType  const * operator[](int c) const         
        { return image_->operator[](c+roiUL_.y) + roiUL_.x; }

    PixelType & operator()(int x, int y)
        { return image_->operator()(x+roiUL_.x,y+roiUL_.y); }

    PixelType const & operator()(int x, int y) const
        { return image_->operator()(x+roiUL_.x,y+roiUL_.y); }

    Iterator upperLeft()
        { return image_->upperLeft()+roiUL_;}
    
    Iterator lowerRight()
        { return upperLeft()+size(); }
        
    ConstIterator upperLeft() const
        { return image_->upperLeft()+roiUL_; }

    ConstIterator lowerRight() const
        { return upperLeft()+size(); }

    iterator begin()
        { return ScanOrderIterator(ScanOrderIteratorPolicy<Iterator>::ROI(upperLeft(), lowerRight())); }
        
    iterator end()
        { return begin() + width()*height(); }
        
    const_iterator begin() const
        { return ConstScanOrderIterator(ScanOrderIteratorPolicy<ConstIterator>::ROI(upperLeft(), lowerRight())); }

    const_iterator end() const
        { return begin() + width()*height(); }

    Accessor accessor() 
        { return Accessor(); }
    
    ConstAccessor accessor() const                    
        { return ConstAccessor(); }
    
  private:
    boost::shared_ptr<InnerImage> image_;
};// end of GrayImage

class SelectBandImage                              
: public SingleBandImage
{
    class ImageHandle
    {
      public:
        virtual ~ImageHandle() {}         /////////das verstehe ich nicht
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
#ifndef NO_COVARIANT_RETURN_TYPES                                           //////AEND_10//////
    typedef GrayImage CloneType;
    typedef SelectBandImage ShallowCopyType;
#else                                                                       //////AEND_10//////
    typedef CloneType ShallowCopyType;
#endif    
  
    template <class IMAGE>                              
    SelectBandImage(IMAGE const & s, int band)
    : imageHandle_(new ImageHandleImpl<typename IMAGE::InnerImage>(s.image_)),    
      band_(band)
    {
        vigra_precondition(band >= 0 && band < s.actualBands(),
            "SelectBandImage(): band out of range.");
        initAdministrationData(getDataPtr(s.image_->begin())+band, s.actualSize(), s.actualBands(), 1);
        if(s.width()*s.height() > 0)
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
    
    /** Entspricht einer deepCopy(), d.h. es wird eine Kopie vom aufrufendem Objekt in Speicher abgelegt
    *   Erzeugt man einen Clone eines Bildes an dem ROI gesetz ist, dann kriegt man ein Bild
    *   mit der Groesse und Pixelinitialisierung der ROI
    */
    virtual CloneType * clone() const           
    {
        CloneType::InnerImage * newimg = new CloneType::InnerImage(size());
        copyImage(srcIterRange(upperLeft(), lowerRight(), accessor()), destImage(*newimg));
        return new CloneType(newimg); 
    }
    
    /** Erzeugt eine flache Kopie vom aufrufendem Objekt, d.h. es wird pointer 
    *   auf das aufrufende Objekt zurueckgeliefert.
    *   Erzeugt man eine shallowCopy eines Bildes an dem ROI gesetz ist, dann ist
    *   die shallowCopy eine Copy des gesamten Bildes und nicht nur der ROI
    *   dabei ist an der shallowCopy genau so die ROI gesetzt.
    */ 
    virtual SelectBandImage * shallowCopy() const         
    {
        return new SelectBandImage(*this); 
    }
    
    virtual void setROI(Diff2D const & ul, Diff2D const & new_size)       
    {
        VariableBandsImage::setROI(ul, new_size);
    }
    
    /** die Methode "setROI(int band)" aendert den selektierten Band.
    * Parameter "band" ist der neue zu selektierender Band des Bildes
    */
    virtual void setROI(int band)                                        
    {
        vigra_precondition(band >= 0 && band < actualBands_,
            "SelectbandImage::setROI(): band out of range.");
        data_ = (data_ - band_ + band);
        band_ = band;
    }
    
    /** Liefert die "Nummer" des selektierten Bandes
    */
    int getSelectedBand()
    { return band_; }
        
 private:
    boost::shared_ptr<ImageHandle> imageHandle_;
    int band_;
};
    
/****************************************************************/

#if 0 // this didn't work with MSVC

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

#endif
////////////////////////////////////////ab hier bis Markierung soll eine Anpassung an den Kompiler sein
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

#define defineArgumentFactories2(Image, BaseImage) \
template <class Accessor> \
inline triple<BaseImage::ConstIterator, BaseImage::ConstIterator, Accessor> \
srcImageRange(Image const & img, Accessor a) \
{ \
    return triple<BaseImage::ConstIterator, BaseImage::ConstIterator, Accessor>( \
            img.upperLeft(), img.lowerRight(), a); \
} \
  \
template <class Accessor> \
inline pair<BaseImage::ConstIterator, Accessor> \
srcImage(Image const & img, Accessor a) \
{ \
    return pair<BaseImage::ConstIterator, Accessor>( \
            img.upperLeft(), a); \
} \
 \
template <class Accessor> \
inline triple<BaseImage::Iterator, BaseImage::Iterator, Accessor> \
destImageRange(Image & img, Accessor a) \
{ \
    return triple<BaseImage::Iterator, BaseImage::Iterator, Accessor>( \
            img.upperLeft(), img.lowerRight(), a); \
} \
 \
template <class Accessor> \
inline pair<BaseImage::Iterator, Accessor> \
destImage(Image & img, Accessor a) \
{ \
    return pair<BaseImage::Iterator, Accessor>( \
            img.upperLeft(), a); \
} \
 \
template <class Accessor> \
inline pair<BaseImage::ConstIterator, Accessor> \
maskImage(Image const & img, Accessor a) \
{ \
    return pair<BaseImage::ConstIterator, Accessor>( \
            img.upperLeft(), a); \
} \
 \
inline triple<BaseImage::ConstIterator, BaseImage::ConstIterator, Image::ConstAccessor> \
srcImageRange(Image const & img) \
{ \
    return triple<BaseImage::ConstIterator, BaseImage::ConstIterator, Image::ConstAccessor>( \
            img.upperLeft(), img.lowerRight(), img.accessor()); \
} \
 \
inline pair<BaseImage::ConstIterator, Image::ConstAccessor> \
srcImage(Image const & img) \
{ \
    return pair<BaseImage::ConstIterator, Image::ConstAccessor>( \
            img.upperLeft(), img.accessor()); \
} \
 \
inline triple<BaseImage::Iterator, BaseImage::Iterator, Image::Accessor> \
destImageRange(Image & img) \
{ \
    return triple<BaseImage::Iterator, BaseImage::Iterator, Image::Accessor>( \
            img.upperLeft(), img.lowerRight(), img.accessor()); \
} \
 \
inline pair<BaseImage::Iterator, Image::Accessor> \
destImage(Image & img) \
{ \
    return pair<BaseImage::Iterator, Image::Accessor>( \
            img.upperLeft(), img.accessor()); \
} \
 \
inline pair<BaseImage::ConstIterator, Image::ConstAccessor> \
maskImage(Image const & img) \
{ \
    return pair<BaseImage::ConstIterator, Image::ConstAccessor>( \
            img.upperLeft(), img.accessor()); \
}

defineArgumentFactories(VariableBandsImage)
defineArgumentFactories2(Vector2Image, FVector2Image)
//defineArgumentFactories2(Vector3Image, FVector3Image)
defineArgumentFactories2(Vector4Image, FVector4Image)
defineArgumentFactories(GrayImage)
defineArgumentFactories(RGBImage)
defineArgumentFactories(SelectBandImage)
defineArgumentFactories(SingleBandImage)
////////////////////////////////////////

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


