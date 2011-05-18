/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
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


#ifndef VIGRA_RANDOM_HXX
#define VIGRA_RANDOM_HXX

#include "mathutil.hxx"
#include "functortraits.hxx"

#include <ctime>

namespace vigra {

enum RandomSeedTag { RandomSeed };

namespace detail {

enum RandomEngineTag { TT800, MT19937 };

template<RandomEngineTag EngineTag>
struct RandomState;

template <RandomEngineTag EngineTag>
void seed(UInt32 theSeed, RandomState<EngineTag> & engine)
{
    engine.state_[0] = theSeed;
    for(UInt32 i=1; i<RandomState<EngineTag>::N; ++i)
    {
        engine.state_[i] = 1812433253U * (engine.state_[i-1] ^ (engine.state_[i-1] >> 30)) + i;
    }
}

template <class Iterator, RandomEngineTag EngineTag>
void seed(Iterator init, UInt32 key_length, RandomState<EngineTag> & engine)
{
    const UInt32 N = RandomState<EngineTag>::N;
    int k = (int)std::max(N, key_length);
    UInt32 i = 1, j = 0;
    Iterator data = init;
    for (; k; --k) 
    {
        engine.state_[i] = (engine.state_[i] ^ ((engine.state_[i-1] ^ (engine.state_[i-1] >> 30)) * 1664525U))
                           + *data + j; /* non linear */
        ++i; ++j; ++data;
        
        if (i >= N) 
        { 
            engine.state_[0] = engine.state_[N-1]; 
            i=1; 
        }
        if (j>=key_length)
        { 
            j=0;
            data = init;
        }
    }

    for (k=N-1; k; --k) 
    {
        engine.state_[i] = (engine.state_[i] ^ ((engine.state_[i-1] ^ (engine.state_[i-1] >> 30)) * 1566083941U))
                           - i; /* non linear */
        ++i;
        if (i>=N) 
        { 
            engine.state_[0] = engine.state_[N-1]; 
            i=1; 
        }
    }

    engine.state_[0] = 0x80000000U; /* MSB is 1; assuring non-zero initial array */ 
}

template <RandomEngineTag EngineTag>
void seed(RandomSeedTag, RandomState<EngineTag> & engine)
{
    static UInt32 globalCount = 0;
    UInt32 init[3] = { (UInt32)time(0), (UInt32)clock(), ++globalCount };
    seed(init, 3, engine);
}


    /* Tempered twister TT800 by M. Matsumoto */
template<>
struct RandomState<TT800>
{
    static const UInt32 N = 25, M = 7;
    
    mutable UInt32 state_[N];
    mutable UInt32 current_;
                   
    RandomState()
    : current_(0)
    {
        UInt32 seeds[N] = { 
	        0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
	        0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
	        0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
	        0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
	        0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
        };
         
        for(UInt32 i=0; i<N; ++i)
            state_[i] = seeds[i];
    }

  protected:  

    UInt32 get() const
    {
        if(current_ == N)
            generateNumbers<void>();
            
        UInt32 y = state_[current_++];
        y ^= (y << 7) & 0x2b5b2500; 
        y ^= (y << 15) & 0xdb8b0000; 
        return y ^ (y >> 16);
    }
    
    template <class DUMMY>
    void generateNumbers() const;

    void seedImpl(RandomSeedTag)
    {
        seed(RandomSeed, *this);
    }

    void seedImpl(UInt32 theSeed)
    {
        seed(theSeed, *this);
    }
    
    template<class Iterator>
    void seedImpl(Iterator init, UInt32 length)
    {
        seed(init, length, *this);
    }
};

template <class DUMMY>
void RandomState<TT800>::generateNumbers() const
{
    UInt32 mag01[2]= { 0x0, 0x8ebfd028 };

    for(UInt32 i=0; i<N-M; ++i)
    {
    	state_[i] = state_[i+M] ^ (state_[i] >> 1) ^ mag01[state_[i] % 2];
    }
    for (UInt32 i=N-M; i<N; ++i) 
    {
        state_[i] = state_[i+(M-N)] ^ (state_[i] >> 1) ^ mag01[state_[i] % 2];
    }
    current_ = 0;
}

    /* Mersenne twister MT19937 by M. Matsumoto */
template<>
struct RandomState<MT19937>
{
    static const UInt32 N = 624, M = 397;
    
    mutable UInt32 state_[N];
    mutable UInt32 current_;
                   
    RandomState()
    : current_(0)
    {
        seed(19650218U, *this);
    }

  protected:  

    UInt32 get() const
    {
        if(current_ == N)
            generateNumbers<void>();
            
        UInt32 x = state_[current_++];
        x ^= (x >> 11);
        x ^= (x << 7) & 0x9D2C5680U;
        x ^= (x << 15) & 0xEFC60000U;
        return x ^ (x >> 18);
    }
    
    template <class DUMMY>
    void generateNumbers() const;

    static UInt32 twiddle(UInt32 u, UInt32 v) 
    {
        return (((u & 0x80000000U) | (v & 0x7FFFFFFFU)) >> 1)
                ^ ((v & 1U) ? 0x9908B0DFU : 0x0U);
    }

    void seedImpl(RandomSeedTag)
    {
        seed(RandomSeed, *this);
        generateNumbers<void>();
    }

    void seedImpl(UInt32 theSeed)
    {
        seed(theSeed, *this);
        generateNumbers<void>();
    }
    
    template<class Iterator>
    void seedImpl(Iterator init, UInt32 length)
    {
        seed(19650218U, *this);
        seed(init, length, *this);
        generateNumbers<void>();
    }
};

template <class DUMMY>
void RandomState<MT19937>::generateNumbers() const
{
    for (unsigned int i = 0; i < (N - M); ++i)
    {
        state_[i] = state_[i + M] ^ twiddle(state_[i], state_[i + 1]);
    }
    for (unsigned int i = N - M; i < (N - 1); ++i)
    {
        state_[i] = state_[i + M - N] ^ twiddle(state_[i], state_[i + 1]);
    }
    state_[N - 1] = state_[M - 1] ^ twiddle(state_[N - 1], state_[0]);
    current_ = 0;
}

} // namespace detail


/** \addtogroup RandomNumberGeneration Random Number Generation

     High-quality random number generators and related functors.
*/
//@{

/** Generic random number generator.

    The actual generator is passed in the template argument <tt>Engine</tt>. Two generators
    are currently available:
    <ul>
    <li> <tt>RandomMT19937</tt>: The state-of-the-art <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">Mersenne Twister</a> with a state length of 2<sup>19937</sup> and very high statistical quality.
    <li> <tt>RandomTT800</tt>: (default) The Tempered Twister, a simpler predecessor of the Mersenne Twister with period length 2<sup>800</sup>.
    </ul>
    
    Both generators have been designed by <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/eindex.html">Makoto Matsumoto</a>. 
    
    <b>Traits defined:</b>
    
    \verbatim FunctorTraits<RandomNumberGenerator<Engine> >::isInitializer \endverbatim
    is true (<tt>VigraTrueType</tt>).
*/
template <class Engine = detail::RandomState<detail::TT800> >
class RandomNumberGenerator
: public Engine
{
    mutable double normalCached_;
    mutable bool normalCachedValid_;
    
  public:
  
        /** Create a new random generator object with standard seed.
            
            Due to standard seeding, the random numbers generated will always be the same. 
            This is useful for debugging.
        */
    RandomNumberGenerator()
    : normalCached_(0.0),
      normalCachedValid_(false)
    {}
  
        /** Create a new random generator object with a random seed.
        
            The seed is obtained from the machines current <tt>clock()</tt> and <tt>time()</tt>
            values.
        
            <b>Usage:</b>
            \code
            RandomNumberGenerator<> rnd = RandomNumberGenerator<>(RandomSeed);
            \endcode
        */
    RandomNumberGenerator(RandomSeedTag)
    : normalCached_(0.0),
      normalCachedValid_(false)
    {
        this->seedImpl(RandomSeed);
    }
  
        /** Create a new random generator object from the given seed.
            
            The same seed will always produce identical random sequences.
        */
    RandomNumberGenerator(UInt32 theSeed)
    : normalCached_(0.0),
      normalCachedValid_(false)
    {
        this->seedImpl(theSeed);
    }

        /** Create a new random generator object from the given seed sequence.
            
            Longer seed sequences lead to better initialization in the sense that the generator's 
            state space is covered much better than is possible with 32-bit seeds alone.
        */
    template<class Iterator>
    RandomNumberGenerator(Iterator init, UInt32 length)
    : normalCached_(0.0),
      normalCachedValid_(false)
    {
        this->seedImpl(init, length);
    }

  
        /** Re-initialize the random generator object with a random seed.
        
            The seed is obtained from the machines current <tt>clock()</tt> and <tt>time()</tt>
            values.
        
            <b>Usage:</b>
            \code
            RandomNumberGenerator<> rnd = ...;
            ...
            rnd.seed(RandomSeed);
            \endcode
        */
    void seed(RandomSeedTag)
    {
        this->seedImpl(RandomSeed);
        normalCachedValid_ = false;
    }

        /** Re-initialize the random generator object from the given seed.
            
            The same seed will always produce identical random sequences.
        */
    void seed(UInt32 theSeed)
    {
        this->seedImpl(theSeed);
        normalCachedValid_ = false;
    }

        /** Re-initialize the random generator object from the given seed sequence.
            
            Longer seed sequences lead to better initialization in the sense that the generator's 
            state space is covered much better than is possible with 32-bit seeds alone.
        */
    template<class Iterator>
    void seed(Iterator init, UInt32 length)
    {
        this->seedImpl(init, length);
        normalCachedValid_ = false;
    }

        /** Return a uniformly distributed integer random number in [0, 2<sup>32</sup>).
            
            That is, 0 &lt;= i &lt; 2<sup>32</sup>. 
        */
    UInt32 operator()() const
    {
        return this->get();
    }

        /** Return a uniformly distributed integer random number in [0, 2<sup>32</sup>).
            
            That is, 0 &lt;= i &lt; 2<sup>32</sup>. 
        */
    UInt32 uniformInt() const
    {
        return this->get();
    }


#if 0 // difficult implementation necessary if low bits are not sufficiently random
        // in [0,beyond)
    UInt32 uniformInt(UInt32 beyond) const
    {
        if(beyond < 2)
            return 0;

        UInt32 factor = factorForUniformInt(beyond);
        UInt32 res = this->get() / factor;

        // Use rejection method to avoid quantization bias.
        // On average, we will need two raw random numbers to generate one.
        while(res >= beyond)
            res = this->get() / factor;
        return res;
    }
#endif /* #if 0 */

        /** Return a uniformly distributed integer random number in [0, <tt>beyond</tt>).
            
            That is, 0 &lt;= i &lt; <tt>beyond</tt>. 
        */
    UInt32 uniformInt(UInt32 beyond) const
    {
        // in [0,beyond) -- simple implementation when low bits are sufficiently random, which is
        // the case for TT800 and MT19937
        if(beyond < 2)
            return 0;

        UInt32 remainder = (NumericTraits<UInt32>::max() - beyond + 1) % beyond;
        UInt32 lastSafeValue = NumericTraits<UInt32>::max() - remainder;
        UInt32 res = this->get();

        // Use rejection method to avoid quantization bias.
        // We will need two raw random numbers in amortized worst case.
        while(res > lastSafeValue)
            res = this->get();
        return res % beyond;
    }
    
        /** Return a uniformly distributed double-precision random number in [0.0, 1.0).
            
            That is, 0.0 &lt;= i &lt; 1.0. All 53-bit bits of the mantissa are random (two 32-bit integers are used to 
            create this number).
        */
    double uniform53() const
    {
	    // make full use of the entire 53-bit mantissa of a double, by Isaku Wada
	    return ( (this->get() >> 5) * 67108864.0 + (this->get() >> 6)) * (1.0/9007199254740992.0); 
    }
    
        /** Return a uniformly distributed double-precision random number in [0.0, 1.0].
            
            That is, 0.0 &lt;= i &lt;= 1.0. This nuber is computed by <tt>uniformInt()</tt> / 2<sup>32</sup>, 
            so it has effectively only 32 random bits. 
        */
    double uniform() const
    {
        return (double)this->get() / 4294967295.0;
    }

        /** Return a uniformly distributed double-precision random number in [lower, upper].
           
            That is, <tt>lower</tt> &lt;= i &lt;= <tt>upper</tt>. This number is computed 
            from <tt>uniform()</tt>, so it has effectively only 32 random bits. 
        */
    double uniform(double lower, double upper) const
    {
        vigra_precondition(lower < upper,
          "RandomNumberGenerator::uniform(): lower bound must be smaller than upper bound."); 
        return uniform() * (upper-lower) + lower;
    }

        /** Return a standard normal variate (Gaussian) random number.
           
            Mean is zero, standard deviation is 1.0. It uses the polar form of the 
            Box-Muller transform.
        */
    double normal() const;
    
        /** Return a normal variate (Gaussian) random number with the given mean and standard deviation.
           
            It uses the polar form of the Box-Muller transform.
        */
    double normal(double mean, double stddev) const
    {
        vigra_precondition(stddev > 0.0,
          "RandomNumberGenerator::normal(): standard deviation must be positive."); 
        return normal()*stddev + mean;
    }
    
        /** Access the global (program-wide) instance of the present random number generator.
        
            Normally, you will create a local generator by one of the constructor calls. But sometimes
            it is useful to have all program parts access the same generator.
        */
    static RandomNumberGenerator & global()
    {
        static RandomNumberGenerator generator;
        return generator;
    }

    static UInt32 factorForUniformInt(UInt32 range)
    {
        return (range > 2147483648U || range == 0)
                     ? 1
                     : 2*(2147483648U / ceilPower2(range));
    }
};

template <class Engine>
double RandomNumberGenerator<Engine>::normal() const
{
    if(normalCachedValid_)
    {
        normalCachedValid_ = false;
        return normalCached_;
    }
    else
    {
        double x1, x2, w;
        do 
        {
             x1 = uniform(-1.0, 1.0);
             x2 = uniform(-1.0, 1.0);
             w = x1 * x1 + x2 * x2;
        } 
        while ( w > 1.0 || w == 0.0);
        
        w = std::sqrt( -2.0 * std::log( w )  / w );

        normalCached_ = x2 * w;
        normalCachedValid_ = true;

        return x1 * w;
    }
}

    /** Shorthand for the TT800 random number generator class.
    */
typedef RandomNumberGenerator<>  RandomTT800; 

    /** Shorthand for the TT800 random number generator class (same as RandomTT800).
    */
typedef RandomNumberGenerator<>  TemperedTwister; 

    /** Shorthand for the MT19937 random number generator class.
    */
typedef RandomNumberGenerator<detail::RandomState<detail::MT19937> > RandomMT19937;

    /** Shorthand for the MT19937 random number generator class (same as RandomMT19937).
    */
typedef RandomNumberGenerator<detail::RandomState<detail::MT19937> > MersenneTwister;

    /** Access the global (program-wide) instance of the TT800 random number generator.
    */
inline RandomTT800   & randomTT800()   { return RandomTT800::global(); }

    /** Access the global (program-wide) instance of the MT19937 random number generator.
    */
inline RandomMT19937 & randomMT19937() { return RandomMT19937::global(); }

template <class Engine>
class FunctorTraits<RandomNumberGenerator<Engine> >
{
  public:
    typedef RandomNumberGenerator<Engine> type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};


/** Functor to create uniformely distributed integer random numbers.

    This functor encapsulates the appropriate functions of the given random number
    <tt>Engine</tt> (usually <tt>RandomTT800</tt> or <tt>RandomMT19937</tt>)
    in an STL-compatible interface. 
    
    <b>Traits defined:</b>
    
    \verbatim FunctorTraits<UniformIntRandomFunctor<Engine> >::isInitializer \endverbatim
    and
    \verbatim FunctorTraits<UniformIntRandomFunctor<Engine> >::isUnaryFunctor \endverbatim
    are true (<tt>VigraTrueType</tt>).
*/
template <class Engine = RandomTT800>
class UniformIntRandomFunctor
{
    UInt32 lower_, difference_, factor_;
    Engine const & generator_;
    bool useLowBits_;

  public:
  
    typedef UInt32 argument_type; ///< STL required functor argument type
    typedef UInt32 result_type; ///< STL required functor result type

        /** Create functor for uniform random integers in the range [0, 2<sup>32</sup>) 
            using the given engine.
            
            That is, the generated numbers satisfy 0 &lt;= i &lt; 2<sup>32</sup>.
        */
    explicit UniformIntRandomFunctor(Engine const & generator = Engine::global() )
    : lower_(0), difference_(0xffffffff), factor_(1),
      generator_(generator),
      useLowBits_(true)
    {}
    
        /** Create functor for uniform random integers in the range [<tt>lower</tt>, <tt>upper</tt>]
            using the given engine.
            
            That is, the generated numbers satisfy <tt>lower</tt> &lt;= i &lt;= <tt>upper</tt>.
            \a useLowBits should be set to <tt>false</tt> when the engine generates
            random numbers whose low bits are significantly less random than the high
            bits. This does not apply to <tt>RandomTT800</tt> and <tt>RandomMT19937</tt>,
            but is necessary for simpler linear congruential generators.
        */
    UniformIntRandomFunctor(UInt32 lower, UInt32 upper, 
                            Engine const & generator = Engine::global(),
                            bool useLowBits = true)
    : lower_(lower), difference_(upper-lower), 
      factor_(Engine::factorForUniformInt(difference_ + 1)),
      generator_(generator),
      useLowBits_(useLowBits)
    {
        vigra_precondition(lower < upper,
          "UniformIntRandomFunctor(): lower bound must be smaller than upper bound."); 
    }
    
        /** Return a random number as specified in the constructor.
        */
    UInt32 operator()() const
    {
        if(difference_ == 0xffffffff) // lower_ is necessarily 0
            return generator_();
        else if(useLowBits_)
            return generator_.uniformInt(difference_+1) + lower_;
        else
        {
            UInt32 res = generator_() / factor_;

            // Use rejection method to avoid quantization bias.
            // On average, we will need two raw random numbers to generate one.
            while(res > difference_)
                res = generator_() / factor_;
            return res + lower_;
        }
    }

        /** Return a uniformly distributed integer random number in the range [0, <tt>beyond</tt>).
        
            That is, 0 &lt;= i &lt; <tt>beyond</tt>. This is a required interface for 
            <tt>std::random_shuffle</tt>. It ignores the limits specified 
            in the constructor and the flag <tt>useLowBits</tt>.
        */
    UInt32 operator()(UInt32 beyond) const
    {
        if(beyond < 2)
            return 0;

        return generator_.uniformInt(beyond);
    }
};

template <class Engine>
class FunctorTraits<UniformIntRandomFunctor<Engine> >
{
  public:
    typedef UniformIntRandomFunctor<Engine> type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraTrueType  isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

/** Functor to create uniformely distributed double-precision random numbers.

    This functor encapsulates the function <tt>uniform()</tt> of the given random number
    <tt>Engine</tt> (usually <tt>RandomTT800</tt> or <tt>RandomMT19937</tt>)
    in an STL-compatible interface. 
    
    <b>Traits defined:</b>
    
    \verbatim FunctorTraits<UniformIntRandomFunctor<Engine> >::isInitializer \endverbatim
    is true (<tt>VigraTrueType</tt>).
*/
template <class Engine = RandomTT800>
class UniformRandomFunctor
{
    double offset_, scale_;
    Engine const & generator_;

  public:
  
    typedef double result_type; ///< STL required functor result type

        /** Create functor for uniform random double-precision numbers in the range [0.0, 1.0] 
            using the given engine.
            
            That is, the generated numbers satisfy 0.0 &lt;= i &lt;= 1.0.
        */
    UniformRandomFunctor(Engine const & generator = Engine::global())
    : offset_(0.0),
      scale_(1.0),
      generator_(generator)
    {}

        /** Create functor for uniform random double-precision numbers in the range [<tt>lower</tt>, <tt>upper</tt>]
            using the given engine.
            
            That is, the generated numbers satisfy <tt>lower</tt> &lt;= i &lt;= <tt>upper</tt>.
        */
    UniformRandomFunctor(double lower, double upper, 
                         Engine & generator = Engine::global())
    : offset_(lower),
      scale_(upper - lower),
      generator_(generator)
    {
        vigra_precondition(lower < upper,
          "UniformRandomFunctor(): lower bound must be smaller than upper bound."); 
    }
    
        /** Return a random number as specified in the constructor.
        */
    double operator()() const
    {
        return generator_.uniform() * scale_ + offset_;
    }
};

template <class Engine>
class FunctorTraits<UniformRandomFunctor<Engine> >
{
  public:
    typedef UniformRandomFunctor<Engine> type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

/** Functor to create normal variate random numbers.

    This functor encapsulates the function <tt>normal()</tt> of the given random number
    <tt>Engine</tt> (usually <tt>RandomTT800</tt> or <tt>RandomMT19937</tt>)
    in an STL-compatible interface. 
    
    <b>Traits defined:</b>
    
    \verbatim FunctorTraits<UniformIntRandomFunctor<Engine> >::isInitializer \endverbatim
    is true (<tt>VigraTrueType</tt>).
*/
template <class Engine = RandomTT800>
class NormalRandomFunctor
{
    double mean_, stddev_;
    Engine const & generator_;

  public:
  
    typedef double result_type; ///< STL required functor result type

        /** Create functor for standard normal random numbers 
            using the given engine.
            
            That is, mean is 0.0 and standard deviation is 1.0.
        */
    NormalRandomFunctor(Engine const & generator = Engine::global())
    : mean_(0.0),
      stddev_(1.0),
      generator_(generator)
    {}

        /** Create functor for normal random numbers with goven mean and standard deviation
            using the given engine.
        */
    NormalRandomFunctor(double mean, double stddev, 
                        Engine & generator = Engine::global())
    : mean_(mean),
      stddev_(stddev),
      generator_(generator)
    {
        vigra_precondition(stddev > 0.0,
          "NormalRandomFunctor(): standard deviation must be positive."); 
    }
    
        /** Return a random number as specified in the constructor.
        */
    double operator()() const
    {
        return generator_.normal() * stddev_ + mean_;
    }

};

template <class Engine>
class FunctorTraits<NormalRandomFunctor<Engine> >
{
  public:
    typedef UniformRandomFunctor<Engine>  type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

//@}

} // namespace vigra 

#endif // VIGRA_RANDOM_HXX
