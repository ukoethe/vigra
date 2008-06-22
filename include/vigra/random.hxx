/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#include <time.h>

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
    UInt32 init[2] = { (UInt32)time(0), (UInt32)clock() };
    seed(init, 2, engine);
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
            generateNumbers();
            
        UInt32 y = state_[current_++];
        y ^= (y << 7) & 0x2b5b2500; 
        y ^= (y << 15) & 0xdb8b0000; 
        return y ^ (y >> 16);
    }
    
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
            generateNumbers();
            
        UInt32 x = state_[current_++];
        x ^= (x >> 11);
        x ^= (x << 7) & 0x9D2C5680U;
        x ^= (x << 15) & 0xEFC60000U;
        return x ^ (x >> 18);
    }
    
    void generateNumbers() const;

    static UInt32 twiddle(UInt32 u, UInt32 v) 
    {
        return (((u & 0x80000000U) | (v & 0x7FFFFFFFU)) >> 1)
                ^ ((v & 1U) ? 0x9908B0DFU : 0x0U);
    }

    void seedImpl(RandomSeedTag)
    {
        seed(RandomSeed, *this);
        generateNumbers();
    }

    void seedImpl(UInt32 theSeed)
    {
        seed(theSeed, *this);
        generateNumbers();
    }
    
    template<class Iterator>
    void seedImpl(Iterator init, UInt32 length)
    {
        seed(19650218U, *this);
        seed(init, length, *this);
        generateNumbers();
    }
};

void RandomState<MT19937>::generateNumbers() const
{
    for (int i = 0; i < (N - M); ++i)
    {
        state_[i] = state_[i + M] ^ twiddle(state_[i], state_[i + 1]);
    }
    for (int i = N - M; i < (N - 1); ++i)
    {
        state_[i] = state_[i + M - N] ^ twiddle(state_[i], state_[i + 1]);
    }
    state_[N - 1] = state_[M - 1] ^ twiddle(state_[N - 1], state_[0]);
    current_ = 0;
}

} // namespace detail

template <class Engine = detail::RandomState<detail::TT800> >
class RandomNumberGenerator
: public Engine
{
    mutable double normalCached_;
    mutable bool normalCachedValid_;
    
  public:
  
    RandomNumberGenerator()
    : normalCached_(0.0),
      normalCachedValid_(false)
    {}

    RandomNumberGenerator(RandomSeedTag)
    : normalCached_(0.0),
      normalCachedValid_(false)
    {
        this->seedImpl(RandomSeed);
    }

    RandomNumberGenerator(UInt32 theSeed)
    : normalCached_(0.0),
      normalCachedValid_(false)
    {
        this->seedImpl(theSeed);
    }

    template<class Iterator>
    RandomNumberGenerator(Iterator init, UInt32 length)
    : normalCached_(0.0),
      normalCachedValid_(false)
    {
        this->seedImpl(init, length);
    }

    void seed(RandomSeedTag)
    {
        this->seedImpl(RandomSeed);
        normalCachedValid_ = false;
    }

    void seed(UInt32 theSeed)
    {
        this->seedImpl(theSeed);
        normalCachedValid_ = false;
    }
    
    template<class Iterator>
    void seed(Iterator init, UInt32 length)
    {
        this->seedImpl(init, length);
        normalCachedValid_ = false;
    }

        // in [0, 2^32)
    UInt32 operator()() const
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

        // in [0,beyond) -- simple implementation when low bits are sufficiently random, which is
        // the case for TT800 and MT19937
    UInt32 uniformInt(UInt32 beyond) const
    {
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
    
        // in [0,1)
    double uniform53() const
    {
	    // make full use of the entire 53-bit mantissa of a double, by Isaku Wada
	    return ( (this->get() >> 5) * 67108864.0 + (this->get() >> 6)) * (1.0/9007199254740992.0); 
    }
    
        // in [0,1]
    double uniform() const
    {
        return (double)this->get() / 4294967295.0;
    }

        // in [lower, upper]
    double uniform(double lower, double upper) const
    {
        vigra_precondition(lower < upper,
          "RandomNumberGenerator::uniform(): lower bound must be smaller than upper bound."); 
        return uniform() * (upper-lower) + lower;
    }
    
    double normal() const;
    
    double normal(double mean, double stddev) const
    {
        vigra_precondition(stddev > 0.0,
          "RandomNumberGenerator::normal(): standard deviation must be positive."); 
        return normal()*stddev + mean;
    }
    
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

typedef RandomNumberGenerator<>  RandomTT800; 
typedef RandomNumberGenerator<detail::RandomState<detail::MT19937> > RandomMT19937;

inline RandomTT800   & randomTT800()   { return RandomTT800::global(); }
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

template <class Engine = RandomTT800>
class UniformIntRandomFunctor
{
    UInt32 lower_, difference_, factor_;
    Engine & generator_;
    bool useLowBits_;

  public:
  
    typedef UInt32 argument_type;
    typedef UInt32 result_type;

    UniformIntRandomFunctor(Engine & generator = Engine::global(),
                            bool useLowBits = false)
    : lower_(0), difference_(0xffffffff), factor_(1),
      generator_(generator),
      useLowBits_(useLowBits)
    {}
    
    UniformIntRandomFunctor(UInt32 lower, UInt32 upper, 
                            Engine & generator = Engine::global(),
                            bool useLowBits = false)
    : lower_(lower), difference_(upper-lower), 
      factor_(Engine::factorForUniformInt(difference_ + 1)),
      generator_(generator),
      useLowBits_(useLowBits)
    {
        vigra_precondition(lower < upper,
          "UniformIntRandomFunctor(): lower bound must be smaller than upper bound."); 
    }
    
    UInt32 operator()() const
    {
        if(difference_ == 0xffffffff) // lower_ is necessarily 0
            return generator_();
        else if(useLowBits_)
            return generator_() % (difference_ + 1) + lower_;
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

        /** Return a uniformly distributed integer random number such that
            <tt>0 <= i < beyond</tt>. This is a required interface for 
            <tt>std::random_shuffle</tt>. It ignores the limits specified 
            in the constructor.
        */
    UInt32 operator()(UInt32 beyond) const
    {
        if(beyond < 2)
            return 0;

        if(useLowBits_)
            return generator_() % beyond;
        else
            return generator_.uniformInt(beyond);
    }
};

template <class Engine>
class FunctorTraits<UniformIntRandomFunctor<Engine> >
{
  public:
    typedef UniformIntRandomFunctor<Engine> type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

template <class Engine = RandomTT800>
class UniformRandomFunctor
{
    double offset_, scale_;
    Engine & generator_;

  public:
  
    typedef double result_type;

    UniformRandomFunctor(Engine & generator = Engine::global())
    : offset_(0.0),
      scale_(1.0),
      generator_(generator)
    {}

    UniformRandomFunctor(double lower, double upper, 
                         Engine & generator = Engine::global())
    : offset_(lower),
      scale_(upper - lower),
      generator_(generator)
    {
        vigra_precondition(lower < upper,
          "UniformRandomFunctor(): lower bound must be smaller than upper bound."); 
    }
    
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

template <class Engine = RandomTT800>
struct NormalRandomFunctor
{
    double mean_, stddev_;
    Engine & generator_;

  public:
  
    typedef double result_type;

    NormalRandomFunctor(Engine & generator = Engine::global())
    : mean_(0.0),
      stddev_(1.0),
      generator_(generator)
    {}

    NormalRandomFunctor(double mean, double stddev, 
                        Engine & generator = Engine::global())
    : mean_(mean),
      stddev_(stddev),
      generator_(generator)
    {
        vigra_precondition(stddev > 0.0,
          "NormalRandomFunctor(): standard deviation must be positive."); 
    }
    
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


} // namespace vigra 

#endif // VIGRA_RANDOM_HXX
