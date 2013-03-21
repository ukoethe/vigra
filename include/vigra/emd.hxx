/************************************************************************/
/*                                                                      */
/*               Copyright 2011-2013 by Ullrich Koethe                  */
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

#ifndef VIGRA_EMD_HXX
#define VIGRA_EMD_HXX


namespace vigra {

template <class FeatureType>
class Signature
{
    // implement this class
    // (used to represent the input histograms or signatures)
};

class EMDFlow
{
    // implement this class
    // (used to return the flow between the histograms, see the paper)
};

class EMDOptions
{
    // implement this class
    // (used to pass options to earthMoverDistance(), 
    //  look at existing options classes in VIGRA to see how this should work)
};

// implement suitable ground distance functors for the most common feature types

/*****************************************************************************/

    /** \brief Compute the earth mover distance between two histograms or signatures.
    
        ADD MORE DOCUMENTATION HERE.
        
        <b> Usage:</b>

        <b>\#include</b> \<vigra/emd.hxx\><br>
        Namespace: vigra
        
        \code
        typedef ... FeatureType;
        Signature<FeatureType> s1, s2;
        ... // fill the signatures with your data
        
            // to compute the distance, pass in the signatures and a functor 
            // that returns ground distances between feature pairs
        double distance = earthMoverDistance(s1, s2, MyGroundDistance());
        
            // optionally, you can also compute the flow between the signatures
        EMDFlow flow;
        double distance = earthMoverDistance(s1, s2, MyGroundDistance(), flow);
        
            // if no ground distance functor is given, the algorithm will use
            // the default ground distance for the given FeatureType
        double distance = earthMoverDistance(s1, s2);
        
            // options can be passed by the associated option object
        double distance = earthMoverDistance(s1, s2,
                                             EMDOptions().setSomeOption());
    */
template<class FeatureType, class GroundDistanceFunctor>
double 
earthMoverDistance(Signature<FeatureType> const & signature1, 
                   Signature<FeatureType> const & signature2, 
                   GroundDistanceFunctor const & groundDistance,
                   EMDFlow & flow,
                   EMDOptions const & options = EMDOptions())
{
}

    // don't compute the flow here
template<class FeatureType, class GroundDistanceFunctor>
void 
earthMoverDistance(Signature<FeatureType> const & signature1, 
                   Signature<FeatureType> const & signature2, 
                   GroundDistanceFunctor const & groundDistance,
                   EMDOptions const & options = EMDOptions())
{
}

    // use the default ground distance for the given FeatureType
    // (the default should be deduced automatically by template matching)
template<class FeatureType>
void 
earthMoverDistance(Signature<FeatureType> const & signature1, 
                   Signature<FeatureType> const & signature2,
                   EMDFlow & flow,
                   EMDOptions const & options = EMDOptions())
{
}

    // likewise, but without computing the flow
template<class FeatureType>
void 
earthMoverDistance(Signature<FeatureType> const & signature1, 
                   Signature<FeatureType> const & signature2,
                   EMDOptions const & options = EMDOptions())
{
}

}

#endif // VIGRA_EMD_HXX

