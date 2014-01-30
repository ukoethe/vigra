/************************************************************************/
/*                                                                      */
/*         Copyright 2002-2003 by Ullrich Koethe, Hans Meine            */
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

#ifndef VIGRA_MULTI_HISTOGRAMM_HXX
#define VIGRA_MULTI_HISTOGRAMM_HXX

#include <vigra/multi_array.hxx> 
#include <vigra/tinyvector.hxx> 
#include <vigra/multi_gridgraph.hxx> 
#include <vigra/multi_convolution.hxx>


namespace vigra{


    template< unsigned int DIM , class T_DATA ,unsigned int CHANNELS, class T_HIST >
    void multi_gaussian_histogram(
        const MultiArrayView<DIM, TinyVector<T_DATA,CHANNELS> > & image,
        const TinyVector<T_DATA,CHANNELS> minVals,
        const TinyVector<T_DATA,CHANNELS> maxVals,
        const size_t bins,
        const float sigma,
        const float sigmaBin,
        MultiArrayView<DIM+2 , T_HIST>    histogram
    ){
        typedef vigra::GridGraph< DIM , boost::undirected_tag> Graph;
        typedef typename Graph::NodeIt graph_scanner;
        typedef typename Graph::Node   Node;
        typedef T_HIST ValueType;
        typedef TinyVector<ValueType,CHANNELS> ChannelsVals;
        typedef typename MultiArrayView<DIM+2 , T_HIST>::difference_type HistCoord;
        const Graph g(image.shape());
        const ChannelsVals nBins(bins);
         // define abreviations for the required iterators


        std::fill(histogram.begin(),histogram.end(),0.0);
        // iterate over all nodes (i.e. pixels)
        for (graph_scanner n(g); n != lemon::INVALID; ++n){
            const Node node(*n);
            ChannelsVals  binIndex = image[node];
            binIndex -=minVals;
            binIndex /=maxVals;
            binIndex *=nBins;



            HistCoord histCoord;
            for(size_t d=0;d<DIM;++d){
                histCoord[d]=node[d];
            }
            for(size_t c=0;c<CHANNELS;++c){
                const float fi = binIndex[c];
                const size_t bi = std::floor(fi+0.5);
                histCoord[DIM]=std::min(bi,static_cast<size_t>(bins-1));
                histCoord[DIM+1]=c;
                histogram[histCoord]+=1.0;
            }
        }

        MultiArray<DIM+2 , T_HIST>    histogramBuffer(histogram);
        Kernel1D<float> gauss,gaussBin;
        gauss.initGaussian(sigma);
        gaussBin.initGaussian(sigmaBin);
        for(size_t c=0;c<CHANNELS;++c){

            // histogram for one channel
            MultiArrayView<DIM+1,T_HIST> histc       = histogram.bindOuter(c);
            MultiArrayView<DIM+1,T_HIST> histcBuffer = histogram.bindOuter(c);

            // convolve along all spatial axis and bin axis
            if(DIM==2){
                convolveMultiArrayOneDimension(srcMultiArrayRange(histc), destMultiArray(histcBuffer), 0, gauss);
                convolveMultiArrayOneDimension(srcMultiArrayRange(histcBuffer), destMultiArray(histc), 1, gauss);
                convolveMultiArrayOneDimension(srcMultiArrayRange(histc), destMultiArray(histcBuffer), 2, gaussBin);
                histc=histcBuffer;
            }


        }


    }


}
//end namespace vigra

#endif //VIGRA_MULTI_HISTOGRAMM_HXX