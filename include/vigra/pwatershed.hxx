#ifndef VIGRA_PWATERSHED_HXX
#define VIGRA_PWATERSHED_HXX

#include "vigra/multi_array.hxx"
#include "vigra/multi_watersheds.hxx"

namespace vigra{


    template<class T>
    void fuseLabelsInplace(
        MultiArrayView<3, T>  & a,
        const MultiArrayView<3, T> & b
    ){
        auto mx = *std::min_element(a.begin(), a.end());
        for(size_t i=0; i<a.size(); ++i){
            a[i] += mx + b[i];
        }
    }

    template<class VALUE_TYPE, class LABEL_TYPE>
    void parallelWatershed3D(
        const MultiArrayView<3, VALUE_TYPE> & evalMap,
        const Int64 w,
        const Int64 m,
        MultiArrayView<3, LABEL_TYPE>       & labels
    ){



        MultiArray<3, LABEL_TYPE>  buffer(labels.shape());

        typedef vigra::TinyVector<Int64, 3> Coord;

        const auto shape = evalMap.shape();

        Coord bCoord(0);
        Coord eCoord(shape);

        

        int left=1;
        int right=2;


        for(Int64 dim=0; dim<1; ++dim){
            for(Int64 x=w; x<shape[dim]; x+=w){

                std::cout<<"x "<<x<<"\n";

                {
                    auto start = std::max(Int64(0), x-w);
                    auto end   = std::min(Int64(shape[dim]),Int64(x+w+1));

                    auto bb = bCoord;
                    auto ee = eCoord;

                    bb[dim] = start;
                    ee[dim] = x;
                    buffer.subarray(bb, ee) = left;
                    bb[dim] = x;
                    ee[dim] = end;
                    buffer.subarray(bb, ee) = right;
                }
                // do the watershed
                
                {
                    auto start = std::max(Int64(0), x-m);
                    auto end   = std::min(Int64(shape[dim]),x+m+1);

                    auto bb = bCoord;
                    auto ee = eCoord;

                    bb[dim] = start;
                    ee[dim] = end;
                    auto subLabels = buffer.subarray(bb, ee);
                    auto subEval   = evalMap.subarray(bb, ee);

                    MultiArray<3, Int64> seeds(subLabels.shape());
                    seeds.bindAt(dim, 0) = left;
                    seeds.bindAt(dim, seeds.shape(dim)-1) = right;


                    watershedsMultiArray(subEval, seeds);
                    subLabels = seeds;
                }
                

                ++left;
                ++right;
            }
            if(dim == 0){
                labels = buffer;
            }
            else{
                fuseLabelsInplace(labels, buffer);
            }
        }
    }




}


#endif /* VIGRA_PWATERSHED_HXX */
