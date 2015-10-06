#ifndef VIGRA_LABEL_DOWNSAMPLE_HXX
#define  VIGRA_LABEL_DOWNSAMPLE_HXX


#include "vigra/multi_resize.hxx"

namespace vigra{

    void downsampleLabels(
        const vigra::MultiArrayView<3, UInt8> & labels, 
        const UInt8 maxLabel,
        vigra::MultiArrayView<3, UInt8> & out 
    ){
        out = 0;

        vigra::MultiArray<3, UInt8> buffer(labels.shape());
        vigra::MultiArray<3, UInt8> bufferSmall(out.shape());

        for(UInt8 l=1; l<=maxLabel; ++l){
            buffer = 0.0f;
            for(size_t i=0; i<labels.size(); ++i){
                if(labels[i] == l){
                    buffer[i] = 1.0;
                }
            }
            resizeMultiArraySplineInterpolation(buffer, bufferSmall);
            for(size_t i=0; i<bufferSmall.size(); ++i){
                if(bufferSmall[i] > 0.5){
                    out[i] = l;
                }
            }
        }
    }

}


#endif /*VIGRA_LABEL_DOWNSAMPLE_HXX*/
