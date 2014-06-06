#ifndef VIGRA_TG_HXX
#define VIGRA_TG_HXX

#include <hash_map>

#include <vigra/multi_array.hxx>
#include <vigra/timing.hxx>

#define TG_DEBUG

/********************************************************/
/*                                                      */
/*           regionVolumeToCrackEdgeVolume              */
/*                                                      */
/********************************************************/

template< class InputValueType, class OutputValueType >
void regionVolumeToCrackEdgeVolume(
    const vigra::MultiArrayView<3, InputValueType>& seg,
    vigra::MultiArrayView<3, OutputValueType>& result
) {
    using vigra::MultiArrayIndex;
    using vigra::Shape3;
    typedef vigra::MultiArray<3, InputValueType> InputArray;
    typedef typename InputArray::difference_type Diff3D;
  
    Shape3 desiredShape(2*seg.shape(0)-1, 2*seg.shape(1)-1, 2*seg.shape(2)-1);
    if( result.shape() != desiredShape ) { 
        std::stringstream ss;
        ss << "result has wrong shape " << result.shape()
           << ", was expecting " << desiredShape;
        throw vigra::PreconditionViolation(ss.str().c_str());
    }
  
    size_t s0 = seg.shape(0);
    size_t s1 = seg.shape(1);
    size_t s2 = seg.shape(2);
    
    //go over all 2-cells, 2 cells have 1 dimension odd, 2 dimensions even
    size_t nCells[3] = {(s2-1)*(s0-1)*(s1-1),
                        s2*( (s0-1)*(s1-1)         ) + (s2-1)*( s1*(s0-1) + s0*(s1-1) ),
                        s2*( s1*(s0-1) + s0*(s1-1) ) + (s2-1)*( s0*s1                 )
    };
   
    #ifdef TG_DEBUG
    size_t n0 = 0;
    size_t n1 = 0;
    size_t n2 = 0;
    #endif
    
    for(MultiArrayIndex i = 0; i<result.shape(0); ++i) {
    for(MultiArrayIndex J = 0; J< ((i%2 != 0) ? 1 : 2); ++J) { //start j first with 0, then with 1
    for(MultiArrayIndex j = J; j<result.shape(1); j+=2) {
    for(MultiArrayIndex k = (i%2 + j%2 == 0)  ? 1 : 0; k<result.shape(2); k+=2) {
        MultiArrayIndex dim = (i%2 != 0) ? 0 : ((j%2) != 0 ? 1 : 2);
        #ifdef TG_DEBUG
        if( (i%2 + j%2 + k%2) != 1 ) {
            std::stringstream ss; ss << "2-cell error at coor " << i << " " << j << " " << k << std::endl;
            throw std::runtime_error(ss.str());
        }
        ++n2;
        #endif
        
        Diff3D d1(i/2,j/2,k/2);
        Diff3D d2(d1);
        d2[dim] += 1;
        if(seg[d1] != seg[d2]) {
            result(i,j,k) = 1;
        }
    }
    }
    }
    }
    
    const MultiArrayIndex nhood2[4][2] = { {1,0}, {0,1}, {-1,0}, {0,-1} };
    
    //go over all 1-cells , 1 cells have 2 dimension odd, 1 dimensions even
    for(MultiArrayIndex i = 0; i<result.shape(0); ++i) {
    for(MultiArrayIndex J = (i%2 ==0) ? 1 : 0; J<2; ++J) { //start j first with 0, then with 1
    for(MultiArrayIndex j = J; j<result.shape(1); j+=2) {
    for(MultiArrayIndex k = (i%2!=0 && j%2!=0) ? 0 : 1; k<result.shape(2); k+=2) {
        #ifdef TG_DEBUG
        if( (i%2 + j%2 + k%2) != 2 ) {
            std::stringstream ss; ss << "1-cell error at coor " << i << " " << j << " " << k << std::endl;
            throw std::runtime_error(ss.str());
        }
        ++n1;
        #endif
    
        MultiArrayIndex dim1 = ((i%2 == 1) ? 0 : ((j%2 == 1) ? 1 : 2));
        MultiArrayIndex dim2 = ((dim1 == 0 && j%2 == 1) ? 1 : 2);
        //std::cout << i << " " << j << " " << k << " ::: " << dim1 << " " << dim2 << std::endl;
       
        for(size_t l=0; l<4; ++l) {
            Diff3D d(i,j,k);
            d[dim1] += nhood2[l][0];
            d[dim2] += nhood2[l][1];
            if( result[d] == 1 ) {
                result(i,j,k) = 1;
                break;
            }
        }
    }
    }
    }
    }
    
    //go over all 0-cells
    const MultiArrayIndex nhood1[6][3] = { {1,0,0}, {0,1,0}, {0,0,1}, {-1,0,0}, {0,-1,0}, {0,0,-1} };
    for(MultiArrayIndex i= 1; i<result.shape(0); i+=2) {
    for(MultiArrayIndex j= 1; j<result.shape(1); j+=2) {
    for(MultiArrayIndex k= 1; k<result.shape(2); k+=2) {
        #ifdef TG_DEBUG
        if( (i%2 + j%2 + k%2) != 3 ) {
            std::stringstream ss; ss << "0-cell error at coor " << i << " " << j << " " << k << std::endl;
            throw std::runtime_error(ss.str());
        }
        ++n0;
        #endif
        
        for(size_t l=0; l<6; ++l) {
            Diff3D d(i,j,k);
            d[0] += nhood1[l][0];
            d[1] += nhood1[l][1];
            d[2] += nhood1[l][2];
            if( result[d] == 1 ) {
                result(i,j,k) = 1;
                break;
            }
        }
    }
    }
    }
    
    #ifdef TG_DEBUG
    if(n0 != nCells[0]) {
        std::stringstream ss; ss << "n0= " << n0 << " <-> " << nCells[0] << std::endl;
        throw std::runtime_error(ss.str());
    }
    if(n1 != nCells[1]) {
        std::stringstream ss; ss << "n1= " << n1 << " <-> " << nCells[1] << std::endl;
        throw std::runtime_error(ss.str());
    }
    if(n2 != nCells[2]) {
        std::stringstream ss; ss << "n2= " << n2 << " <-> " << nCells[2] << std::endl;
    }
    #endif
}
    
    
/********************************************************/
/*                                                      */
/*           tg2surface                                 */
/*                                                      */
/********************************************************/

template<class TgValueType, class SegValueType, class CoordinateValueType>
void vectorialMaxDistanceBruteForce(
    const vigra::MultiArrayView<3, TgValueType>& tg,
    const vigra::MultiArrayView<3, SegValueType>& seg,
    vigra::MultiArrayView<3, vigra::TinyVector<CoordinateValueType, 3> >& maxDist
) {
    std::cout << "tg2surface" << std::endl;
    
    if(seg.shape() != maxDist.shape()) {
        std::stringstream err; err << "shape mismatch: ";
        err << "seg.shape() = " << seg.shape() << " vs. " << maxDist.shape();
        throw std::runtime_error(err.str());
    }
    if(tg.shape() != vigra::Shape3(2*seg.shape(0)-1, 2*seg.shape(1)-1, 2*seg.shape(2)-1)) {
        throw std::runtime_error("tg has wrong shape");
    }
    
    using vigra::MultiArrayIndex;
    typedef vigra::TinyVector<CoordinateValueType, 3> CoordinateType;
    typedef __gnu_cxx::hash_map<SegValueType, std::vector<CoordinateType> > Surfaces;
    
    Surfaces surfaces;
    
    const MultiArrayIndex nhood26[26][3] = {
        { 0, 0, 1}, /* 0 */
        { 0, 0,-1}, /* 1 */
        { 0, 1, 0}, /* 2 */
        { 0,-1, 0}, /* 3 */
        { 0, 1, 1}, /* 4 */
        { 0, 1,-1}, /* 5 */
        { 0,-1, 1}, /* 6 */
        { 0,-1,-1}, /* 7 */
        { 1, 0, 0}, /* 8 */
        {-1, 0, 0}, /* 9 */
        { 1, 0, 1}, /* 10 */
        { 1, 0,-1}, /* 11 */
        {-1, 0, 1}, /* 12 */
        {-1, 0,-1}, /* 13 */
        { 1, 1, 0}, /* 14 */
        { 1,-1, 0}, /* 15 */
        {-1, 1, 0}, /* 16 */
        {-1,-1, 0}, /* 17 */
        { 1, 1, 1}, /* 18 */
        { 1, 1,-1}, /* 19 */
        { 1,-1, 1}, /* 20 */
        {-1, 1, 1}, /* 21 */
        { 1,-1,-1}, /* 22 */
        {-1,-1, 1}, /* 23 */
        {-1, 1,-1}, /* 24 */
        {-1,-1,-1} /* 25 */
    };
    
    USETICTOC;
    TIC
    for(MultiArrayIndex k=0; k<tg.shape(2); k+=2) {
    for(MultiArrayIndex j=0; j<tg.shape(1); j+=2) {
    for(MultiArrayIndex i=0; i<tg.shape(0); i+=2) {
        const SegValueType& label = seg(i/2,j/2,k/2);
       
        for(MultiArrayIndex l=0; l<26; ++l) {
            const MultiArrayIndex ii = i+nhood26[l][0];
            const MultiArrayIndex jj = j+nhood26[l][1];
            const MultiArrayIndex kk = k+nhood26[l][2];
            
            if(ii>=0 && jj>=0 && kk>=0 && ii<tg.shape(0) && jj<tg.shape(1) && kk < tg.shape(2)) {
                if(tg(i,j,k) != tg(ii,jj,kk)) {
                    surfaces[label].push_back( CoordinateType(ii,jj,kk) );
                }
            }
        }
    }
    }
    }
    std::cout << "  building surfaces took " <<  TOCS << std::endl;
    
    for(typename Surfaces::iterator it=surfaces.begin();
        it!=surfaces.end(); ++it)
    {
        std::vector<CoordinateType>& v = it->second;
        std::sort(v.begin(), v.end());
        typename std::vector<CoordinateType>::iterator e = std::unique(v.begin(),  v.end());
        v.erase(e, v.end());
    }
   
    double dist = 0.0;
    double theMaxDist;
    size_t theMaxCoord = -1;
    CoordinateType currentCoor;
    for(MultiArrayIndex k=0; k<tg.shape(2); k+=2) {
        TIC
        std::cout << "  " << k << "/" << tg.shape(2) << " ... " << std::flush;
    for(MultiArrayIndex j=0; j<tg.shape(1); j+=2) {
    for(MultiArrayIndex i=0; i<tg.shape(0); i+=2) {
        const SegValueType& label = seg(i/2,j/2,k/2);
        if(label == 0) {
            continue;
        }
        const std::vector<CoordinateType>& surface = surfaces[label];

        currentCoor = CoordinateType(i,j,k);
        theMaxDist  = std::numeric_limits<double>::min();
       
        size_t idx=0;
        for(typename std::vector<CoordinateType>::const_iterator surfaceCoor=surface.begin();
            surfaceCoor!=surface.end(); ++surfaceCoor, ++idx)
        {
            dist = (currentCoor - *surfaceCoor).squaredMagnitude();
            if(dist > theMaxDist) {
                theMaxDist  = dist;
                theMaxCoord = idx;
            }
        }
        maxDist(i/2,j/2,k/2) = (surface[theMaxCoord]-currentCoor)/2.0;
    }
    }
        std::cout << TOCS << std::endl;
    }
}
    
#endif