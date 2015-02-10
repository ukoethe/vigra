#ifndef VIGRA_ILASTIKTOOLS_CARVING_HXX
#define VIGRA_ILASTIKTOOLS_CARVING_HXX

#include "../adjacency_list_graph.hxx"
#include "../timing.hxx"
#include "../multi_gridgraph.hxx"
#include "../timing.hxx"

#include <omp.h>

namespace vigra{



    template< unsigned int DIM, class LABELS>
    class GridRag
    : public AdjacencyListGraph
    {
    public:
        typedef GridGraph<DIM, boost::undirected_tag>  GridGraphType;
        typedef LABELS LabelType;
        typedef TinyVector<MultiArrayIndex, DIM>  Shape;
        typedef TinyVector<MultiArrayIndex,   1>  Shape1;
        GridRag() : AdjacencyListGraph(){

        }

        int findEdgeFromIds(const LabelType lu, const LabelType lv){
            const Edge e  = this->findEdge(this->nodeFromId(lu), this->nodeFromId(lv));
            return this->id(e);
        }

        void assignLabels(const MultiArrayView<DIM, LABELS> & labels){



            labelView_ = labels;

            LABELS minLabel, maxLabel;
            labelView_.minmax(&minLabel, &maxLabel);

            if(minLabel!=1){
                throw std::runtime_error("Labels need to start at 1 !");
            }

            this->assignNodeRange(1, maxLabel+1);

            const Shape shape = labelView_.shape();

            if(DIM == 2){
                for(size_t y=0; y<shape[1]; ++y)
                for(size_t x=0; x<shape[0]; ++x){
                    const LabelType l  = labelView_(x, y);
                    if(x+1 < shape[0] )
                        maybeAddEdge(l, labelView_(x+1, y));
                    if(y+1 < shape[1])
                        maybeAddEdge(l, labelView_(x, y+1));
                }
            }
            else if(DIM==3){
                for(size_t z=0; z<shape[2]; ++z)
                for(size_t y=0; y<shape[1]; ++y)
                for(size_t x=0; x<shape[0]; ++x){
                    const LabelType l  = labelView_(x, y, z);
                    if(x+1 < shape[0] )
                        maybeAddEdge(l, labelView_(x+1, y, z));
                    if(y+1 < shape[1])
                        maybeAddEdge(l, labelView_(x, y+1, z));
                    if(z+1 < shape[2])
                        maybeAddEdge(l, labelView_(x, y, z+1));
                }
            }
            else{
                throw std::runtime_error("currently only 2D and 3D");
            }
        }

        void assignLabelsFromSerialization(
            const MultiArrayView<DIM, LABELS> & labels,
            const MultiArrayView<1, LABELS> & serialization
        ){
            labelView_ = labels;
            this->deserialize(serialization.begin(), serialization.end());
        }

        template<class WEIGHTS_IN, class WEIGHTS_OUT>
        void accumulateEdgeFeatures(
            const MultiArrayView<DIM, WEIGHTS_IN> & featuresIn,
            MultiArrayView<1, WEIGHTS_OUT > & featuresOut
        ){

 
            const Shape shape = labelView_.shape();
            MultiArray<1, UInt32>   counting(Shape1(this->edgeNum()));



            // initiaize output with zeros
            featuresOut = WEIGHTS_OUT(0);
            
            omp_lock_t * edgeLocks = new omp_lock_t[this->edgeNum()];
            #pragma omp parallel for
            for(size_t i=0; i<this->edgeNum();++i){
                omp_init_lock(&(edgeLocks[i]));
            }


            if(DIM==2){
                //do the accumulation
                for(size_t y=0; y<shape[1]; ++y)
                for(size_t x=0; x<shape[0]; ++x){
                    const LabelType lu  = labelView_(x, y);
                    if(x+1 < shape[0]){
                        const LabelType lv = labelView_(x+1, y);
                        if(lu!=lv){
                            const int eid = findEdgeFromIds(lu, lv);
                            counting[eid]+=2;
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x,y));
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x+1,y));
                        }
                    }

                    if(y+1 < shape[1]){
                        const LabelType lv = labelView_(x, y+1);
                        if(lu!=lv){
                            const int eid = findEdgeFromIds(lu, lv);
                            counting[eid]+=2;
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x,y));
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x,y+1));
                        }
                    }
                }
            }
            else if(DIM==3){


                //do the accumulation
                #pragma omp parallel for
                for(size_t z=0; z<shape[2]; ++z)
                for(size_t y=0; y<shape[1]; ++y)
                for(size_t x=0; x<shape[0]; ++x){
                    const LabelType lu  = labelView_(x, y, z);
                    if(x+1 < shape[0]){
                        const LabelType lv = labelView_(x+1, y, z);
                        if(lu!=lv){
                            const int eid = findEdgeFromIds(lu, lv);
                            omp_set_lock(&(edgeLocks[eid]));
                            counting[eid]+=2;
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x,y,z));
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x+1,y,z));
                            omp_unset_lock(&(edgeLocks[eid]));
                        }
                    }
                    if(y+1 < shape[1]){
                        const LabelType lv = labelView_(x, y+1, z);
                        if(lu!=lv){
                            const int eid = findEdgeFromIds(lu, lv);
                            omp_set_lock(&(edgeLocks[eid]));
                            counting[eid]+=2;
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x,y,z));
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x,y+1,z));
                            omp_unset_lock(&(edgeLocks[eid]));
                        }
                    }
                    if(z+1 < shape[2]){
                        const LabelType lv = labelView_(x, y, z+1);
                        if(lu!=lv){
                            const int eid = findEdgeFromIds(lu, lv);
                            omp_set_lock(&(edgeLocks[eid]));
                            counting[eid]+=2;
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x,y,z));
                            featuresOut[eid]+=static_cast<WEIGHTS_OUT>(featuresIn(x,y,z+1));
                            omp_unset_lock(&(edgeLocks[eid]));
                        }
                    }
                }
            }

            #pragma omp parallel for
            for(size_t i=0; i<this->edgeNum();++i){
                omp_destroy_lock(&(edgeLocks[i]));
            }
            delete[] edgeLocks;


            // normalize
            #pragma omp parallel for
            for(size_t i=0; i<this->edgeNum();++i){
                featuresOut[i]/=counting[i];
            }

        }

        const vigra::MultiArrayView< DIM, LABELS> & labels()const{
            return labelView_;
        }
    private:
        void maybeAddEdge(const LabelType lu, const LabelType lv){
            if(lu != lv){
                this->addEdge( this->nodeFromId(lu),this->nodeFromId(lv));
            }
        }


        vigra::MultiArrayView< DIM, LABELS> labelView_;
    };

    template<class T>
    struct GridSegmentorEdgeMap{

        typedef T Value;
        typedef T & Reference;
        typedef const T & ConstReference;
        GridSegmentorEdgeMap(MultiArrayView<1, T> & values)
        : values_(values){
            
        }
        template<class K>
        Reference operator[](const K key){
            return values_[key.id()];
        }

        template<class K>
        ConstReference operator[](const K key)const{
            return values_[key.id()];
        }

        MultiArrayView<1, T> values_;
    };

    template<class T>
    struct GridSegmentorNodeMap{

        typedef T Value;
        typedef T & Reference;
        typedef const T & ConstReference;
        GridSegmentorNodeMap(MultiArrayView<1, T> & values)
        : values_(values){

        }
        template<class K>
        Reference operator[](const K key){
            return values_[key.id()];
        }

        template<class K>
        ConstReference operator[](const K key)const{
            return values_[key.id()];
        }

        MultiArrayView<1, T> & values_;
    };

    template<  
        unsigned int DIM, 
        class LABELS, 
        class VALUE_TYPE
    >
    class GridSegmentor{

    public:
        typedef GridRag<DIM, LABELS> Graph;
        typedef TinyVector<MultiArrayIndex, 1>   Shape1;
        typedef TinyVector<MultiArrayIndex, DIM> ShapeN;

        typedef MultiArrayView<DIM, LABELS> LabelView;
        typedef typename LabelView::const_iterator  LabelViewIter;
        GridSegmentor()
        :   graph_(),
            edgeWeights_(),
            nodeSeeds_()
        {
        }


        template<class WEIGHTS_IN>
        void preprocessing(const MultiArrayView< DIM, LABELS> & labels,
                      const MultiArrayView< DIM,  WEIGHTS_IN> & weightArray
        )
        {

            USETICTOC;

            std::cout<<"get RAG\n";
            TIC;
            // get the RAG
            graph_.assignLabels(labels);
            TOC;



            // now we have the RAG  we can 
            // reshape weights and seeds
            // and resultSegmentation
            edgeWeights_.reshape(Shape1(graph_.edgeNum()));
            nodeSeeds_.reshape(Shape1(graph_.maxNodeId()+ 1));
            resultSegmentation_.reshape(Shape1(graph_.maxNodeId()+ 1));

            std::cout<<"get edge weights\n";
            TIC;
            // accumulate the edge weights
            graph_.accumulateEdgeFeatures(weightArray, edgeWeights_);
            TOC;
        }

        void preprocessingFromSerialization(
            const MultiArrayView< DIM, LABELS> & labels,
            const MultiArrayView< 1,  LABELS> & serialization,
            const MultiArrayView< 1,  VALUE_TYPE> & edgeWeights,
            const MultiArrayView< 1,  UInt8> & nodeSeeds,
            const MultiArrayView< 1,  UInt8> & resultSegmentation
        )
        {

            USETICTOC;

            std::cout<<"get RAG from serialization\n";
            TIC;
            // get the RAG
            graph_.assignLabelsFromSerialization(labels, serialization);
            TOC;


            // assign weights and seeds
            // and resultSegmentation
            edgeWeights_ = edgeWeights;
            nodeSeeds_ = nodeSeeds;
            resultSegmentation_ = resultSegmentation;

        }

        void clearSeeds(){
            #pragma omp parallel for
            for(size_t i=0; i<this->nodeNum();++i){
                nodeSeeds_[i]=0;
            }
        }


        void run(float bias, float noBiasBelow){

            GridSegmentorNodeMap<UInt8> nodeSeeds(nodeSeeds_);
            GridSegmentorNodeMap<UInt8> resultSegmentation(resultSegmentation_);
            GridSegmentorEdgeMap<VALUE_TYPE> edgeWeights(edgeWeights_);
            std::cout<<"run with bias "<<bias<<" no noBiasBelow "<<noBiasBelow<<"\n";
            carvingSegmentation(graph_, edgeWeights, nodeSeeds, 1, bias,noBiasBelow, resultSegmentation);
        }

        template<class PIXEL_LABELS>
        void addLabels(
            const MultiArrayView<DIM, PIXEL_LABELS> & brushStroke,
            const TinyVector<MultiArrayIndex, DIM> roiBegin,
            const TinyVector<MultiArrayIndex, DIM> roiEnd,
            const PIXEL_LABELS maxValidLabel
        ){
            const LabelView & labels = graph_.labels();
            const LabelView roiLabels = labels.subarray(roiBegin, roiEnd);

            typedef typename MultiArrayView<DIM, PIXEL_LABELS>::const_iterator  BrushIter;

            if(roiLabels.shape() != brushStroke.shape()){
                std::cout<<"ROI LABELS shape  "<<roiLabels.shape()<<"\n";
                std::cout<<"brushStroke shape "<<brushStroke.shape()<<"\n";
                throw std::runtime_error("wrong shapes");
            }
            LabelViewIter labelIter(roiLabels.begin());
            LabelViewIter labelIterEnd(roiLabels.end());
            BrushIter brushIter(brushStroke.begin());

            for(; labelIter<labelIterEnd; ++labelIter,++brushIter){
                const int brushLabel = int(*brushIter);
                const LABELS nodeId = *labelIter;

                //std::cout<<"brush label "<< int(brushLabel)<<"\n";
                if(brushLabel == 0 ){

                }
                else if(brushLabel == 1 || brushLabel == 2 ){
                    nodeSeeds_[nodeId] = brushLabel;
                }
                else{
                    nodeSeeds_[nodeId] = 0;
                }
            }
        }


        template<class PIXEL_LABELS>
        void getSegmentation(
            const TinyVector<MultiArrayIndex, DIM> roiBegin,
            const TinyVector<MultiArrayIndex, DIM> roiEnd,
            MultiArrayView<DIM, PIXEL_LABELS> & segmentation
        )const{
            const LabelView & labels = graph_.labels();
            const LabelView roiLabels = labels.subarray(roiBegin, roiEnd);

            typedef typename MultiArrayView<DIM, PIXEL_LABELS>::iterator  SegIter;

            LabelViewIter labelIter(roiLabels.begin());
            LabelViewIter labelIterEnd(roiLabels.end());
            SegIter segIter(segmentation.begin());

            for(; labelIter<labelIterEnd; ++labelIter,++segIter){
                const LABELS nodeId = *labelIter;
                *segIter = resultSegmentation_[nodeId];
            }

        }

        void getSuperVoxelSeg(
            MultiArrayView<1, UInt8> & segmentation
        )const{
            std::copy(resultSegmentation_.begin(), resultSegmentation_.end(), segmentation.begin());
        }

        void getSuperVoxelSeeds(
            MultiArrayView<1, UInt8> & seeds
        )const{
            std::copy(nodeSeeds_.begin(), nodeSeeds_.end(), seeds.begin());
        }

        const GridRag<DIM, LABELS> & graph()const{
            return graph_;
        }

        GridRag<DIM, LABELS> & graph(){
            return graph_;
        }

        size_t nodeNum() const{
            return graph_.nodeNum();
        }
        size_t edgeNum()const{
            return graph_.edgeNum();
        }
        size_t maxNodeId()const{
            return graph_.maxNodeId();
        }
        size_t maxEdgeId()const{
            return graph_.maxEdgeId();
        }

        void setSeeds(
            const MultiArray<1 , UInt8> & fgSeeds,
            const MultiArray<1 , UInt8> & bgSeeds
        ){
            nodeSeeds_ = 0;
            for(size_t i=0; i<fgSeeds.size(); ++i){
                nodeSeeds_[fgSeeds[i]] = 2;
            }
            for(size_t i=0; i<fgSeeds.size(); ++i){
                nodeSeeds_[bgSeeds[i]] = 1;
            }
        }

        const MultiArray<1 , VALUE_TYPE> & edgeWeights()const{
            return edgeWeights_;
        }
        const MultiArray<1 , UInt8>      & nodeSeeds()const{
            return nodeSeeds_;
        }
        const MultiArray<1 , UInt8>      & resultSegmentation()const{
            return resultSegmentation_;
        }
    private:
        GridRag<DIM, LABELS> graph_;
        MultiArray<1 , VALUE_TYPE> edgeWeights_;
        MultiArray<1 , UInt8>      nodeSeeds_;
        MultiArray<1 , UInt8>      resultSegmentation_;

    };

}


#endif /*VIGRA_ILASTIKTOOLS_CARVING_HXX*/
