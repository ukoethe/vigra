#ifndef VIGRA_ILASTIKTOOLS_CARVING_HXX
#define VIGRA_ILASTIKTOOLS_CARVING_HXX

#include "../adjacency_list_graph.hxx"
#include "../timing.hxx"
#include "../multi_gridgraph.hxx"

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
            const Edge e  = this->findEdge(this->nodeFromId(lu), this->nodeFromId(lu));
            return this->id(e);
        }

        void assignLabels(const MultiArrayView<DIM, LABELS> & labels){
            labelView_ = labels;

            LABELS minLabel, maxLabel;
            labelView_.minmax(&minLabel, &maxLabel);

            if(minLabel!=1){
                throw std::runtime_error("Labels need to start at 1 !");
            }


            USETICTOC;


            std::cout<<"add nodes\n";
            TIC;
            this->assignNodeRange(1, maxLabel+1);
            TOC;

     


            
            std::cout<<"add edges\n";
            TIC;
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
            TOC;
        }

        template<class WEIGHTS_IN, class WEIGHTS_OUT>
        void accumulateEdgeFeatures(
            const MultiArrayView<DIM, WEIGHTS_IN> & featuresIn,
            const MultiArrayView<1, typename vigra::NumericTraits<WEIGHTS_IN>::RealPromote > & featuresOut
        ){
            typedef typename vigra::NumericTraits<WEIGHTS_IN>::RealPromote RealType;
            const Shape shape = labelView_.shape();
            MultiArray<1, UInt32>   counting(Shape1(this->edgeNum()));



            // initiaize output with zeros
            featuresOut = RealType(0);

            //do the accumulation
            for(size_t z=0; z<shape[2]; ++z)
            for(size_t y=0; y<shape[1]; ++y)
            for(size_t x=0; x<shape[0]; ++x){
                const LabelType lu  = labelView_(x, y, z);
                if(x+1 < shape[0]){
                    const LabelType lv = labelView_(x+1, y, z);
                    if(lu!=lv){
                        const int eid = findEdgeFromIds(lu, lv);
                        counting[eid]+=2;
                        featuresOut[eid]+=static_cast<RealType>(featuresIn(x,y,z));
                        featuresOut[eid]+=static_cast<RealType>(featuresIn(x+1,y,z));

                    }
                }

                if(y+1 < shape[1]){
                    const LabelType lv = labelView_(x, y+1, z);
                    if(lu!=lv){
                        const int eid = findEdgeFromIds(lu, lv);
                        counting[eid]+=2;
                        featuresOut[eid]+=static_cast<RealType>(featuresIn(x,y,z));
                        featuresOut[eid]+=static_cast<RealType>(featuresIn(x,y+1,z));

                    }
                }
                    
                if(z+1 < shape[2]){
                    const LabelType lv = labelView_(x, y, z+1);
                    if(lu!=lv){
                        const int eid = findEdgeFromIds(lu, lv);
                        counting[eid]+=2;
                        featuresOut[eid]+=static_cast<RealType>(featuresIn(x,y,z));
                        featuresOut[eid]+=static_cast<RealType>(featuresIn(x,y,z+1));
                    }
                }
            }
        }

    private:
        void maybeAddEdge(const LabelType lu, const LabelType lv){
            if(lu != lv){
                this->addEdge( this->nodeFromId(lu),this->nodeFromId(lv));
            }
        }


        vigra::MultiArrayView< DIM, LABELS> labelView_;
	};




}


#endif /*VIGRA_ILASTIKTOOLS_CARVING_HXX*/