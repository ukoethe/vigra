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
        GridRag() : AdjacencyListGraph(){

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

     

            if(false){
                GridGraphType gridGraph(labelView_.shape());
                for(typename  GridGraphType::EdgeIt e(gridGraph); e!=lemon::INVALID; ++e){
                    const typename GridGraphType::Edge edge(*e);
                    const LabelType lu = labelView_[gridGraph.u(edge)];
                    const LabelType lv = labelView_[gridGraph.v(edge)];
                    if(  lu!=lv   ){
                        this->addEdge( this->nodeFromId(lu),this->nodeFromId(lv));
                    }
                }
            }
            else{


                omp_lock_t  nodeLocks;
                omp_init_lock(&nodeLocks);
                


                std::cout<<"add edges\n";
                TIC;
                const Shape shape = labelView_.shape();

                //#pragma omp parallel for
                for(size_t z=0; z<shape[2]; ++z){
                    for(size_t y=0; y<shape[1]; ++y){
                        for(size_t x=0; x<shape[0]; ++x){

                            const LabelType l  = labelView_(x, y, z);
                            

                            if(x+1 < shape[0] ){
                                const LabelType ol  = labelView_(x+1, y, z);
                                if(l != ol){
                                    //omp_set_lock(&nodeLocks);
                                    this->addEdge( this->nodeFromId(l),this->nodeFromId(ol));
                                    //omp_unset_lock(&nodeLocks);
                                }
                            }
                            if(y+1 < shape[1]){
                                const LabelType ol  = labelView_(x, y+1, z);
                                if(l != ol){
                                    //omp_set_lock(&nodeLocks);
                                    this->addEdge( this->nodeFromId(l),this->nodeFromId(ol));
                                    //omp_unset_lock(&nodeLocks);
                                }
                            }
                            if(z+1 < shape[2]){
                                const LabelType ol  = labelView_(x, y, z+1);
                                if(l != ol){
                                    //omp_set_lock(&nodeLocks);
                                    this->addEdge( this->nodeFromId(l),this->nodeFromId(ol));
                                    //omp_unset_lock(&nodeLocks);
                                }
                            }
                           
                        }
                    }
                }
                TOC;



                //omp_destroy_lock(nodeLocks);
                //TOC;


            }


        }


    private:
        vigra::MultiArrayView< DIM, LABELS> labelView_;
	};

}


#endif /*VIGRA_ILASTIKTOOLS_CARVING_HXX*/