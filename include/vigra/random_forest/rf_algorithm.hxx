/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by Rahul Nair                             */
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
#define VIGRA_RF_ALGORTIHM_HXX

#include <vector>
#include "splices.hxx"
#include <queue>
#include <fstream>
namespace vigra
{
 
namespace rf
{
/** This namespace contains all algorithms developed for feature 
 * selection
 *
 */
namespace algorithms
{

namespace detail
{
    /** create a MultiArray containing only columns supplied between iterators
        b and e
    */
    template<class OrigMultiArray,
             class Iter,
             class DestMultiArray>
    void choose(OrigMultiArray     const & in,
                Iter               const & b,
                Iter               const & e,
                DestMultiArray        & out)
    {
        int columnCount = std::distance(b, e);
        int rowCount     = in.shape(0);
        out.reshape(MultiArrayShape<2>::type(rowCount, columnCount));
        int ii = 0;
        for(Iter iter = b; iter != e; ++iter, ++ii)
        {
            columnVector(out, ii) = columnVector(in, *iter);
        }
    }
}



/** Standard random forest Errorrate callback functor
 *
 * returns the random forest error estimate when invoked. 
 */
class RFErrorCallback
{
    RandomForestOptions options;
    
    public:
    /** Default constructor
     *
     * optionally supply options to the random forest classifier
     * \sa RandomForestOptions
     */
    RFErrorCallback(RandomForestOptions opt = RandomForestOptions())
    : options(opt)
    {}

    /** returns the RF OOB error estimate given features and 
     * labels
     */
    template<class Feature_t, class Response_t>
    double operator() (Feature_t const & features,
                       Response_t const & response)
    {
        RandomForest<>             rf(options);
        visitors::OOB_Error        oob;
        rf.learn(features, 
                 response, 
                 visitors::create_visitor(oob ));
        return oob.oob_breiman;
    }
};


/** Structure to hold Variable Selection results
 */
class VariableSelectionResult
{
    bool initialized;

  public:
    VariableSelectionResult()
    : initialized(false)
    {}

    typedef std::vector<int> FeatureList_t;
    typedef std::vector<double> ErrorList_t;
    typedef FeatureList_t::iterator Pivot_t;

    Pivot_t pivot;

    /** list of features. 
     */
    FeatureList_t selected;
    
    /** vector of size (number of features)
     *
     * the i-th entry encodes the error rate obtained
     * while using features [0 - i](including i) 
     *
     * if the i-th entry is -1 then no error rate was obtained
     * this may happen if more than one feature is added to the
     * selected list in one step of the algorithm.
     *
     * during initialisation error[m+n-1] is always filled
     */
    ErrorList_t errors;
    

    /** errorrate using no features
     */
    double no_features;

    template<class FeatureT, 
             class ResponseT, 
             class Iter,
             class ErrorRateCallBack>
    bool init(FeatureT const & all_features,
              ResponseT const & response,
              Iter b,
              Iter e,
              ErrorRateCallBack errorcallback)
    {
        bool ret_ = init(all_features, response, errorcallback); 
        if(!ret_)
            return false;
        vigra_precondition(std::distance(b, e) == (std::ptrdiff_t)selected.size(),
                           "Number of features in ranking != number of features matrix");
        std::copy(b, e, selected.begin());
        return true;
    }
    
    template<class FeatureT, 
             class ResponseT, 
             class Iter>
    bool init(FeatureT const & all_features,
              ResponseT const & response,
              Iter b,
              Iter e)
    {
        RFErrorCallback ecallback;
        return init(all_features, response, b, e, ecallback);
    }


    template<class FeatureT, 
             class ResponseT>
    bool init(FeatureT const & all_features,
              ResponseT const & response)
    {
        return init(all_features, response, RFErrorCallback());
    }
    /**initialization routine. Will be called only once in the lifetime
     * of a VariableSelectionResult. Subsequent calls will not reinitialize
     * member variables.
     *
     * This is intended, to allow continuing variable selection at a point 
     * stopped in an earlier iteration. 
     *
     * returns true if initialization was successful and false if 
     * the object was already initialized before.
     */
    template<class FeatureT, 
             class ResponseT,
             class ErrorRateCallBack>
    bool init(FeatureT const & all_features,
              ResponseT const & response,
              ErrorRateCallBack errorcallback)
    {
        if(initialized)
        {
            return false;
        }
        // calculate error with all features
        selected.resize(all_features.shape(1), 0);
        for(unsigned int ii = 0; ii < selected.size(); ++ii)
            selected[ii] = ii;
        errors.resize(all_features.shape(1), -1);
        errors.back() = errorcallback(all_features, response);

        // calculate error rate if no features are chosen 
        // corresponds to max(prior probability) of the classes
        std::map<typename ResponseT::value_type, int>     res_map;
        std::vector<int>                                 cts;
        int                                             counter = 0;
        for(int ii = 0; ii < response.shape(0); ++ii)
        {
            if(res_map.find(response(ii, 0)) == res_map.end())
            {
                res_map[response(ii, 0)] = counter;
                ++counter;
                cts.push_back(0);
            }
            cts[res_map[response(ii,0)]] +=1;
        }
        no_features = double(*(std::max_element(cts.begin(),
                                                 cts.end())))
                    / double(response.shape(0));

        /*init not_selected vector;
        not_selected.resize(all_features.shape(1), 0);
        for(int ii = 0; ii < not_selected.size(); ++ii)
        {
            not_selected[ii] = ii;
        }
        initialized = true;
        */
        pivot = selected.begin();
        return true;
    }
};


    
/** Perform forward selection
 *
 * \param features    IN:     n x p matrix containing n instances with p attributes/features
 *                             used in the variable selection algorithm
 * \param response  IN:     n x 1 matrix containing the corresponding response
 * \param result    IN/OUT: VariableSelectionResult struct which will contain the results
 *                             of the algorithm. 
 *                             Features between result.selected.begin() and result.pivot will
 *                             be left untouched.
 *                             \sa VariableSelectionResult
 * \param errorcallback
 *                     IN, OPTIONAL: 
 *                             Functor that returns the error rate given a set of 
 *                             features and labels. Default is the RandomForest OOB Error.
 *
 * Forward selection subsequently chooses the next feature that decreases the Error rate most.
 *
 * usage:
 * \code
 *         MultiArray<2, double>     features = createSomeFeatures();
 *         MultiArray<2, int>        labels   = createCorrespondingLabels();
 *         VariableSelectionResult  result;
 *         forward_selection(features, labels, result);
 * \endcode
 * To use forward selection but ensure that a specific feature e.g. feature 5 is always 
 * included one would do the following
 *
 * \code
 *         VariableSelectionResult result;
 *         result.init(features, labels);
 *         std::swap(result.selected[0], result.selected[5]);
 *         result.setPivot(1);
 *         forward_selection(features, labels, result);
 * \endcode
 *
 * \sa VariableSelectionResult
 *
 */                    
template<class FeatureT, class ResponseT, class ErrorRateCallBack>
void forward_selection(FeatureT          const & features,
                       ResponseT          const & response,
                       VariableSelectionResult & result,
                       ErrorRateCallBack          errorcallback)
{
    VariableSelectionResult::FeatureList_t & selected         = result.selected;
    VariableSelectionResult::ErrorList_t &     errors            = result.errors;
    VariableSelectionResult::Pivot_t       & pivot            = result.pivot;    
    int featureCount = features.shape(1);
    // initialize result struct if in use for the first time
    if(!result.init(features, response, errorcallback))
    {
        //result is being reused just ensure that the number of features is
        //the same.
        vigra_precondition((int)selected.size() == featureCount,
                           "forward_selection(): Number of features in Feature "
                           "matrix and number of features in previously used "
                           "result struct mismatch!");
    }
    

    int not_selected_size = std::distance(pivot, selected.end());
    while(not_selected_size > 1)
    {
        std::vector<int> current_errors;
        VariableSelectionResult::Pivot_t next = pivot;
        for(int ii = 0; ii < not_selected_size; ++ii, ++next)
        {
            std::swap(*pivot, *next);
            MultiArray<2, double> cur_feats;
            detail::choose( features, 
                            selected.begin(), 
                            pivot+1, 
                            cur_feats);
            double error = errorcallback(cur_feats, response);
            current_errors.push_back(error);
            std::swap(*pivot, *next);
        }
        int pos = std::distance(current_errors.begin(),
                                std::min_element(current_errors.begin(),
                                                   current_errors.end()));
        next = pivot;
        std::advance(next, pos);
        std::swap(*pivot, *next);
        errors[std::distance(selected.begin(), pivot)] = current_errors[pos];
        ++pivot;
        not_selected_size = std::distance(pivot, selected.end());
    }
}
template<class FeatureT, class ResponseT>
void forward_selection(FeatureT          const & features,
                       ResponseT          const & response,
                       VariableSelectionResult & result)
{
    forward_selection(features, response, result, RFErrorCallback());
}


/** Perform backward elimination
 *
 * \param features    IN:     n x p matrix containing n instances with p attributes/features
 *                             used in the variable selection algorithm
 * \param response  IN:     n x 1 matrix containing the corresponding response
 * \param result    IN/OUT: VariableSelectionResult struct which will contain the results
 *                             of the algorithm. 
 *                             Features between result.pivot and result.selected.end() will
 *                             be left untouched.
 *                             \sa VariableSelectionResult
 * \param errorcallback
 *                     IN, OPTIONAL: 
 *                             Functor that returns the error rate given a set of 
 *                             features and labels. Default is the RandomForest OOB Error.
 *
 * Backward elimination subsequently eliminates features that have the least influence
 * on the error rate
 *
 * usage:
 * \code
 *         MultiArray<2, double>     features = createSomeFeatures();
 *         MultiArray<2, int>        labels   = createCorrespondingLabels();
 *         VariableSelectionResult  result;
 *         backward_elimination(features, labels, result);
 * \endcode
 * To use backward elimination but ensure that a specific feature e.g. feature 5 is always 
 * excluded one would do the following:
 *
 * \code
 *         VariableSelectionResult result;
 *         result.init(features, labels);
 *         std::swap(result.selected[result.selected.size()-1], result.selected[5]);
 *         result.setPivot(result.selected.size()-1);
 *         backward_elimination(features, labels, result);
 * \endcode
 *
 * \sa VariableSelectionResult
 *
 */                    
template<class FeatureT, class ResponseT, class ErrorRateCallBack>
void backward_elimination(FeatureT              const & features,
                             ResponseT         const & response,
                          VariableSelectionResult & result,
                          ErrorRateCallBack         errorcallback)
{
    int featureCount = features.shape(1);
    VariableSelectionResult::FeatureList_t & selected         = result.selected;
    VariableSelectionResult::ErrorList_t &     errors            = result.errors;
    VariableSelectionResult::Pivot_t       & pivot            = result.pivot;    
    
    // initialize result struct if in use for the first time
    if(!result.init(features, response, errorcallback))
    {
        //result is being reused just ensure that the number of features is
        //the same.
        vigra_precondition((int)selected.size() == featureCount,
                           "backward_elimination(): Number of features in Feature "
                           "matrix and number of features in previously used "
                           "result struct mismatch!");
    }
    pivot = selected.end() - 1;    

    int selected_size = std::distance(selected.begin(), pivot);
    while(selected_size > 1)
    {
        VariableSelectionResult::Pivot_t next = selected.begin();
        std::vector<int> current_errors;
        for(int ii = 0; ii < selected_size; ++ii, ++next)
        {
            std::swap(*pivot, *next);
            MultiArray<2, double> cur_feats;
            detail::choose( features, 
                            selected.begin(), 
                            pivot, 
                            cur_feats);
            double error = errorcallback(cur_feats, response);
            current_errors.push_back(error);
            std::swap(*pivot, *next);
        }
        int pos = std::distance(current_errors.begin(),
                                std::max_element(current_errors.begin(),
                                                   current_errors.end()));
        next = selected.begin();
        std::advance(next, pos);
        std::swap(*pivot, *next);
//        std::cerr << std::distance(selected.begin(), pivot) << " " << pos << " " << current_errors.size() << " " << errors.size() << std::endl;
        errors[std::distance(selected.begin(), pivot)] = current_errors[pos];
        selected_size = std::distance(selected.begin(), pivot);
        --pivot;
    }
}

template<class FeatureT, class ResponseT>
void backward_elimination(FeatureT              const & features,
                             ResponseT         const & response,
                          VariableSelectionResult & result)
{
    backward_elimination(features, response, result, RFErrorCallback());
}

/** Perform rank selection using a predefined ranking
 *
 * \param features    IN:     n x p matrix containing n instances with p attributes/features
 *                             used in the variable selection algorithm
 * \param response  IN:     n x 1 matrix containing the corresponding response
 * \param result    IN/OUT: VariableSelectionResult struct which will contain the results
 *                             of the algorithm. The struct should be initialized with the
 *                             predefined ranking.
 *                         
 *                             \sa VariableSelectionResult
 * \param errorcallback
 *                     IN, OPTIONAL: 
 *                             Functor that returns the error rate given a set of 
 *                             features and labels. Default is the RandomForest OOB Error.
 *
 * Often some variable importance, score measure is used to create the ordering in which
 * variables have to be selected. This method takes such a ranking and calculates the 
 * corresponding error rates. 
 *
 * usage:
 * \code
 *         MultiArray<2, double>     features = createSomeFeatures();
 *         MultiArray<2, int>        labels   = createCorrespondingLabels();
 *         std::vector<int>        ranking  = createRanking(features);
 *         VariableSelectionResult  result;
 *         result.init(features, labels, ranking.begin(), ranking.end());
 *         backward_elimination(features, labels, result);
 * \endcode
 *
 * \sa VariableSelectionResult
 *
 */                    
template<class FeatureT, class ResponseT, class ErrorRateCallBack>
void rank_selection      (FeatureT              const & features,
                             ResponseT         const & response,
                          VariableSelectionResult & result,
                          ErrorRateCallBack         errorcallback)
{
    VariableSelectionResult::FeatureList_t & selected         = result.selected;
    VariableSelectionResult::ErrorList_t &     errors            = result.errors;
    VariableSelectionResult::Pivot_t       & iter            = result.pivot;
    int featureCount = features.shape(1);
    // initialize result struct if in use for the first time
    if(!result.init(features, response, errorcallback))
    {
        //result is being reused just ensure that the number of features is
        //the same.
        vigra_precondition((int)selected.size() == featureCount,
                           "forward_selection(): Number of features in Feature "
                           "matrix and number of features in previously used "
                           "result struct mismatch!");
    }
    
	int ii = 0;
    for(; iter != selected.end(); ++iter)
    {
//        std::cerr << ii<< std::endl;
        ++ii;
        MultiArray<2, double> cur_feats;
        detail::choose( features, 
                        selected.begin(), 
                        iter, 
                        cur_feats);
        double error = errorcallback(cur_feats, response);
        errors[std::distance(selected.begin(), iter)] = error;

    }
}

template<class FeatureT, class ResponseT>
void rank_selection      (FeatureT              const & features,
                             ResponseT         const & response,
                          VariableSelectionResult & result)
{
    rank_selection(features, response, result, RFErrorCallback());
}



enum ClusterLeafTypes{c_Leaf = 95, c_Node = 99};

/* View of a Node in the hierarchical clustering 
 * class 
 * For internal use only - 
 * \sa NodeBase
 */
class ClusterNode
: public NodeBase
{
    public:

    typedef NodeBase BT;

        /**constructors **/
    ClusterNode():NodeBase(){}
    ClusterNode(    int                      nCol,
                    BT::T_Container_type    &   topology,
                    BT::P_Container_type    &   split_param)
                :   BT(nCol + 5, 5,topology, split_param)
    {
        status() = 0; 
        BT::column_data()[0] = nCol;
        if(nCol == 1)
            BT::typeID() = c_Leaf;
        else
            BT::typeID() = c_Node;
    }

    ClusterNode(           BT::T_Container_type  const  &   topology,
                    BT::P_Container_type  const  &   split_param,
                    int                  n             )
                :   NodeBase(5 , 5,topology, split_param, n)
    {
        //TODO : is there a more elegant way to do this?
        BT::topology_size_ += BT::column_data()[0];
    }

    ClusterNode( BT & node_)
        :   BT(5, 5, node_) 
    {
        //TODO : is there a more elegant way to do this?
        BT::topology_size_ += BT::column_data()[0];
        BT::parameter_size_ += 0;
    }
    int index()
    {
        return static_cast<int>(BT::parameters_begin()[1]);
    }
    void set_index(int in)
    {
        BT::parameters_begin()[1] = in;
    }
    double& mean()
    {
        return BT::parameters_begin()[2];
    }
    double& stdev()
    {
        return BT::parameters_begin()[3];
    }
    double& status()
    {
        return BT::parameters_begin()[4];
    }
};

/** Stackentry class for HClustering class
 */
struct HC_Entry
{
    int parent;
    int level;
    int addr; 
    bool infm;
    HC_Entry(int p, int l, int a, bool in)
        : parent(p), level(l), addr(a), infm(in)
    {}
};


/** Hierarchical Clustering class. 
 * Performs single linkage clustering
 * \code
 *         Matrix<double> distance = get_distance_matrix();
 *      linkage.cluster(distance);
 *      // Draw clustering tree.
 *      Draw<double, int> draw(features, labels, "linkagetree.graph");
 *      linkage.breadth_first_traversal(draw);
 * \endcode
 * \sa ClusterImportanceVisitor
 *
 * once the clustering has taken place. Information queries can be made
 * using the breadth_first_traversal() method and iterate() method
 *
 */
class HClustering
{
public:
    typedef MultiArrayShape<2>::type Shp;
    ArrayVector<int>         topology_;
    ArrayVector<double>     parameters_;
    int                     begin_addr;

    // Calculates the distance between two 
    double dist_func(double a, double b)
    {
        return std::min(a, b); 
    }

    /** Visit each node with a Functor 
     * in creation order (should be depth first)
     */
    template<class Functor>
    void iterate(Functor & tester)
    {

        std::vector<int> stack; 
        stack.push_back(begin_addr); 
        while(!stack.empty())
        {
            ClusterNode node(topology_, parameters_, stack.back());
            stack.pop_back();
            if(!tester(node))
            {
                if(node.columns_size() != 1)
                {
                    stack.push_back(node.child(0));
                    stack.push_back(node.child(1));
                }
            }
        }
    }

    /** Perform breadth first traversal of hierarchical cluster tree
     */
    template<class Functor>
    void breadth_first_traversal(Functor & tester)
    {

        std::queue<HC_Entry> queue; 
        int level = 0;
        int parent = -1;
        int addr   = -1;
        bool infm  = false;
        queue.push(HC_Entry(parent,level,begin_addr, infm)); 
        while(!queue.empty())
        {
            level  = queue.front().level;
            parent = queue.front().parent;
            addr   = queue.front().addr;
            infm   = queue.front().infm;
            ClusterNode node(topology_, parameters_, queue.front().addr);
            ClusterNode parnt;
            if(parent != -1)
            {
                parnt = ClusterNode(topology_, parameters_, parent); 
            }
            queue.pop();
            bool istrue = tester(node, level, parnt, infm);
            if(node.columns_size() != 1)
            {
                queue.push(HC_Entry(addr, level +1,node.child(0),istrue));
                queue.push(HC_Entry(addr, level +1,node.child(1),istrue));
            }
        }
    }
    /**save to HDF5 - defunct - has to be updated to new HDF5 interface
     */
#if HasHDF5        
    void save(std::string file, std::string prefix)
    {

        vigra::writeHDF5(file.c_str(), (prefix + "topology").c_str(), 
                               MultiArrayView<2, int>(
                                    Shp(topology_.size(),1),
                                    topology_.data()));
        vigra::writeHDF5(file.c_str(), (prefix + "parameters").c_str(), 
                               MultiArrayView<2, double>(
                                    Shp(parameters_.size(), 1),
                                    parameters_.data()));
        vigra::writeHDF5(file.c_str(), (prefix + "begin_addr").c_str(), 
                               MultiArrayView<2, int>(Shp(1,1), &begin_addr));

    }
#endif                               

    /**Perform single linkage clustering
     * \param distance distance matrix used. \sa CorrelationVisitor
     */
    template<class T, class C>
    void cluster(MultiArrayView<2, T, C> distance)
    {
        MultiArray<2, T> dist(distance); 
        std::vector<std::pair<int, int> > addr; 
        typedef std::pair<int, int>  Entry;
        int index = 0;
        for(int ii = 0; ii < distance.shape(0); ++ii)
        {
            addr.push_back(std::make_pair(topology_.size(), ii));
            ClusterNode leaf(1, topology_, parameters_);
            leaf.set_index(index);
            ++index;
            leaf.columns_begin()[0] = ii;
        }

        while(addr.size() != 1)
        {
            //find the two nodes with the smallest distance
            int ii_min = 0;
            int jj_min = 1;
            double min_dist = dist((addr.begin()+ii_min)->second, 
                              (addr.begin()+jj_min)->second);
            for(unsigned int ii = 0; ii < addr.size(); ++ii)
            {
                for(unsigned int jj = ii+1; jj < addr.size(); ++jj)
                {
                    if(  dist((addr.begin()+ii_min)->second, 
                              (addr.begin()+jj_min)->second)
                       > dist((addr.begin()+ii)->second, 
                              (addr.begin()+jj)->second))
                    {
                        min_dist = dist((addr.begin()+ii)->second, 
                              (addr.begin()+jj)->second);
                        ii_min = ii; 
                        jj_min = jj;
                    }
                }
            }

            //merge two nodes
            int col_size = 0;
            // The problem is that creating a new node invalidates the iterators stored
            // in firstChild and secondChild.
            {
                ClusterNode firstChild(topology_, 
                                       parameters_, 
                                       (addr.begin() +ii_min)->first);
                ClusterNode secondChild(topology_, 
                                       parameters_, 
                                       (addr.begin() +jj_min)->first);
                col_size = firstChild.columns_size() + secondChild.columns_size();
            }
            int cur_addr = topology_.size();
            begin_addr = cur_addr;
//            std::cerr << col_size << std::endl;
            ClusterNode parent(col_size,
                               topology_,
                               parameters_); 
            ClusterNode firstChild(topology_, 
                                   parameters_, 
                                   (addr.begin() +ii_min)->first);
            ClusterNode secondChild(topology_, 
                                   parameters_, 
                                   (addr.begin() +jj_min)->first);
            parent.parameters_begin()[0] = min_dist;
            parent.set_index(index);
            ++index;
            std::merge(firstChild.columns_begin(), firstChild.columns_end(),
                       secondChild.columns_begin(),secondChild.columns_end(),
                       parent.columns_begin());
            //merge nodes in addr
            int to_keep;
            int to_desc;
            int ii_keep;
            if(*parent.columns_begin() ==  *firstChild.columns_begin())
            {
                parent.child(0) = (addr.begin()+ii_min)->first;
                parent.child(1) = (addr.begin()+jj_min)->first;
                (addr.begin()+ii_min)->first = cur_addr;
                ii_keep = ii_min;
                to_keep = (addr.begin()+ii_min)->second;
                to_desc = (addr.begin()+jj_min)->second;
                addr.erase(addr.begin()+jj_min);
            }
            else
            {
                parent.child(1) = (addr.begin()+ii_min)->first;
                parent.child(0) = (addr.begin()+jj_min)->first;
                (addr.begin()+jj_min)->first = cur_addr;
                ii_keep = jj_min;
                to_keep = (addr.begin()+jj_min)->second;
                to_desc = (addr.begin()+ii_min)->second;
                addr.erase(addr.begin()+ii_min);
            }
            //update distances;
            
            for(int jj = 0 ; jj < (int)addr.size(); ++jj)
            {
                if(jj == ii_keep)
                    continue;
                double bla = dist_func(
                                  dist(to_desc, (addr.begin()+jj)->second),
                                  dist((addr.begin()+ii_keep)->second,
                                        (addr.begin()+jj)->second));

                dist((addr.begin()+ii_keep)->second,
                     (addr.begin()+jj)->second) = bla;
                dist((addr.begin()+jj)->second,
                     (addr.begin()+ii_keep)->second) = bla;
            }
        }
    }

};


/** Normalize the status value in the HClustering tree (HClustering Visitor)
 */
class NormalizeStatus
{
public:
    double n;
    /** Constructor
     * \param m normalize status() by m
     */
    NormalizeStatus(double m)
        :n(m)
    {}
    template<class Node>
    bool operator()(Node& node)
    {
        node.status()/=n;
        return false;
    }
};


/** Perform Permutation importance on HClustering clusters
 * (See visit_after_tree() method of visitors::VariableImportance to 
 * see the basic idea. (Just that we apply the permutation not only to
 * variables but also to clusters))
 */
template<class Iter, class DT>
class PermuteCluster
{
public:
    typedef MultiArrayShape<2>::type Shp;
    Matrix<double> tmp_mem_;
    MultiArrayView<2, double> perm_imp;
    MultiArrayView<2, double> orig_imp;
    Matrix<double> feats_;
    Matrix<int>    labels_;
    const int      nPerm;
    DT const &           dt;
    int index;
    int oob_size;

    template<class Feat_T, class Label_T>
    PermuteCluster(Iter  a, 
                   Iter  b,
                   Feat_T const & feats,
                   Label_T const & labls, 
                   MultiArrayView<2, double> p_imp, 
                   MultiArrayView<2, double> o_imp, 
                   int np,
                   DT const  & dt_)
        :tmp_mem_(_spl(a, b).size(), feats.shape(1)),
         perm_imp(p_imp),
         orig_imp(o_imp),
         feats_(_spl(a,b).size(), feats.shape(1)),
         labels_(_spl(a,b).size(),1),
         nPerm(np),
         dt(dt_),
         index(0),
         oob_size(b-a)
    {
        copy_splice(_spl(a,b),
                    _spl(feats.shape(1)),
                    feats,
                    feats_);
        copy_splice(_spl(a,b),
                    _spl(labls.shape(1)),
                    labls,
                    labels_);
    }

    template<class Node>
    bool operator()(Node& node)
    {
        tmp_mem_ = feats_;
        RandomMT19937 random;
        int class_count = perm_imp.shape(1) - 1;
        //permute columns together
        for(int kk = 0; kk < nPerm; ++kk)
        {
            tmp_mem_ = feats_;
            for(int ii = 0; ii < rowCount(feats_); ++ii)
            {
                int index = random.uniformInt(rowCount(feats_) - ii) +ii;
                for(int jj = 0; jj < node.columns_size(); ++jj)
                {
                    if(node.columns_begin()[jj] != feats_.shape(1))
                        tmp_mem_(ii, node.columns_begin()[jj]) 
                            = tmp_mem_(index, node.columns_begin()[jj]);
                }
            }
            
            for(int ii = 0; ii < rowCount(tmp_mem_); ++ii)
            {
                if(dt
                        .predictLabel(rowVector(tmp_mem_, ii)) 
                    ==  labels_(ii, 0))
                {
                    //per class
                    ++perm_imp(index,labels_(ii, 0));
                    //total
                    ++perm_imp(index, class_count);
                }
            }
        }
        double node_status  = perm_imp(index, class_count);
        node_status /= nPerm;
        node_status -= orig_imp(0, class_count);
        node_status *= -1;
        node_status /= oob_size;
        node.status() += node_status;
        ++index;
         
        return false;
    }
};

/** Convert ClusteringTree into a list (HClustering visitor)
 */
class GetClusterVariables
{
public:
    /** NumberOfClusters x NumberOfVariables MultiArrayView containing
     * in each row the variable belonging to a cluster
     */
    MultiArrayView<2, int>    variables;
    int index;
    GetClusterVariables(MultiArrayView<2, int> vars)
        :variables(vars), index(0)
    {}
#if HasHDF5
    void save(std::string file, std::string prefix)
    {
        vigra::writeHDF5(file.c_str(), (prefix + "_variables").c_str(), 
                               variables);
    }
#endif

    template<class Node>
    bool operator()(Node& node)
    {
        for(int ii = 0; ii < node.columns_size(); ++ii)
            variables(index, ii) = node.columns_begin()[ii];
        ++index;
        return false;
    }
};
/** corrects the status fields of a linkage Clustering (HClustering Visitor)
 *  
 *  such that status(currentNode) = min(status(parent), status(currentNode))
 *  \sa cluster_permutation_importance()
 */
class CorrectStatus
{
public:
    template<class Nde>
    bool operator()(Nde & cur, int level, Nde parent, bool infm)
    {
        if(parent.hasData_)
            cur.status() = std::min(parent.status(), cur.status());
        return true;
    }
};


/** draw current linkage Clustering (HClustering Visitor)
 *
 * create a graphviz .dot file
 * usage:
 * \code
 *         Matrix<double> distance = get_distance_matrix();
 *      linkage.cluster(distance);
 *      Draw<double, int> draw(features, labels, "linkagetree.graph");
 *      linkage.breadth_first_traversal(draw);
 * \endcode 
 */
template<class T1,
         class T2, 
         class C1 = UnstridedArrayTag,
         class C2 = UnstridedArrayTag> 
class Draw
{
public:
    typedef MultiArrayShape<2>::type Shp;
    MultiArrayView<2, T1, C1> const &   features_;
    MultiArrayView<2, T2, C2> const &   labels_;
    std::ofstream graphviz;


    Draw(MultiArrayView<2, T1, C1> const & features, 
         MultiArrayView<2, T2, C2> const& labels,
         std::string const  gz)
        :features_(features), labels_(labels), 
        graphviz(gz.c_str(), std::ios::out)
    {
        graphviz << "digraph G\n{\n node [shape=\"record\"]";
    }
    ~Draw()
    {
        graphviz << "\n}\n";
        graphviz.close();
    }

    template<class Nde>
    bool operator()(Nde & cur, int level, Nde parent, bool infm)
    {
        graphviz << "node" << cur.index() << " [style=\"filled\"][label = \" #Feats: "<< cur.columns_size() << "\\n";
        graphviz << " status: " << cur.status() << "\\n";
        for(int kk = 0; kk < cur.columns_size(); ++kk)
        {
                graphviz  << cur.columns_begin()[kk] << " ";
                if(kk % 15 == 14)
                    graphviz << "\\n";
        }
        graphviz << "\"] [color = \"" <<cur.status() << " 1.000 1.000\"];\n";
        if(parent.hasData_)
        graphviz << "\"node" << parent.index() << "\" -> \"node" << cur.index() <<"\";\n";
        return true;
    }
};

/** calculate Cluster based permutation importance while learning. (RandomForestVisitor)
 */
class ClusterImportanceVisitor : public visitors::VisitorBase
{
    public:

    /** List of variables as produced by GetClusterVariables
     */
    MultiArray<2, int>          variables;
    /** Corresponding importance measures
     */
    MultiArray<2, double>       cluster_importance_;
    /** Corresponding error
     */
    MultiArray<2, double>       cluster_stdev_;
    int                         repetition_count_;
    bool                        in_place_;
    HClustering            &    clustering;


#ifdef HasHDF5
    void save(std::string filename, std::string prefix)
    {
        std::string prefix1 = "cluster_importance_" + prefix;
        writeHDF5(filename.c_str(), 
                        prefix1.c_str(), 
                        cluster_importance_);
        prefix1 = "vars_" + prefix;
        writeHDF5(filename.c_str(), 
                        prefix1.c_str(), 
                        variables);
    }
#endif

    ClusterImportanceVisitor(HClustering & clst, int rep_cnt = 10) 
    :   repetition_count_(rep_cnt), clustering(clst)

    {}

    /** Allocate enough memory 
     */
    template<class RF, class PR>
    void visit_at_beginning(RF const & rf, PR const & pr)
    {
        Int32 const  class_count = rf.ext_param_.class_count_;
        Int32 const  column_count = rf.ext_param_.column_count_+1;
        cluster_importance_
            .reshape(MultiArrayShape<2>::type(2*column_count-1, 
                                                class_count+1));
        cluster_stdev_
            .reshape(MultiArrayShape<2>::type(2*column_count-1, 
                                                class_count+1));
        variables
            .reshape(MultiArrayShape<2>::type(2*column_count-1, 
                                                column_count), -1);
        GetClusterVariables gcv(variables);
        clustering.iterate(gcv);
        
    }

    /**compute permutation based var imp. 
     * (Only an Array of size oob_sample_count x 1 is created.
     *  - apposed to oob_sample_count x feature_count in the other method.
     * 
     * \sa FieldProxy
     */
    template<class RF, class PR, class SM, class ST>
    void after_tree_ip_impl(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        typedef MultiArrayShape<2>::type Shp_t;
        Int32                   column_count = rf.ext_param_.column_count_ +1;
        Int32                   class_count  = rf.ext_param_.class_count_;  
        
        // remove the const cast on the features (yep , I know what I am 
        // doing here.) data is not destroyed.
        typename PR::Feature_t & features 
            = const_cast<typename PR::Feature_t &>(pr.features());

        //find the oob indices of current tree. 
        ArrayVector<Int32>      oob_indices;
        ArrayVector<Int32>::iterator
                                iter;
        
        if(rf.ext_param_.actual_msample_ < pr.features().shape(0)- 10000)
        {
            ArrayVector<int> cts(2, 0);
            ArrayVector<Int32> indices(pr.features().shape(0));
            for(int ii = 0; ii < pr.features().shape(0); ++ii)
               indices.push_back(ii); 
            std::random_shuffle(indices.begin(), indices.end());
            for(int ii = 0; ii < rf.ext_param_.row_count_; ++ii)
            {
                if(!sm.is_used()[indices[ii]] && cts[pr.response()(indices[ii], 0)] < 3000)
                {
                    oob_indices.push_back(indices[ii]);
                    ++cts[pr.response()(indices[ii], 0)];
                }
            }
        }
        else
        {
            for(int ii = 0; ii < rf.ext_param_.row_count_; ++ii)
                if(!sm.is_used()[ii])
                    oob_indices.push_back(ii);
        }

        // Random foo
        RandomMT19937           random(RandomSeed);
        UniformIntRandomFunctor<RandomMT19937>  
                                randint(random);

        //make some space for the results
        MultiArray<2, double>
                    oob_right(Shp_t(1, class_count + 1)); 
        
        // get the oob success rate with the original samples
        for(iter = oob_indices.begin(); 
            iter != oob_indices.end(); 
            ++iter)
        {
            if(rf.tree(index)
                    .predictLabel(rowVector(features, *iter)) 
                ==  pr.response()(*iter, 0))
            {
                //per class
                ++oob_right[pr.response()(*iter,0)];
                //total
                ++oob_right[class_count];
            }
        }
        
        MultiArray<2, double>
                    perm_oob_right (Shp_t(2* column_count-1, class_count + 1)); 
        
        PermuteCluster<ArrayVector<Int32>::iterator,typename RF::DecisionTree_t>
            pc(oob_indices.begin(), oob_indices.end(), 
                            pr.features(),
                            pr.response(),
                            perm_oob_right,
                            oob_right,
                            repetition_count_,
                            rf.tree(index));
        clustering.iterate(pc);

        perm_oob_right  /=  repetition_count_;
        for(int ii = 0; ii < rowCount(perm_oob_right); ++ii)
            rowVector(perm_oob_right, ii) -= oob_right;

        perm_oob_right       *= -1;
        perm_oob_right       /= oob_indices.size();
        cluster_importance_  += perm_oob_right;
    }

    /** calculate permutation based impurity after every tree has been 
     * learned  default behaviour is that this happens out of place.
     * If you have very big data sets and want to avoid copying of data 
     * set the in_place_ flag to true. 
     */
    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {    
            after_tree_ip_impl(rf, pr, sm, st, index);
    }

    /** Normalise variable importance after the number of trees is known.
     */
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
        NormalizeStatus nrm(rf.tree_count());
        clustering.iterate(nrm);
        cluster_importance_ /= rf.trees_.size();
    }
};

/** Perform hierarchical clustering of variables and assess importance of clusters
 *
 * \param features    IN:     n x p matrix containing n instances with p attributes/features
 *                             used in the variable selection algorithm
 * \param response  IN:     n x 1 matrix containing the corresponding response
 * \param linkage    OUT:    Hierarchical grouping of variables.
 * \param distance  OUT:    distance matrix used for creating the linkage
 *
 * Performs Hierarchical clustering of variables. And calculates the permutation importance 
 * measures of each of the clusters. Use the Draw functor to create human readable output
 * The cluster-permutation importance measure corresponds to the normal permutation importance
 * measure with all columns corresponding to a cluster permuted. 
 * The importance measure for each cluster is stored as the status() field of each clusternode
 * \sa HClustering
 *
 * usage:
 * \code
 *         MultiArray<2, double>     features = createSomeFeatures();
 *         MultiArray<2, int>        labels   = createCorrespondingLabels();
 *         HClustering                linkage;
 *         MultiArray<2, double>    distance;
 *         cluster_permutation_importance(features, labels, linkage, distance)
 *        // create graphviz output
 *
 *      Draw<double, int> draw(features, labels, "linkagetree.graph");
 *      linkage.breadth_first_traversal(draw);
 *
 * \endcode
 *
 *
 */                    
template<class FeatureT, class ResponseT>
void cluster_permutation_importance(FeatureT              const & features,
                                         ResponseT         const &     response,
                                    HClustering               & linkage,
                                    MultiArray<2, double>      & distance)
{

        RandomForestOptions opt;
        opt.tree_count(100);
        if(features.shape(0) > 40000)
            opt.samples_per_tree(20000).use_stratification(RF_EQUAL);


        vigra::RandomForest<int> RF(opt); 
        visitors::RandomForestProgressVisitor             progress;
        visitors::CorrelationVisitor                     missc;
        RF.learn(features, response,
                 create_visitor(missc, progress));
        distance = missc.distance;
        /*
           missc.save(exp_dir + dset.name() + "_result.h5", dset.name()+"MACH");
           */


        // Produce linkage
        linkage.cluster(distance);
        
        //linkage.save(exp_dir + dset.name() + "_result.h5", "_linkage_CC/");
        vigra::RandomForest<int> RF2(opt); 
        ClusterImportanceVisitor          ci(linkage);
        RF2.learn(features, 
                  response,
                  create_visitor(progress, ci));
        
        
        CorrectStatus cs;
        linkage.breadth_first_traversal(cs);

        //ci.save(exp_dir + dset.name() + "_result.h5", dset.name());
        //Draw<double, int> draw(dset.features(), dset.response(), exp_dir+ dset.name() + ".graph");
        //linkage.breadth_first_traversal(draw);

}

    
template<class FeatureT, class ResponseT>
void cluster_permutation_importance(FeatureT              const & features,
                                         ResponseT         const &     response,
                                    HClustering               & linkage)
{
    MultiArray<2, double> distance;
    cluster_permutation_importance(features, response, linkage, distance);
}
}//namespace algorithms
}//namespace rf
}//namespace vigra
