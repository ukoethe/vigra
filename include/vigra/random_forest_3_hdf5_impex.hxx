/************************************************************************/
/*                                                                      */
/*    Copyright 2009,2014, 2015 by Sven Peter, Philip Schill,           */
/*                                 Rahul Nair and Ullrich Koethe        */
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

#ifndef VIGRA_RF3_IMPEX_HDF5_HXX
#define VIGRA_RF3_IMPEX_HDF5_HXX

#include <string>
#include <sstream>
#include <iomanip>
#include <stack>

#include "config.hxx"
#include "random_forest_3/random_forest.hxx"
#include "random_forest_3/random_forest_common.hxx"
#include "random_forest_3/random_forest_visitors.hxx"
#include "hdf5impex.hxx"

namespace vigra 
{
namespace rf3
{

// needs to be in sync with random_forest_hdf5_impex for backwards compatibility
static const char *const rf_hdf5_ext_param     = "_ext_param";
static const char *const rf_hdf5_options       = "_options";
static const char *const rf_hdf5_topology      = "topology";
static const char *const rf_hdf5_parameters    = "parameters";
static const char *const rf_hdf5_tree          = "Tree_";
static const char *const rf_hdf5_version_group = ".";
static const char *const rf_hdf5_version_tag   = "vigra_random_forest_version";
static const double      rf_hdf5_version       =  0.1;

// keep in sync with include/vigra/random_forest/rf_nodeproxy.hxx
enum NodeTags
{
    rf_UnFilledNode        = 42,
    rf_AllColumns          = 0x00000000,
    rf_ToBePrunedTag       = 0x80000000,
    rf_LeafNodeTag         = 0x40000000,

    rf_i_ThresholdNode     = 0,
    rf_i_HyperplaneNode    = 1,
    rf_i_HypersphereNode   = 2,
    rf_e_ConstProbNode     = 0 | rf_LeafNodeTag,
    rf_e_LogRegProbNode    = 1 | rf_LeafNodeTag
};

static const unsigned int rf_tag_mask = 0xf0000000;
static const unsigned int rf_type_mask = 0x00000003;
static const unsigned int rf_zero_mask = 0xffffffff & ~rf_tag_mask & ~rf_type_mask;

namespace detail
{
    inline std::string get_cwd(HDF5File & h5context)
    {
        return h5context.get_absolute_path(h5context.pwd());
    }
}

template <typename FEATURES, typename LABELS>
typename DefaultRF<FEATURES, LABELS>::type
random_forest_import_HDF5(HDF5File & h5ctx, std::string const & pathname = "")
{
    typedef typename DefaultRF<FEATURES, LABELS>::type RF;
    typedef typename RF::Graph Graph;
    typedef typename RF::Node Node;
    typedef typename RF::SplitTests SplitTest;
    typedef typename LABELS::value_type LabelType;
    typedef typename RF::AccInputType AccInputType;
    typedef typename AccInputType::value_type AccValueType;

    std::string cwd;

    if (pathname.size()) {
        cwd = detail::get_cwd(h5ctx);
        h5ctx.cd(pathname);
    }

    if (h5ctx.existsAttribute(rf_hdf5_version_group, rf_hdf5_version_tag)) {
        double version;
        h5ctx.readAttribute(rf_hdf5_version_group, rf_hdf5_version_tag, version);
        vigra_precondition(version <= rf_hdf5_version, "random_forest_import_HDF5(): unexpected file format version.");
    }

    // Read ext params.
    size_t actual_mtry;
    size_t num_instances;
    size_t num_features;
    size_t num_classes;
    size_t msample;
    int is_weighted_int;
    MultiArray<1, LabelType> distinct_labels_marray;
    MultiArray<1, double> class_weights_marray;

    h5ctx.cd(rf_hdf5_ext_param);
    h5ctx.read("actual_msample_", msample);
    h5ctx.read("actual_mtry_", actual_mtry);
    h5ctx.read("class_count_", num_classes);
    h5ctx.readAndResize("class_weights_", class_weights_marray);
    h5ctx.read("column_count_", num_features);
    h5ctx.read("is_weighted_", is_weighted_int);
    h5ctx.readAndResize("labels", distinct_labels_marray);
    h5ctx.read("row_count_", num_instances);
    h5ctx.cd_up();

    bool is_weighted = is_weighted_int == 1 ? true : false;

    // Read options.
    size_t min_num_instances;
    int mtry;
    int mtry_switch_int;
    int bootstrap_sampling_int;
    int tree_count;
    h5ctx.cd(rf_hdf5_options);
    h5ctx.read("min_split_node_size_", min_num_instances);
    h5ctx.read("mtry_", mtry);
    h5ctx.read("mtry_switch_", mtry_switch_int);
    h5ctx.read("sample_with_replacement_", bootstrap_sampling_int);
    h5ctx.read("tree_count_", tree_count);
    h5ctx.cd_up();

    RandomForestOptionTags mtry_switch = (RandomForestOptionTags)mtry_switch_int;
    bool bootstrap_sampling = bootstrap_sampling_int == 1 ? true : false;

    std::vector<LabelType> const distinct_labels(distinct_labels_marray.begin(), distinct_labels_marray.end());
    std::vector<double> const class_weights(class_weights_marray.begin(), class_weights_marray.end());

    auto const pspec = ProblemSpec<LabelType>()
                               .num_features(num_features)
                               .num_instances(num_instances)
                               .num_classes(num_classes)
                               .distinct_classes(distinct_labels)
                               .actual_mtry(actual_mtry)
                               .actual_msample(msample);

    auto options = RandomForestOptions()
                            .min_num_instances(min_num_instances)
                            .bootstrap_sampling(bootstrap_sampling)
                            .tree_count(tree_count);
    options.features_per_node_switch_ = mtry_switch;
    options.features_per_node_ = mtry;
    if (is_weighted)
        options.class_weights(class_weights);

    Graph gr;
    typename RF::template NodeMap<SplitTest>::type split_tests;
    typename RF::template NodeMap<AccInputType>::type leaf_responses;

    auto const groups = h5ctx.ls();
    for (auto const & groupname : groups) {
        if (groupname.substr(0, std::char_traits<char>::length(rf_hdf5_tree)).compare(rf_hdf5_tree) != 0) {
            continue;
        }

        MultiArray<1, unsigned int> topology;
        MultiArray<1, double> parameters;
        h5ctx.cd(groupname);
        h5ctx.readAndResize(rf_hdf5_topology, topology);
        h5ctx.readAndResize(rf_hdf5_parameters, parameters);
        h5ctx.cd_up();

        vigra_precondition(topology[0] == num_features, "random_forest_import_HDF5(): number of features mismatch.");
        vigra_precondition(topology[1] == num_classes, "random_forest_import_HDF5(): number of classes mismatch.");

        Node const n = gr.addNode();

        std::queue<std::pair<unsigned int, Node> > q;
        q.emplace(2, n);
        while (!q.empty()) {
            auto const el = q.front();

            unsigned int const index = el.first;
            Node const parent = el.second;

            vigra_precondition((topology[index] & rf_zero_mask) == 0, "random_forest_import_HDF5(): unexpected node type: type & zero_mask > 0");

            if (topology[index] & rf_LeafNodeTag) {
                unsigned int const probs_start = topology[index+1] + 1;

                vigra_precondition((topology[index] & rf_tag_mask) == rf_LeafNodeTag, "random_forest_import_HDF5(): unexpected node type: additional tags in leaf node");

                std::vector<AccValueType> node_response;

                for (unsigned int i = 0; i < num_classes; ++i) {
                    node_response.push_back(parameters[probs_start + i]);
                }

                leaf_responses.insert(parent, node_response);

            } else {
                vigra_precondition(topology[index] == rf_i_ThresholdNode, "random_forest_import_HDF5(): unexpected node type.");

                Node const left = gr.addNode();
                Node const right = gr.addNode();

                gr.addArc(parent, left);
                gr.addArc(parent, right);

                split_tests.insert(parent, SplitTest(topology[index+4], parameters[topology[index+1]+1]));

                q.push(std::make_pair(topology[index+2], left));
                q.push(std::make_pair(topology[index+3], right));
            }

            q.pop();
        }
    }

    if (cwd.size()) {
        h5ctx.cd(cwd);
    }

    RF rf(gr, split_tests, leaf_responses, pspec);
    rf.options_ = options;
    return rf;
}

namespace detail
{
    class PaddedNumberString
    {
    public:

        PaddedNumberString(int n)
        {
            ss_ << (n-1);
            width_ = ss_.str().size();
        }

        std::string operator()(int k) const
        {
            ss_.str("");
            ss_ << std::setw(width_) << std::setfill('0') << k;
            return ss_.str();
        }

    private:

        mutable std::ostringstream ss_;
        unsigned int width_;
    };
}

template <typename RF>
void random_forest_export_HDF5(
        RF const & rf,
        HDF5File & h5context,
        std::string const & pathname = ""
){
    typedef typename RF::LabelType LabelType;
    typedef typename RF::Node Node;

    std::string cwd;
    if (pathname.size()) {
        cwd = detail::get_cwd(h5context);
        h5context.cd_mk(pathname);
    }

    // version attribute
    h5context.writeAttribute(rf_hdf5_version_group, rf_hdf5_version_tag,
                             rf_hdf5_version);


    auto const & p = rf.problem_spec_;
    auto const & opts = rf.options_;
    MultiArray<1, LabelType> distinct_classes(Shape1(p.distinct_classes_.size()), p.distinct_classes_.data());
    MultiArray<1, double> class_weights(Shape1(p.num_classes_), 1.0);
    int is_weighted = 0;
    if (opts.class_weights_.size() > 0)
    {
        is_weighted = 1;
        for (size_t i = 0; i < opts.class_weights_.size(); ++i)
            class_weights(i) = opts.class_weights_[i];
    }

    // Save external parameters.
    h5context.cd_mk(rf_hdf5_ext_param);
    h5context.write("column_count_", p.num_features_);
    h5context.write("row_count_", p.num_instances_);
    h5context.write("class_count_", p.num_classes_);
    h5context.write("actual_mtry_", p.actual_mtry_);
    h5context.write("actual_msample_", p.actual_msample_);
    h5context.write("labels", distinct_classes);
    h5context.write("is_weighted_", is_weighted);
    h5context.write("class_weights_", class_weights);
    h5context.write("precision_", 0.0);
    h5context.write("problem_type_", 1.0);
    h5context.write("response_size_", 1.0);
    h5context.write("used_", 1.0);
    h5context.cd_up();

    // Save the options.
    h5context.cd_mk(rf_hdf5_options);
    h5context.write("min_split_node_size_", opts.min_num_instances_);
    h5context.write("mtry_", opts.features_per_node_);
    h5context.write("mtry_func_", 0.0);
    h5context.write("mtry_switch_", opts.features_per_node_switch_);
    h5context.write("predict_weighted_", 0.0);
    h5context.write("prepare_online_learning_", 0.0);
    h5context.write("sample_with_replacement_", opts.bootstrap_sampling_ ? 1 : 0);
    h5context.write("stratification_method_", 3.0);
    h5context.write("training_set_calc_switch_", 1.0);
    h5context.write("training_set_func_", 0.0);
    h5context.write("training_set_proportion_", 1.0);
    h5context.write("training_set_size_", 0.0);
    h5context.write("tree_count_", opts.tree_count_);
    h5context.cd_up();

    // Save the trees.
    detail::PaddedNumberString tree_number(rf.num_trees());
    for (size_t i = 0; i < rf.num_trees(); ++i)
    {
        // Create the topology and parameters arrays.
        std::vector<UInt32> topology;
        std::vector<double> parameters;
        topology.push_back(p.num_features_);
        topology.push_back(p.num_classes_);

        auto const & probs = rf.node_responses_;
        auto const & splits = rf.split_tests_;
        auto const & gr = rf.graph_;
        auto const root = gr.getRoot(i);

        // Write the tree nodes using a depth-first search.
        // When a node is created, the indices of the child nodes are unknown.
        // Therefore, they have to be updated once the child nodes are created.
        // The stack holds the node and the topology-index that must be updated.
        std::stack<std::pair<Node, std::ptrdiff_t> > stack;
        stack.emplace(root, -1);
        while (!stack.empty())
        {
            auto const n = stack.top().first; // the node descriptor
            auto const i = stack.top().second; // index from the parent node that must be updated
            stack.pop();

            // Update the index in the parent node.
            if (i != -1)
                topology[i] = topology.size();

            if (gr.numChildren(n) == 0)
            {
                // The node is a leaf.
                // Topology: leaf node tag, index of weight in parameters array.
                // Parameters: node weight, class probabilities.
                topology.push_back(rf_LeafNodeTag);
                topology.push_back(parameters.size());
                auto const & prob = probs.at(n);
                auto const weight = std::accumulate(prob.begin(), prob.end(), 0.0);
                parameters.push_back(weight);
                parameters.insert(parameters.end(), prob.begin(), prob.end());
            }
            else
            {
                // The node is an inner node.
                // Topology: threshold tag, index of weight in parameters array, index of left child, index of right child, split dimension.
                // Parameters: node weight, split value.
                topology.push_back(rf_i_ThresholdNode);
                topology.push_back(parameters.size());
                topology.push_back(-1); // index of left children (currently unknown, will be updated when the child node is taken from the stack)
                topology.push_back(-1); // index of right children (see above)
                topology.push_back(splits.at(n).dim_);
                parameters.push_back(1.0); // inner nodes have the weight 1.
                parameters.push_back(splits.at(n).val_);
                
                // Place the children on the stack.
                stack.emplace(gr.getChild(n, 0), topology.size()-3);
                stack.emplace(gr.getChild(n, 1), topology.size()-2);
            }
        }

        // Convert the vectors to multi arrays.
        MultiArray<1, UInt32> topo(Shape1(topology.size()), topology.data());
        MultiArray<1, double> para(Shape1(parameters.size()), parameters.data());

        auto const name = rf_hdf5_tree + tree_number(i);
        h5context.cd_mk(name);
        h5context.write(rf_hdf5_topology, topo);
        h5context.write(rf_hdf5_parameters, para);
        h5context.cd_up();
    }

    if (pathname.size())
        h5context.cd(cwd);
}



} // namespace rf3
} // namespace vigra

#endif // VIGRA_NEW_RANDOM_FOREST_IMPEX_HDF5_HXX
