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

#ifndef VIGRA_NEW_RANDOM_FOREST_IMPEX_HDF5_HXX
#define VIGRA_NEW_RANDOM_FOREST_IMPEX_HDF5_HXX

#include <string>
#include <sstream>
#include <iomanip>

#include "config.hxx"
#include "random_forest_new/random_forest.hxx"
#include "random_forest_new/random_forest_common.hxx"
#include "random_forest_new/random_forest_visitors.hxx"
#include "hdf5impex.hxx"

namespace vigra 
{

// needs to be in sync with random_forest_hdf5_impex for backwards compatibility
static const char *const rf_hdf5_ext_param     = "_ext_param";
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

    size_t mtry;
    size_t num_instances;
    size_t num_features;
    size_t num_classes;
    MultiArray<1, LabelType> distinct_labels_marray;

    h5ctx.cd(rf_hdf5_ext_param);
    h5ctx.read("column_count_", num_features);
    h5ctx.read("row_count_", num_instances);
    h5ctx.read("class_count_", num_classes);
    h5ctx.readAndResize("labels", distinct_labels_marray);
    h5ctx.read("actual_mtry_", mtry);
    h5ctx.cd_up();

    std::vector<LabelType> const distinct_labels(distinct_labels_marray.begin(), distinct_labels_marray.end());

    auto const pspec = ProblemSpecNew<LabelType>()
                               .num_features(num_features)
                               .num_instances(num_instances)
                               .num_classes(num_classes)
                               .distinct_classes(distinct_labels)
                               .actual_mtry(mtry);

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

    return RF(gr, split_tests, leaf_responses, pspec);
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
void rf_export_HDF5(
        RF const & rf,
        HDF5File & h5context,
        std::string const & pathname = ""
){
    std::string cwd;
    if (pathname.size()) {
        cwd = detail::get_cwd(h5context);
        h5context.cd_mk(pathname);
    }

    // version attribute
    h5context.writeAttribute(rf_hdf5_version_group, rf_hdf5_version_tag,
                             rf_hdf5_version);

    // Save external parameters.
    auto const & p = rf.problem_spec_;
    h5context.cd_mk(rf_hdf5_ext_param);
    h5context.write("column_count_", p.num_features_);
    h5context.write("row_count_", p.num_instances_);
    h5context.write("class_count_", p.num_classes_);
    h5context.write("labels", p.distinct_classes_);
    h5context.write("actual_mtry_", p.actual_mtry_);
    h5context.cd_up();

    // Save the trees.
    detail::PaddedNumberString tree_number(rf.num_trees());
    for (size_t i = 0; i < rf.num_trees(); ++i)
    {
        auto const name = rf_hdf5_tree + tree_number(i);
        h5context.cd_mk(name);

        // TODO: write topology
        // TODO: write parameters

        h5context.cd_up();
    }

    if (pathname.size())
        h5context.cd(cwd);
}


} // namespace vigra

#endif // VIGRA_NEW_RANDOM_FOREST_IMPEX_HDF5_HXX
