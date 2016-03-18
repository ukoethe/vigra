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

#include "config.hxx"
#include "random_forest_new/random_forest.hxx"
#include "random_forest_new/random_forest_common.hxx"
#include "random_forest_new/random_forest_visitors.hxx"
#include "hdf5impex.hxx"
#include <string>

namespace vigra 
{

// needs to be in sync with random_forest_hdf5_impex for backwards compatibility
static const char *const rf_hdf5_options       = "_options";
static const char *const rf_hdf5_ext_param     = "_ext_param";
static const char *const rf_hdf5_labels        = "labels";
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
        cwd = h5ctx.get_absolute_path(h5ctx.pwd());
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
    size_t num_trees;
    MultiArray<1, LabelType> distinct_labels_marray;

    h5ctx.read("/_ext_param/actual_mtry_", mtry);
    h5ctx.read("/_ext_param/row_count_", num_instances);
    h5ctx.read("/_ext_param/column_count_", num_features);
    h5ctx.read("/_ext_param/class_count_", num_classes);
    h5ctx.read("/_options/tree_count_", num_trees);
    h5ctx.readAndResize("/_ext_param/labels", distinct_labels_marray);

    std::vector<LabelType> distinct_labels(distinct_labels_marray.begin(), distinct_labels_marray.end());

    auto const pspec = ProblemSpecNew<LabelType>()
                               .num_instances(num_instances)
                               .actual_mtry(mtry)
                               .num_features(num_features)
                               .distinct_classes(distinct_labels);

    Graph gr;
    typename RF::template NodeMap<SplitTest>::type split_tests;
    typename RF::template NodeMap<AccInputType>::type leaf_responses;

    auto groups = h5ctx.ls();
    for (auto const & groupname : groups) {
        if (groupname.substr(0, std::char_traits<char>::length(rf_hdf5_tree)).compare(rf_hdf5_tree) != 0) {
            continue;
        }

        MultiArray<1, unsigned int> topology;
        MultiArray<1, double> parameters;
        h5ctx.cd(groupname);
        h5ctx.readAndResize("topology", topology);
        h5ctx.readAndResize("parameters", parameters);
        h5ctx.cd_up();

        vigra_precondition(topology[0] == num_features, "random_forest_import_HDF5(): number of features mismatch.");
        vigra_precondition(topology[1] == num_classes, "random_forest_import_HDF5(): number of classes mismatch.");

        Node n = gr.addNode();
        unsigned int index = 2;

        std::queue<std::pair<unsigned int, Node> > q;
        q.emplace(2, n);
        while (!q.empty()) {
            auto el = q.front();

            index = el.first;
            Node parent = el.second;

            vigra_precondition((topology[index] & rf_zero_mask) == 0, "random_forest_import_HDF5(): unexpected node type: type & zero_mask > 0");

            if (topology[index] & rf_LeafNodeTag) {
                unsigned int probs_start = topology[index+1] + 1;

                vigra_precondition((topology[index] & rf_tag_mask) == rf_LeafNodeTag, "random_forest_import_HDF5(): unexpected node type: additional tags in leaf node");

                std::vector<AccValueType> node_response;

                for (unsigned int i = 0; i < num_classes; ++i) {
                    node_response.push_back(parameters[probs_start + i]);
                }

                leaf_responses.insert(parent, node_response);

            } else {
                vigra_precondition(topology[index] == rf_i_ThresholdNode, "random_forest_import_HDF5(): unexpected node type.");

                Node left = gr.addNode();
                Node right = gr.addNode();

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


} // namespace vigra

#endif // VIGRA_NEW_RANDOM_FOREST_IMPEX_HDF5_HXX
