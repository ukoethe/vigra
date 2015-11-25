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
static const unsigned int rf_LeafNodeTag    = 0x40000000;
static const unsigned int rf_ThresholdNode  = 0;



template <typename FEATURES, typename LABELS>
typename DefaultRF<FEATURES, LABELS>::type
random_forest_import_HDF5(HDF5File &h5ctx, const std::string &pathname = "")
{
	typedef typename DefaultRF<FEATURES, LABELS>::type RF;
	typedef typename RF::Graph Graph;
	typedef typename RF::Node Node;
	typedef typename RF::SplitTests SplitTest;
	typedef typename LABELS::value_type LabelsType;
	typedef typename RF::AccInputType::value_type AccValueType;

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
	MultiArray<1, LabelsType> distinct_labels_marray;
	std::vector<LabelsType> distinct_labels;

	h5ctx.read("/_ext_param/actual_mtry_", mtry);
	h5ctx.read("/_ext_param/row_count_", num_instances);
	h5ctx.read("/_ext_param/column_count_", num_features);
	h5ctx.read("/_ext_param/class_count_", num_classes);
	h5ctx.read("/_options/tree_count_", num_trees);
	h5ctx.readAndResize("/_ext_param/labels", distinct_labels_marray);

	for (auto i = distinct_labels_marray.begin(); i != distinct_labels_marray.end(); ++i) {
		distinct_labels.push_back(*i);
	}

	auto const pspec = ProblemSpecNew<typename LABELS::value_type>().num_instances(num_instances).actual_mtry(mtry).num_features(num_features).distinct_classes(distinct_labels);

	Graph gr;
	PropertyMap<Node, SplitTest, RF::container_tag> split_tests;
	PropertyMap<Node, typename RF::AccInputType, RF::container_tag> leaf_responses;

	auto groups = h5ctx.ls();
	for (auto i = groups.begin(); i != groups.end(); ++i) {
		std::string groupname = *i;
		if (groupname.substr(0, std::char_traits<char>::length(rf_hdf5_tree)).compare(rf_hdf5_tree) != 0) {
			continue;
		}

		MultiArray<1, unsigned int> topology;
		MultiArray<1, double> parameters;
		h5ctx.cd(*i);
		h5ctx.readAndResize("topology", topology);
		h5ctx.readAndResize("parameters", parameters);
		h5ctx.cd_up();

		vigra_precondition(topology[0] == num_features, "random_forest_import_HDF5(): number of features mismatch.");
		vigra_precondition(topology[1] == num_classes, "random_forest_import_HDF5(): number of classes mismatch.");

		Node n = gr.addNode();
		unsigned int index = 2;

		typedef std::pair<unsigned int, Node> QueuePair;
		std::queue<QueuePair> q;
		q.push(std::make_pair(2, n));

		while (!q.empty()) {
			QueuePair el = q.front();

			index = el.first;
			Node parent = el.second;

			if (topology[index] & rf_LeafNodeTag) {
				unsigned int probs_start = topology[index+1] + 1;

				std::vector<AccValueType> node_response;

				for (unsigned int i = 0; i < num_classes; ++i) {
					node_response.push_back(parameters[probs_start + i]);
				}

				leaf_responses.insert(parent, node_response);

			} else {
				vigra_precondition(topology[index] == rf_ThresholdNode, "random_forest_import_HDF5(): unexpected node type.");

				Node left = gr.addNode();
				Node right = gr.addNode();

				gr.addArc(parent, left);
				gr.addArc(parent, right);

				split_tests.insert(parent, SplitTest(topology[index+4], parameters[topology[index+1]]));

				q.push(std::make_pair(topology[index+2], left));
				q.push(std::make_pair(topology[index+3], right));
			}

			q.pop();
		}
	}


	RF rf = RF(gr, split_tests, leaf_responses, pspec);
	return rf;
}


} // namespace vigra

#endif // VIGRA_NEW_RANDOM_FOREST_IMPEX_HDF5_HXX
