/************************************************************************/
/*                                                                      */
/*       Copyright 2009 by Rahul Nair and  Ullrich Koethe               */
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

#ifndef VIGRA_RANDOM_FOREST_IMPEX_HDF5_HXX
#define VIGRA_RANDOM_FOREST_IMPEX_HDF5_HXX

#include "config.hxx"
#include "random_forest.hxx"
#include "hdf5impex.hxx"
#include <string>

namespace vigra 
{

static const char *const rf_hdf5_options    = "_options";
static const char *const rf_hdf5_ext_param  = "_ext_param";
static const char *const rf_hdf5_labels     = "labels";
static const char *const rf_hdf5_topology   = "topology";
static const char *const rf_hdf5_parameters = "parameters";

namespace detail
{

VIGRA_EXPORT void options_import_HDF5(HDF5File &, RandomForestOptions &,
                                      const std::string &);

VIGRA_EXPORT void options_export_HDF5(HDF5File &, const RandomForestOptions &,
                                      const std::string &);

VIGRA_EXPORT void dt_import_HDF5(HDF5File &, detail::DecisionTree &,
                                 const std::string &);

VIGRA_EXPORT void dt_export_HDF5(HDF5File &, const detail::DecisionTree &,
                                 const std::string &);

template<class X>
void rf_import_HDF5_to_map(HDF5File & h5context, X & param,
                           const char *const ignored_label = 0)
{
    // read a map containing all the double fields
    typedef typename X::map_type map_type;
    typedef std::pair<typename map_type::iterator, bool> inserter_type;
    typedef typename map_type::value_type value_type;
    typedef typename map_type::mapped_type mapped_type;

    map_type serialized_param;
    bool ignored_seen = ignored_label == 0;

    std::vector<std::string> names = h5context.ls();
    std::vector<std::string>::const_iterator j;
    for (j = names.begin(); j != names.end(); ++j)
    {
        if (ignored_label && *j == ignored_label)
        {
            ignored_seen = true;
            continue;
        }
        // get sort of an iterator to a new empty array vector in the map ...
        inserter_type new_array
            = serialized_param.insert(value_type(*j, mapped_type()));
        // ... and read the data into that place.
        h5context.readAndResize(*j, (*(new_array.first)).second);
    }
    vigra_precondition(ignored_seen, "rf_import_HDF5_to_map(): "
                                                         "labels are missing.");
    param.make_from_map(serialized_param);
}

template<class T>
void problemspec_import_HDF5(HDF5File & h5context, ProblemSpec<T> & param,
                             const std::string & name)
{
    h5context.cd(name);
    rf_import_HDF5_to_map(h5context, param, rf_hdf5_labels);
    // load_class_labels
    ArrayVector<T> labels;
    h5context.readAndResize(rf_hdf5_labels, labels);
    param.classes_(labels.begin(), labels.end());
    h5context.cd_up();
}

template<class X>
void rf_export_map_to_HDF5(HDF5File & h5context, const X & param)
{
    typedef typename X::map_type map_type;
    map_type serialized_param;
    // get a map containing all the double fields
    param.make_map(serialized_param);
    typename map_type::const_iterator j;
    for (j = serialized_param.begin(); j != serialized_param.end(); ++j)
        h5context.write(j->first, j->second);
}

template<class T>
void problemspec_export_HDF5(HDF5File & h5context, ProblemSpec<T> const & param,
                             const std::string & name)
{
    h5context.cd_mk(name);
    rf_export_map_to_HDF5(h5context, param);
    h5context.write(rf_hdf5_labels, param.classes);
    h5context.cd_up();
}

struct padded_number_string_data;
class VIGRA_EXPORT padded_number_string
{
private:
    padded_number_string_data* padded_number;
protected:
    padded_number_string(const padded_number_string &);
    void operator=(const padded_number_string &);
public:
    padded_number_string(int n);
    std::string operator()(int k) const;
    ~padded_number_string();
};

} // namespace detail

template<class T, class Tag>
void rf_export_HDF5(const RandomForest<T, Tag> & rf,
                    HDF5File & h5context,
                    const std::string & pathname = "")
{
    if (pathname.size())
        h5context.cd_mk(pathname);

    // save serialized options
    detail::options_export_HDF5(h5context, rf.options(), rf_hdf5_options);
    // save external parameters
    detail::problemspec_export_HDF5(h5context, rf.ext_param(),
                                    rf_hdf5_ext_param);
    // save trees
    int tree_count = rf.options_.tree_count_;
    detail::padded_number_string tree_number(tree_count);
    for (int i = 0; i < tree_count; ++i)
        detail::dt_export_HDF5(h5context, rf.tree(i), "Tree_" + tree_number(i));

    if (pathname.size())
        h5context.cd_up();
}

template<class T, class Tag>
void rf_export_HDF5(const RandomForest<T, Tag> & rf,
                    const std::string & filename, 
                    const std::string & pathname = "")
{
    HDF5File h5context(filename , HDF5File::Open);
    rf_export_HDF5(rf, h5context, pathname);
}

template<class T, class Tag>
bool rf_import_HDF5(RandomForest<T, Tag> & rf,
                    HDF5File & h5context,
                    const std::string & pathname = "")
{
    if (pathname.size())
        h5context.cd(pathname);
    // get serialized options
    detail::options_import_HDF5(h5context, rf.options_, rf_hdf5_options);
    // get external parameters
    detail::problemspec_import_HDF5(h5context, rf.ext_param_,
                                    rf_hdf5_ext_param);
    // get all groups in base path
    // no check for the "Tree_" prefix...
    std::vector<std::string> names = h5context.ls();
    std::vector<std::string>::const_iterator j;
    for (j = names.begin(); j != names.end(); ++j)
    {
        if ((*j->rbegin() == '/') && (*j->begin() != '_')) // skip the above
        {
            rf.trees_.push_back(detail::DecisionTree(rf.ext_param_));
            detail::dt_import_HDF5(h5context, rf.trees_.back(), *j);
        }
    }
    if (pathname.size())
        h5context.cd_up();
    return true;
}

template<class T, class Tag>
bool rf_import_HDF5(RandomForest<T, Tag> & rf,
                    const std::string & filename, 
                    const std::string & pathname = "")
{
    HDF5File h5context(filename, HDF5File::Open);
    return rf_import_HDF5(rf, h5context, pathname);
}

} // namespace vigra

#endif // VIGRA_RANDOM_FOREST_HDF5_IMPEX_HXX
