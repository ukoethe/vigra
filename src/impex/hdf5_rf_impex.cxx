/************************************************************************/
/*                                                                      */
/*             Copyright 2009-2010 by Ullrich Koethe                    */
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

#ifdef HasHDF5

#include "vigra/random_forest_hdf5_impex.hxx"
#include "vigra/multi_array.hxx"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

namespace vigra {

namespace detail {

struct padded_number_string_data : public std::ostringstream
{
    unsigned w;
};
padded_number_string::padded_number_string(int n)
    : padded_number(new padded_number_string_data())
{
    (*padded_number) << (n - 1);
    padded_number->w = padded_number->str().size();
}
std::string padded_number_string::operator()(int k) const
{
    padded_number->str("");
    (*padded_number) << std::setw(padded_number->w) << std::setfill('0') << k;
    return padded_number->str();
}
padded_number_string::~padded_number_string()
{
    delete padded_number;
}

void options_import_HDF5(HDF5File & h5context,
                        RandomForestOptions & opt,
                        const std::string & name)
{
    h5context.cd(name);
    rf_import_HDF5_to_map(h5context, opt);
    h5context.cd_up();
}

void options_export_HDF5(HDF5File & h5context, RandomForestOptions const & opt,
                         const std::string & name)
{
    h5context.cd_mk(name);
    rf_export_map_to_HDF5(h5context, opt);
    h5context.cd_up();
}

void dt_import_HDF5(HDF5File & h5context, detail::DecisionTree & tree,
                            const std::string & name)
{
    // check if ext_param was written(?) and read it if not
    if (tree.ext_param_.actual_msample_ == 0)
    {
        problemspec_import_HDF5(h5context, tree.ext_param_, rf_hdf5_ext_param);
        tree.classCount_ = tree.ext_param_.class_count_;
    }
    h5context.cd(name);
    // read topology
    h5context.readAndResize(rf_hdf5_topology, tree.topology_);
    // read parameters
    h5context.readAndResize(rf_hdf5_parameters, tree.parameters_);
    h5context.cd_up();
}

void dt_export_HDF5(HDF5File & h5context,
                    detail::DecisionTree const & tree,
                    const std::string & name)
{
    // make the folder for the tree.
    h5context.cd_mk(name);
    // write down topology
    h5context.write(rf_hdf5_topology, tree.topology_);
    // write down parameters
    h5context.write(rf_hdf5_parameters, tree.parameters_);
    h5context.cd_up();
}

} // namespace detail
} // namespace vigra

#endif // HasHDF5
