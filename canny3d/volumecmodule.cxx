#include <cassert>
#include <limits>
#include <valarray>
#include <fstream>
#include <math.h>
#include <algorithm>

#include <boost/python.hpp>
#include <boost/ref.hpp>
#include "boost/tuple/tuple.hpp"

#include <numpy/arrayobject.h>
#include "vigranumpyimpex.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/multi_impex.hxx"
#include "vigra/pythonimage.hxx"
#include "vigra/pythonutil.hxx"  

#include "splineVolumeView.hxx"
#include "gaussian_volumes.hxx"
#include "commontypes.hxx"


//FIXME: difference_type arg
#define SPLINE_VIEW_MEBMER_FUNCTION_XYZ(view,function)			\
  volume_value_type (view::*function##xyz)(real_difference_type const &d) const \
    = &view::function;


template <class VolumeView>
struct ConstTraits   
{
  static VolumeView *
  makeView(volume_input_type const & vol)
  {  return new VolumeView(vol); }
};

template <>
struct ConstTraits<GaussPlane>
{
  static GaussPlane *
  makeView1(real_difference_type d, 
       real_difference_value_type trans, 
       real_difference_value_type sigma, 
       volume_value_type p, volume_value_type q,
       difference_type shape) 
  { std::cerr << d << std::endl;
    return new GaussPlane(d,trans,sigma,p,q,shape); }

  static GaussPlane *
  makeView2(real_difference_type d,  
       real_difference_value_type sigma, 
       volume_value_type p, volume_value_type q,
       difference_type shape) 
  { return new GaussPlane(d,0.0,sigma,p,q,shape); }


};

template <>
struct ConstTraits<GaussSphere>
{
  static GaussSphere *
  makeView1(real_type radius, real_type sigma, 
       real_type p, real_type q, 
       difference_type shape,
       real_difference_type trans)  
  { return new GaussSphere(radius,sigma,p,q,shape,trans); }
  
  static GaussSphere *
  makeView2(real_type radius, real_type sigma, 
       real_type p, real_type q, 
       difference_type shape)  
  { return new GaussSphere(radius,sigma,p,q,shape); }

};

// template <>
// struct ConstTraits<vigra::G2VolumeView<GaussSphere> >
// {
//   static vigra::G2VolumeView<GaussSphere> *
//   makeView(real_type radius, real_type sigma, 
// 	   real_type p, real_type q, 
// 	   difference_type shape)  
//   { return new vigra::G2VolumeView<GaussSphere>
//       (GaussSphere(radius,sigma,p,q,shape)); }
// };


template <>
struct ConstTraits<GaussCylinder>
{
  static GaussCylinder *
  makeView1(real_type radius, real_type sigma, 
        real_type p, real_type q, 
        difference_type shape,
        real_type a,real_type b,
       real_difference_value_type trans)  
  { return new GaussCylinder(radius,sigma,p,q,shape,a,b,trans); }

  static GaussCylinder *
  makeView2(real_type radius, real_type sigma, 
        real_type p, real_type q, 
        difference_type shape)  
  { return new GaussCylinder(radius,sigma,p,q,shape); }
};

template <class SplineView>
SplineView *
pySplineView(volume_input_type const & vol)
{
  return new SplineView(vol);
}

#define GETITEM(name,function)						\
  template <class VolumeView, class VectorPos, class RetType>		\
  RetType getitem_##name(VolumeView const & self, VectorPos const & d)	\
  {									\
    if(d.size() != 3)							\
      {									\
    PyErr_SetString(PyExc_IndexError, "index must have size 3.");	\
    throw_error_already_set();					\
      }									\
    return self.function(d);						\
  }									\
                                    
GETITEM(callf,operator())
  GETITEM(calld,d)
  GETITEM(callg2,g2)


template <class VolumeView>	
object shape(VolumeView const & self)					
{ return make_tuple(self.shape()[0],self.shape()[1],self.shape()[2]); }
  
  
#define GRID_FUNCS(name,function)					\
  template <class VolumeView, class ArrayType>				\
  PyObject *  grid_##name(VolumeView const & self, ArrayType const & grid) \
  {									\
    if(grid.shape()[0] != 3)						\
      {									\
    PyErr_SetString(PyExc_IndexError, "only grid with dim=3 supported."); \
    throw_error_already_set();					\
      }									\
                                        \
    int w = grid.shape()[1];						\
    int h = grid.shape()[2];						\
    int d = grid.shape()[3];						\
                                        \
    PyObject* pyret =							\
      vigra::createNumpyArray<int,3>(difference_type(w,h,d),NPY_FLOAT64); \
    vigra::MultiArrayView<3, double> ret(difference_type(w,h,d), (double*) PyArray_DATA(pyret)); \
    for(unsigned int z = 0; z < d; z++)					\
      for(unsigned int y = 0; y < h; y++)				\
    for(unsigned int x = 0; x < w; x++) \ 
  ret(x,y,z) = self.function(grid(0,x,y,z),grid(1,x,y,z),grid(2,x,y,z)); \
                                    \
return pyret;								\
}									\

GRID_FUNCS(callf,operator())
  GRID_FUNCS(callg2,g2)
  GRID_FUNCS(calldx,dx)
  GRID_FUNCS(calldy,dy)
  GRID_FUNCS(calldz,dz)
  
template <class VolumeView>
class_<VolumeView> &
defVolumeView(char const * name)
{
  volume_value_type (VolumeView::*fxyz)(difference_value_type,
                    difference_value_type,
                    difference_value_type) const              
    
    = &VolumeView::operator();
  
//   SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dx);
//   SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dy);
//   SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dz);
  
  
  static class_<VolumeView> theclass(name, no_init);
  theclass
    .def("__init__", make_constructor(&ConstTraits<VolumeView>::makeView1))
    .def("__init__", make_constructor(&ConstTraits<VolumeView>::makeView2))

    //.def("__call__", fxyz)
    .def("__getitem__", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_callf)
    .def("__getitem__", (PyObject* (*)(VolumeView const &, UMArray4Float64 const &))&grid_callf)
    //     .def("__getitem__", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callf)
    //     .def("__getitem__", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callf)
    .def("__getitem__", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callf)
//     .def("d",  (real_difference_type  (*)(VolumeView const &, real_difference_type const &))&getitem_calld)
//     .def("g2", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_callg2)
//     .def("g2", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callg2)
//     .def("dx", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_calldx)
//     .def("dx", dxxyz)							
//     .def("dy", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_calldy)
//     .def("dy", dyxyz)							
//     .def("dz", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_calldz)
//     .def("dz", dzxyz)						       
    //.def("isInside", &VolumeView::isInside)
    .def("size", &shape<VolumeView>)						       
    .add_property("shape",&shape<VolumeView>)
    ;
  return theclass;
}

real_type planeDistance(GaussPlane const & plane, 
            real_difference_type const & v)
{
  return plane.distance(v);
}

object planeDistance(GaussPlane const & plane, 
            UMArray2FLOAT64 const & arr)
{
  if(arr.shape()[0] != 3)
    {
      PyErr_SetString(PyExc_ValueError, 
              " must be (3,...)");
      throw_error_already_set();
    }  

  list ret;
  unsigned int size = arr.shape()[1];
  for(unsigned int i = 0; i < size; i++) {
    real_difference_type 
      p(arr(0,i),arr(1,i),arr(2,i));
    ret.append(plane.distance(p));
  }
  return ret ;
}

real_type cylDistance(GaussCylinder const & cyl, 
              real_difference_type const & v)
{
  return cyl.distance(v);
}

object cylDistance(GaussCylinder const & cyl, 
           UMArray2FLOAT64 const & arr)
{
  if(arr.shape()[0] != 3)
    {
      PyErr_SetString(PyExc_ValueError, 
              " must be (3,...)");
      throw_error_already_set();
    }  

  list ret;
  unsigned int size = arr.shape()[1];
  for(unsigned int i = 0; i < size; i++) {
    real_difference_type 
      p(arr(0,i),arr(1,i),arr(2,i));
    ret.append(cyl.distance(p));
  }
  return ret ;
}


template <class VolumeView>
class_<VolumeView> &
defSplineVolumeView(char const * name)
{
  volume_value_type (VolumeView::*fxyz)(difference_value_type,
                    difference_value_type,
                    difference_value_type) const              
    
    = &VolumeView::operator();
  
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dx);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dy);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dz);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dxx);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dyy);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dzz);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dxy);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dxz);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dyz);
  
  static class_<VolumeView> theclass(name, no_init);
  theclass
    .def("__init__", make_constructor(&ConstTraits<VolumeView>::makeView))

    //.def("__call__", fxyz)
    .def("__getitem__", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_callf)
    .def("__getitem__", (PyObject* (*)(VolumeView const &, UMArray4Float64 const &))&grid_callf)
    //     .def("__getitem__", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callf)
    //     .def("__getitem__", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callf)
    .def("__getitem__", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callf)
    .def("d",  (real_difference_type  (*)(VolumeView const &, real_difference_type const &))&getitem_calld)
    .def("g2", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_callg2)
    .def("g2", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callg2)
    .def("dx", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_calldx)
    .def("dx", dxxyz)							
    .def("dy", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_calldy)
    .def("dy", dyxyz)							
    .def("dz", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_calldz)
    .def("dz", dzxyz)						       
    .def("dxx", dxxxyz)							
    .def("dyy", dyyxyz)							
    .def("dzz", dzzxyz)							
    .def("dxy", dxyxyz)							
    .def("dxz", dxzxyz)							
    .def("dyz", dyzxyz)							
    //.def("isInside", &VolumeView::isInside)						       
    .def("size", &shape<VolumeView>)						       
    .add_property("shape",&shape<VolumeView>)
    ;
  return theclass;
}

template <class VolumeView>
class_<VolumeView> &
defTestVolumeView(char const * name)
{
  volume_value_type (VolumeView::*fxyz)(difference_value_type,
                    difference_value_type,
                    difference_value_type) const   
    = &VolumeView::operator();
  
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dx);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dy);
  SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dz);
//   SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dxx);
//   SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dyy);
//   SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dzz);
//   SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dxy);
//   SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dxz);
//   SPLINE_VIEW_MEBMER_FUNCTION_XYZ(VolumeView,dyz);
  
  static class_<VolumeView> theclass(name, no_init);
  theclass
    .def("__init__", make_constructor(&ConstTraits<VolumeView>::makeView))
    
    //.def("__call__", fxyz)
    .def("__getitem__", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_callf)
    .def("__getitem__", (PyObject* (*)(VolumeView const &, UMArray4Float64 const &))&grid_callf)
    //     .def("__getitem__", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callf)
    //     .def("__getitem__", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callf)
    .def("__getitem__", (volume_value_type (*)(VolumeView const &, real_difference_type const &))&getitem_callf)
    .def("d",  (real_difference_type  (*)(VolumeView const &, real_difference_type const &))&getitem_calld)
    .def("dx", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_calldx)
    .def("dx", dxxyz)							
    .def("dy", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_calldy)
    .def("dy", dyxyz)							
    .def("dz", (PyObject* (*)(VolumeView const &, UMArray4Int32 const &))&grid_calldz)
    .def("dz", dzxyz)						       
    //.def("isInside", &VolumeView::isInside)						       
    .def("size", &shape<VolumeView>)						    .add_property("shape",&shape<VolumeView>)
    ;
  return theclass;
}


template <class Diff>
struct DiffFromPython
{
  typedef Diff Type;
  typedef typename Type::value_type value_type;
  
  DiffFromPython()
  {
    converter::registry::insert(&convertible, &construct, type_id<Type>());
  }
  
  static void* convertible(PyObject* obj)
  {
    bool iscon = (((PyTuple_Check(obj) || PyList_Check(obj)) &&
           PySequence_Size(obj) == 3 && 
           extract<value_type>
           (PySequence_Fast_GET_ITEM(obj, 0)).check() &&
           extract<value_type>
           (PySequence_Fast_GET_ITEM(obj, 1)).check() &&
           extract<value_type>
           (PySequence_Fast_GET_ITEM(obj, 2)).check()) ||
          (PyArray_Check(obj) && 
           PyArray_NDIM(obj) == 1 && 
           PyArray_ITEMSIZE(obj) == sizeof(value_type) &&
           PyArray_SIZE(obj) == 3)
          );
    
    if(iscon)
      return obj;
//     if(PyTuple_Check(obj) || PyList_Check(obj))
//        std::cerr <<   extract<value_type>
//        (PySequence_Fast_GET_ITEM(obj, 0)).check() << std::endl;
//          std::cerr << "convertible called." << std::endl;
//          std::cerr << "check:" << PyArray_Check(obj) << std::endl;
//          std::cerr << "behd :" << PyArray_ISBEHAVED_RO(obj) << std::endl;
//          std::cerr << "dim  :" << PyArray_NDIM(obj)   << std::endl;
//          std::cerr << "size :" << (sizeof(value_type) == PyArray_ITEMSIZE(obj)) << std::endl;
     
    return 0;
  }

  static void construct(PyObject* obj, converter::rvalue_from_python_stage1_data* data)
  {
    value_type x,y,z;
    if(PyTuple_Check(obj) || PyList_Check(obj)) {
      x = extract<value_type>(PySequence_Fast_GET_ITEM(obj, 0))();
      y = extract<value_type>(PySequence_Fast_GET_ITEM(obj, 1))();
      z = extract<value_type>(PySequence_Fast_GET_ITEM(obj, 2))();
    } else {
      value_type * data = (value_type*) PyArray_DATA(obj);
      x = data[0];
      y = data[1];
      z = data[2];
    }
    void* const storage = ((converter::rvalue_from_python_storage<Type>*)data)->storage.bytes;
    new (storage) Type(x,y,z);
    data->convertible = storage;
  }
};

PyObject * 
importVolume(const char * filename) {
  
  vigra::VolumeImportInfo info(filename);

  if(!info.isGrayscale())
    {
      PyErr_SetString(PyExc_ValueError, 
              "volume type not supported");
    }  

  std::cerr << "Loading volume '" << info.name() << "' "
        << info.shape() << ", "
        << (info.isColor() ? "color" : "grayscale") << "...";

  PyObject* pyret = 
    vigra::createNumpyArray<int,3>(info.shape(),NPY_FLOAT64);
  vigra::MultiArrayView<3, double> 
    ret(info.shape(), (double*) PyArray_DATA(pyret)); 
  
  vigra::MultiArray<3, double> tmp(info.shape()); 
  vigra::importVolume(info,tmp);
  vigra::copyMultiArray(srcMultiArrayRange(tmp), destMultiArray(ret));

  std::cerr << "done" << std::endl;

  return pyret;
}

void defVolume()
{
  import_array();
  exportVigraNumpyImpex();
  
  defSplineVolumeView<vigra::SplineVolumeView<2, volume_value_type> >("SplineVolumeView2");
  defSplineVolumeView<vigra::SplineVolumeView<3, volume_value_type> >("SplineVolumeView3");
  defSplineVolumeView<vigra::SplineVolumeView<4, volume_value_type> >("SplineVolumeView4");
  defSplineVolumeView<vigra::SplineVolumeView<5, volume_value_type> >("SplineVolumeView5");
  
  defVolumeView<GaussPlane>("GaussianPlane");
  defVolumeView<GaussCylinder>("GaussianCylinder");
  defVolumeView<GaussSphere>("GaussianSphere");

  //defTestVolumeView<vigra::G2VolumeView<GaussSphere> >("MagGaussianSphere");

  DiffFromPython<difference_type>();
  DiffFromPython<real_difference_type>();
  
  def("importVolume", &importVolume);

  def("planeDistance", (object (*)(GaussPlane const &, UMArray2FLOAT64 const &))&planeDistance);
  def("planeDistance", (real_type (*)(GaussPlane const &, real_difference_type const &))&planeDistance);

  def("cylinderDistance", (object (*)(GaussCylinder const &, UMArray2FLOAT64 const &))&cylDistance);
  def("cylinderDistance", (real_type (*)(GaussCylinder const &, real_difference_type const &))&cylDistance);

}  
