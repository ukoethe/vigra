#include <cassert>
#include <limits>
#include <valarray>
#include <fstream>
#include <math.h>
#include <algorithm>

#include <boost/python.hpp>
#include <boost/ref.hpp>
#include "boost/tuple/tuple.hpp"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/tuple.hpp>

#include <numpy/arrayobject.h>
#include "vigranumpyimpex.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/multi_impex.hxx"
#include "vigra/pythonimage.hxx"
#include "vigra/pythonutil.hxx"

#include <gsl/gsl_multifit.h>

#include "splineVolumeView.hxx"
#include "gaussian_volumes.hxx"
#include "optimize.hxx"
#include "commontypes.hxx"
#include "thinning.hxx"
#include "more_thuente.hxx"

typedef vigra::Surfel<difference_type, real_type> Surfel;
typedef std::vector<Surfel> Surfels;
typedef std::vector<real_difference_type> Points;

typedef vigra::MultiArray<3, real_difference_type> gradient_volume_type;
typedef vigra::MultiArray<1, real_difference_type> real_diff_array;


template<typename VOLUME>
PyObject * 
gaussianGradient(VOLUME const & vol, double sigma) {

  vigra::TinyVector<int, 4> shape(3,vol.shape()[0],vol.shape()[1],vol.shape()[2]);
  PyObject* pyret = 
    vigra::createNumpyArray<int,4>(shape,NPY_FLOAT64);
  vigra::MultiArrayView<3, real_difference_type> ret(vol.shape(), (real_difference_type*) PyArray_DATA(pyret)); 

  vigra::gaussianGradientMultiArray(vigra::srcMultiArrayRange(vol), vigra::destMultiArray(ret), 
                    sigma); 
  return pyret;
}

template<typename VOLUME>
PyObject * 
gradientMagnitude(VOLUME const & vol, double sigma) {

  gradient_volume_type gvol(vol.shape());
  
  vigra::gaussianGradientMultiArray(vigra::srcMultiArrayRange(vol), vigra::destMultiArray(gvol), 
                    sigma); 

  PyObject* pyret = 
    vigra::createNumpyArray<int,3>(vol.shape(),NPY_FLOAT64);
  vigra::MultiArrayView<3, double> ret(vol.shape(), (double*) PyArray_DATA(pyret)); 

  vigra::transformMultiArray(vigra::srcMultiArrayRange(gvol), vigra::destMultiArray(ret), 
                 vigra::VectorNormFunctor<real_difference_type>());

  return pyret;
}

template<typename VOLUME>
PyObject * 
gaussianSmoothing(VOLUME const & vol, double sigma) 
{
  PyObject* pyret = 
    vigra::createNumpyArray<int,3>(vol.shape(),NPY_FLOAT64);
  vigra::MultiArrayView<3, double> ret(vol.shape(), (double*) PyArray_DATA(pyret)); 
  

  vigra::gaussianSmoothMultiArray
    (vigra::srcMultiArrayRange(vol), 
     vigra::destMultiArray(ret),sigma);

  return pyret;
}

PyObject * 
createSphere(double r, double sigma, difference_type shape, real_difference_type o) 
{
  PyObject* pyret = 
    vigra::createNumpyArray<int,3>(shape,NPY_FLOAT64);
  vigra::MultiArrayView<3, double> ret(shape, (double*) PyArray_DATA(pyret)); 
  vigra::MultiArray<3, double> vol(shape); 

  real_difference_type t1,t2,t3;
//   real_difference_type o
//     ((shape-vigra::NumericTraits<difference_type>::one())/2.0);

  t1[0]=t1[1]=t1[2]=0.0;
  for(;t1[2]<shape[2];t1[2]++)
    for(t2=t1;t2[1]<shape[1];t2[1]++)
      for(t3=t2;t3[0]<shape[0];t3[0]++)
    if(vigra::norm(t3-o)<=r)
      vol[t3] = 1; 
    else
      vol[t3] = 0; 

  vigra::gaussianSmoothMultiArray
    (vigra::srcMultiArrayRange(vol), 
     vigra::destMultiArray(ret),sigma);

  return pyret;
}

PyObject * 
createCylinder(double r, double sigma, 
           difference_type shape, real_difference_type o,
           real_type a,real_type b) 
{
  PyObject* pyret = 
    vigra::createNumpyArray<int,3>(shape,NPY_FLOAT64);
  vigra::MultiArrayView<3, double> ret(shape, (double*) PyArray_DATA(pyret)); 
  vigra::MultiArray<3, double> vol(shape); 

  real_difference_type t1,t2,t3,p;
  vigra::TinyVector<real_type, 2> oz(o[0],o[1]); 
//     ((shape-vigra::NumericTraits<difference_type>::one())/2.0);

  t1[0]=t1[1]=t1[2]=0.0;
  for(;t1[2]<shape[2];t1[2]++)
    for(t2=t1;t2[1]<shape[1];t2[1]++)
      for(t3=t2;t3[0]<shape[0];t3[0]++) {
    p = t3-o;
    p[0] = p[0]*cos(b)+p[2]*sin(b);
    p[1] = p[1]*cos(a)-p[2]*cos(b)*sin(a)+p[0]*sin(b)*sin(a);
    p[2] = p[2]*cos(b)*cos(a)-p[0]*cos(a)*sin(b)+p[1]*sin(a);

    if(sqrt(p[0]*p[0]+p[1]*p[1])<=r)
      vol[t3] = 1; 
    else
      vol[t3] = 0; 
      }

  vigra::gaussianSmoothMultiArray
    (vigra::srcMultiArrayRange(vol), 
     vigra::destMultiArray(ret),sigma);

  return pyret;
}


template<typename GRAD_VOLUME, typename MAG_VOLUME, typename SURFEL_ARRAY>
void findSurfels(GRAD_VOLUME const & grd, MAG_VOLUME const & mgn, 
         real_type threshold, SURFEL_ARRAY & surfels)
{
  if(threshold <= 0)
    {
      PyErr_SetString(PyExc_ValueError, 
              "threshold must be positive");
      throw_error_already_set();
    }  

  difference_type gshape(grd.shape()[1],grd.shape()[2],grd.shape()[3]);
  if(grd.shape()[0] != 3)
    {
      PyErr_SetString(PyExc_ValueError, 
              "gradient volume shape mismatch, must be (3,...)");
      throw_error_already_set();
    }  
  if(gshape != mgn.shape())
    {
      PyErr_SetString(PyExc_ValueError, 
              "gradient grdume and magnitude volume shape mismatch");
      throw_error_already_set();
    }  

  vigra::MultiArrayView<3, real_difference_type> 
    gvol(gshape, (real_difference_type*) grd.data()); 

  real_type t = 0.5 / VIGRA_CSTD::sin(M_PI/8.0);
  real_type sqrt3 = vigra::sqrt(3.);
  real_difference_type zf(0.5);
  difference_type t1,t2,t3,o;
  
  
  t1[0]=t1[1]=t1[2]=1;
  for (;t1[2]<gshape[2]-1;t1[2]++)
    for (t2=t1;t2[1]<gshape[1]-1;t2[1]++)
      for (t3=t2;t3[0]<gshape[0]-1;t3[0]++) {
    real_type mag = mgn[t3];
    if (mag < threshold)
      continue;
    
    o=vigra::floor(gvol[t3]/mag * t + zf);
    //o = gvol[t3]; o=sqrt3 * 1./vigra::norm(o);
    real_type m1=mgn[t3-o];
    real_type m2=mgn[t3];
    real_type m3=mgn[t3+o];
// 	std::cerr << m1 << std::endl;
// 	std::cerr << m2 << std::endl;
// 	std::cerr << m3 << std::endl;
    if (m1<m2 && m2>=m3) {
      surfels.push_back(Surfel(t3,m2)); //fix orientation
    }
      }
}

template<typename OPTIMIZER, typename FUNCTOR, typename MULTIARRAY,
     typename RECORDER>
object 
optimizeImpl(OPTIMIZER & opt, FUNCTOR const & func, MULTIARRAY const & points, 
         RECORDER & rec, bool debug)
{
  list ret;
  list iii;
  list tra;
  
  int terrc  = 0;
  int titerc = 0;
  int psize  = points.shape()[1];
  for(unsigned int i = 0; i < psize; ++i) 
    {  
      if(debug)
    rec.clear();

      int cerrCount = 0;
      real_difference_type 
    p(points(0,i),points(1,i),points(2,i));

       
      real_type xx,yy,zz;
      int iterc=0;
      real_type fretb = func(p);
      try {
    iii.append(make_tuple(p[0],p[1],p[2]));
    //p -= func.dir(p);
    opt(p,iterc);
    real_type freta = func(p);
    //if(freta >= fretb)
    //npoints->push_back(p); //CHANGE
    ret.append(make_tuple(p[0],p[1],p[2]));
    if(debug)
      tra.append(rec.ret);
      }
      catch(std::exception & e)
    {
      //std::cerr << e.what() << std::endl; // print message
      terrc++;
    }
      catch(std::string & e)
    {
      //std::cerr << e << std::endl; // print message
      terrc++;
    }
      titerc += iterc;
    }
  if(debug)
    return make_tuple(iii,ret,tra);
  return ret;
}


struct IterationStepsRecorder
{
  list ret;

  template<typename VALUETYPE>
  void operator()(VALUETYPE const & p) 
  {
    ret.append(make_tuple(p[0],p[1],p[2]));
  }

  void clear() {ret = list();}
};

template<typename VOLUME, typename MULTIARRAY>
object 
cgBrend(VOLUME const & vol, MULTIARRAY const & points, int maxiter, bool debug)
{
  typedef vigra::InvertFunctor<VOLUME> InvertFunctor;

  if(points.shape()[0] != 3)
    {
      PyErr_SetString(PyExc_ValueError, 
              "point array shape mismatch (should be (3,x))");
      throw_error_already_set();
    }  
  
  InvertFunctor func(vol);     
  IterationStepsRecorder * rec = NULL;
  if(debug)
    rec = new IterationStepsRecorder();
    
  vigra::CG<vigra::LineSearchBrend<InvertFunctor>,
    IterationStepsRecorder > cg(func, maxiter, rec);
  
  return optimizeImpl(cg, func, points, *rec, debug);
}

template<typename VOLUME, typename MULTIARRAY>
object 
cgMore(VOLUME const & vol, MULTIARRAY const & points, int maxiter, bool debug)
{
  typedef vigra::InvertFunctor<VOLUME> InvertFunctor;

  if(points.shape()[0] != 3)
    {
      PyErr_SetString(PyExc_ValueError, 
              "point array shape mismatch (should be (3,x))");
      throw_error_already_set();
    }  
  
  InvertFunctor func(vol);     
  IterationStepsRecorder * rec = NULL;
  if(debug)
    rec = new IterationStepsRecorder();
    
  vigra::CG<vigra::LineSearchMore<InvertFunctor>,
    IterationStepsRecorder > cg(func, maxiter, rec);
  
  return optimizeImpl(cg, func, points, *rec, debug);
}

template<typename VOLUME, typename MULTIARRAY>
object 
mcgMore(VOLUME const & vol, VOLUME const & mgn,
    MULTIARRAY const & points, int maxiter, bool debug)
{
  typedef vigra::InvertFunctor<VOLUME> InvertFunctor;

  if(points.shape()[0] != 3)
    {
      PyErr_SetString(PyExc_ValueError, 
              "point array shape mismatch (should be (3,x))");
      throw_error_already_set();
    }  

  vigra::TestFunctor<VOLUME> func(vol, mgn);
  IterationStepsRecorder * rec = NULL;
  if(debug)
    rec = new IterationStepsRecorder();

  vigra::CG<vigra::LineSearchMore
    <vigra::TestFunctor<VOLUME> > ,
    IterationStepsRecorder> 
    cg(func, maxiter, rec);
      
  return optimizeImpl(cg, func, points, *rec, debug);
}

template<typename VOLUME, typename MULTIARRAY>
object 
mcgBrend(VOLUME const & vol, VOLUME const & mgn,
     MULTIARRAY const & points, int maxiter, bool debug)
{
  typedef vigra::InvertFunctor<VOLUME> InvertFunctor;

  if(points.shape()[0] != 3)
    {
      PyErr_SetString(PyExc_ValueError, 
              "point array shape mismatch (should be (3,x))");
      throw_error_already_set();
    }  

  vigra::TestFunctor<VOLUME> func(vol, mgn);
  IterationStepsRecorder * rec = NULL;
  if(debug)
    rec = new IterationStepsRecorder();

  vigra::CG<vigra::LineSearchBrend
    <vigra::TestFunctor<VOLUME> > ,
    IterationStepsRecorder> 
    cg(func, maxiter, rec);
      
  return optimizeImpl(cg, func, points, *rec, debug);
}

template<typename VOLUME, typename MULTIARRAY>
object 
mcgBfgs(VOLUME const & vol, VOLUME const & mgn,
     MULTIARRAY const & points, int maxiter, bool debug)
{
  typedef vigra::InvertFunctor<VOLUME> InvertFunctor;

  if(points.shape()[0] != 3)
    {
      PyErr_SetString(PyExc_ValueError, 
              "point array shape mismatch (should be (3,x))");
      throw_error_already_set();
    }  

  vigra::TestFunctor<VOLUME> func(vol, mgn);
  IterationStepsRecorder * rec = NULL;
  if(debug)
    rec = new IterationStepsRecorder();

  vigra::BFGS<vigra::LineSearchBrend
    <vigra::TestFunctor<VOLUME> > >
    bfgs(func, maxiter);
      
  return optimizeImpl(bfgs, func, points, *rec, false);
}


template<typename GRAD_VOLUME, typename MGN_VOLUME, 
     typename MULTIARRAY>
object 
parabolaFits(GRAD_VOLUME const & grd, MGN_VOLUME const & mgn,
         MULTIARRAY const & points,  bool debug)
{
  if(points.shape()[0] != 3)
    {
      PyErr_SetString(PyExc_ValueError, 
              "point array shape mismatch (should be (3,x))");
      throw_error_already_set();
    }  

  
//     p = vigra.Vector3(p[0],p[1],p[2])
//     if d is None:
//         d = vol.d(p)

  vigra::SplineVolumeView<1, real_type> smv(mgn);

  real_type t = 0.5 / VIGRA_CSTD::sin(M_PI/8.0);
  real_difference_type zf(0.5);

  list ret;
  int n = 3;
  int x[3] = {-1,0,1};
  
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;
  X = gsl_matrix_alloc (n, 3);
  y = gsl_vector_alloc (n);
  w = gsl_vector_alloc (n);  
  c = gsl_vector_alloc (3);
  cov = gsl_matrix_alloc (3, 3);
  
  for(unsigned int j = 0; j < points.shape()[1]; j++) {
    difference_type p;
    p[0] = points(0,j);
    p[1] = points(1,j);
    p[2] = points(2,j);

    real_difference_type d = grd[p];
    difference_type di = vigra::floor(d/mgn[p] * t + zf);
    
    d = d/vigra::norm(d);
    d = d*vigra::dot(di,d);
    //real_type f = vigra::sqrt(d[0]+d[1]+[2]);
    
    for (unsigned int i = 0; i < n; i++)
      {           
    gsl_matrix_set (X, i, 0, 1.0);
    gsl_matrix_set (X, i, 1, x[i]);
    gsl_matrix_set (X, i, 2, x[i]*x[i]);
      
    gsl_vector_set (y, i, smv[p + x[i]*d]);
      }
 
    gsl_multifit_linear_workspace * work 
      = gsl_multifit_linear_alloc (n, 3);
    gsl_multifit_linear (X, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);
    
#define C(i) (gsl_vector_get(c,(i)))
    real_difference_type pr = p + d*(-C(1)/(2.*C(2)));

    //vigra_precondition(C(2)<0.0,"parabola invalid value");

    ret.append(make_tuple(pr[0],pr[1],pr[2]));
  }
  
  gsl_matrix_free (X);
  gsl_vector_free (y);
  gsl_vector_free (w);
  gsl_vector_free (c);
  gsl_matrix_free (cov);

  return ret;
}


void thinning(Surfels const & surfels, Surfels  & surfels_thin, 
          difference_type shape,
          bool erode_vol, bool csurfs, bool ccurvs, 
          bool osurfs, bool ocurvs, bool spoints)
{
//   if(points.shape()[0] != 3)
//     {
//       PyErr_SetString(PyExc_ValueError, 
// 		      "point array shape mismatch (should be (3,x))");
//       throw_error_already_set();
//     }  

//   if(points.shape()[1] != priorities.shape()[0])
//     {
//       PyErr_SetString(PyExc_ValueError, 
// 		      "number of priorities must match number of points");
//       throw_error_already_set();
//     }  

  //int size = points.shape()[1];
//   vigra::ArrayVector<Surfel> surfels(size);
//   for(unsigned int i = 0; i < size; i++){
//     surfels[i] = Surfel(difference_type(points(0,i),
// 					points(1,i),
// 					points(2,i)), 
// 			priorities[i]);
//   }

  //std::cerr << points.shape() << std::endl; 
  
  //Surfels  surfels_thin;// = new Surfels();  
  vigra::surfelThinning(surfels,surfels_thin,shape,
            erode_vol, csurfs, ccurvs, osurfs, ocurvs, spoints);

  //  std::cerr << surfels_thin.size() << std::endl;

//   list ret;
//   for(unsigned int i = 0; i < surfels_thin.size(); i++){
//     difference_type & p = surfels_thin[i].p;
//       ret.append(make_tuple(p[0],p[1],p[2]));
//   }
//   return ret;
}

template<typename VOLUMEVIEW>
real_type brentLineSearch(real_difference_type p, real_difference_type u, VOLUMEVIEW const & vol)
{
  typedef vigra::InvertFunctor<VOLUMEVIEW> Functor;
  Functor func(vol);
  vigra::LineSearchBrend<Functor> brent(func);
  return brent(p,u);
}

template<typename VOLUMEVIEW>
real_type moreLineSearch(real_difference_type p, real_difference_type u, VOLUMEVIEW const & vol)
{
  typedef vigra::InvertFunctor<VOLUMEVIEW> Functor;
  Functor func(vol);
  vigra::LineSearchMore<Functor> more(func);
  return more(p,u);
}



template<typename VOLUME>
object canny(VOLUME const & vol, real_type sigma, 
         real_type thresh, const char * optimizer,
         bool thin, bool print,
         int maxiter,
         int dist, double shift,
         bool erode_vol, bool csurfs, bool ccurvs, 
         bool osurfs, bool ocurvs, bool spoints)
{
//   vigra::VolumeImportInfo info(filename);
//   if(!info.isGrayscale())
//     {
//       PyErr_SetString(PyExc_ValueError, 
// 		      "volume type not supported");
//     }  
//   if(print)
//     std::cerr << "Loading volume '" << info.name() << "' "
// 	      << info.shape() << ", "
// 	      << (info.isColor() ? "color" : "grayscale") << "...";
//   vigra::MultiArray<3, real_type> vol(info.shape()); 
//   vigra::importVolume(info,vol);
//   if(print) std::cerr << "done" << std::endl;

  bool debug = false;

  if(print) std::cerr << "computing gaussian gradient ...";
  vigra::MultiArray<3, real_difference_type> grd(vol.shape()); 
  vigra::gaussianGradientMultiArray(vigra::srcMultiArrayRange(vol), 
                    vigra::destMultiArray(grd), 
                    sigma);
  //dirty FIX
  vigra::TinyVector<int, 4> shape(3,
                  vol.shape()[0],
                  vol.shape()[1],
                  vol.shape()[2]);
  UMArray4FLOAT64 tmp(shape, (real_type*) grd.data());
  if(print) std::cerr << "done" << std::endl;

  if(print) std::cerr << "computing gaussian gradient magnitude ...";
  vigra::MultiArray<3, real_type> mgn(vol.shape()); 
  vigra::transformMultiArray(vigra::srcMultiArrayRange(grd), 
                 vigra::destMultiArray(mgn), 
                 vigra::VectorNormFunctor<real_difference_type>());
  if(print) std::cerr << "done" << std::endl;

  vigra::SplineVolumeView<5, real_type> smv(mgn);

  Surfels * surfels = new Surfels();
  findSurfels(tmp,mgn,thresh,*surfels);
  //findSurfels(tmp,smv,thresh,*surfels);
  if(thin) {
    Surfels * surfels_thin = new Surfels();
    if(print) std::cerr << "apply thinning ...";
    thinning(*surfels, *surfels_thin, vol.shape(),
         erode_vol, csurfs, ccurvs, osurfs, ocurvs, spoints);
    delete surfels;
    surfels = surfels_thin;
    if(print) std::cerr << "done" << std::endl;
    if(print) std::cerr << "found " << surfels->size() 
            << " surfels"  << std::endl;
  }

  //remove unness copying
  if(print) std::cerr << "copying ...";
  unsigned int size = surfels->size();
  typedef vigra::MultiArray<2, real_type> Array2D;
  Array2D ptmp(Array2D::difference_type(3,size));
  real_difference_type o = 
    (vol.shape()-vigra::NumericTraits<real_difference_type>::one())/2.;
  int s = 0;
  real_type min = 100000000;
  o += real_difference_type(shift,shift,shift);
  for(unsigned int i = 0; i < size; i++) {
    real_difference_type p =  (*surfels)[i].p;
    real_difference_type ps = p-o;
    if(norm(ps) < dist) { //Filter ptmp
      ptmp(0,s) = p[0];
      ptmp(1,s) = p[1];
      ptmp(2,s) = p[2];
      s++;
    }
  }
  vigra::MultiArrayView<2, real_type> 
    points(Array2D::difference_type(3,s),ptmp.data());
  if(print) std::cerr << "done" << std::endl;

  vigra::MultiArray<3, double> sol(vol.shape()); 
  vigra::gaussianSmoothMultiArray
    (vigra::srcMultiArrayRange(vol), 
     vigra::destMultiArray(sol),sigma);
  vigra::SplineVolumeView<5, real_type> svv(sol);

  if (!strcmp(optimizer,"mcgMore")) // is mcgMore1Iter
    {
      if(print) std::cerr << "apply subvoxel refinement";

      return 
    mcgMore(svv, smv, points, maxiter, debug);
      if(print) std::cerr << "done" << std::endl;
    }
  else if (!strcmp(optimizer,"mcgBrend")) // is mcgMore1Iter
    {
      if(print) std::cerr << "apply subvoxel refinement";

      return 
    mcgBrend(svv, smv, points, maxiter, debug);
      if(print) std::cerr << "done" << std::endl;
    }
  else if (!strcmp(optimizer,"mcgBfgs")) // is mcgMore1Iter
    {
      if(print) std::cerr << "apply subvoxel refinement";

      return 
    mcgBfgs(svv, smv, points, maxiter, false);
      if(print) std::cerr << "done" << std::endl;
    }
  else if (!strcmp(optimizer,"brend")) // is mcgMore1Iter
    {
      if(print) std::cerr << "apply subvoxel refinement";
      return 
    mcgBrend(svv, smv, points, 1, debug);
      if(print) std::cerr << "done" << std::endl;
    }
  else if (!strcmp(optimizer,"cgMore")) // is mcgMore1Iter
    {
      if(print) std::cerr << "apply subvoxel refinement";
      return 
    cgMore(smv, points, maxiter, debug);
      if(print) std::cerr << "done" << std::endl;
    }
  else if (!strcmp(optimizer,"cgBrend")) // is mcgMore1Iter
    {
      if(print) std::cerr << "apply subvoxel refinement";
      return 
    cgBrend(smv, points, maxiter, debug);
      if(print) std::cerr << "done" << std::endl;
    }
  else if (!strcmp(optimizer,"parabolafits")) // is mcgMore1Iter
    {
      if(print) std::cerr << "apply subvoxel refinement";
      return 
    parabolaFits(grd, mgn, points, debug);
      if(print) std::cerr << "done" << std::endl;
    }
  else
    {
      std::cerr << optimizer << std::endl;
      PyErr_SetString(PyExc_ValueError, 
              "optimizing flag unknown");
      throw_error_already_set();
    }    
}

template<typename SURFEL>
struct surfel_to_tuple
{
  static PyObject* convert(SURFEL const& x)
  { return incref( make_tuple(x.p[0],x.p[1],x.p[2] ).ptr() ); }
};



void defCanny()
{
  import_array();
  exportVigraNumpyImpex();
    
  def("gaussianSmoothing", &gaussianSmoothing<volume_input_type>);
  def("gradientMagnitude", &gradientMagnitude<volume_input_type>);
  def("gaussianGradient", &gaussianGradient<volume_input_type>);
  def("findSurfels", &findSurfels<UMArray4FLOAT64, volume_input_type,Surfels>);
  def("thinning", &thinning,
      (arg("erode_vol")=true,arg("csurfs")=true,arg("ccurvs")=true,
       arg("osurfs")=true,arg("ocurvs")=true,arg("spoints")=true));
  def("thinning", &thinning,
      (arg("erode_vol")=true,arg("csurfs")=true,arg("ccurvs")=true,
       arg("osurfs")=true,arg("ocurvs")=true,arg("spoints")=true));
  def("cgBrend", &cgBrend<SplineVolumeView2, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("cgBrend", &cgBrend<SplineVolumeView5, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("cgBrend", &cgBrend<vigra::G2VolumeView<GaussSphere>, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("cgMore", &cgMore<SplineVolumeView2, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("cgMore", &cgMore<SplineVolumeView5, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("cgMore", &cgMore<vigra::G2VolumeView<GaussSphere>, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("mcgMore", &mcgMore<SplineVolumeView2, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("mcgMore", &mcgMore<SplineVolumeView5, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("mcgBrend", &mcgMore<SplineVolumeView5, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("mcgBfgs", &mcgBfgs<SplineVolumeView5, UMArray2FLOAT64>,
      (arg("maxiter")=100,arg("debug")=false));
  def("brentLineSearch", &brentLineSearch<SplineVolumeView2>);
  def("brentLineSearch", &brentLineSearch<SplineVolumeView5>);
  def("moreLineSearch", &moreLineSearch<SplineVolumeView2>);
  def("moreLineSearch", &moreLineSearch<SplineVolumeView5>);
  def("moreLineSearch", &moreLineSearch<vigra::G2VolumeView<GaussSphere> >);

  def("createSphere", &createSphere);
  def("createCylinder", &createCylinder);


  def("canny", &canny<volume_input_type>,
      (arg("thin")=true,arg("print")=true,
       arg("maxiter")=100,arg("dist")=10000,arg("shift")=0.0,
       arg("erode_vol")=true,arg("csurfs")=true,arg("ccurvs")=true,
       arg("osurfs")=true,arg("ocurvs")=true,arg("spoints")=true));

  to_python_converter<Surfel, surfel_to_tuple<Surfel> >();

  class_<Surfels >("Surfels")
    .def(vector_indexing_suite<Surfels,true >())
    ;
}
