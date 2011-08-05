#ifndef COMMONTYPES_HXX
#define COMMONTYPES_HXX

#include "vigra/numerictraits.hxx"
#include <vigra/multi_array.hxx>

typedef double    real_type;
typedef real_type volume_value_type;
typedef real_type difference_value_type;
typedef real_type real_difference_value_type;

typedef vigra::TinyVector<volume_value_type,3> real_difference_type;
typedef vigra::TinyVector<int,3> difference_type;
typedef vigra::TinyVector<long int,3> long_difference_typeqq;

typedef vigra::MultiArrayView<3, volume_value_type,  vigra::UnstridedArrayTag> volume_input_type;
typedef vigra::MultiArray<3, real_difference_type> gradient_volume_type;
typedef vigra::MultiArrayView<3, volume_value_type,  vigra::UnstridedArrayTag> real_volume_type;

typedef vigra::SplineVolumeView<2, volume_value_type> SplineVolumeView2;
typedef vigra::SplineVolumeView<5, volume_value_type> SplineVolumeView5;

typedef vigra::GaussianPlane<volume_value_type>  GaussPlane;
typedef vigra::GaussianSphere<volume_value_type>  GaussSphere;
typedef vigra::GaussianCylinder<volume_value_type>  GaussCylinder;

#endif
