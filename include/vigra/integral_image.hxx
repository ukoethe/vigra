#ifndef INTEGRALIMAGE_HXX
#define INTEGRALIMAGE_HXX

#include <vigra/multi_array.hxx>
// #include <vigra/numpy_array.hxx>
// #include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_pointoperators.hxx>
#include <vigra/utilities.hxx>
#include <vigra/functorexpression.hxx>

namespace vigra {

template <class T1, class S1, class T2, class S2, class FUNCTOR>
void 
integralImage(MultiArrayView<2, T1, S1> const & image, 
              MultiArrayView<2, T2, S2> intimage,
              FUNCTOR const & functor)
{
    int width = image.shape()[1];
    int height = image.shape()[0];
    
    T2 s = 0;
    for (int i=0; i<width; ++i){
        s += functor(image(0, i));
        intimage(0, i)=s;
    }
    for (int i=1; i<height; ++i){
        s=0;
        for (int j=0; j<width; ++j){
            s += functor(image(i, j));
            intimage(i, j) = s+intimage(i-1, j);
        }
    }
}

template <class T1, class S1, class T2, class S2>
inline void 
integralImage(MultiArrayView<2, T1, S1> const & image, 
              MultiArrayView<2, T2, S2> intimage)
{
    using vigra::functor::Identity;
    integralImage(image, intimage, Identity());
}

template <class T1, class S1, class T2, class S2, class FUNCTOR>
inline void 
integralImage(MultiArrayView<3, Multiband<T1>, S1> const & image, 
              MultiArrayView<3, Multiband<T2>, S2> intimage,
              FUNCTOR const & functor)
{
    for(int c=0; c<image.shape(2); ++c)
        integralImage(image.bindOuter(c), intimage.bindOuter(c), functor);
}

template <class T1, class S1, class T2, class S2>
inline void 
integralImage(MultiArrayView<3, Multiband<T1>, S1> const & image, 
              MultiArrayView<3, Multiband<T2>, S2> intimage)
{
    using vigra::functor::Identity;
    integralImage(image, intimage, Identity());
}

template <class T1, class S1, class T2, class S2>
inline void 
integralImage2(MultiArrayView<2, T1, S1> const & image, 
               MultiArrayView<2, T2, S2> intimage)
{
    using namespace vigra::functor;
    integralImage(image, intimage, sq(Arg1()));
}

template <class T1, class S1, class T2, class S2>
inline void 
integralImage2(MultiArrayView<3, Multiband<T1>, S1> const & image, 
               MultiArrayView<3, Multiband<T2>, S2> intimage)
{
    for(int c=0; c<image.shape(2); ++c)
        integralImage2(image.bindOuter(c), intimage.bindOuter(c));
}

template <class T1, class S1, class T2, class S2>
void 
integralVolume(MultiArrayView<3, T1, S1> const & volume, 
               MultiArrayView<3, T2, S2> intvolume)
{
    //compute the integral volume of the the given volume
    //for multichannel volumes, compute channelwise
    //For more info, see [Ke, Suktankar, Hebert, "Efficient Visual Event Detection using Volumetric Features"]
    int nx = volume.shape()[0];
    int ny = volume.shape()[1];
    int nz = volume.shape()[2];

    //this vector will store s2(x-1, y) for all values of y
    std::vector<T2> s2_temp(ny, 0);
    
    
    T2 s1 = 0;
    T2 s2 = 0;
    for (int iy=0; iy<ny; ++iy){
        s2_temp[iy] = 0;
    }
    for (int ix=0; ix<nx; ++ix){
        s1 = 0;
        for (int iy=0; iy<ny; ++iy){
            s1 += volume(ix, iy, 0);
            s2 = s2_temp[iy] + s1;
            s2_temp[iy] = s2;
            intvolume(ix, iy, 0) = s2;
        }
    }
    
    for (int iz=1; iz<nz; ++iz){
        for (int iy=0; iy<ny; ++iy){
            s2_temp[iy] = 0;
        }
        for (int ix=0; ix<nx; ++ix){
            s1 = 0;
            for (int iy=0; iy<ny; ++iy){
                s1 += volume(ix, iy, iz);
                s2 = s2_temp[iy] + s1;
                s2_temp[iy] = s2;
                intvolume(ix, iy, iz) = s2 + intvolume(ix, iy, iz-1);
            }
        }
    }
}

template <class T1, class S1, class T2, class S2>
inline void 
integralVolume(MultiArrayView<4, Multiband<T1>, S1> const & volume, 
               MultiArrayView<4, Multiband<T2>, S2> intvolume)
{
    for(int c=0; c<volume.shape(3); ++c)
        integralVolume(volume.bindOuter(c), intvolume.bindOuter(c));
}

template <class T1, class S1, class T2, class S2>
void 
integralVolume2(MultiArrayView<3, T1, S1> const & volume, 
                MultiArrayView<3, T2, S2> intvolume)
{
    //compute the integral volume of the given volume squared (needed for variance calculation)
    //for multichannel volumes, compute channelwise
    //For more info, see [Ke, Suktankar, Hebert, "Efficient Visual Event Detection using Volumetric Features"]
    int nx = volume.shape()[0];
    int ny = volume.shape()[1];
    int nz = volume.shape()[2];

    //this vector will store s2(x-1, y) for all values of y
    std::vector<T2> s2_temp(ny, 0);
    
    
    T2 s1 = 0;
    T2 s2 = 0;
    for (int iy=0; iy<ny; ++iy){
        s2_temp[iy] = 0;
    }
    for (int ix=0; ix<nx; ++ix){
        s1 = 0;
        for (int iy=0; iy<ny; ++iy){
            s1 += volume(ix, iy, 0)*volume(ix, iy, 0);
            s2 = s2_temp[iy] + s1;
            s2_temp[iy] = s2;
            intvolume(ix, iy, 0) = s2;
        }
    }
    
    for (int iz=1; iz<nz; ++iz){
        for (int iy=0; iy<ny; ++iy){
            s2_temp[iy] = 0;
        }
        for (int ix=0; ix<nx; ++ix){
            s1 = 0;
            for (int iy=0; iy<ny; ++iy){
                s1 += volume(ix, iy, iz)*volume(ix, iy, iz);
                s2 = s2_temp[iy] + s1;
                s2_temp[iy] = s2;
                intvolume(ix, iy, iz) = s2 + intvolume(ix, iy, iz-1);
            }
        }
    }
}

template <class T, class S1, class S2>
void integralVolume2(MultiArrayView<4, T, S1>& volume, MultiArrayView<4, T, S2>& intvolume)
{
    //compute the integral volume of the given volume squared (needed for variance calculation)
    //for multichannel volumes, compute channelwise
    //For more info, see [Ke, Suktankar, Hebert, "Efficient Visual Event Detection using Volumetric Features"]
    int nx = volume.shape()[0];
    int ny = volume.shape()[1];
    int nz = volume.shape()[2];
    int nc = volume.shape()[3];

    //this vector will store s2(x-1, y) for all values of y
    std::vector<T> s2_temp(ny, 0);
    
    
    for (int ic=0; ic<nc; ++ic) {
        T s1 = 0;
        T s2 = 0;
        for (int iy=0; iy<ny; ++iy){
            s2_temp[iy] = 0;
        }
        for (int ix=0; ix<nx; ++ix){
            s1 = 0;
            for (int iy=0; iy<ny; ++iy){
                s1 += volume(ix, iy, 0, ic)*volume(ix, iy, 0, ic);
                s2 = s2_temp[iy] + s1;
                s2_temp[iy] = s2;
                intvolume(ix, iy, 0, ic) = s2;
            }
        }
        
        for (int iz=1; iz<nz; ++iz){
            for (int iy=0; iy<ny; ++iy){
                s2_temp[iy] = 0;
            }
            for (int ix=0; ix<nx; ++ix){
                s1 = 0;
                for (int iy=0; iy<ny; ++iy){
                    s1 += volume(ix, iy, iz, ic)*volume(ix, iy, iz, ic);
                    s2 = s2_temp[iy] + s1;
                    s2_temp[iy] = s2;
                    intvolume(ix, iy, iz, ic) = s2 + intvolume(ix, iy, iz-1, ic);
                }
            }
        }
    }
}

} // namespace vigra
    
#endif
