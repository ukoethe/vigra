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
    int width = image.shape(0);
    int height = image.shape(1);
    
    T2 s = T2();
    for (int x=0; x<width; ++x)
    {
        s += functor(image(x, 0));
        intimage(x, 0) = s;
    }
    for (int y=1; y<height; ++y)
    {
        s = T2();
        for (int x=0; x<width; ++x)
        {
            s += functor(image(x, y));
            intimage(x, y) = s + intimage(x, y-1);
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

template <class T1, class S1, class T2, class S2, class FUNCTOR>
void 
integralVolume(MultiArrayView<3, T1, S1> const & volume, 
               MultiArrayView<3, T2, S2> intvolume,
               FUNCTOR const & functor)
{
    int nx = volume.shape()[0];
    int ny = volume.shape()[1];
    int nz = volume.shape()[2];

    //this vector will store s2(x, y-1) for all values of x
    MultiArray<1, T2> s2_temp(nx);
    
    T2 s1 = T2();
    T2 s2 = T2();

    for (int iy=0; iy<ny; ++iy)
    {
        s1 = T2();
        for (int ix=0; ix<nx; ++ix)
        {
            s1 += functor(volume(ix, iy, 0));
            s2 = s2_temp(ix) + s1;
            s2_temp(ix) = s2;
            intvolume(ix, iy, 0) = s2;
        }
    }
    
    for (int iz=1; iz<nz; ++iz)
    {    
        s2_temp = T2();
        
        for (int iy=0; iy<ny; ++iy)
        {
            s1 = T2();
            for (int ix=0; ix<nx; ++ix)
            {
                s1 += functor(volume(ix, iy, iz));
                s2 = s2_temp(ix) + s1;
                s2_temp(ix) = s2;
                intvolume(ix, iy, iz) = s2 + intvolume(ix, iy, iz-1);
            }
        }
    }
}

template <class T1, class S1, class T2, class S2>
inline void 
integralVolume(MultiArrayView<3, T1, S1> const & volume, 
               MultiArrayView<3, T2, S2> intvolume)
{
    using vigra::functor::Identity;
    integralVolume(volume, intvolume, Identity());
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
inline void 
integralVolume2(MultiArrayView<3, T1, S1> const & volume, 
                MultiArrayView<3, T2, S2> intvolume)
{
    using namespace vigra::functor;
    integralVolume(volume, intvolume, sq(Arg1()));
}

template <class T1, class S1, class T2, class S2>
inline void 
integralVolume2(MultiArrayView<4, Multiband<T1>, S1> const & volume,
                MultiArrayView<4, Multiband<T2>, S2> intvolume)
{
    for(int c=0; c<volume.shape(3); ++c)
        integralVolume2(volume.bindOuter(c), intvolume.bindOuter(c));
}

} // namespace vigra
    
#endif
