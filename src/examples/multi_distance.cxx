#include <iostream>
#include <iomanip>
#include <cmath>

#include "vigra/vectorial_distance.hxx"
#include "vigra/vectorial_boundary_distance.hxx"
#include "vigra/multi_distance.hxx"
#include "vigra/labelimage.hxx"
#include "vigra/labelvolume.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/impex.hxx"
#include "vigra/tg.hxx"
#include "vigra/hdf5impex.hxx"
#include "vigra/random.hxx"
#include "vigra/timing.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/eigensystem.hxx"

#include "vis_vectorial_dist.hxx"

typedef vigra::MultiArray<2,int> IntImage;
typedef vigra::MultiArray<2,vigra::TinyVector<int, 2> > IntVectorImage;
typedef vigra::MultiArray<2,vigra::TinyVector<double, 2> > DoubleVectorImage;
using vigra::Shape3;
using vigra::Shape2;
using vigra::MultiArrayIndex;

void example_vectorialDist(const IntImage& image, const std::string& fname = "") {
    DoubleVectorImage distImage(Shape2(image.shape()));

    //vigra::ArrayVector<double> pitch(2); pitch[0] = 1.0; pitch[1] = 1.5;
    separableMultiVectorialDist(srcMultiArrayRange(image), destMultiArray(distImage), true);

    for(MultiArrayIndex i=0; i<image.shape(0); ++i) {
        for(MultiArrayIndex j=0; j<image.shape(1); ++j) {
            std::cout << "(";
            for(MultiArrayIndex k=0; k<2; ++k) {
                std::cout << std::setw(3) << distImage(i,j)[k];
                if(k<2-1) {
                    std::cout << " ";
                }
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
    }
    
    if(fname.size() > 0) {
        VectorialDistancePlotter plotter(image, distImage);
        plotter.plot(fname+"1.png"); 
        plotter.setModeToUnitvectors();
        plotter.plot(fname+"2.png"); 
    }
}

//void example_vectorialBoundaryDist(const IntImage& image, const std::string& saveAs = "") {

//    std::cout << "Input segmentation:" << std::endl;
//    for(MultiArrayIndex i=0; i<image.shape(0); ++i) {
//        for(MultiArrayIndex j=0; j<image.shape(1); ++j) {
//            std::cout << std::setw(3) << image(i,j);
//            std::cout << " ";
//        }
//        std::cout << std::endl;
//    }
    
//    std::pair<DoubleVectorImage, vigra::MultiArray<2, unsigned char> > p;
    
//    p = vigra::separableMultiVectorialBoundaryDist<2, int>(image);
    
//    DoubleVectorImage distImage = p.first;
//    vigra::MultiArray<2, unsigned char> cornersImage = p.second;
    
//    for(MultiArrayIndex i=0; i<image.shape(0); ++i) {
//        for(MultiArrayIndex j=0; j<image.shape(1); ++j) {
//            std::cout << "(";
//            for(MultiArrayIndex k=0; k<2; ++k) {
//                std::cout << std::setprecision(3) << std::setw(4) << distImage(i,j)[k];
//                if(k<2-1) {
//                    std::cout << " ";
//                }
//            }
//            std::cout << ") ";
//        }
//        std::cout << std::endl;
//    }
    
//    if(saveAs.size() > 0) {
//        VectorialDistancePlotter plotter(image, distImage);
//        plotter.setCornersImage(cornersImage);
//        plotter.plot(saveAs);
//    }
//}

void example_euclideanDist(const IntImage& image) {
    IntImage result(image.shape());
    separableMultiDistSquared(srcMultiArrayRange(image), destMultiArray(result), true); 
    for(MultiArrayIndex i=0; i<image.shape(0); ++i) {
        for(MultiArrayIndex j=0; j<image.shape(1); ++j) {
            std::cout << std::setw(3) << result(i,j);
            std::cout << " ";
        }
        std::cout << std::endl;
    }
}

IntImage tg(const IntImage& seg) {
    using vigra::MultiArrayIndex;
    
    IntImage a(Shape2(2*seg.shape(0)-1, 2*seg.shape(1)-1));
    
    //go over all 2-cells
    for(MultiArrayIndex i=0; i<a.shape(0); ++i) {
        if( (i%2) == 0) {
            for(MultiArrayIndex j=1; j<a.shape(1); j+=2) {
                if( seg(i/2, (j-1)/2) != seg(i/2, (j+1)/2) ) {
                    a(i,j) = 1;
                }
            }
        }
        else {
            for(MultiArrayIndex j=0; j<a.shape(1); j+=2) {
                if( seg((i-1)/2, j/2) != seg((i+1)/2, j/2) ) {
                    a(i,j) = 1;
                }
            }
        }
    }
    
    return a;
}

struct Divider {
    inline vigra::TinyVector<double, 3> operator()(const vigra::TinyVector<double, 3>& a,  double b) const {
        return a / static_cast<double>(b);
    }
};

void example_vectorialMaxDistance() {
    using namespace vigra;
    
    ImageImportInfo iii("seg_gimp.png");
    MultiArray<2, unsigned int> seg(Shape2(iii.width(), iii.height()));
    importImage(iii, destImage(seg));
    
    labelImageWithBackground(srcImageRange(seg), destImage(seg), false, 0); 
    
    MultiArray<3, unsigned int> seg3D(Shape3(seg.shape(0), seg.shape(1), 3));
    seg3D.bind<2>(0) = seg;
    seg3D.bind<2>(1) = seg;
    seg3D.bind<2>(2) = seg;
    
    MultiArray<3, unsigned char> tg(Shape3(2*seg3D.shape(0)-1, 2*seg3D.shape(1)-1, 2*seg3D.shape(2)-1));
    
    USETICTOC
    TIC
    regionVolumeToCrackEdgeVolume(seg3D, tg);
    std::cout << "regionVolumeToCrackEdgeVolume took " << TOCS << " seconds" << std::endl;
    
    MultiArray<3, TinyVector<double, 3> > distImage(seg3D.shape());
    if(true) {
        vectorialMaxDistanceBruteForce(tg, seg3D, distImage);
        {
            HDF5File f("test_maxDist.h5", HDF5File::New);
            f.write("vectorialMaxDistance", distImage);
            std::cout << "wrote test_maxDist.h5" << std::endl;
        }
    }
    else {
        HDF5File f(std::string("test_maxDist.h5"), HDF5File::OpenReadOnly);
        f.readAndResize("vectorialMaxDistance", distImage);
    }
    
    std::vector<QRgb> cmap = VectorialDistancePlotter::generateRandomColormap();
   
    for(int z=0; z<seg3D.shape(2); ++z) {
        const MultiArrayView<2, TinyVector<double, 3>, StridedArrayTag> slice = distImage.bind<2>(z);
        const MultiArrayView<2, unsigned int> segSlice = seg3D.bind<2>(z);
        MultiArray<2, TinyVector<double, 2> > vectors(slice.shape());
        for(int i=0; i<slice.shape(0); ++i) {
        for(int j=0; j<slice.shape(1); ++j) {
            vectors(i,j)[0] = slice(i,j)[0];
            vectors(i,j)[1] = slice(i,j)[1];
        }
        }
        std::stringstream fname; fname << "tttt_" << std::setw(5) << std::setfill('0') << z << ".png";
        MultiArray<2, int> segSliceInt(segSlice);
        VectorialDistancePlotter plotter(segSliceInt, vectors);
        plotter.setColormap(cmap);
        plotter.plot(fname.str());
        
        plotter.setModeToUnitvectors();
        std::stringstream fname2; fname2 << "uuuu_" << std::setw(5) << std::setfill('0') << z << ".png";
        plotter.plot(fname2.str());
        
    }
    
}

void example2d_seg_small_tiff() {
    using namespace vigra;
    
    ImageImportInfo iii("seg_small.tiff");
    IntImage seg(Shape2(iii.width(), iii.height()));
    importImage(iii, destImage(seg));
    
    IntImage crackEdgeImage(Shape2(2*seg.shape(0)-1, 2*seg.shape(1)-1));
    regionImageToCrackEdgeImage(srcImageRange(seg), destImage(crackEdgeImage), false);
    using namespace functor;
    transformImage(srcImageRange(crackEdgeImage), destImage(crackEdgeImage), 
                            ifThenElse(Arg1() == Param(0), Param(1.0), Param(0.0))); 

    DoubleVectorImage distImage(Shape2(crackEdgeImage.shape(0), crackEdgeImage.shape(1)));
    separableMultiVectorialDist(srcMultiArrayRange(crackEdgeImage), destMultiArray(distImage), true);
    
    DoubleVectorImage distImage2(Shape2(seg.shape(0), seg.shape(1)));
    for(MultiArrayIndex i=0; i<crackEdgeImage.shape(0); i+=2) {
        for(MultiArrayIndex j=0; j<crackEdgeImage.shape(1); j+=2) {
            distImage2(i/2, j/2) = distImage(i,j)/2.0;
        }
    }
    
    for(MultiArrayIndex i=0; i<distImage2.shape(0); ++i) {
    for(MultiArrayIndex j=0; j<distImage2.shape(1); ++j) {
        for(MultiArrayIndex k=0; k<2; ++k) {
            if(isnan(distImage2(i,j)[k])) {
                throw std::runtime_error("oh no nan");
            }
        }
        double d = distImage2(i,j).squaredMagnitude();
        if(d == 0) {
            std::stringstream ss; ss << "at (" << i << "," << j << ") = " << distImage2(i,j) << std::endl;
            throw std::runtime_error(ss.str());
        }
    }
    }
    
    VectorialDistancePlotter plotter(seg, distImage2);
    plotter.plot("/tmp/D.png");
    
    example_vectorialDist(crackEdgeImage, "/tmp/A");
//    example_vectorialBoundaryDist(seg, "/tmp/B");
}

int main() {
    using namespace vigra;
    using namespace vigra::functor;

    //example_vectorialMaxDistance();
    //example_euclideanDist(image);
    //example2d_seg_small_tiff()
    
    //return 0;
    
    //HDF5File f(std::string("seg3d_small.h5"), HDF5File::OpenReadOnly);
    HDF5File f(std::string("cylinder/20_cylinder.h5"), HDF5File::OpenReadOnly);
    MultiArray<3, unsigned int> seg3D;
    f.readAndResize("seg", seg3D);
    
    typedef MultiArrayShape<3>::type Diff3;
    std::cout << seg3D.shape() << "xxxxx" << std::endl;
    //seg3D = seg3D.subarray(Diff3(0,0,0), Diff3(10,150,100));
    MultiArray<3, unsigned int> seg3D_copy(seg3D);
    labelVolumeWithBackground(srcMultiArrayRange(seg3D_copy), destMultiArray(seg3D), NeighborCode3DSix(), 0);
    
    MultiArray<3, unsigned char> tg(Shape3(2*seg3D.shape(0)-1, 2*seg3D.shape(1)-1, 2*seg3D.shape(2)-1));
    
    //compute topological grid
    USETICTOC
    TIC
    regionVolumeToCrackEdgeVolume(seg3D, tg);
    std::cout << "regionVolumeToCrackEdgeVolume took " << TOCS << " seconds" << std::endl;
    
    //NEW
    /*
    MultiArray<3, TinyVector<double, 3> > maxVec(tg.shape());
    vectorialMaxDistanceBruteForce(tg, seg3D, maxVec);
    {
        HDF5File f("maxDist.h5", HDF5File::New);
        f.write("vectorialMaxDistance", maxVec);
        std::cout << "wrote maxDist.h5" << std::endl;
    }
    */
    //NEW
    
    //compute vectorial distance transform (3D)
    MultiArray<3, TinyVector<double,3> > distImage(Shape3(seg3D.shape()));
    {
        TIC
        MultiArray<3, TinyVector<double,3> > distImageTG(Shape3(tg.shape()));
        separableMultiVectorialDist(srcMultiArrayRange(tg), destMultiArray(distImageTG), true);
    
        //project back to original shape
        for(MultiArrayIndex i=0; i<tg.shape(0); i+=2) {
            for(MultiArrayIndex j=0; j<tg.shape(1); j+=2) {
                for(MultiArrayIndex k=0; k<tg.shape(2); k+=2) {
                    distImage(i/2, j/2, k/2) = distImageTG(i,j,k)/2.0;
                }
            }
        }
        std::cout << "computing vectorial distance took " << TOCS << " seconds" << std::endl;
    }
    {
        HDF5File f("vectorialDist.h5", HDF5File::New);
        f.write("dist", distImage); std::cout << "wrote vectorialDist.h5" << std::endl;
        f.write("seg", seg3D);
    }
    
    //split distImage into vector length and vector direction
    TIC
    MultiArray<3, double> magnitude(distImage.shape());
    MultiArray<3, TinyVector<double, 3> > directions(distImage.shape());
    transformMultiArray(srcMultiArrayRange(distImage), destMultiArray(magnitude), norm(Arg1()));
    combineTwoMultiArrays(srcMultiArrayRange(distImage), srcMultiArray(magnitude),
                            destMultiArray(directions), Divider() );
    std::cout << "split into magnitude/direction took " << TOCS << " seconds" << std::endl;
    
    MultiArray<3, TinyVector<double, 3> > gradients[3];
    for(int i=0; i<3; ++i) {
        TIC
        gradients[i].reshape(distImage.shape());
        gaussianGradientMultiArray(srcMultiArrayRange(directions.bindElementChannel(i)), destMultiArray(gradients[i]), 3.0);
        std::cout << "gaussianGradientMultiArray (channel=" << i << ") took " << TOCS << " seconds" << std::endl;
    }
    
    {
        HDF5File f("vectorialDist.h5", HDF5File::Open);
        for(int i=0; i<3; ++i) {
            std::stringstream ss; ss << "gradient-" << i;
            f.write(ss.str(), gradients[i]); 
        }
    }
    
    /// CURL ///
    
    MultiArray<3, TinyVector<double, 3> > curl(Shape3(seg3D.shape(0), seg3D.shape(1), seg3D.shape(2)));
    for(MultiArrayIndex I=0; I<tg.shape(0); I+=2) {
    for(MultiArrayIndex J=0; J<tg.shape(1); J+=2) {
    for(MultiArrayIndex K=0; K<tg.shape(2); K+=2) {
        curl(I/2,J/2,K/2) = TinyVector<double, 3>(
            gradients[2](I/2,J/2,K/2)[1] /*dF_z/dy */ - gradients[1](I/2,J/2,K/2)[2] /*dF_y/dz */,
            gradients[0](I/2,J/2,K/2)[2] /*dF_x/dz */ - gradients[2](I/2,J/2,K/2)[0] /*dF_z/dx */,
            gradients[1](I/2,J/2,K/2)[0] /*dF_y/dx */ - gradients[0](I/2,J/2,K/2)[1] /*dF_x/dy */
        );
    }
    }
    }
    {
        HDF5File f("vectorialDist.h5", HDF5File::Open);
        f.write("curl", curl); std::cout << "wrote curl.h5" << std::endl;
    }
    
    /// DIV ///
    MultiArray<3, double> div(Shape3(seg3D.shape(0), seg3D.shape(1), seg3D.shape(2)));
    for(MultiArrayIndex I=0; I<tg.shape(0); I+=2) {
    for(MultiArrayIndex J=0; J<tg.shape(1); J+=2) {
    for(MultiArrayIndex K=0; K<tg.shape(2); K+=2) {
        div(I/2,J/2,K/2) = gradients[0](I/2,J/2,K/2)[0] + gradients[1](I/2,J/2,K/2)[1] + gradients[2](I/2,J/2,K/2)[2]; 
    }
    }
    }
    {
        HDF5File f("vectorialDist.h5", HDF5File::Open);
        f.write("div", div); std::cout << "wrote div.h5" << std::endl;
    }
    
    
    MultiArray<4, double> ew(Shape4(magnitude.shape(0), magnitude.shape(1), magnitude.shape(2), 3));
    MultiArray<5, double> ev(Shape5(magnitude.shape(0), magnitude.shape(1), magnitude.shape(2), 3, 3));
    
    TIC
    MultiArray<2, double> m(Shape2(3,3));
    for(MultiArrayIndex K=0; K<magnitude.shape(2); ++K) {
    for(MultiArrayIndex J=0; J<magnitude.shape(1); ++J) {
    for(MultiArrayIndex I=0; I<magnitude.shape(0); ++I) {
        // A^T * A
        const double a = gradients[0](I,J,K)[0];
        const double b = gradients[0](I,J,K)[1];
        const double c = gradients[0](I,J,K)[2];
        const double d = gradients[1](I,J,K)[0];
        const double e = gradients[1](I,J,K)[1];
        const double f = gradients[1](I,J,K)[2];
        const double g = gradients[2](I,J,K)[0];
        const double h = gradients[2](I,J,K)[1];
        const double i = gradients[2](I,J,K)[2];
        
        // A * A^T
        /*
        const double a = gradients[0](I,J,K)[0];
        const double b = gradients[1](I,J,K)[0];
        const double c = gradients[2](I,J,K)[0];
        const double d = gradients[0](I,J,K)[1];
        const double e = gradients[1](I,J,K)[1];
        const double f = gradients[2](I,J,K)[1];
        const double g = gradients[0](I,J,K)[2];
        const double h = gradients[1](I,J,K)[2];
        const double i = gradients[2](I,J,K)[2];
        */
        
        m(0,0) = sq(a)+sq(b)+sq(c);
        m(0,1) = a*d+b*e+c*f;
        m(0,2) = a*g+b*h+c*i;
        m(1,1) = sq(d)+sq(e)+sq(f);
        m(1,2) = d*g+e*h+f*i;
        m(2,2) = sq(g)+sq(h)+sq(i);
        
        m(1,0) = m(0,1);
        m(2,0) = m(0,2);
        m(2,1) = m(1,2);

        MultiArray<2, double> ew_col(Shape2(3,1));
        
        linalg::symmetricEigensystem(m, ew_col, m);
        
        for(MultiArrayIndex u=0; u<3; ++u) {
            ew(I,J,K,u) = ew_col(u,0);
        }
        ev.bind<0>(I).bind<0>(J).bind<0>(K) = m;
    }
    }
    }
    std::cout << "computing ev/ew took " << TOCS << " seconds" << std::endl;
    
    {
        HDF5File f("vectorialDist.h5", HDF5File::Open);
        f.write("ew", ew);
        f.write("ev", ev);
        f.write("magnitude", magnitude);
    }
    
    std::vector<QRgb> cmap = VectorialDistancePlotter::generateRandomColormap();
    
    for(int z=0; z<seg3D.shape(2); ++z) {
        const MultiArrayView<2, TinyVector<double, 3>, StridedArrayTag> slice = distImage.bind<2>(z);
        const MultiArrayView<2, unsigned int> segSlice = seg3D.bind<2>(z);
        MultiArray<2, TinyVector<double, 2> > vectors(slice.shape());
        for(int i=0; i<slice.shape(0); ++i) {
        for(int j=0; j<slice.shape(1); ++j) {
            vectors(i,j)[0] = slice(i,j)[0];
            vectors(i,j)[1] = slice(i,j)[1];
        }
        }
        std::stringstream fname; fname << std::setw(5) << std::setfill('0') << z << ".png";
        MultiArray<2, int> segSliceInt(segSlice);
        VectorialDistancePlotter plotter(segSliceInt, vectors);
        plotter.setColormap(cmap);
        plotter.plot(fname.str());
        
    }
       
}
