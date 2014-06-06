#ifndef VIS_VECTORIAL_DIST__HXX
#define VIS_VECTORIAL_DIST__HXX

#include <vigra/multi_array.hxx>
#include <QColor>

class VectorialDistancePlotter {
    enum Mode {
        ModeUnitvectors,
        ModeVectors
    };
    
    public:
    VectorialDistancePlotter(
        const vigra::MultiArrayView<2, int>& img,
        const vigra::MultiArrayView<2, vigra::TinyVector<double, 2> >& vectors
    );
    
    static std::vector<QRgb> generateRandomColormap();
    
    void setCornersImage(const vigra::MultiArrayView<2, unsigned char>& corners);
    
    void setModeToUnitvectors();
    void setModeToVectors();
    
    void randomizeColormap();
    void setColormap(const std::vector<QRgb>& colormap) { colormap_ = colormap; }
   
    void plot(const std::string& fname);
    
    private:
    vigra::MultiArrayView<2, int> img_;
    vigra::MultiArrayView<2, vigra::TinyVector<double, 2> > vectors_;
    vigra::MultiArrayView<2, unsigned char> corners_;
    Mode mode_; 
    std::vector<QRgb> colormap_;
};

#endif /* VIS_VECTORIAL_DIST__HXX */