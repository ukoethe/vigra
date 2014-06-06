#include "vis_vectorial_dist.hxx"

#include <iostream>

#include <vigra/linear_algebra.hxx>
#include <vigra/random.hxx>

#include <QPainter>
#include <QRectF>
#include <QLineF>
#include <QColor>
#include <qdebug.h>


typedef vigra::TinyVector<double, 2> Vec;

void drawArrow(QPainter& painter, const Vec& start, const Vec& offset, const QColor& color, double size=1.0) {
    double arrowH = size;   //height of arrow tip
    double arrowW = 0.5*size;   //width of arrow tip
    double tailR  = size/4.0; //radius of tail circle
    
    QPen pen = QPen(color); 
    painter.setPen(pen);
    painter.setBrush(color);
   
    Vec end = start+offset;
    
    QLineF line(QPointF(start[0], start[1]), QPointF(end[0], end[1]));

    if(line.length() == 0.0) {
        return;
    }

    Vec b = end-start;
    Vec B = b / std::sqrt(b.squaredMagnitude()); //unit vector in arrow direction
    Vec Bp = Vec(-B[1], B[0]);           //unit vector perpendicular to arrow direction
    //d,e: points defining arrow head
    Vec c = end - arrowH*B;
    Vec d = c+arrowW/2.0*Bp;
    Vec e = c-arrowW/2.0*Bp;
   
    //draw arrow line
    painter.drawLine(line);
    
    //draw arrow head
    QPolygonF poly;
    poly << QPointF(d[0],d[1]) << QPointF(e[0],e[1]) << QPointF(end[0], end[1]);
    painter.drawPolygon(poly, Qt::WindingFill);
   
    //draw arrow tail
    painter.drawEllipse(QRectF(line.p1().x()-tailR, line.p1().y()-tailR, 2.0*tailR, 2.0*tailR));
}

VectorialDistancePlotter::VectorialDistancePlotter(
    const vigra::MultiArrayView<2, int>& img,
    const vigra::MultiArrayView<2, vigra::TinyVector<double, 2> >& vectors
) : img_(img), vectors_(vectors), mode_(ModeVectors)
{
    randomizeColormap();
}
    
void VectorialDistancePlotter::setCornersImage(
    const vigra::MultiArrayView<2, unsigned char>& corners
) {
    corners_ = corners;
}
    
void VectorialDistancePlotter::setModeToUnitvectors() {
    mode_ = ModeUnitvectors;
}
void VectorialDistancePlotter::setModeToVectors() {
    mode_ = ModeVectors;
}

std::vector<QRgb> VectorialDistancePlotter::generateRandomColormap()
{
    std::vector<QRgb> cmap(256);
    cmap[0] =  QColor(255,255,255).rgb();
    cmap[1] = QColor(128,128,128).rgb();
    for(int i=2; i<256; ++i) {
        cmap[i] = QColor(
            vigra::RandomTT800::global().uniformInt(256),
            vigra::RandomTT800::global().uniformInt(256),
            vigra::RandomTT800::global().uniformInt(256)
        ).rgb();
    }
    return cmap;
}

void VectorialDistancePlotter::randomizeColormap() {
    colormap_ = generateRandomColormap();
}

void VectorialDistancePlotter::plot(const std::string& fname) {
    int margin = 1;
    
    vigra_precondition(img_.shape(0) == vectors_.shape(0) && img_.shape(1) == vectors_.shape(1),
                       "img_, vectors_ disagree on shape");
    
    std::cout << "exporting vectorial dist image " << fname << " for image of size " << img_.shape(0) << "x" << img_.shape(1) << std::endl;
    
    using vigra::MultiArrayIndex;
    
    double N = 30.0;
   
    QImage canvas(N*(2*margin+img_.shape(0)), N*(2*margin+img_.shape(1)), QImage::Format_ARGB32);
    QPainter p(&canvas);
    p.setRenderHints(QPainter::Antialiasing);
   
    canvas.fill(Qt::white);
    if(true) {
        QImage qimg(img_.shape(0), img_.shape(1), QImage::Format_ARGB32);
        for(MultiArrayIndex i=0; i<img_.shape(0); ++i) {
            for(MultiArrayIndex j=0; j<img_.shape(1); ++j) {
                qimg.setPixel(i,j, colormap_[img_(i,j)%colormap_.size()]);
            }
        }
        p.drawImage(QRectF(N*margin,N*margin,N*qimg.width(),N*qimg.height()), qimg);
    }
    
    //draw grid
    p.setPen(QPen(Qt::blue));
    
    for(MultiArrayIndex i=0; i<img_.shape(0)+1; ++i) {
        p.drawLine(QPointF(N*(i+margin),N*margin), QPointF(N*(i+margin), N*(img_.shape(1)+margin)));
    }
    for(MultiArrayIndex j=0; j<img_.shape(1)+1; ++j) {
        p.drawLine(QPointF(N*margin,N*(j+margin)), QPointF(N*(img_.shape(0)+margin), N*(j+margin)));
    }
   
    p.setPen(QPen(Qt::green, 3));
    for(MultiArrayIndex i=0; i<img_.shape(0); ++i) {
        for(MultiArrayIndex j=0; j<img_.shape(1); ++j) {
            if(i+1 < img_.shape(0) && img_(i,j) != img_(i+1,j)) {
                p.drawLine(QPointF(N*(i+1+margin),N*(j+margin) ), QPointF(N*(i+1+margin), N*(j+1+margin)));
            }
            if(j+1 < img_.shape(1) && img_(i,j) != img_(i,j+1)) {
                p.drawLine(QPointF(N*(i+margin),N*(j+1+margin)), QPointF(N*(i+1+margin), N*(j+1+margin)));
            }
        }
    }
    
    if(corners_.size() > 0) {
        for(MultiArrayIndex i=1; i<corners_.shape(0); ++i) {
            for(MultiArrayIndex j=1; j<corners_.shape(1); ++j) {
                if(corners_(i-1,j-1) == 1) {
                    double sz = N/8;
                    p.drawEllipse( QRectF( N*(i-1+margin)-sz, N*(j-1+margin)-sz, 2*sz, 2*sz) );
                }
            }
        }
    }

    QColor arrowColor(255,0,0);
    if(mode_ == ModeVectors) {
        for(MultiArrayIndex i=0; i<img_.shape(0); ++i) {
            for(MultiArrayIndex j=0; j<img_.shape(1); ++j) {
                const double I = vectors_(i,j)[0];
                const double J = vectors_(i,j)[1];
                drawArrow(p, Vec(N*(i+0.5+margin),N*(j+0.5+margin)), Vec(N*I,N*J), arrowColor, N/4.0);
            }
        }
    }
    else {
        for(MultiArrayIndex i=0; i<img_.shape(0); ++i) {
            for(MultiArrayIndex j=0; j<img_.shape(1); ++j) {
                Vec v = vectors_(i,j);
                double mag = v.squaredMagnitude();
                if(mag == 0.0) {
                    continue;
                }
                v /= std::sqrt(mag);
                v *= 0.666;
                drawArrow(p, Vec(N*(i+0.5+margin-0.5*v[0]),N*(j+0.5+margin-0.5*v[1])), N*v, arrowColor, N/4.0);
            }
        }
    }
    canvas.save(QString::fromStdString(fname));
}

