#include <vigra/tinyvector.hxx>
#include <vigra/linear_algebra.hxx>
#include <cmath>
#include <vigra/hdf5impex.hxx>
#include <vigra/random.hxx>

using namespace vigra;
typedef vigra::TinyVector<double, 3> Vec;

class Cylinder {
    public:
    Cylinder(const Vec& B, const Vec& N, double L, double R) : b(B), n(N), l(L), r(R) {
        n /= sqrt(squaredNorm(n));
    }
    Cylinder(const Cylinder& o) : b(o.b), n(o.n), l(o.l), r(o.r) {}
    Vec b;
    Vec n;
    double l;
    double r;

    Vec stop() const {
        return b+l*n;
    }
};

int main() {
    int shape[] = {100,200,300};

    RandomNumberGenerator<RandomMT19937> rng = RandomNumberGenerator<RandomMT19937>::global();

    std::vector<Cylinder> cylinders;

    Cylinder c1(Cylinder(Vec(0,0,0), Vec(1,2,1) , 100, 5));
    cylinders.push_back(c1);

    Cylinder c2(c1.stop(), Vec(1,4,1), 100, 5);
    cylinders.push_back(c2);

    Cylinder c3(c1.stop(), Vec(1,1,1), 100, 5);
    cylinders.push_back(c3);

    typedef MultiArrayShape<3>::type VolShape;
    MultiArray<3, double> vol(VolShape(shape[0], shape[1], shape[2]));

    double noise = 30.0;

    std::vector<double> dists(cylinders.size());

    for(int i=0; i<shape[0]; ++i) {
    for(int j=0; j<shape[1]; ++j) {
    for(int k=0; k<shape[2]; ++k) {
        Vec p(i,j,k);

        for(int c=0; c<cylinders.size(); ++c) {
            const Cylinder& cyl = cylinders[c];
            double t = vigra::dot(cyl.n, p-cyl.b);

            Vec x = p- ( cyl.b+t*cyl.n );
            double y = sqrt(vigra::squaredNorm(x))-cyl.r;

            double D;
            if(t < 0 || t > cyl.l) {
                double d = (t<0) ? std::fabs(t) : t-d;
                D = sqrt(sq(y) + sq(d));
            }
            else {
                D = y;
            }

            dists[c] = (D < 0) ? 0 : D;
        }

        vol(i,j,k) = *std::min_element(dists.begin(), dists.end());

        /*
        if(vol(i,j,k) > 0 && vol(i,j,k) < noise) {
            vol(i,j,k) += std::max(0.0, (rng.uniform()-0.5)*noise);
        }
        */

    }
    }
    }

    vigra::HDF5File f("cylinder.h5", vigra::HDF5File::Open);
    f.write("seg", vol);
}
