/************************************************************************/
/*                                                                      */
/*                 Copyright 2014 by Benjamin Seppke                    */      
/*                                                                      */
/************************************************************************/

#include <iostream>

#include <vigra/unittest.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/convolution.hxx>
#include <vigra/impex.hxx>

#include <vigra/affine_registration.hxx>
#include <vigra/projective_registration.hxx>
#include <vigra/polynomial_registration.hxx>
#include <vigra/rbf_registration.hxx>

using namespace vigra;

static int pointdata[] = {
/*  x_1     y_1     x_2     y_2 */
    362  ,  379  ,  783  ,  658,
    362  ,  250  ,  786  ,  582,
    378  ,  250  ,  797  ,  581,
    470  ,  250  ,  859  ,  581,
    504  ,  251  ,  879  ,  583,
    674  ,  249  ,  991  ,  582,
    858  ,  255  , 1121  ,  581,
    857  ,  284  , 1121  ,  605,
    913  ,  285  , 1161  ,  606,
    913  ,  311  , 1161  ,  620,
    913  ,  376  , 1159  ,  657,
    824  ,  376  , 1093  ,  657,
    785  ,  377  , 1067  ,  657,
    646  ,  377  ,  972  ,  658,
    535  ,  379  ,  896  ,  659,
    535  ,  393  ,  896  ,  668,
    377  ,  392  ,  793  ,  667,
    377  ,  367  ,  793  ,  655,
    535  ,  367  ,  896  ,  655,
    377  ,  334  ,  794  ,  638,
    377  ,  281  ,  796  ,  607,
    504  ,  281  ,  878  ,  607,
    503  ,  334  ,  877  ,  638,
    674  ,  277  ,  991  ,  602,
    754  ,  285  , 1049  ,  611,
    672  ,  365  ,  990  ,  654,
    659  ,  352  ,  981  ,  648,
    659  ,  339  ,  981  ,  642,
    654  ,  339  ,  978  ,  642,
    654  ,  352  ,  978  ,  648,
    738  ,  353  , 1036  ,  645,
    790  ,  363  , 1071  ,  648,
    774  ,  365  , 1060  ,  651,
    806  ,  364  , 1081  ,  649,
    807  ,  316  , 1083  ,  623,
    790  ,  316  , 1071  ,  623,
    790  ,  330  , 1071  ,  631,
    806  ,  331  , 1081  ,  631,
    822  ,  350  , 1093  ,  644,
    835  ,  350  , 1102  ,  644,
    835  ,  364  , 1102  ,  651,
    902  ,  364  , 1152  ,  651,
    902  ,  350  , 1152  ,  643,
    889  ,  350  , 1144  ,  643,
    618  ,  103  ,  960  ,  448,
    505  ,  104  ,  883  ,  448,
    505  ,   53  ,  883  ,  392,
    452  ,   10  ,  853  ,  339,
    273  ,   10  ,  732  ,  340,
    269  ,  212  ,  724  ,  545,
    503  ,  215  ,  879  ,  545,
    271  ,   83  ,  728  ,  428,
    504  ,   84  ,  883  ,  428,
     17  ,  424  ,  558  ,  681,
     19  ,  378  ,  559  ,  662,
    155  ,  381  ,  645  ,  663,
    174  ,  382  ,  659  ,  663,
    317  ,  385  ,  750  ,  663,
    317  ,  430  ,  750  ,  682,
    174  ,  428  ,  658  ,  682,
    154  ,  427  ,  644  ,  682,
    758  ,  255  , 1052  ,  582,
    461  ,  368  ,  848  ,  656,
    461  ,  392  ,  848  ,  668,
    461  ,  487  ,  848  ,  702,
    540  ,  492  ,  897  ,  704
};

static int point_count = 66;
static int pointdata_size = point_count*4;
static double test_epsilon = 1.0e-5;
static double rbf_test_epsilon = 1.0e-2;

static std::vector<TinyVector<double,2> > srcPoints()
{
    std::vector<TinyVector<double,2> >  result(point_count);
    
    for (int i=0; i<pointdata_size; i+=4)
    {
        result[i/4][0] = pointdata[i  ];
        result[i/4][1] = pointdata[i+1];
    }
    
    return result;
}

static std::vector<TinyVector<double,2> > destPoints()
{
    std::vector<TinyVector<double,2> >  result(point_count);
    
    for (int i=0; i<pointdata_size; i+=4)
    {
        result[i/4][0] = pointdata[i+2];
        result[i/4][1] = pointdata[i+3];
    }  
    
    return result;
}

void printMatrix(const Matrix<double> & m)
{    
    for (int i=0; i<m.rowCount(); ++i)
    {
        for (int j=0; j<m.columnCount(); ++j)
        {
            printf("m(%d,%d) = %10.15f;\n", i , j , m(i,j));
        }        
    }
}

void shouldEqualToleranceMatrices(const Matrix<double> & m1, const Matrix<double> & m2, double eps=test_epsilon)
{
    should(m1.rowCount() == m2.rowCount());
    should(m1.columnCount() == m2.columnCount());
    
    for (int i=0; i<m1.rowCount(); ++i)
    {
        for (int j=0; j<m1.columnCount(); ++j)
        {
            shouldEqualTolerance(m1(i,j), m2(i,j), eps);
        }        
    }
}

struct EstimateGlobalRotationTranslationTest
{
    
    unsigned int size;
    unsigned int fill_size_h;
    unsigned int fill_size_v;
    
    BImage s_img;
    BImage d_img;
    
    EstimateGlobalRotationTranslationTest()
    : size(200),
      fill_size_h(20),
      fill_size_v(50)
    {
        s_img.resize(size,size);
        d_img.resize(size,size);
        
        for(int y=(size-fill_size_v)/2; y<(size+fill_size_v)/2; ++y)
            for(int x=(size-fill_size_h)/2; x<(size+fill_size_h)/2; ++x)
                s_img(x,y) = 255;
        
        gaussianSmoothing(srcImageRange(s_img),destImage(s_img), 1.0);
    }
    void testInit()
    {            
        TinyVector<double,2> trans(-25,30);
        double theta = 24.75;
        
        Matrix<double> t_m = translationMatrix2D(trans);
        Matrix<double> r_m = rotationMatrix2DDegrees(theta, TinyVector<double,2>(size,size)/2.0);
        
        //Normaly, we would write t*r. since the matrix is reversed to transform I2 based on tranformed I1,
        // we have to change our point of view and reverse the order here to.
        Matrix<double> known_transformation =  r_m * t_m;
        std::cerr << "Known transformation: " << known_transformation << "\n"; 
        
        Matrix<double> estimated_transformation;
        
        affineWarpImage(SplineImageView<2,double>(srcImageRange(s_img)), destImageRange(d_img), known_transformation);
        
        exportImage(srcImageRange(s_img), ImageExportInfo("img_quad_1.png"));
        exportImage(srcImageRange(d_img), ImageExportInfo("img_quad_2.png"));
        
        double corr_rot, corr_trans;
        
        estimateGlobalRotationTranslation(srcImageRange(s_img), 
        								  destImageRange(d_img), 
        								  estimated_transformation, 
        								  corr_rot, 
        								  corr_trans);
        
        std::cerr << "Estimated transformation: \n" << estimated_transformation << "(rot-certainty: " << corr_rot << ", trans-certainty: " << corr_trans << ")\n\n\n";
        
        d_img=0;
        affineWarpImage(SplineImageView<2,double>(srcImageRange(s_img)), destImageRange(d_img), estimated_transformation);
        exportImage(srcImageRange(d_img), ImageExportInfo("img_quad_corr_1.png"));
        
        shouldEqualToleranceMatrices(known_transformation, estimated_transformation, 0.1);
        should(corr_rot   > 0.95);
        should(corr_trans > 0.99);
        
    }
};
struct EstimateGlobalRotationTranslationRealImageTest
{    
    BImage s_img;
    BImage d_img;
    
    double rotation;
    TinyVector<double,2> translation;
    
    Diff2D border;
    
    EstimateGlobalRotationTranslationRealImageTest()
    : rotation(24.25),
      translation(10,-30),
      border(400,200)
    {        
		ImageImportInfo info1("nuernberg-1995.png");
		s_img.resize(info1.width(), info1.height());
		d_img.resize(info1.width(), info1.height());
		importImage(info1, destImage(s_img));
    }
    
    void testInit()
    {
        Matrix<double> r_m = rotationMatrix2DDegrees(rotation, TinyVector<double,2>(s_img.width(), s_img.height())/2.0);
        Matrix<double> t_m = translationMatrix2D(translation);
        
        //Normaly, we would write t*r. since the matrix is reversed to transform I2 based on transformed I1,
        // we have to change our point of view and reverse the order here to.
        Matrix<double> known_transformation =  r_m * t_m;
        std::cerr << "Known transformation: \n" << known_transformation << "\n"; 
        
        affineWarpImage(SplineImageView<2,double>(srcImageRange(s_img)), destImageRange(d_img), known_transformation);
        
        exportImage(srcIterRange(s_img.upperLeft()+border, s_img.lowerRight()-border, s_img.accessor()),
                           ImageExportInfo("img_nb_1.png"));
        exportImage(srcIterRange(d_img.upperLeft()+border, d_img.lowerRight()-border, d_img.accessor()),
                           ImageExportInfo("img_nb_2.png"));
        
        Matrix<double> estimated_transformation;
        
        double corr_rot, corr_trans;
        
        
        estimateGlobalRotationTranslation(srcImageRange(s_img),
                                          srcImageRange(d_img),
                                          estimated_transformation,
                                          corr_rot, 
                                          corr_trans,
                                          border);
        
        std::cerr << "Estimated transformation: \n" << estimated_transformation << "(rot-certainty: " << corr_rot << ", trans-certainty: " << corr_trans << ")\n\n\n";
        
        d_img=0;
        affineWarpImage(SplineImageView<2,double>(srcImageRange(s_img)),
                        destImageRange(d_img), 
                        estimated_transformation);
        exportImage(srcIterRange(d_img.upperLeft()+border, d_img.lowerRight()-border, d_img.accessor()),
                    ImageExportInfo("img_nb_corr_1.png"));
        
        shouldEqualToleranceMatrices(known_transformation, estimated_transformation, 0.1);
        should(corr_rot   > 0.85);
        should(corr_trans > 0.95);
    }
};

struct EstimateGlobalRotationTranslationTestSuite
: public test_suite
{
    EstimateGlobalRotationTranslationTestSuite()
    : test_suite("EstimateGlobalRotationTranslationTestSuite")
    {
        add( testCase( &EstimateGlobalRotationTranslationTest::testInit));
        add( testCase( &EstimateGlobalRotationTranslationRealImageTest::testInit));
    }
};

struct ProjectiveIdentityTest
{
    std::vector<TinyVector<double,2> > s_points;
    
    ProjectiveIdentityTest()
    : s_points(srcPoints())
    {
    }
    
    void testInit()
    {
        /**
         * First test: If point sets are equal -> identity matrix should be the result!
         */
        Matrix<double> identity = projectiveMatrix2DFromCorrespondingPoints(s_points.begin(), s_points.end(), s_points.begin());
        
        shouldEqualToleranceMatrices(identity, linalg::identityMatrix<double>(3), test_epsilon);
    }
    
};

struct ProjectiveRegistrationTest
{
    BImage s_img;
    BImage d_img;
    
    std::vector<TinyVector<double,2> > s_points;
    std::vector<TinyVector<double,2> > d_points;
    
    ProjectiveRegistrationTest()
    : s_points(srcPoints()),
      d_points(destPoints())
    {
		ImageImportInfo info1("nuernberg-1991.png");
		s_img.resize(info1.width(), info1.height());
		importImage(info1, destImage(s_img));
		
		ImageImportInfo info2("nuernberg-1995.png");
		d_img.resize(info2.width(), info2.height());
		importImage(info2, destImage(d_img));
    }
    
    void testInit()
    {        
        /**
         * Test with well-known point sets and a known result matrix
         */
        Matrix<double> proj = projectiveMatrix2DFromCorrespondingPoints(s_points.begin(), s_points.end(), d_points.begin());
        
        /** 
         * Estimated result:
         *
         *  0.7908564596 -0.2771243619 -246.1575689662 
         * -0.0156283750  0.6189887443 -195.9746959403 
         * -0.0000258457 -0.0006847562    1.0000000000 
         */
        
        Matrix<double> reference(3,3);
        reference(0,0) =  0.7908564596; reference(0,1) = -0.2771243619; reference(0,2) = -246.1575689662; 
        reference(1,0) = -0.0156283750; reference(1,1) =  0.6189887443; reference(1,2) = -195.9746959403; 
        reference(2,0) = -0.0000258457; reference(2,1) = -0.0006847562; reference(2,2) =    1.000000; 
        
        shouldEqualToleranceMatrices(proj, reference, test_epsilon);
        
        /**
         * visual interpretation by means of the warped image:
         */
        
         BImage temp = d_img;
         projectiveWarpImage(SplineImageView<2,unsigned char>(srcImageRange(s_img)), destImageRange(temp), proj);
         exportImage(srcImageRange(temp), ImageExportInfo("res-proj.png"));
         
    }
    
};
struct ProjectiveRegistrationTestSuite
: public test_suite
{
    ProjectiveRegistrationTestSuite()
    : test_suite("ProjectiveRegistrationTestSuite")
    {
        add( testCase( &ProjectiveIdentityTest::testInit));
        add( testCase( &ProjectiveRegistrationTest::testInit));
   }
};
struct PolynomialIdentityTest
{
    std::vector<TinyVector<double,2> > s_points;
    
    PolynomialIdentityTest()
    : s_points(srcPoints())
    {
    }
    
    void testInit()
    {
        /**
         * First test: If point sets are equal -> identity w.r.t polynom matrix representation should be the result!
         */
        Matrix<double> identity = polynomialMatrix2DFromCorrespondingPoints<3>(s_points.begin(), s_points.end(), s_points.begin());
        
        /** 
         * Estimated result:
         *
         * Simple polygon: x -> (0 + 1*x + 0*y + 0*x^2 + 0*x*y + 0*y^2 + 0*x^3 + 0*x^2*y + 0*x*y^2 + 0*y^3)
         *                 y -> (0 + 0*x + 1*y + 0*x^2 + 0*x*y + 0*y^2 + 0*x^3 + 0*x^2*y + 0*x*y^2 + 0*y^3)
         *
         * In matrix notation: [0.00, 0.00]
         *                     [1.00, 0.00] x^1 , y^0  
         *                     [0.00, 1.00] x^0 , y^1 
         *                     [0.00, 0.00] x^2 , y^0  
         *                     [0.00, 0.00] x^1 , y^1 
         *                     [0.00, 0.00] x^0 , y^2  
         *                     [0.00, 0.00] x^3 , y^0  
         *                     [0.00, 0.00] x^2 , y^1  
         *                     [0.00, 0.00] x^1 , y^2  
         *                     [0.00, 0.00] x^0 , y^3 
         */
         
        Matrix<double> reference(10,2, 0.0);
        reference(1,0) = 1.0; reference(2,1) = 1.0;
        shouldEqualToleranceMatrices(identity, reference, test_epsilon);
    }
};
        
        
struct PolynomialRegistrationTest
{
    BImage s_img;
    BImage d_img;
    
    std::vector<TinyVector<double,2> > s_points;
    std::vector<TinyVector<double,2> > d_points;
    
    PolynomialRegistrationTest()
    : s_points(srcPoints()),
      d_points(destPoints())
    {
        ImageImportInfo info1("nuernberg-1991.png");
        s_img.resize(info1.width(), info1.height());
        importImage(info1, destImage(s_img));
        
        ImageImportInfo info2("nuernberg-1995.png");
        d_img.resize(info2.width(), info2.height());
        importImage(info2, destImage(d_img));
    }
    
    void testDegree0()
    {                
        /**
         * Test with well-known point sets and a known result matrix
         */
        Matrix<double> poly = polynomialMatrix2DFromCorrespondingPoints<0>(s_points.begin(), s_points.end(), d_points.begin());
        
        /** 
         * Estimated result:
         */        
        Matrix<double> reference(1,2);
        reference(0,0) =  571.1969696970; reference(0,1) = 313.2727272727; 
        
        shouldEqualToleranceMatrices(poly, reference, test_epsilon);
    }
    
    void testDegree1()
    {                
        /**
         * Test with well-known point sets and a known result matrix
         */
        Matrix<double> poly = polynomialMatrix2DFromCorrespondingPoints<1>(s_points.begin(), s_points.end(), d_points.begin());
        
        /** 
         * Estimated result:
         */        
        Matrix<double> reference(3,2);
        reference(0,0) =  -824.6829238728; reference(0,1) = -431.5519142745; 
        reference(1,0) =     1.4856474242; reference(1,1) =   -0.0376429179;
        reference(2,0) =     0.0350960498; reference(2,1) =    1.2727692988;  
        
        shouldEqualToleranceMatrices(poly, reference, test_epsilon);
        
        /**
         * visual interpretation by means of the warped image:
         */
        BImage temp = d_img;
        polynomialWarpImage<1>(SplineImageView<2,unsigned char>(srcImageRange(s_img)), destImageRange(temp), poly);
        exportImage(srcImageRange(temp), ImageExportInfo("res-poly<1>.png"));
    }
    
    void testDegree2()
    {                
        /**
         * Test with well-known point sets and a known result matrix
         */
        Matrix<double> poly = polynomialMatrix2DFromCorrespondingPoints<2>(s_points.begin(), s_points.end(), d_points.begin());
        
        /** 
         * Estimated result:
         */        
        Matrix<double> reference(6,2);
        reference(0,0) =  -929.9603797537; reference(0,1) = 62.7580017829; 
        reference(1,0) =     1.7368027425; reference(1,1) =  0.1279186082;
        reference(2,0) =     0.0013088970; reference(2,1) = -1.0443071694; 
        reference(3,0) =    -0.0001356183; reference(3,1) = -0.0000165442;
        reference(4,0) =    -0.0000115962; reference(4,1) = -0.0001563045; 
        reference(5,0) =     0.0000516651; reference(5,1) =  0.0022927547; 
        
        shouldEqualToleranceMatrices(poly, reference, test_epsilon);
        
        /**
         * visual interpretation by means of the warped image:
         */
        BImage temp = d_img;
        polynomialWarpImage<2>(SplineImageView<2,unsigned char>(srcImageRange(s_img)), destImageRange(temp), poly);
        exportImage(srcImageRange(temp), ImageExportInfo("res-poly<2>.png"));
    }
    
    void testDegree3()
    {                
        /**
         * Test with well-known point sets and a known result matrix
         */
        Matrix<double> poly = polynomialMatrix2DFromCorrespondingPoints<3>(s_points.begin(), s_points.end(), d_points.begin());
        
        /** 
         * Estimated result:
         */        
        Matrix<double> reference(10,2);
        reference(0,0) =  -694.265682659954450; reference(0,1) = -762.451820047407978; 
        reference(1,0) =     0.931229115560906; reference(1,1) =   -0.268572962775846;
        reference(2,0) =     0.048740058262766; reference(2,1) =    4.944397113671227; 
        reference(3,0) =     0.000541326855916; reference(3,1) =   -0.000071716238508;
        reference(4,0) =     0.000742558256144; reference(4,1) =    0.001026085249118; 
        reference(5,0) =    -0.000725321990759; reference(5,1) =   -0.010416513068178; 
        reference(6,0) =    -0.000000240139823; reference(6,1) =    0.000000135932974; 
        reference(7,0) =    -0.000000079399077; reference(7,1) =   -0.000000426260228;
        reference(8,0) =    -0.000000502113506; reference(8,1) =   -0.000000232882327; 
        reference(9,0) =     0.000000746973881; reference(9,1) =    0.000008095034260;  
        
        shouldEqualToleranceMatrices(poly, reference, test_epsilon);
        
        /**
         * visual interpretation by means of the warped image:
         */
        BImage temp = d_img;
        polynomialWarpImage<3>(SplineImageView<2,unsigned char>(srcImageRange(s_img)), destImageRange(temp), poly);
        exportImage(srcImageRange(temp), ImageExportInfo("res-poly<3>.png"));
    }
};

struct PolynomialRegistrationTestSuite
: public test_suite
{
    PolynomialRegistrationTestSuite()
    : test_suite("PolynomialRegistrationTestSuite")
    {
        add( testCase( &PolynomialIdentityTest::testInit));
        add( testCase( &PolynomialRegistrationTest::testDegree0));
        add( testCase( &PolynomialRegistrationTest::testDegree1));
        add( testCase( &PolynomialRegistrationTest::testDegree2));
        add( testCase( &PolynomialRegistrationTest::testDegree3));
    }
};


template<class RBF>
struct RBFNameTraits
{
    static std::string name() { return "unknown"; }
};

template<>
struct RBFNameTraits<ThinPlateSplineFunctor>
{
    static std::string name() { return "tps"; }
};

template<int N>
struct RBFNameTraits<DistancePowerFunctor<N> >
{
    static std::string name() 
    {
        char num_string[16];
        sprintf ( num_string, "%d", N );
        return std::string("dist<") + std::string(num_string) + std::string(">"); 
    }
};

        
        
template<class RadialBasisFunctor>
struct RadialBasisRegistrationTest
{
    BImage s_img;
    BImage d_img;
    
    std::vector<TinyVector<double,2> > s_points;
    std::vector<TinyVector<double,2> > d_points;
    
    RadialBasisRegistrationTest()
    : s_points(srcPoints()),
    d_points(destPoints())
    {
        ImageImportInfo info1("nuernberg-1991.png");
        s_img.resize(info1.width(), info1.height());
        importImage(info1, destImage(s_img));
        
        ImageImportInfo info2("nuernberg-1995.png");
        d_img.resize(info2.width(), info2.height());
        importImage(info2, destImage(d_img));
    }   
    
    void testIdentity()
    {
        /**
         * First test: If point sets are equal -> identity matrix should be the result!
         */
        RadialBasisFunctor rbf;
        Matrix<double> identity_weight_matrix = rbfMatrix2DFromCorrespondingPoints(s_points.begin(), s_points.end(), s_points.begin(),rbf);
        
        //Reference
        Matrix<double> m(69, 2, 0.0);
        
        m(67,0) = 1;  
        m(68,1) = 1;
        shouldEqualToleranceMatrices(identity_weight_matrix, m, rbf_test_epsilon);
    }
    
    void testEssential()
    {
        /**
         * First test: If point sets are equal -> identity matrix should be the result!
         */
        RadialBasisFunctor rbf;
        Matrix<double> weight_matrix = rbfMatrix2DFromCorrespondingPoints(s_points.begin(), s_points.end(), d_points.begin(),rbf);
        
        for(int j=0; j< d_points.size(); j++)
        {
            double x = d_points[j][0];
            double y = d_points[j][1];
            //Affine part		
            double	sx = weight_matrix(point_count,0)+weight_matrix(point_count+1,0)*x+ weight_matrix(point_count+2,0)*y,
                    sy = weight_matrix(point_count,1)+weight_matrix(point_count+1,1)*x+ weight_matrix(point_count+2,1)*y;
            
            //RBS part
            for(int i=0; i<d_points.size(); i++)
            {
                double weight = rbf(d_points[i], d_points[j]);
                sx += weight_matrix(i,0)*weight;
                sy += weight_matrix(i,1)*weight;
            }
            shouldEqualTolerance(sx, s_points[j][0], rbf_test_epsilon);
            shouldEqualTolerance(sy, s_points[j][1], rbf_test_epsilon);
        }
        
        /**
         * visual interpretation by means of the warped image: 
         */
        BImage temp = d_img;
        rbfWarpImage(SplineImageView<2,unsigned char>(srcImageRange(s_img)), 
                                        destImageRange(temp), 
                                        d_points.begin(), d_points.end(),
                                        weight_matrix,
                                        rbf);
        std::string filename = std::string("res-rbf(") + RBFNameTraits<RadialBasisFunctor>::name() + std::string(").png");
        exportImage(srcImageRange(temp), ImageExportInfo(filename.c_str()));
    }
};

struct ThinPlateSplineRegistrationTest
{
    BImage s_img;
    BImage d_img;
    
    std::vector<TinyVector<double,2> > s_points;
    std::vector<TinyVector<double,2> > d_points;
    
    ThinPlateSplineRegistrationTest()
    : s_points(srcPoints()),
      d_points(destPoints())
    {
        ImageImportInfo info1("nuernberg-1991.png");
        s_img.resize(info1.width(), info1.height());
        importImage(info1, destImage(s_img));
        
        ImageImportInfo info2("nuernberg-1995.png");
        d_img.resize(info2.width(), info2.height());
        importImage(info2, destImage(d_img));
    }   
    
    void testInit()
    {
        /**
         * First test: If point sets are equal -> identity matrix should be the result!
         */
        ThinPlateSplineFunctor rbf;
        Matrix<double> weight_matrix = rbfMatrix2DFromCorrespondingPoints(s_points.begin(), s_points.end(), d_points.begin(),rbf);
        
        //Reference
        Matrix<double> m(69,2);
        m(0,0)  = -0.001377498384002;        m(0,1)  =  0.009226600866451;
        m(1,0)  =  0.000032503339089;        m(1,1)  = -0.001509720626847;
        m(2,0)  =  0.000349849805277;        m(2,1)  =  0.001938153880867;
        m(3,0)  = -0.001340286960980;        m(3,1)  =  0.000001064271251;
        m(4,0)  =  0.000539325657797;        m(4,1)  = -0.000543378211795;
        m(5,0)  =  0.000348484938355;        m(5,1)  = -0.000915966809432;
        m(6,0)  =  0.000554791916605;        m(6,1)  =  0.001775753335912;
        m(7,0)  = -0.000365505447921;        m(7,1)  = -0.001756188922282;
        m(8,0)  =  0.000099526739683;        m(8,1)  = -0.003328881346579;
        m(9,0)  = -0.000855945840757;        m(9,1)  =  0.003122855949267;
        m(10,0) = -0.000273792030977;        m(10,1) =  0.000917143372878;
        m(11,0) =  0.002653140382189;        m(11,1) = -0.000731506595669;
        m(12,0) = -0.000069457355061;        m(12,1) = -0.003435613347967;
        m(13,0) =  0.001058065031450;        m(13,1) =  0.004231105333029;
        m(14,0) = -0.000424334717220;        m(14,1) =  0.052130535965185;
        m(15,0) = -0.000981691205091;        m(15,1) = -0.020823555463346;
        m(16,0) = -0.000678829351729;        m(16,1) = -0.003871154682503;
        m(17,0) =  0.001657429743661;        m(17,1) = -0.004684537034900;
        m(18,0) =  0.001836431431159;        m(18,1) = -0.035442938589403;
        m(19,0) = -0.000251036604155;        m(19,1) =  0.000423307680013;
        m(20,0) = -0.000458224650186;        m(20,1) = -0.001421408760879;
        m(21,0) =  0.000713493348626;        m(21,1) = -0.000687257870581;
        m(22,0) = -0.001296157341485;        m(22,1) =  0.001729551998941;
        m(23,0) =  0.000179155018211;        m(23,1) =  0.002047551397835;
        m(24,0) = -0.001276131259978;        m(24,1) = -0.006769539699348;
        m(25,0) = -0.001669855759675;        m(25,1) = -0.001840405300212;
        m(26,0) =  0.006625210595690;        m(26,1) =  0.004544925562541;
        m(27,0) =  0.005061616020780;        m(27,1) = -0.000998799955481;
        m(28,0) = -0.004699886019629;        m(28,1) = -0.001953236466070;
        m(29,0) = -0.006707966667593;        m(29,1) = -0.005335617405465;
        m(30,0) = -0.000008995171240;        m(30,1) =  0.003085439002388;
        m(31,0) = -0.003332042704579;        m(31,1) =  0.007086300850810;
        m(32,0) =  0.000931177328601;        m(32,1) = -0.002990278471355;
        m(33,0) =  0.002247575545516;        m(33,1) =  0.003224790578197;
        m(34,0) = -0.003940778852971;        m(34,1) = -0.002643673196841;
        m(35,0) =  0.003378921135544;        m(35,1) =  0.005651350720849;
        m(36,0) = -0.003256131814242;        m(36,1) = -0.004754825881629;
        m(37,0) =  0.006689872663251;        m(37,1) =  0.004849334797544;
        m(38,0) = -0.005304870007691;        m(38,1) = -0.008107337815761;
        m(39,0) =  0.004746451338148;        m(39,1) =  0.000065750506455;
        m(40,0) = -0.002867905561383;        m(40,1) =  0.003350233719098;
        m(41,0) = -0.001059664485314;        m(41,1) = -0.002596464294288;
        m(42,0) =  0.004792819318783;        m(42,1) = -0.000662922218361;
        m(43,0) = -0.003831224776921;        m(43,1) =  0.002107270860370;
        m(44,0) = -0.000170035532265;        m(44,1) =  0.000071254611424;
        m(45,0) =  0.000293692366501;        m(45,1) = -0.000242570567639;
        m(46,0) =  0.000812835564245;        m(46,1) =  0.000188878966033;
        m(47,0) = -0.000365508959215;        m(47,1) =  0.000248625368764;
        m(48,0) =  0.000059556845286;        m(48,1) =  0.000423880323634;
        m(49,0) =  0.000006617506705;        m(49,1) =  0.000937127030016;
        m(50,0) =  0.000260640568960;        m(50,1) =  0.001410021894571;
        m(51,0) =  0.000022519101131;        m(51,1) = -0.000395009938888;
        m(52,0) = -0.000695392842208;        m(52,1) = -0.000444357574640;
        m(53,0) = -0.000497156317044;        m(53,1) =  0.001898946517854;
        m(54,0) =  0.000233316447714;        m(54,1) = -0.002215381428407;
        m(55,0) =  0.001607979140179;        m(55,1) = -0.001191652089368;
        m(56,0) = -0.001806711090928;        m(56,1) = -0.000974487074949;
        m(57,0) =  0.001122008742798;        m(57,1) = -0.003154750122584;
        m(58,0) =  0.000000977433459;        m(58,1) =  0.001396272823941;
        m(59,0) = -0.000271382198722;        m(59,1) =  0.001355857463705;
        m(60,0) =  0.000417816858966;        m(60,1) =  0.000776812795073;
        m(61,0) = -0.000050555955893;        m(61,1) =  0.002676273718959;
        m(62,0) =  0.001117169088679;        m(62,1) =  0.003193743999999;
        m(63,0) = -0.000317593888323;        m(63,1) = -0.004971039750199;
        m(64,0) = -0.000544110734374;        m(64,1) =  0.002762020120260;
        m(65,0) =  0.000616392901421;        m(65,1) =  0.002574136813928;
        m(66,0) =  0.000539075187236;        m(66,1) =  0.000035460307110;
        m(67,0) =  1.033807852169952;        m(67,1) = -0.676702412782481;
        m(68,0) = -0.260446157273618;        m(68,1) =  1.130619201887931;
        
        shouldEqualToleranceMatrices(weight_matrix, m, test_epsilon);
    }
};


struct RadialBasisRegistrationTestSuite
: public test_suite
{
    RadialBasisRegistrationTestSuite()
    : test_suite("RadialBasisRegistrationTestSuite")
    {
        //TPS warping
        add( testCase( &RadialBasisRegistrationTest<ThinPlateSplineFunctor>::testIdentity));
        add( testCase( &RadialBasisRegistrationTest<ThinPlateSplineFunctor>::testEssential));
        add( testCase( &ThinPlateSplineRegistrationTest::testInit));
        //DistancePowerFunctor warping
        add( testCase( &RadialBasisRegistrationTest<DistancePowerFunctor<1> >::testIdentity));
        add( testCase( &RadialBasisRegistrationTest<DistancePowerFunctor<1> >::testEssential));
        add( testCase( &RadialBasisRegistrationTest<DistancePowerFunctor<3> >::testIdentity));
        add( testCase( &RadialBasisRegistrationTest<DistancePowerFunctor<3> >::testEssential));        
    }
};



struct RegistrationTestCollection
: public test_suite
{
    RegistrationTestCollection()
    : test_suite("RegistrationTestCollection")
    {
        add( new EstimateGlobalRotationTranslationTestSuite);
        add( new ProjectiveRegistrationTestSuite);
        add( new PolynomialRegistrationTestSuite);
        add( new RadialBasisRegistrationTestSuite);
   }
};

int main(int argc, char ** argv)
{
    RegistrationTestCollection test;
 
    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

