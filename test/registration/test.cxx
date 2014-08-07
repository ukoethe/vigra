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

struct RegistrationTestCollection
: public test_suite
{
    RegistrationTestCollection()
    : test_suite("RegistrationTestCollection")
    {
        add( new EstimateGlobalRotationTranslationTestSuite);
        add( new ProjectiveRegistrationTestSuite);
   }
};

int main(int argc, char ** argv)
{
    RegistrationTestCollection test;
 
    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

