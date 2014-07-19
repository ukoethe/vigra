#include <functional>
#include "vigra/unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/openmp_vigra.h"
#include "vigra/multi_array.hxx"
#include "vigra/error.hxx"

using namespace std;
using namespace vigra;

//#define PARA_CHECK;

class ParallelTestFunctor
{
	int ID;
	ParallelTestFunctor(ParallelTestFunctor const&) :
			ID(omp_get_thread_num()) {
	}

	template<class T>
	T operator()(T const& arg1, T const& arg2) const {
		return arg1 + arg2 + ID;
	}
};

class CountIterationFunctor
{
private:
	unsigned int iteration_num;

public:
	CountIterationFunctor()
	:iteration_num(0)
	{ }

	void operator()()
	{
	//#pragma omp atomic update
	#pragma omp critical
		iteration_num++;
	}
	unsigned int getIterationNum()
	{
		return iteration_num;
	}
};

struct OpenMPWrapperTest {
	typedef vigra::DImage Image;

	OpenMPWrapperTest() :
			img(1000, 100000), mask(Shape2(3, 3)) {
		Image::Accessor acc = img.accessor();
		Image::ScanOrderIterator i = img.begin();

		for (int j = 1; i < img.end(); ++i, ++j) {
			acc.set(j, i);
		}

		mask.init(1);
		mask.begin()[0] = 0;
		mask.begin()[8] = 0;
	}

	void copyImageTest() {
		Image img1(1000, 100000);

		CountIterationFunctor functor;

		vigra::omp::copyImage(srcImageRange(img), destImage(img1), functor);

		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		should(functor.getIterationNum() == img.height());

#ifndef PARA_CHECK
		for ( ; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
#else
		int onetheard = (img1.width() * img1.height()) / 3;
#pragma omp parallel sections
		{
#pragma omp section
			{
				i1 = img1.begin();
				i = img.begin();
				std::cout << "First id = " << omp_get_thread_num() << std::endl;
				for ( ; i != img.begin() + onetheard; ++i, ++i1)
				{
					should(acc(i) == acc(i1));
				}
			}
#pragma omp section
			{
				i1 = img1.begin() + onetheard;
				i = img.begin() + onetheard;
				std::cout << "Second id = " << omp_get_thread_num() << std::endl;
				for ( ; i != img.begin()+ 2*onetheard; ++i, ++i1)
				{
					should(acc(i) == acc(i1));
				}
			}
#pragma omp section
			{
				i1 = img1.begin() + 2 * onetheard;
				i = img.begin() + 2 * onetheard;
				std::cout << "Third id = " << omp_get_thread_num() << std::endl;
				for ( ; i != img.end(); ++i, ++i1)
				{
					should(acc(i) == acc(i1));
				}
			}
		} //omp sections
#endif //PARA_CHECK
	}

	void combineTwoImagesTest()
	{
		Image img1(1000, 100000);

		CountIterationFunctor iter_count;

		std::plus<Image::value_type> add;

		vigra::omp::combineTwoImages(srcImageRange(img), srcImage(img), destImage(img1), add, iter_count);

		using namespace functor; //TODO: For what?

		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(2.0*acc(i) == acc(i1));
		}

		should( iter_count.getIterationNum() == img.height());
	}

	void combineTwoImagesIfTest()
	{
	}

	void combineThreeImagesTest()
	{
		//TODO: combineThreeImages does not run in parallel
				//vigra::omp::combineThreeImages(View(img), View(img), View(img1), View(img3), Arg1() + Arg2() + Arg3());
				//vigra_precondition(View(img).shape() == View(img1).shape() && View(img).shape() == View(img3).shape() && View(img).shape() == View(img3).shape(),"combineThreeImages(): shape mismatch between inputs and/or output.");
				//vigra::omp::combineThreeImages(srcImageRange(img), srcImage(img), srcImage(img1), destImage(img3), Arg1() + Arg2() + Arg3());
	}

	void combineThreeImagesIfTest()
	{

	}


	void copyTest()
	{
		Image source(3,3);
		source = 9.9;
		Image destination(3,3);
		destination = 1.1;

		Image::Iterator src_upperleft = source.upperLeft();
		Image::Iterator src_lowerright = source.lowerRight();
		Image::Accessor src_acc = source.accessor();

		Image::Iterator dest_upperleft = destination.upperLeft();
		Image::Accessor dest_acc = destination.accessor();

		//vigra::omp::copyImage(src_upperleft,src_lowerright,src_acc,dest_upperleft, dest_acc);

		Image::ScanOrderIterator src_iterator = source.begin();
		Image::ScanOrderIterator dest_iterator = destination.begin();
		Image::Accessor acc_1 = source.accessor();

		for(; src_iterator != source.end(); ++src_iterator, ++dest_iterator)
		{
			should(acc_1(src_iterator) == acc_1(dest_iterator));
		}
		std::cout << "Done copyTest. Max threads: " << omp_get_max_threads() << std::endl;

		//TODO: Print source, destination
		src_iterator = source.begin();
		for(; src_iterator != source.end(); ++src_iterator)
		{
			std::cout << acc_1(src_iterator) << " " << std::endl;
		}

		dest_iterator = destination.begin();
		std::cout << "\nDestination\n" << std::endl;
		for(; dest_iterator != destination.end(); ++dest_iterator)
		{
			std::cout << acc_1(dest_iterator)<< " " << std::endl;
		}

	}

	void transformImageTest()
	{
		Image source(3,3);
		source = 9;
		Image destination(3,3);
		destination = 1.1;

		Image::Iterator src_upperleft = source.upperLeft();
		Image::Iterator src_lowerright = source.lowerRight();
		Image::Accessor src_acc = source.accessor();

		Image::Iterator dest_upperleft = destination.upperLeft();
		Image::Accessor dest_acc = destination.accessor();

		vigra::omp::transformImage(src_upperleft, src_lowerright, src_acc, dest_upperleft, dest_acc, (double(*)(double))&std::sqrt);

		Image::ScanOrderIterator src_iterator = source.begin();
		Image::ScanOrderIterator dest_iterator = destination.begin();
		Image::Accessor acc_1 = source.accessor();
		Image::Accessor acc_2 = source.accessor();

		for(; src_iterator != source.end(); ++src_iterator, ++dest_iterator)
		{
			should(std::sqrt(acc_1(src_iterator)) == acc_2(dest_iterator));
		}
		std::cout << "Done Transform Test. Max threads: " << omp_get_max_threads() << std::endl;

		//TODO: Print source, destination
		src_iterator = source.begin();
		for(; src_iterator != source.end(); ++src_iterator)
		{
			std::cout << acc_1(src_iterator) << " " << std::endl;
		}

		dest_iterator = destination.begin();
		std::cout << "\nDestination\n" << std::endl;
		for(; dest_iterator != destination.end(); ++dest_iterator)
		{
			std::cout << acc_2(dest_iterator)<< " " << std::endl;
		}
	}
	void transformImageIfTest()
	{
		Image source(3,3);
		source = 9;
		Image destination(3,3);
		destination = 1.1;

		vigra::MultiArray<2, unsigned char> mask(Shape2(2,2)); // New api inside old api :D
		mask.init(1);
		mask.begin()[0] = 0;
		mask.begin()[1] = 0;

		Image::Iterator src_upperleft = source.upperLeft();
		Image::Iterator src_lowerright = source.lowerRight();
		Image::Accessor src_acc = source.accessor();

		//Image::Iterator mask_upperleft = maskImage(mask).;
		//Image::Accessor mask_acc = maskImage(mask).second;

		Image::Iterator dest_upperleft = destination.upperLeft();
		Image::Accessor dest_acc = destination.accessor();

		//vigra::omp::transformImageIf(src_upperleft, src_lowerright, src_acc, mask_upperleft, mask_acc,dest_upperleft, dest_acc, (double(*)(double))&std::sqrt);
		vigra::omp::transformImageIf(srcImageRange(source), maskImage(mask), destImage(destination), (double(*)(double))&std::sqrt);

		Image::ScanOrderIterator src_iterator = source.begin();
		Image::ScanOrderIterator dest_iterator = destination.begin();
		Image::Accessor acc_1 = source.accessor();
		Image::Accessor acc_2 = source.accessor();

//		for(; src_iterator != source.end(); ++src_iterator, ++dest_iterator)
//		{
//			should(std::sqrt(acc_1(src_iterator)) == acc_2(dest_iterator));
//		}
//		std::cout << "Done Transform Test If. Max threads: " << omp_get_max_threads() << std::endl;

		//TODO: Print source, destination
		src_iterator = source.begin();
		for(; src_iterator != source.end(); ++src_iterator)
		{
			std::cout << acc_1(src_iterator) << " " << std::endl;
		}

		dest_iterator = destination.begin();
		std::cout << "\nDestination\n" << std::endl;
		for(; dest_iterator != destination.end(); ++dest_iterator)
		{
			std::cout << acc_2(dest_iterator)<< " " << std::endl;
		}
	}

	Image img;
	vigra::MultiArray<2, unsigned char> mask;
};

struct OpenMPWrapperTestSuite: public vigra::test_suite {
	OpenMPWrapperTestSuite() :
			vigra::test_suite("OpenMPWrapperTestSuite") {

		add( testCase( &OpenMPWrapperTest::copyImageTest));
		add( testCase( &OpenMPWrapperTest::combineTwoImagesTest));
		add( testCase( &OpenMPWrapperTest::transformImageTest));
	}
};

int main(int argc, char ** argv) {
	OpenMPWrapperTestSuite test;

	int failed = test.run(vigra::testsToBeExecuted(argc, argv));

	std::cout << test.report() << std::endl;

	return (failed != 0);
}
