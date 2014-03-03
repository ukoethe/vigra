#include <iostream>
#include <vigra/matrix.hxx>
#include <vigra/regression.hxx>
#include <vigra/quadprog.hxx>

double A_data[] = {
     1, -3,  2,
    -3, 10, -5,
     2, -5,  6
};

double b_data[] = {
     27, 
    -78, 
     64
};

int main()
{
    using namespace vigra;
    using namespace vigra::linalg;
    
    Matrix<double> A(Shape2(3,3), A_data);
    Matrix<double> b(Shape2(3,1), b_data);
    Matrix<double> x(Shape2(3,1));
    
     // minimize (A*x-b)^2  s.t. x>=0 using least angle regression (LARS algorithm)
    nonnegativeLeastSquares(A, b, x);

    std::cout << "solution LARS: ";
    for(int k=0; k<3; ++k)
        std::cout << x(k,0) << ", ";
    std::cout << "\n";
    
    // expected output:
    // solution LARS: 18.4493, 0, 4.50725,
    
    Matrix<double> eye(identityMatrix<double>(3)),
                   zeros(Shape2(3,1)),
                   empty,
                   U =  transpose(A)*A,
                   v = -transpose(A)*b;
    x = 0;
    
     // minimize (A*x-b)^2  s.t. x>=0 using the Goldfarb-Idnani algorithm
    quadraticProgramming(U, v, empty, empty, eye, zeros, x);
    
    std::cout << "solution Goldfarb-Idnani: ";
    for(int k=0; k<3; ++k)
        std::cout << x(k,0) << ", ";
    std::cout << "\n";
    
    // expected output:
    // solution Goldfarb-Idnani: 18.4493, 0, 4.50725,
}
