#include <iostream>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"


using namespace std;
using namespace Eigen;


/**
* Evaluate a polynomial defined by coeffs at x.
*/
double polyeval(
    const VectorXd coeffs,
    const double x
) {
    double result = 0.0;

    const unsigned int coeffsSize = coeffs.size();

    for (unsigned int i = 0; i < coeffsSize; ++i) {
        result += coeffs[i] * pow(x, i);
    }

    return result;
}


/**
* Fit a polynomial.
* Code adapted from https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
*/
VectorXd polyfit(
    const VectorXd xvals,
    const VectorXd yvals,
    const unsigned int order
) {
    const unsigned int xvalsSize = xvals.size();

    assert(xvalsSize == yvals.size());
    assert(order >= 1 && order <= xvalsSize - 1);

    // Create matrix A:
    MatrixXd A(xvalsSize, order + 1);

    // Initialise first column to 1s:
    for (unsigned int i = 0; i < xvalsSize; ++i) {
        A(i, 0) = 1.0;
    }
    
    // Calculate the remaining columns:
    for (unsigned int j = 0; j < xvalsSize; ++j) {
        const double current_xval = xvals(j);

        for (unsigned int i = 0; i < order; ++i) {
            A(j, i + 1) = A(j, i) * current_xval;
        }
    }

    return A.householderQr().solve(yvals);
}


int main() {
    VectorXd xvals(6);
    VectorXd yvals(6);

    // x waypoint coordinates:
    xvals << 9.261977, -2.06803, -19.6663, -36.868, -51.6263, -66.3482;

    // y waypoint coordinates:
    yvals << 5.17, -2.25, -15.306, -29.46, -42.85, -57.6116;

    // Fit a third order polynomial to the (x, y) coordinates:
    const VectorXd coeffs = polyfit(xvals, yvals, 3);

    // Evaluate the previous polynomial in [0, 1, ..., 20] 
    VectorXd result(21);

    for (unsigned int x = 0; x <= 20; ++x) {
        result(x) = polyeval(coeffs, x);
    }

    // TEST:

    VectorXd expected(21);
    expected <<
        -0.905562,
        -0.226606,
        0.447594,
        1.11706,
        1.7818,
        2.44185,
        3.09723,
        3.74794,
        4.39402,
        5.03548,
        5.67235,
        6.30463,
        6.93236,
        7.55555,
        8.17423,
        8.7884,
        9.3981,
        10.0033,
        10.6041,
        11.2005,
        11.7925;

    if ((result - expected).norm() < 1e-4) {
        cout << endl << "Great, output is:" << endl << endl << result << endl << endl;
    } else {
        cout << endl << "Ops, output is:" << endl << endl << result << endl << endl << "...but it should be:" << endl << endl << expected << endl << endl;
    }
}