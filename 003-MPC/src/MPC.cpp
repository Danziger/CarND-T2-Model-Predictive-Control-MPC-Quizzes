#include "MPC.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "matplotlibcpp.h"


namespace plt = matplotlibcpp;
using CppAD::AD;


// CONSTANTS:

// We want T to be 2 seconds, so we could use dt = 0.05 seconds and N = 20;
const double dt = 0.05;
const size_t N = 20;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius:
const double Lf = 2.67;

// NOTE: Feel free to play around with this or do something completely different:
const double V_REF = 40;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should establish where one variable
// starts and another ends to make our lifes easier:
const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;


class FG_eval {

public:

    // Coefficients of the fitted polynomial:
    Eigen::VectorXd coeffs_;
    
    FG_eval(Eigen::VectorXd coeffs) {
        coeffs_ = coeffs;
    }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator() (
        ADvector& fg,
        const ADvector& vars
    ) {
        // `fg` is a vector containing the cost and constraints.
        // `vars` is a vector containing the variable values (state & actuators).

        CppAD::AD<double> cost = 0;

        // COSTS:
        // Define the cost related to the reference state and anything you think may be beneficial:

        // Reference State Cost:

        for (size_t t = 0; t < N; ++t) {
            cost += CppAD::pow(vars[cte_start + t], 2);
            cost += CppAD::pow(vars[epsi_start + t], 2);
            cost += CppAD::pow(vars[v_start + t] - V_REF, 2);
        }

        // Control Cost (minimize the use of actuators):

        for (size_t t = 0; t < N - 1; ++t) {
            cost += CppAD::pow(vars[delta_start + t ], 2);
            cost += CppAD::pow(vars[a_start + t], 2);
        }

        // Control Cost (minimize the value gap between sequential actuations to achieve temporal smoothness):

        for (size_t t = 0; t < N - 2; ++t) {
            cost += CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
            cost += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
        }

        // The cost is stored is the first element of fg, fg[0]:
        fg[0] = cost;

        // CONSTRAINTS:
        // In this section you'll setup the model constraints.

        // Initial constraints:
        // We add 1 to each of the starting indices due to cost being located at index 0 of fg:
        fg[1 + x_start] = vars[x_start];
        fg[1 + y_start] = vars[y_start];
        fg[1 + psi_start] = vars[psi_start];
        fg[1 + v_start] = vars[v_start];
        fg[1 + cte_start] = vars[cte_start];
        fg[1 + epsi_start] = vars[epsi_start];

        // The rest of the constraints:
        for (size_t t = 1; t < N; ++t) {
            // State at t + 1:
            AD<double> x1 = vars[x_start + t];
            AD<double> y1 = vars[y_start + t];
            AD<double> psi1 = vars[psi_start + t];
            AD<double> v1 = vars[v_start + t];
            AD<double> cte1 = vars[cte_start + t];
            AD<double> epsi1 = vars[epsi_start + t];

            // State at t:
            AD<double> x0 = vars[x_start + t - 1];
            AD<double> y0 = vars[y_start + t - 1];
            AD<double> psi0 = vars[psi_start + t - 1];
            AD<double> v0 = vars[v_start + t - 1];
            AD<double> cte0 = vars[cte_start + t - 1];
            AD<double> epsi0 = vars[epsi_start + t - 1];

            // Actuators at t:
            AD<double> delta0 = vars[delta_start + t - 1];
            AD<double> a0 = vars[a_start + t - 1];

            // The idea here is to constraint this value to be 0.
            // NOTE: The use of `AD<double>` and use of `CppAD`!
            // This is also CppAD can compute derivatives and pass these to the solver.

            AD<double> f0 = coeffs_[0] + coeffs_[1] * x0;
            AD<double> psides0 = CppAD::atan(coeffs_[1]);

            fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
            fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
            fg[1 + psi_start + t] = psi1 - (psi0 + v0/Lf * delta0 * dt);
            fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
            fg[1 + cte_start + t] = cte1 - (f0 - y0 + v0 * CppAD::sin(epsi0) * dt);
            fg[1 + epsi_start + t] = epsi1 - (psi0 - psides0 + v0/Lf * delta0 * dt);
        }
    }
};


// MPC CLASS DEFINITION:

MPC::MPC() {}

MPC::~MPC() {}

vector<double> MPC::Solve(
    const Eigen::VectorXd x0,
    const Eigen::VectorXd coeffs
) {
    typedef CPPAD_TESTVECTOR(double) Dvector;

    const double x = x0[0];
    const double y = x0[1];
    const double psi = x0[2];
    const double v = x0[3];
    const double cte = x0[4];
    const double epsi = x0[5];

    // Number of independent variables (includes both states and inputs):
    // For example, if the state is a 4 element vector, the actuators is a 2 element vector
    // and there are 10 timesteps. The number of variables is: 4 * 10 + 2 * 9

    // N timesteps = N - 1 actuations:
    size_t n_vars = N * 6 + (N - 1) * 2;

    // Number of constraints:
    size_t n_constraints = N * 6;

    // Initial value of the independent variables.
    Dvector vars(n_vars);

    // Should be 0 except for the initial values:
    for (size_t i = 0; i < n_vars; ++i) {
        vars[i] = 0.0;
    }

    // Set the initial variable values:
    vars[x_start] = x;
    vars[y_start] = y;
    vars[psi_start] = psi;
    vars[v_start] = v;
    vars[cte_start] = cte;
    vars[epsi_start] = epsi;

    // Lower and upper limits for x:
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);

    // Set all non-actuators upper and lowerlimits
    // to the max negative and positive values:
    for (size_t i = 0; i < delta_start; ++i) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    // The upper and lower limits of delta are set to -25 and 25
    // degrees (values in radians):
    // NOTE: Feel free to change this to something else.
    for (size_t i = delta_start; i < a_start; ++i) {
        vars_lowerbound[i] = -0.436332;
        vars_upperbound[i] = 0.436332;
    }

    // Acceleration/decceleration upper and lower limits:
    // NOTE: Feel free to change this to something else.
    for (size_t i = a_start; i < n_vars; ++i) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }

    // Lower and upper limits for constraints:
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);

    // Should be 0 except for the initial values:
    for (size_t i = 0; i < n_constraints; ++i) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }

    // Set the initial variable values:

    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[psi_start] = psi;
    constraints_lowerbound[v_start] = v;
    constraints_lowerbound[cte_start] = cte;
    constraints_lowerbound[epsi_start] = epsi;

    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[psi_start] = psi;
    constraints_upperbound[v_start] = v;
    constraints_upperbound[cte_start] = cte;
    constraints_upperbound[epsi_start] = epsi;

    // Object that computes objective and constraints:
    FG_eval fg_eval(coeffs);

    // Options:
    // NOTE: You don't have to worry about these options.

    string options;

    // Uncomment this if you'd like more print information
    options += "Integer print_level  0\n";

    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER. If you
    // can uncomment 1 of these and see if it makes a difference or not but
    // if you uncomment both the computation time should go up in orders of
    // magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";

    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    // options += "Numeric max_cpu_time          0.5\n";

    // Solve the problem:

    CppAD::ipopt::solve_result<Dvector> solution;

    CppAD::ipopt::solve<Dvector, FG_eval>(
        options,
        vars,
        vars_lowerbound,
        vars_upperbound,
        constraints_lowerbound,
        constraints_upperbound,
        fg_eval,
        solution
    );

    // Check some of the solution values:

    bool ok = solution.status == CppAD::ipopt::solve_result<Dvector>::success;
    auto cost = solution.obj_value;

    if (ok) {
        cout << "  COST = " << cost << endl;
    } else {
        cout << "There was an error calculating the solution." << endl << endl;
    }    

    // {...} is shorthand for creating a vector:

    return {
        solution.x[x_start + 1],
        solution.x[y_start + 1],
        solution.x[psi_start + 1],
        solution.x[v_start + 1],
        solution.x[cte_start + 1],
        solution.x[epsi_start + 1],
        solution.x[delta_start],
        solution.x[a_start]
    };
}


// HELPER FUNCTIONS:

/**
* Evaluate a polynomial defined by coeffs at x.
*/
double polyeval(const Eigen::VectorXd coeffs, const double x) {
    double result = 0.0;

    for (int i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * pow(x, i);
    }

    return result;
}

/**
* Fit a polynomial.
* Adapted from https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
*/
Eigen::VectorXd polyfit(
    const Eigen::VectorXd xvals,
    const Eigen::VectorXd yvals,
    const unsigned int order
) {
    const unsigned int xvalsSize = xvals.size();

    assert(xvalsSize == yvals.size());
    assert(order >= 1 && order <= xvalsSize - 1);

    Eigen::MatrixXd A(xvalsSize, order + 1);

    for (Eigen::Index j = 0; j < xvalsSize; ++j) {
        A(j, 0) = 1.0;

        for (size_t i = 0; i < order; ++i) {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }
    
    return A.householderQr().solve(yvals);
}

int main() {
    MPC mpc;

    Eigen::VectorXd ptsx(2);
    Eigen::VectorXd ptsy(2);
    ptsx << -100, 100;
    ptsy << -1, -1;

    // Fit a polynomial to the above x and y coordinates:
    const auto coeffs = polyfit(ptsx, ptsy, 1);

    // NOTE: free feel to play around with these
    const unsigned int iters = 50;
    const double x = -1;
    const double y = 10;
    const double psi = 0;
    const double v = 10;

    // Calculate the cross track error:
    const double cte = polyeval(coeffs, x) - y;

    // Calculate the orientation error:
    // Due to the sign starting at 0, the orientation error is -f'(x).
    // derivative of coeffs[0] + coeffs[1] * x -> coeffs[1]
    const double epsi = psi - atan(coeffs[1]);

    Eigen::VectorXd state(6);
    state << x, y, psi, v, cte, epsi;

    vector<double> x_vals = {state[0]};
    vector<double> y_vals = {state[1]};
    vector<double> psi_vals = {state[2]};
    vector<double> v_vals = {state[3]};
    vector<double> cte_vals = {state[4]};
    vector<double> epsi_vals = {state[5]};
    vector<double> delta_vals = {};
    vector<double> a_vals = {};

    for (size_t i = 0; i < iters; ++i) {
        cout << endl << " ITERATION " << i << ":" << endl;
        cout << "────────────────────────────────" << endl << endl;

        auto vars = mpc.Solve(state, coeffs);

        x_vals.push_back(vars[0]);
        y_vals.push_back(vars[1]);
        psi_vals.push_back(vars[2]);
        v_vals.push_back(vars[3]);
        cte_vals.push_back(vars[4]);
        epsi_vals.push_back(vars[5]);

        delta_vals.push_back(vars[6]);
        a_vals.push_back(vars[7]);

        state << vars[0], vars[1], vars[2], vars[3], vars[4], vars[5];

        // LOG:

        cout << "     X = " << vars[0] << endl;
        cout << "     Y = " << vars[1] << endl;
        cout << "   PSI = " << vars[2] << endl;
        cout << "     V = " << vars[3] << endl;
        cout << "   CTE = " << vars[4] << endl;
        cout << "  EPSI = " << vars[5] << endl;
        cout << " DELTA = " << vars[6] << endl;
        cout << "     A = " << vars[7] << endl;
        cout << endl;
    }

    // PLOT:
    // NOTE: Feel free to play around with this, it's useful for debugging!

    plt::subplot(3, 1, 1);
    plt::title("CTE");
    plt::plot(cte_vals);

    plt::subplot(3, 1, 2);
    plt::title("Delta (Radians)");
    plt::plot(delta_vals);

    plt::subplot(3, 1, 3);
    plt::title("Velocity");
    plt::plot(v_vals);

    plt::show();
}