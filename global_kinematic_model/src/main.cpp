#include <math.h>
#include <iostream>
#include "Eigen-3.3/Eigen/Core"


#define LF 2


using namespace std;
using namespace Eigen;


/**
* Helper functions
*/

double deg2rad(const double x) {
    return x * M_PI / 180;
}

double rad2deg(const double x) {
    return x * 180 / M_PI;
}


/**
* Global kinematic model.
*
* @return the next state.
*/
Eigen::VectorXd globalKinematic(
    const VectorXd state,
    const VectorXd actuators,
    const double dt
) {
    // NOTE: state is [x, y, psi, v]:

    const double x = state(0);
    const double y = state(1);
    const double p = state(2);
    const double v = state(3);

    // NOTE: actuators is [delta, a]:

    const double d = actuators(0);
    const double a = actuators(1);

    VectorXd next_state(state.size());

    next_state <<
        x + v * cos(p) * dt,
        y + v * sin(p) * dt,
        p + v * d / LF * dt,
        v + a * dt;

    return next_state;
}


int main() {
    VectorXd state(4); // [x, y, psi, v]
    state << 0, 0, deg2rad(45), 1;
    
    VectorXd actuators(2); // [delta, v]
    actuators << deg2rad(5), 1;

    // TEST:

    VectorXd expected(4);
    expected << 0.212132, 0.212132, 0.798488, 1.3;

    auto next_state = globalKinematic(state, actuators, 0.3);

    if ((next_state - expected).norm() < 1e-6) {
        cout << endl << "Great, output is:" << endl << endl << next_state << endl << endl;
    } else {
        cout << endl << "Ops, output is:" << endl << endl << next_state << endl << endl << "...but it should be:" << endl << endl << expected << endl << endl;
    }
}
