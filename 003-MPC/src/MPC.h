#ifndef MPC_H
#define MPC_H


#include <vector>
#include "Eigen-3.3/Eigen/Core"


using namespace std;


class MPC {

public:
    MPC();

    virtual ~MPC();

    /**
    * Solve the model given an initial state.
    * @return the next state and actuations as a vector.
    */
    vector<double> Solve(const Eigen::VectorXd x0, const Eigen::VectorXd coeffs);
};


#endif /* MPC_H */
