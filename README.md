CarND · T2 · Model Predictive Control Quizzes
=============================================

[![Udacity - Self-Driving Car NanoDegree](https://s3.amazonaws.com/udacity-sdc/github/shield-carnd.svg)](http://www.udacity.com/drive)


Project Overview
----------------

[The original starter code and solutions can be found here](https://github.com/udacity/CarND-MPC-Quizzes).

Quizzes for *Vehicle Models* and *Model Predictive Control* sections.

1. [Global Kinematic Model Quiz](./001-Global-Kinematic-Model):
   
   Implement the *Global Kinematic Model*.
   
2. [Polynomial Fitting Quiz](./002-Polyfit):

   Fit and evaluate polynomials.

3. [Mind The Line Quiz](./003-MPC):

   Implement MPC and minimize cross track and orientation errors to a straight line trajectory.
   
   See [this document](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/install_Ipopt_CppAD.md) for setup tips for executing the plotting code in the ```MPC.cpp``` solution file.


Dependencies
------------

- [`cmake >= 3.5`](https://cmake.org/install/)
- `make >= 4.1` (Linux / [Mac](https://developer.apple.com/xcode/features/)), [`3.81` (Windows)](http://gnuwin32.sourceforge.net/packages/make.htm)
- `gcc/g++ >= 5.4` (Linux / [Mac](https://developer.apple.com/xcode/features/)), [`MinGW` (Windows)](http://www.mingw.org/)

The *Global Kinematic Quiz* and *Polynomial Fitting* quizzes have all the dependencies in repo.

For the *MPC* quiz, you'll have to install Ipopt and CppAD. Please, refer to [this document](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/install_Ipopt_CppAD.md) for installation instructions.


Build
-----

To do a quiz:

1. Go to quiz directory.
2. Create a build directory and navigate to it: `mkdir build && cd build`
3. Compile the project: `cmake .. && make`.
4. Run it: `./global_kinematic_model`, `./polyfit` or `./mpc`.
