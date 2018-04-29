#! /bin/bash

# update
sudo apt-get update

# gfortran dependency
sudo apt-get install gfortran

# get unzip
sudo apt-get install unzip

# get, unzip and remove Ipopt
wget https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.9.zip && unzip Ipopt-3.12.9.zip && rm Ipopt-3.12.9.zip

# install Ipopt
sudo ./install_ipopt.sh Ipopt-3.12.9/

# CppAD
sudo apt-get install cppad

# Gnuplot
sudo apt-get install gnuplot

# python and matplotlib
sudo apt-get install python-matplotlib
sudo apt-get install python2.7-dev
