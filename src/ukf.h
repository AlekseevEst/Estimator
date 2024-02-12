#pragma once
#include <array>
#include <map>
#include <fstream>
#include <iostream>
// #include <math.h>
#include <cmath>
#include <vector>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

struct Cv
{
    MatrixXd Xe;
    MatrixXd Pe;
};

struct Hcv
{
    MatrixXd Ze;
    MatrixXd Se;
    MatrixXd Pxz;
};

struct Predict
{
    MatrixXd Xe;
    MatrixXd Pe;
    MatrixXd Se;

    MatrixXd Ze;
    MatrixXd K;
};

struct Correct
{
    MatrixXd X;
    MatrixXd P;
};

class unscent_filter
{

public:

    unscent_filter();
    ~unscent_filter();

    MatrixXd predictUkf(const MatrixXd X);
    MatrixXd correctUkf(const MatrixXd Z);


private:
    
};
