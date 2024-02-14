#pragma once
#include <Eigen/Dense>

using MatrixType = Eigen::MatrixXd;
// #include "utils.h"
    struct Cv
{
    MatrixType Xe;
    MatrixType Pe;
};

struct Hcv
{
    MatrixType Ze;
    MatrixType Se;
};

struct Predict
{
    MatrixType Xe;
    MatrixType Pe;
    MatrixType Se;

    MatrixType Ze;
    MatrixType K;
};

struct Correct
{
    MatrixType X;
    MatrixType P;
};