#pragma once

#include <array>
#include <fstream>
#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include "structs.h"
#include "utils.h"
using namespace std;
using namespace Eigen;

template <class M = MatrixXd> class ExtendedKalmanFilter
{

public:

    ExtendedKalmanFilter();
    ~ExtendedKalmanFilter();

    M predictEkf(const M X);
    M correctEkf(const M Z);

private:

    Predict predictStruct;
    Correct correctStruct;

    M P = M::Zero(6, 6);
    M R_sph_rad = M::Zero(3, 3);
    M R_sph_deg = M::Zero(3, 3);

    double sa;        // СКО ускорения
    double dispRgn_R; // дисперс. задается в рад.
    double dispAz_R_rad;
    double dispUm_R_rad;
};
