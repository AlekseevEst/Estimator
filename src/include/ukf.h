#pragma once

#include <array>
#include <map>
#include <iostream>
#include <cmath>
#include <vector>
#include <boost/math/distributions/chi_squared.hpp>
#include "utils.h"
#include "structs.h"
using namespace std;
using namespace Eigen;

// Хорошо, но давай чуть переделаем. Повысим читаемость.
// template <class M = MatrixXd> class unscent_filter {....}
// В классе не должно фигурировать тип матрицы, везде

template <class M = MatrixXd> class unscent_filter
{

public:
    unscent_filter();
    ~unscent_filter();

    M predictUkf(const M X);
    M correctUkf(const M Z);

private:

    Predict predictStruct;
    Correct correctStruct;

    M P = M::Zero(6, 6);
    M R_sph_rad = M::Zero(3, 3);;
    M R_sph_deg = M::Zero(3, 3);;

    double dispRgn_R; // дисперс. задается в рад.
    double dispAz_R_rad;
    double dispUm_R_rad;

    double sa;

    std::map<double, double> w;
    std::map<double, M> Xue;
};