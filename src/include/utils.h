#pragma once

#include "Eigen/Dense"
#include <iostream>
using namespace std;
using namespace Eigen;

template <class M = MatrixXd> class Utils
{
public:

    static M rot_Z(const double &val);
    static M rot_Y(const double &val);
    static std::pair<M, M> sph2cartcov(const M &sphCov, const double &r, const double &az, const double &el);
    static M do_cart_P(std::pair<M, M> cartCov);

private:
};