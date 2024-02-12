#pragma once
#include "Eigen/Dense"
#include <iostream>
using namespace std;
using namespace Eigen;

namespace utils
{

    Matrix3d rot_Z(const double &val)
    {
        Matrix3d R = Matrix3d::Identity();
        double ca = cos(val * (M_PI/180.0));
        double sa = sin(val * (M_PI/180.0));
        R(0, 0) = ca;
        R(0, 1) = -sa;
        R(1, 0) = sa;
        R(1, 1) = ca;
        return R;
    }

    Matrix3d rot_Y(const double &val)
    {
        Matrix3d R = Matrix3d::Identity();
        double ca = cos(val * (M_PI/180.0));
        double sa = sin(val * (M_PI/180.0));
        R(0, 0) = ca;
        R(0, 2) = sa;
        R(2, 0) = -sa;
        R(2, 2) = ca;
        return R;
    }
 
    std::pair<MatrixXd, MatrixXd> sph2cartcov(const MatrixXd &sphCov, const double &r, const double &az, const double &el)
    {
        double rngSig = sqrt(sphCov(0, 0));
        double azSig = sqrt(sphCov(1, 1));
        double elSig = sqrt(sphCov(2, 2));

        Matrix3d Rpos;
        Rpos << pow(rngSig, 2.0), 0.0, 0.0,
                 0.0, pow(r * cos(el * (M_PI/180.0)) * azSig * (M_PI / 180.0), 2.0), 0.0,
                 0.0, 0.0, pow(r * elSig * (M_PI / 180.0), 2.0);


        MatrixXd rot = rot_Z(az) * rot_Y(el).transpose();
        MatrixXd posCov = rot * Rpos * rot.transpose();
        MatrixXd velCov = MatrixXd::Zero(3, 3);
        velCov.diagonal() << 100.0, 100.0, 100.0;
        return std::make_pair(posCov, velCov);
    }

    MatrixXd do_cart_P(std::pair<MatrixXd, MatrixXd> cartCov)
    {
        MatrixXd posCov = cartCov.first;
        MatrixXd velCov = cartCov.second;
        MatrixXd Hp(3, 6);
        Hp << 1, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 1, 0;

        MatrixXd Hv(3, 6);
        Hv << 0, 1, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1;

        MatrixXd P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv;

        // cout << "P ="<< P << endl; 

        return P;
    }
}