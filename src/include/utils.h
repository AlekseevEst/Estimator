#pragma once
#include "Eigen/Dense"
#include "structs.h"
template <class M> 
class Utils
{
public:

    static M rot_Z(const double &val);
    static M rot_Y(const double &val);
    static std::pair<M, M> sph2cartcov(const M &sphCov, const double &r, const double &az, const double &el);
    static M do_cart_P(std::pair<M, M> cartCov);
    static M doMatrixNoiseProc_Q(M procVar, double T);
    static Measurement make_Z0(const M &X);
    static M RSph2Rcart(const M &R);

private:
};

template <class M>
M Utils<M>::rot_Z(const double &val)
{
    M R(3, 3);
    R << 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;

    double ca = cos(val * (M_PI / 180.0));
    double sa = sin(val * (M_PI / 180.0));
    R(0, 0) = ca;
    R(0, 1) = -sa;
    R(1, 0) = sa;
    R(1, 1) = ca;
    return R;
}

template <class M>
M Utils<M>::rot_Y(const double &val)
{
    M R(3, 3);

    R << 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;

    double ca = cos(val * (M_PI / 180.0));
    double sa = sin(val * (M_PI / 180.0));
    R(0, 0) = ca;
    R(0, 2) = sa;
    R(2, 0) = -sa;
    R(2, 2) = ca;
    return R;
}
template <class M>
std::pair<M, M> Utils<M>::sph2cartcov(const M &sphCov, const double &r, const double &az, const double &el)
{
    double rngSig = sqrt(sphCov(0, 0));
    double azSig = sqrt(sphCov(1, 1));
    double elSig = sqrt(sphCov(2, 2));

    M Rpos (3,3);
    Rpos << pow(rngSig, 2.0), 0.0, 0.0,
        0.0, pow(r * cos(el * (M_PI / 180.0)) * azSig * (M_PI / 180.0), 2.0), 0.0,
        0.0, 0.0, pow(r * elSig * (M_PI / 180.0), 2.0);

    M rot = rot_Z(az) * rot_Y(el).transpose();
    M posCov = rot * Rpos * rot.transpose();
    M velCov = M::Zero(3, 3);
    velCov.diagonal() << 100.0, 100.0, 100.0;
    return std::make_pair(posCov, velCov);
}

template <class M>
M Utils<M>::do_cart_P(std::pair<M, M> cartCov)
{
    M posCov = cartCov.first;
    M velCov = cartCov.second;
    M Hp(3, 6);
    Hp << 1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 1, 0;

    M Hv(3, 6);
    Hv << 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 1;

    M P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv;

    return P;
}

template <class M>
M Utils<M>::doMatrixNoiseProc_Q(M Q, double T)
{
    M G (6,3);
    G << (T * T) / 2.0, 0.0, 0.0,
        T, 0.0, 0.0,
        0.0, (T * T) / 2.0, 0.0,
        0.0, T, 0.0,
        0.0, 0.0, (T * T) / 2.0,
        0.0, 0.0, T;
    
    M Qp = (G * Q) * G.transpose();

    return Qp;
}

template <class M>
Measurement Utils<M>::make_Z0(const M &X)
{
    Measurement MeasZ0;
    MeasZ0.r_meas = sqrt(pow(X(0, 0), 2) + pow(X(2, 0), 2) + pow(X(4, 0), 2));
    MeasZ0.az_meas = atan2(X(2, 0), X(0, 0)) * (180 / M_PI);
    MeasZ0.um_meas = atan2(X(4, 0), sqrt(pow(X(0, 0), 2) + pow(X(2, 0), 2))) * (180 / M_PI);
    return MeasZ0;
}


enum class SphPos{

    POS_RANGE = 0,
    POS_AZIM,
    POS_ELEV 
    
};

template <class M>
M Utils<M>::RSph2Rcart(const M &R)
{
    M R_sph_deg(R.rows(), R.cols());
    double dispersRgn_R = R(static_cast<int>(SphPos::POS_RANGE), static_cast<int>(SphPos::POS_RANGE));
    double dispersAz_R_rad = R(static_cast<int>(SphPos::POS_AZIM), static_cast<int>(SphPos::POS_AZIM));
    double dispersUm_R_rad = R(static_cast<int>(SphPos::POS_ELEV), static_cast<int>(SphPos::POS_ELEV));

    R_sph_deg << (dispersRgn_R), 0.0, 0.0, // известная ковариационная матрица ошибок измерении (СКО измерении).// ОШИБКИ ДОЛЖНЫ БЫТЬ ИЗВЕСТНЫМИ(ОШИБКИ ДАТЧИКОВ)
        0.0, (dispersAz_R_rad * (180 / M_PI)), 0.0,
        0.0, 0.0, (dispersUm_R_rad * (180 / M_PI));

    return R_sph_deg;
}
