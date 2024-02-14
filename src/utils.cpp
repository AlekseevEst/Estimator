#include "utils.h"

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

    // cout << "P ="<< P << endl;

    return P;
}
template class Utils<MatrixXd>;