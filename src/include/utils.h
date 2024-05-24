#pragma once
#include "Eigen/Dense"
#include <iostream>
#include "structs.h"

#define ENUM_TO_INT(x) static_cast<int>(x)
#define PRINTM(x) std::cerr << #x << std::endl << x << std::endl<< __FILE__ << ":" << __LINE__ << std::endl << std::endl
template <class M> 
class Utils
{
public:

    static M rot_Z(const double &val);
    static M rot_Y(const double &val);
    static std::pair<M, M> sph2cartcov(const M &sphCov, const double &r, const double &az, const double &el);
    static M do_cart_P0(std::pair<M, M> cartCov, int numOfParameters);
    static M doMatrixNoiseProc_Q(M procVar, double T, int size);
    static Measurement make_Z0(const M &X);
    // static M RsphRad2RsphDeg(const M &R);
    static double ComputeAngleDifference(double angle1, double angle2);
    static M sqrtMat(const M& P);
    static bool CheckingConditionsMat(const M& P);

private:
};

enum class SizeMat
{   
    ROW0 = 0,
    ROW1 = 1,
    ROW3 = 3,
    ROW4 = 4,
    ROW6 = 6,
    ROW7 = 7,
    ROW9 = 9,
    COL0 = 0,
    COL1 = 1,
    COL3 = 3,
    COL4 = 4,
    COL6 = 6,
    COL7 = 7,
    COL9 = 9
};
enum class MeasPositionMat
{
    RANGE = 0,
    AZ = 1,
    EL = 2
};
enum class CoordPositionMat
{
    X = 0,
    VX =1,
    Y = 2,
    VY =3,
    Z = 4,
    VZ =5,
    W = 6,

    X_CA = 0,
    Y_CA = 3,
    Z_CA = 6,
    VX_CA = 1,
    VY_CA = 4,
    VZ_CA = 7
};

enum class SphPos{

    POS_RANGE = 0,
    POS_AZIM,
    POS_ELEV, 
    POS_VR
    
};

template <class M>
M Utils<M>::rot_Z(const double &val)
{

    M R(ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL3));
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

    M R(ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL3));

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

    M Rpos (ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL3));
    Rpos << pow(rngSig, 2.0), 0.0, 0.0,
        0.0, pow(r * cos(el * (M_PI / 180.0)) * azSig * (M_PI / 180.0), 2.0), 0.0,
        0.0, 0.0, pow(r * elSig * (M_PI / 180.0), 2.0);

    M rot = rot_Z(az) * rot_Y(el).transpose();
    M posCov = rot * Rpos * rot.transpose();
    M velCov = M::Zero(ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL3));
    
    if (sphCov.rows()== ENUM_TO_INT(SizeMat::ROW4) && sphCov.cols()== ENUM_TO_INT(SizeMat::COL4))
    {
        double rrSig = sqrt(sphCov(3,3));
        double crossVelSig = pow(100,2);
        M Rvel (ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL3));     
        Rvel << pow(rrSig,2), 0.0, 0.0,
                0.0, pow(crossVelSig,2), 0.0,
                0.0, 0.0, pow(crossVelSig,2);
        velCov = rot * Rvel * rot.transpose(); 
    }
    else
    {
        velCov.diagonal() << pow(200,2), pow(200,2), pow(200,2);
    }
    return std::make_pair(posCov, velCov);
}

template <class M>
M Utils<M>::do_cart_P0(std::pair<M, M> cartCov, int numOfParameters)
{
    M posCov = cartCov.first;
    M velCov = cartCov.second;

    if (numOfParameters == ENUM_TO_INT(SizeMat::ROW9))
    {
        M Hp(ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL9));
        Hp << 1, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 1, 0, 0;

        M Hv(ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL9));
        Hv << 0, 1, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1, 0;

        M Ha(ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL9));
        Ha << 0, 0, 1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1;

        M AccelerationCov = M::Zero(ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL3));
        AccelerationCov.diagonal() << pow(50,2), pow(50,2), pow(50,2);
        M P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv + Ha.transpose() * AccelerationCov * Ha;
        return P;
    }
        M Hp(ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL6));
        Hp << 1, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 1, 0;

        M Hv(ENUM_TO_INT(SizeMat::ROW3), ENUM_TO_INT(SizeMat::COL6));
        Hv << 0, 1, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1;

        M P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv;

    if (numOfParameters == ENUM_TO_INT(SizeMat::ROW7))
    {
        M PAngleVel = M::Zero(ENUM_TO_INT(SizeMat::ROW7), ENUM_TO_INT(SizeMat::COL7));
        PAngleVel.block(ENUM_TO_INT(SizeMat::ROW0),ENUM_TO_INT(SizeMat::COL0),ENUM_TO_INT(SizeMat::ROW6),ENUM_TO_INT(SizeMat::COL6)) = P;
        PAngleVel(ENUM_TO_INT(SizeMat::ROW6),ENUM_TO_INT(SizeMat::COL6)) = pow(22,2);
        return PAngleVel;
    }
    return P;
}

template <class M>
M Utils<M>::doMatrixNoiseProc_Q(M Q, double T, int size)
{
    if (size == ENUM_TO_INT(SizeMat::ROW6))
    {
        M G(ENUM_TO_INT(SizeMat::ROW6), ENUM_TO_INT(SizeMat::COL3));
        G << (T * T) / 2.0,          0.0,            0.0,
                    T,               0.0,            0.0,
                   0.0,         (T * T) / 2.0,       0.0,
                   0.0,               T,             0.0,
                   0.0,              0.0,       (T * T) / 2.0,
                   0.0,              0.0,             T;
        M Qp = G * Q * G.transpose();
        return Qp;
    }
    if (size == ENUM_TO_INT(SizeMat::ROW7))
    {
    M G(ENUM_TO_INT(SizeMat::ROW7), ENUM_TO_INT(SizeMat::COL4));
    G << (T * T) / 2.0,      0.0,               0.0,          0.0,
               T,            0.0,               0.0,          0.0,
              0.0,      (T * T) / 2.0,          0.0,          0.0,
              0.0,            T,                0.0,          0.0,
              0.0,           0.0,          (T * T) / 2.0,     0.0,
              0.0,           0.0,                T,           0.0,
              0.0,           0.0,               0.0,          1.0;

    M Qp = G * Q * G.transpose();
    return Qp;
    }

    M G(ENUM_TO_INT(SizeMat::ROW9), ENUM_TO_INT(SizeMat::COL3));
    G << (T * T) / 2.0,      0.0,               0.0,
               T,            0.0,               0.0,
              1.0,           0.0,               0.0,         
              0.0,      (T * T) / 2.0,          0.0,         
              0.0,            T,                0.0,
              0.0,           1.0,               0.0,         
              0.0,           0.0,          (T * T) / 2.0,    
              0.0,           0.0,                T,          
              0.0,           0.0,               1.0;

    M Qp = G * Q * G.transpose();
    // PRINTM(Qp);
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


// template <class M>
// M Utils<M>::RsphRad2RsphDeg(const M &R)
// {
//     M R_sph_deg(R.rows(), R.cols());
//     double dispersRgn_R = R(ENUM_TO_INT(SphPos::POS_RANGE), ENUM_TO_INT(SphPos::POS_RANGE));
//     double dispersAz_R_rad = R(ENUM_TO_INT(SphPos::POS_AZIM), ENUM_TO_INT(SphPos::POS_AZIM));
//     double dispersUm_R_rad = R(ENUM_TO_INT(SphPos::POS_ELEV), ENUM_TO_INT(SphPos::POS_ELEV));

    
//         if (R.rows()== ENUM_TO_INT(SizeMat::ROW3))
//         R_sph_deg << dispersRgn_R, 0.0, 0.0, 
//                     0.0, (dispersAz_R_rad * (180 / M_PI)), 0.0,
//                     0.0, 0.0, (dispersUm_R_rad * (180 / M_PI));

//         else
//         {
//                 double dispersVr = R(ENUM_TO_INT(SphPos::POS_VR), ENUM_TO_INT(SphPos::POS_VR));
//                 R_sph_deg <<    dispersRgn_R, 0.0, 0.0, 0.0,
//                                 0.0, (dispersAz_R_rad * (180 / M_PI)), 0.0, 0.0,
//                                 0.0, 0.0, (dispersUm_R_rad * (180 / M_PI)), 0.0,
//                                 0.0, 0.0, 0.0, dispersVr;
//         }

//     return R_sph_deg;
// }


template <class M>
double Utils<M>::ComputeAngleDifference(double angle1, double angle2)
{

    double diff = std::arg(std::complex<double>(cos(angle1 - angle2), sin((angle1 - angle2))));
    return diff;
}

template <class M>
bool Utils<M>::CheckingConditionsMat(const M &P)
{
    if ((P.transpose().isApprox(P, 1e-8)) && (P.llt().info() == Eigen::Success) && (P.determinant() != 0))
        return true;
    else
        return false;
}
template <class M>
M Utils<M>::sqrtMat(const M& P)
{
     Eigen::LLT<M>lltofP(P);
    if (lltofP.info() != Eigen::Success)
    {
        throw std::runtime_error("cholesky decomposition ERROR");
    }
    M L = lltofP.matrixL();
    return L;
}