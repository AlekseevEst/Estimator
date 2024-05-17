
#pragma once
#include "Eigen/Dense"
#include "structs.h"
#include "utils.h"
#include "sigma_points.h"
template <class M>
struct UnscentedKalmanFilterMath
{

    UnscentedKalmanFilterMath(const M &measNoiseMat, const M &procNoise, Points &p);
    M make_P_cart(const M& P, const M& X);
    M doSigmaVectors(const M& X, const M& P);
    M doExtrapolatedStateVector(const M &Xue);
    M doCovMatExtrapolatedStateVector(const M &Xue, const M& Xe, double t);
    M doExtrapolatedMeasVector(const M &Zue);
    M doCovMatExtrapolatedMeasVector(const M &Zue,const M& Ze);
    M calcGainFilter(const M &Xue, const M &Xe, const M &Zue, const M &Ze, const M &Se);
    M correctState(const M& Xe,const M& Z, const M& Ze, const M& K);
    M correctCov(const M& Pe, const M& K, const M& Se);
    

    private:

    M R_sph_deg;
    M Q;
    size_t n;
    Points points;
    SigmaPoints<M> sigmaPoints;
    
};

template <class M>
UnscentedKalmanFilterMath<M>::UnscentedKalmanFilterMath(const M &measNoiseMat, const M &procNoise, Points &p)
{   
    SigmaPoints<M> sigmaPoints;
    R_sph_deg = measNoiseMat;
    Q = procNoise; 
    points = p;
}

template <class M>
M UnscentedKalmanFilterMath<M>::make_P_cart(const M& P, const M& X)
{   
    if (P.isZero())
    {
        // M R_sph_deg = Utils<M>::RsphRad2RsphDeg(R_sph_rad);
        Measurement measZ0 = Utils<M>::make_Z0(X);
        int numOfParameters = X.rows();
        M P0 = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(R_sph_deg, measZ0.r_meas, measZ0.az_meas, measZ0.um_meas),numOfParameters);
        PRINTM(P0);
        return P0;
    }
    // PRINTM(P);
    return P;
}

template <class M>
M UnscentedKalmanFilterMath<M>::doSigmaVectors(const M& X, const M& P)
{
    //----------СОЗДАЕМ Xu СИГМА-ВЕКТОРОВ------------------
    M Xu = sigmaPoints.compute_sigma_points(X, P, points);
    sigmaPoints.compute_weights(points);
    return Xu;
}

template <class M>
M UnscentedKalmanFilterMath<M>::doExtrapolatedStateVector(const M &Xue)
{
    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ----------

    M Xe = M::Zero(Xue.rows(), 1);

    for (int i = 0; i < Xue.cols(); i++)
    {
        // PRINTM(sigmaPoints.Wm[i]);
        Xe = Xe + sigmaPoints.Wm[i] * Xue.col(i);
    }
    // PRINTM(Xe);
    return Xe;
    
}

template <class M>
M UnscentedKalmanFilterMath<M>::doCovMatExtrapolatedStateVector(const M &Xue, const M& Xe, double t)
{
    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ
    M Pe = M::Zero(Xue.rows(), Xue.rows());
    
    for (int i = 0; i < Xue.cols(); i++)
    {
        M dX = Xue.col(i) - Xe;
        Pe = Pe + sigmaPoints.Wc[i] * (dX * dX.transpose());

    }
    
    Pe = Pe + Utils<M>::doMatrixNoiseProc_Q(Q, t, Xe.rows());
    // PRINTM(Pe);
    
    return Pe;
}


template <class M>
M UnscentedKalmanFilterMath<M>::doExtrapolatedMeasVector(const M &Zue)
{
    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СФЕР.  ------------------
    M Ze = M::Zero(Zue.rows(),1);
    for (int i = 0; i < Zue.cols(); i++)
    {
        Ze = Ze + sigmaPoints.Wm[i] * Zue.col(i);
    }
    // PRINTM(Ze);
    return Ze;
}

template <class M>
M UnscentedKalmanFilterMath<M>::doCovMatExtrapolatedMeasVector(const M &Zue, const M &Ze)
{
    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА ИЗМЕРЕНИИ СФЕР.
    M Pzz = M::Zero(Zue.rows(), Zue.rows());
    for (int i = 0; i < Zue.cols(); i++)
    {
        M v;
        v = Zue.col(i) - Ze;
        Pzz = Pzz + sigmaPoints.Wc[i] * (v * v.transpose());
    }
    M Se = Pzz + R_sph_deg;
    // PRINTM(Se);
    return Se;
}

//----------------------------------------------------------------------

template <class M>
M UnscentedKalmanFilterMath<M>::calcGainFilter(const M &Xue, const M &Xe, const M &Zue, const M &Ze, const M &Se)
{
    M Pxz = M::Zero(Xue.rows(), Zue.rows());

    for (int i = 0; i < Zue.cols(); i++)
    {
        // PRINTM(sigmaPoints.Wc[i]);
        M dX = Xue.col(i) - Xe;
        // PRINTM(dX);
        M v = Zue.col(i) - Ze;
        // PRINTM(v);
        Pxz = Pxz + sigmaPoints.Wc[i] * dX * v.transpose();
       
    }
    // PRINTM(Pxz);
    M gainKalman = Pxz * Se.inverse();
    // PRINTM(gainKalman);
    return gainKalman;
}

template <class M>
M UnscentedKalmanFilterMath<M>::correctState(const M &Xe, const M &Z, const M &Ze, const M &K)

{
    M X = Xe + K * (Z - Ze);
    // PRINTM(X);
    return X;
}

template <class M>
M UnscentedKalmanFilterMath<M>::correctCov(const M &Pe, const M &K, const M &Se)
{
    M P = Pe - K * Se * K.transpose();
    // PRINTM(P);
    if (Utils<M>::CheckingConditionsMat(P)) // проверка на симметричность, положительно определённость и не вырожденность
        return P;
    else
        throw std::runtime_error("СheckingСonditionsMat ERROR");
}