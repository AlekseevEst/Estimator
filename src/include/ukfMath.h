
#pragma once
#include "Eigen/Dense"
#include "structs.h"
#include "utils.h"

template <class M>
struct UnscentedKalmanFilterMath
{

    UnscentedKalmanFilterMath(const M &measNoiseMatRadian, const M &procNoise, double &k);
    M make_P_cart(const M& P, const M& X);
    M sqrt_matrix_p(const M& P);
    M doSigmaVectors(const M& X, const M& U);
    std::vector<double> calculationVectorWeights();
    M doExtrapolatedStateVector(const M &Xue, std::vector<double> w);
    M doCovMatExtrapolatedStateVector(const M &Xue, const M& Xe, std::vector<double> w, double t);
    
    M doExtrapolatedMeasVector(const M &Zue, std::vector<double> w);
    M doCovMatExtrapolatedMeasVector(const M &Zue,const M& Ze, std::vector<double> w);
    M calcGainFilter(const M &Xue, const M &Xe, const M &Zue, const M &Ze, const M &Se, std::vector<double> w);
    M correctState(const M& Xe,const M& Z, const M& Ze, const M& K);
    M correctCov(const M& Pe, const M& K, const M& Se);

    private:

    M R_sph_rad;
    M Q;
    size_t n;
    double kappa;
};

template <class M>
UnscentedKalmanFilterMath<M>::UnscentedKalmanFilterMath(const M &measNoiseMatRadian, const M &procNoise, double &k)
{   
    R_sph_rad = measNoiseMatRadian;
    Q = procNoise;
    kappa = k;
}

template <class M>
M UnscentedKalmanFilterMath<M>::make_P_cart(const M& P, const M& X)
{   
    if (P.isZero())
    {
        M R_sph_deg = Utils<M>::RsphRad2RsphDeg(R_sph_rad);
        Measurement measZ0 = Utils<M>::make_Z0(X);

        M P0 = Utils<M>::do_cart_P(Utils<M>::sph2cartcov(R_sph_deg, measZ0.r_meas, measZ0.az_meas, measZ0.um_meas));
        n = P0.cols();
        return P0;
    }
    return P;
}

template <class M>
M UnscentedKalmanFilterMath<M>::sqrt_matrix_p(const M& P)
{
// -----------ВЗЯТИЕ МАТРИЧНОГО КОРНЯ----------------

    M L = Utils<M>::CholeskyLowerTriangularTransposition(P);
    M U = sqrt(n + kappa) * L; // Масшатбирующий коэффициент умноженный на Матричный корень
    return U;
}

template <class M>
M UnscentedKalmanFilterMath<M>::doSigmaVectors(const M& X, const M& U)
{
    //----------СОЗДАЕМ Xu СИГМА-ВЕКТОРОВ------------------

    M Xu(X.rows(),2*n+1);
    // Первый компонент
    Xu.col(0) = X; // в качестве первого сигма вектора берется текущий вектор состояния.

    // Второй компонент. В качестве n/2 берется сумма среднего и некоторого отклонения U.col(i)
    for (size_t i = 0; i < n; i++)
    {   
        Xu.col(i + 1) = X + U.col(i);
    }
    // Третий компонент. В качестве n/2 берется разность среднего и некоторого отклонения U.col(i)
    for (size_t i = 0; i < n; i++)
    {  
        Xu.col(i + n + 1) = X - U.col(i);

    }
    
    return Xu;
}
template <class M>
std::vector<double> UnscentedKalmanFilterMath<M>::calculationVectorWeights()
{
    // ------------РАСЧЕТ ВЕСОВ ВЕКТОРОВ--------------

    std::vector<double> w(2*n+1);

    w[0] = kappa / (n + kappa);

    for (size_t i = 0; i < 2 * n ; i++)
    {
        w[i + 1] = 0.5 / (n + kappa);
    }

    return w;
}

template <class M>
M UnscentedKalmanFilterMath<M>::doExtrapolatedStateVector(const M &Xue, std::vector<double> w)
{
    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ----------

    M Xe = M::Zero(Xue.rows(), 1);

    for (int i = 0; i < Xue.cols(); i++)
    {
        Xe = Xe + w[i] * Xue.col(i);
    }
    return Xe;
    
}

template <class M>
M UnscentedKalmanFilterMath<M>::doCovMatExtrapolatedStateVector(const M &Xue, const M& Xe, std::vector<double> w, double t)
{

    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ
    M Pe = M::Zero(Xue.rows(), Xue.rows());
    
    for (int i = 0; i < Xue.cols(); i++)
    {
        M dX = Xue.col(i) - Xe;
        Pe = Pe + w[i] * (dX * dX.transpose());
    }

    Pe = Pe + Utils<M>::doMatrixNoiseProc_Q(Q, t);
    return Pe;
}

template <class M>
M UnscentedKalmanFilterMath<M>::doExtrapolatedMeasVector(const M &Zue, std::vector<double> w)
{
    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СФЕР.  ------------------
    M Ze = M::Zero(Zue.rows(),1);
    for (int i = 0; i < Zue.cols(); i++)
    {
        Ze = Ze + w[i] * Zue.col(i);
    }
    return Ze;
}

template <class M>
M UnscentedKalmanFilterMath<M>::doCovMatExtrapolatedMeasVector(const M &Zue, const M &Ze, std::vector<double> w)
{
    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА ИЗМЕРЕНИИ СФЕР.
    M Pzz = M::Zero(Zue.rows(), Zue.rows());
    for (int i = 0; i < Zue.cols(); i++)
    {
        M v;
        v = Zue.col(i) - Ze;
        Pzz = Pzz + w[i] * (v * v.transpose());
    }
    M Se = Pzz + R_sph_rad;
    return Se;
}

//----------------------------------------------------------------------

template <class M>
M UnscentedKalmanFilterMath<M>::calcGainFilter(const M &Xue, const M &Xe, const M &Zue, const M &Ze, const M &Se, std::vector<double> w)
{
    M Pxz = M::Zero(Xue.rows(), Zue.rows());

    for (int i = 0; i < Zue.cols(); i++)
    {
        M dX = Xue.col(i) - Xe;
        M v = Zue.col(i) - Ze;
        Pxz = Pxz + w[i] * dX * v.transpose();
    }

    M gainKalman = Pxz * Se.inverse();
    return gainKalman;
}

template <class M>
M UnscentedKalmanFilterMath<M>::correctState(const M &Xe, const M &Z, const M &Ze, const M &K)

{
    M X = Xe + K * (Z - Ze);
    return X;
}

template <class M>
M UnscentedKalmanFilterMath<M>::correctCov(const M &Pe, const M &K, const M &Se)
{
    M P = Pe - (K * Se) * K.transpose();

    if (Utils<M>::СheckingСonditionsMat(P)) // проверка на симметричность, положительно определённость и не вырожденность
        return P;
    else
        throw std::runtime_error("СheckingСonditionsMat ERROR");
}