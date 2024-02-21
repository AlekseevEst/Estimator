
#pragma once
#include "Eigen/Dense"
#include "structs.h"
#include "utils.h"

template <class M>
struct UnscentedKalmanFilterMath
{

    M P;
    M R_sph_rad;
    M Q;
    size_t n;
    double kappa;
    double t;
    Predict<M> predictStruct;
    Correct<M> correctStruct;

    UnscentedKalmanFilterMath(const M &covMat, const M &measNoiseMatRadian, const M &procNoise, double T, double &k);
    M make_P_cart(const M &P, const M &X);
    M sqrt_matrix_p(const M &covMat, const M &X);
    std::vector<M> doSigmaVectors(const M &X,const M &U);
    std::vector<M> calculationVectorWeights();
    void doExtrapolatedStateVector(std::vector<M> &Xue, std::vector<M> w);
    void doCovMatExtrapolatedStateVector(const std::vector<M> &Xue, std::vector<M> w, const M &Q);
    std::vector<M> doSigmaVectorsSph(const std::vector<M> &Xue);
    void doExtrapolatedStateVectorSph(const std::vector<M> &Zue, std::vector<M> w);
    void doCovMatExtrapolatedStateVectorSph(const std::vector<M> &Zue, std::vector<M> w);
    void calcGainFilter(const std::vector<M> &Xue, const std::vector<M> &Zue, std::vector<M> w);
    M correctState(const M& Z);
    M correctCov();

};

template <class M>
UnscentedKalmanFilterMath<M>::UnscentedKalmanFilterMath(const M &covMat,const M &measNoiseMatRadian, const M &procNoise, double T, double &k)
{
    P = covMat;
    R_sph_rad = measNoiseMatRadian;
    Q = procNoise;
    n = covMat.cols();
    kappa = k;
    t = T;
}

template <class M>
M UnscentedKalmanFilterMath<M>::make_P_cart(const M &P, const M &X)
{

    M cov = P;
    if (cov.isZero())
    {
        M R_sph_deg = Utils<M>::RSph2Rcart(R_sph_rad);
        Measurement measZ0 = Utils<M>::make_Z0(X);
        cov = Utils<M>::do_cart_P(Utils<M>::sph2cartcov(R_sph_deg, measZ0.r_meas, measZ0.az_meas, measZ0.um_meas));
        return cov;
    }
    return cov;

}


template <class M>
M UnscentedKalmanFilterMath<M>::sqrt_matrix_p(const M &covMat, const M &X)
{
// -----------ВЗЯТИЕ МАТРИЧНОГО КОРНЯ----------------

    Eigen::LLT<M>lltofP(covMat);
    if (lltofP.info() != Eigen::Success)
    {
        std::cout << " cholesky decomposition ERROR "; // необходима проверка матриц. чтобы не было вырождения
    }
    
    M L = lltofP.matrixL();
    M U = sqrt(n + kappa) * L; // Масшатбирующий коэффициент умноженный на Матричный корень

    return U;
}

template <class M>
std::vector<M> UnscentedKalmanFilterMath<M>::doSigmaVectors(const M &X, const M &U)
{
    //----------СОЗДАЕМ Xu СИГМА-ВЕКТОРОВ------------------
    std::vector<M> Xu(n);

    // Первый компонент
    Xu[0] = X; // в качестве первого сигма вектора берется текущий вектор состояния.
    // cout << "\nXu[0] = " << Xu[0] << endl;

    // Второй компонент. В качестве n/2 берется сумма среднего и некоторого отклонения U.col(i)
    for (int i = 0; i < n; i++)
    {
        Xu[i + 1] = X + U.col(i);
    }
    // Третий компонент. В качестве n/2 берется разность среднего и некоторого отклонения U.col(i)
    for (double i = 0; i < n; i++)
    {
        Xu[i + n + 1] = X - U.col(i);
    }
    return Xu;
}

template <class M>
std::vector<M> UnscentedKalmanFilterMath<M>::calculationVectorWeights()
{
    // ------------РАСЧЕТ ВЕСОВ ВЕКТОРОВ--------------

    std::vector<M> w(n);

    w[0] = kappa / (n + kappa);

    for (int i = 0; i < 2 * n; i++)
    {
        w[i + 1] = 1.0 / (2.0 * (n + kappa));
    }

    return w;
}

template <class M>
void UnscentedKalmanFilterMath<M>::doExtrapolatedStateVector(std::vector<M> &Xue, std::vector<M> w)
{
    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ----------

    predictStruct.Xe = M::Zero(Xue[0].rows(), Xue[0].cols());

    for (int i = 0; i < Xue.size(); i++)
    {
        predictStruct.Xe = predictStruct.Xe + w[i] * Xue[i];
    }
}

template <class M>
void UnscentedKalmanFilterMath<M>::doCovMatExtrapolatedStateVector(const std::vector<M> &Xue, std::vector<M> w , const M &Q)
{

    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ

    predictStruct.Pe = M::Zero(P.rows(), P.cols());

    for (int i = 0; i < Xue.size(); i++)
    {
        M dX = Xue[i] - predictStruct.Xe;
        predictStruct.Pe = predictStruct.Pe + w[i] * dX * dX.transpose();
    }
    predictStruct.Pe = predictStruct.Pe + Utils<M>::doMatrixProcNoise_Q(Q, t);
}

template <class M>
std::vector<M> UnscentedKalmanFilterMath<M>::doSigmaVectorsSph(const std::vector<M> &Xue)
{
    //----------ЭКСТРАПОЛИРОВАНЫЕ СИГМА-ВЕКТОРА ИЗМЕРЕНИЙ ПО НЕЛИНЕЙНЫМ ФУНКЦИЯМ------------------

    std::vector<M> Zue(Xue.size());

    for (int i = 0; i < Xue.size(); i++)
    {
        M zTmp(3, 1);
        zTmp << sqrt(pow(Xue[i](0, 0), 2) + pow(Xue[i](2, 0), 2) + pow(Xue[i](4, 0), 2)), atan2(Xue[i](2, 0), Xue[i](0, 0)), atan2(Xue[i](4, 0), sqrt(pow(Xue[i](0, 0), 2) + pow(Xue[i](2, 0), 2)));
        Zue[i] = zTmp;
    }
    return Zue;
}

template <class M>
void UnscentedKalmanFilterMath<M>::doExtrapolatedStateVectorSph(const std::vector<M> &Zue, std::vector<M> w)
{
    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СФЕР.  ------------------
    predictStruct.Ze = M::Zero(3, 1);
    for (int i = 0; i < Zue.size(); i++)
    {
        predictStruct.Ze = predictStruct.Ze + (w[i] * Zue[i]);
    }
}

template <class M>
void UnscentedKalmanFilterMath<M>::doCovMatExtrapolatedStateVectorSph(const std::vector<M> &Zue, std::vector<M> w)
{
    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СФЕР.
    M Pzz = M::Zero(Zue[0].rows(), Zue[0].rows());
    for (int i = 0; i < Zue.size(); i++)
    {
        M v;
        v = Zue[i] - predictStruct.Ze;
        Pzz = Pzz + w[i] * v * v.transpose();
    }
    predictStruct.Se = Pzz + R_sph_rad;
}

//----------------------------------------------------------------------

template <class M>
void UnscentedKalmanFilterMath<M>::calcGainFilter(const std::vector<M> &Xue, const std::vector<M> &Zue, std::vector<M> w)
{
    M Pxz = M::Zero(Xue[0].rows(), Zue[0].rows());

    for (int i = 0; i < Zue.size(); i++)
    {
        M dX = Xue[i] - predictStruct.Xe;
        M v = Zue[i] - predictStruct.Ze;
        Pxz = Pxz + w[i] * dX * v.transpose();
    }

    predictStruct.K = Pxz * predictStruct.Se.inverse();
}

template <class M>
M UnscentedKalmanFilterMath<M>::correctState(const M& Z)
{   
    correctStruct.X = predictStruct.Xe + predictStruct.K * (Z - predictStruct.Ze); 
    return correctStruct.X;
}

template <class M>
M UnscentedKalmanFilterMath<M>::correctCov()
{   
    correctStruct.P = predictStruct.Pe - (predictStruct.K * predictStruct.Se) * predictStruct.K.transpose();
    P = correctStruct.P;
    return correctStruct.P;
}