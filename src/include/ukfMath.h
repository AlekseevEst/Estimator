
#pragma once
#include "Eigen/Dense"
#include "structs.h"
#include "utils.h"

template <class M>
struct UnscentedKalmanFilterMath
{
    M X;
    M P;
    M R_sph_rad;
    M Q;
    M U;
    size_t n;
    double kappa;
    double t;
    Predict<M> predictStruct;
    Correct<M> correctStruct;

    UnscentedKalmanFilterMath(const M &state, const M &measNoiseMatRadian, const M &procNoise, double T, double &k);
    void make_P_cart();
    void sqrt_matrix_p();
    std::vector<M> doSigmaVectors();
    std::vector<double> calculationVectorWeights();
    void doExtrapolatedStateVector(std::vector<M> &Xue, std::vector<double> w);
    void doCovMatExtrapolatedStateVector(const std::vector<M> &Xue, std::vector<double> w);
    
    void doExtrapolatedStateVectorSph(const std::vector<M> &Zue, std::vector<double> w);
    void doCovMatExtrapolatedStateVectorSph(const std::vector<M> &Zue, std::vector<double> w);
    void calcGainFilter(const std::vector<M> &Xue, const std::vector<M> &Zue, std::vector<double> w);
    M correctState(const M& Z);
    M correctCov();

};

template <class M>
UnscentedKalmanFilterMath<M>::UnscentedKalmanFilterMath(const M &state, const M &measNoiseMatRadian, const M &procNoise, double T, double &k)
{   
    X = state;
    P.resize(X.rows(),X.rows());
    P.setZero();
    R_sph_rad = measNoiseMatRadian;
    Q = procNoise;
    n = P.cols();
    kappa = k;
    t = T;
}

template <class M>
void UnscentedKalmanFilterMath<M>::make_P_cart()
{  
    
    if (P.isZero())
    {
        M R_sph_deg = Utils<M>::RsphRad2RsphDeg(R_sph_rad);
        Measurement measZ0 = Utils<M>::make_Z0(X);
        P = Utils<M>::do_cart_P(Utils<M>::sph2cartcov(R_sph_deg, measZ0.r_meas, measZ0.az_meas, measZ0.um_meas));
    }

}

template <class M>
void UnscentedKalmanFilterMath<M>::sqrt_matrix_p()
{
// -----------ВЗЯТИЕ МАТРИЧНОГО КОРНЯ----------------

    M L = Utils<M>::CholeskyLowerTriangularTransposition(P);
    U = sqrt(n + kappa) * L; // Масшатбирующий коэффициент умноженный на Матричный корень
}

template <class M>
std::vector<M> UnscentedKalmanFilterMath<M>::doSigmaVectors()
{
    //----------СОЗДАЕМ Xu СИГМА-ВЕКТОРОВ------------------
    std::vector<M> Xu(2*n+1);
    // Первый компонент
    Xu[0] = X; // в качестве первого сигма вектора берется текущий вектор состояния.

    // Второй компонент. В качестве n/2 берется сумма среднего и некоторого отклонения U.col(i)
    for (int i = 0; i < n; i++)
    {   
        Xu[i + 1] = X + U.col(i);
    }
    // Третий компонент. В качестве n/2 берется разность среднего и некоторого отклонения U.col(i)
    for (int i = 0; i < n; i++)
    {  
        Xu[i + n + 1] = X - U.col(i);

    }
    
    return Xu;
}
template <class M>
std::vector<double> UnscentedKalmanFilterMath<M>::calculationVectorWeights()
{
    // ------------РАСЧЕТ ВЕСОВ ВЕКТОРОВ--------------

    std::vector<double> w(2*n+1);

    w[0] = kappa / (n + kappa);

    for (int i = 0; i < 2 * n ; i++)
    {
        w[i + 1] = 0.5 / (n + kappa);
    }

    return w;
}

template <class M>
void UnscentedKalmanFilterMath<M>::doExtrapolatedStateVector(std::vector<M> &Xue, std::vector<double> w)
{
    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ----------

    predictStruct.Xe = M::Zero(Xue[0].rows(), Xue[0].cols());

    for (int i = 0; i < Xue.size(); i++)
    {
        predictStruct.Xe = predictStruct.Xe + w[i] * Xue[i];
    }
    
}

template <class M>
void UnscentedKalmanFilterMath<M>::doCovMatExtrapolatedStateVector(const std::vector<M> &Xue, std::vector<double> w)
{

    //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ
    predictStruct.Pe = M::Zero(P.rows(), P.cols());
    
    for (int i = 0; i < Xue.size(); i++)
    {
        M dX = Xue[i] - predictStruct.Xe;
        predictStruct.Pe = predictStruct.Pe + w[i] * (dX * dX.transpose());
    }

    predictStruct.Pe = predictStruct.Pe + Utils<M>::doMatrixNoiseProc_Q(Q, t);

}

template <class M>
void UnscentedKalmanFilterMath<M>::doExtrapolatedStateVectorSph(const std::vector<M> &Zue, std::vector<double> w)
{
    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СФЕР.  ------------------

    predictStruct.Ze = M::Zero(Zue[0].rows(),Zue[0].cols());
    for (int i = 0; i < Zue.size(); i++)
    {
        predictStruct.Ze = predictStruct.Ze + w[i] * Zue[i];
    }
}

template <class M>
void UnscentedKalmanFilterMath<M>::doCovMatExtrapolatedStateVectorSph(const std::vector<M> &Zue, std::vector<double> w)
{
    //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СФЕР.
    M Pzz = M::Zero(Zue[0].rows(), Zue[0].rows());
    for (int i = 0; i < Zue.size(); i++)
    {
        M v;
        v = Zue[i] - predictStruct.Ze;
        Pzz = Pzz + w[i] * (v * v.transpose());
    }
    predictStruct.Se = Pzz + R_sph_rad;
}

//----------------------------------------------------------------------

template <class M>
void UnscentedKalmanFilterMath<M>::calcGainFilter(const std::vector<M> &Xue, const std::vector<M> &Zue, std::vector<double> w)
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
    X = correctStruct.X; 
    return correctStruct.X;
}

template <class M>
M UnscentedKalmanFilterMath<M>::correctCov()
{   
    correctStruct.P = predictStruct.Pe - (predictStruct.K * predictStruct.Se) * predictStruct.K.transpose();

    if(Utils<M>::СheckingСonditionsMat(correctStruct.P)) // проверка на симметричность, положительно определённость и не вырожденность
    {
        P = correctStruct.P;
    }
        else {
        std::cout<< "\n ERROR: the matrix does not meet the requirements."<<correctStruct.P<<std::endl;
    }
     
    return correctStruct.P;
    
}