
#pragma once
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "structs.h"
#include "utils.h"

template <class M>
struct UnscentedKalmanFilterMath
{

    void compute_weights(ParamSigmaPoints &paramSigmaPoints, double &lamda, double &c, const int &dim_x, std::vector<double> &Wc, std::vector<double> &Wm)
    {
        lamda = pow(paramSigmaPoints.alpha, 2) * (dim_x + paramSigmaPoints.kappa) - dim_x;
        c = 0.5 / (dim_x + lamda);
        Wc.assign(2 * dim_x + 1, c);
        Wm.assign(2 * dim_x + 1, c);
        Wc[0] = lamda / (dim_x + lamda) + (1 - pow(paramSigmaPoints.alpha, 2) + paramSigmaPoints.beta);
        Wm[0] = lamda / (dim_x + lamda);
    }

    void compute_sigma_points(const M &X, const M &P, const double &lamda, M &Xu, M &U)
    {
        U = sqrt(lamda + X.rows()) * Utils<M>::sqrtMatSpectral(P);

        Xu.col(0) = X;

        for (long int i = 0; i < X.rows(); i++)
        {
            Xu.col(i + 1) = X + U.col(i);
        }
        for (long int i = 0; i < X.rows(); i++)
        {
            Xu.col(i + X.rows() + 1) = X - U.col(i);
        }
    }

    void doExtrapolatedStateVector(const M &Xue, M &Xe, const std::vector<double>& Wm)
    {
        //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ----------
        Xe.setZero();
        
        
        for (int i = 0; i < Xue.cols(); i++)
        {
            Xe += Wm[i] * Xue.col(i);
        }
    }

    void doCovMatExtrapolatedStateVector(const M &Xue, const M &Xe, const M &G, const M &Q, M &Pe, const std::vector<double>& Wc)
    {
        //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ
        Pe.setZero();
        for (int i = 0; i < Xue.cols(); i++)
        {
          
            Pe += Wc[i] * ((Xue.col(i) - Xe) * (Xue.col(i) - Xe).transpose());
        }

        Pe += G * Q * G.transpose();
   
    }

    void doExtrapolatedMeasVector(const M &Zue, M &Ze, const std::vector<double>& Wm)
    {
        //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СФЕР.  ------------------
        Ze.setZero();
        for (int i = 0; i < Zue.cols(); i++)
        {
            Ze += Wm[i] * Zue.col(i);
        }
    }

    void doCovMatExtrapolatedMeasVector(const M &Zue, const M &Ze, M &Pzz, const std::vector<double>& Wc)
    {
        //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА ИЗМЕРЕНИИ СФЕР.
        Pzz.setZero();
        for (int i = 0; i < Zue.cols(); i++)
        {

            Pzz += Wc[i] * ((Zue.col(i) - Ze) * (Zue.col(i) - Ze).transpose());
        }
    }

    void doCovMatInnovation(const M &Pzz, const M &R_sph_deg, M &Se)
    {

        Se = Pzz + R_sph_deg; // innovation covariance

    }

    //----------------------------------------------------------------------

    void calcGainFilter(const M &Xue, const M &Xe, const M &Zue, const M &Ze, const M &Se, M &Pxz, M &gainKalman, const std::vector<double>& Wc)
    {
        Pxz.setZero();

        for (int i = 0; i < Zue.cols(); i++)
        {

            Pxz += Wc[i] * (Xue.col(i) - Xe) * (Zue.col(i) - Ze).transpose();
        }
      
        gainKalman = Pxz * Se.inverse();
        
    }

    void correctState(const M &Xe, const M &Z, const M &Ze, const M &K, M &X)

    {
        X = Xe + K * (Z - Ze);
     
    }

    void correctCov(const M &Pe, const M &K, const M &Se, M &P)
    {
        P = Pe - K * Se * K.transpose();
        if (!Utils<M>::CheckingConditionsMat(P)) // проверка на симметричность, положительно определённость и не вырожденность
            throw std::runtime_error("СheckingСonditionsMat ERROR");
    }

private:

};