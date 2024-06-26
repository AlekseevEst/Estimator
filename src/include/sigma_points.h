#pragma once
#include "utils.h"
#include <vector>


template <class M>
struct SigmaPoints
{
public:
    std::vector<double> Wc, Wm;
    M compute_sigma_points(const M& X, const M& P, Points& points);
    void compute_weights(Points & points);
    
    double lamda;
    size_t n;
private:
    
};

template <class M>
M SigmaPoints<M>::compute_sigma_points (const M& X, const M& P, Points & points)
{   
    n = X.rows();
    lamda = pow(points.alpha,2) * (n + points.kappa) - n;
    M U = sqrt(lamda + n) * Utils<M>::sqrtMat(P);

    M Xu (n, 2*n+1);
    // Первый компонент
    Xu.col(0) = X; // в качестве первого сигма вектора берется текущий вектор состояния.
    // PRINTM(Xu.col(0));
    // Второй компонент. В качестве n/2 берется сумма среднего и некоторого отклонения U.col(i)
    for (size_t i = 0; i < n; i++)
    {   
        Xu.col(i + 1) = X + U.col(i);
        // PRINTM(Xu.col(i + 1));
    }
    // Третий компонент. В качестве n/2 берется разность среднего и некоторого отклонения U.col(i)
    for (size_t i = 0; i < n; i++)
    {  
        Xu.col(i + n + 1) = X - U.col(i);
        // PRINTM(Xu.col(i + n + 1));
    }
    
    return Xu;

}

template <class M>
void SigmaPoints<M>::compute_weights(Points& points)
{
    lamda = pow(points.alpha,2) * (n + points.kappa) - n;
    double c = 0.5/(n + lamda);
    Wc.assign(2 * n + 1, c);
    Wm.assign(2 * n + 1, c);
    Wc[0] = lamda / (n + lamda) + (1 - pow(points.alpha,2) + points.beta);
    Wm[0] = lamda / (n + lamda);
    // PRINTM(Wc[0]);
    // PRINTM(Wm[0]);

}