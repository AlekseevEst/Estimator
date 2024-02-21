// #pragma once

// #include <array>
// #include <fstream>
// #include <iostream>
// #include <math.h>
// #include <cmath>
// #include <iostream>
// #include <vector>
// #include "structs.h"
// #include "utils.h"

// template <class M> class ExtendedKalmanFilter
// {

// public:

//     ExtendedKalmanFilter();
//     ~ExtendedKalmanFilter();

//     M predictEkf(const M& X);
//     M correctEkf(const M& Z);

// private:

//     Predict<M> predictStruct;
//     Correct<M> correctStruct;

//     M P = M::Zero(6, 6);
//     M R_sph_rad = M::Zero(3, 3);
//     M R_sph_deg = M::Zero(3, 3);

//     double sa;        // СКО ускорения
//     double dispRgn_R; // дисперс. задается в рад.
//     double dispAz_R_rad;
//     double dispUm_R_rad;
// };

// template <class M>
// ExtendedKalmanFilter<M>::ExtendedKalmanFilter()
// {
//     R_sph_rad << (dispRgn_R), 0.0, 0.0, // известная ковариационная матрица ошибок измерении (СКО измерении).// ОШИБКИ ДОЛЖНЫ БЫТЬ ИЗВЕСТНЫМИ(ОШИБКИ ДАТЧИКОВ)
//     0.0, (dispAz_R_rad), 0.0,
//     0.0, 0.0, (dispUm_R_rad);

//     R_sph_deg << (dispRgn_R), 0.0, 0.0, // известная ковариационная матрица ошибок измерении (СКО измерении).// ОШИБКИ ДОЛЖНЫ БЫТЬ ИЗВЕСТНЫМИ(ОШИБКИ ДАТЧИКОВ)
//     0.0, (dispAz_R_rad * (180/M_PI)), 0.0,
//     0.0, 0.0, (dispUm_R_rad * (180/M_PI));

//     dispRgn_R = 1.0; // дисперс. в рад.
//     dispAz_R_rad = 1e-4; //
//     dispUm_R_rad = 1e-4; //
    
//     sa = 0.5; 
// }
// template <class M>
// ExtendedKalmanFilter<M>::~ExtendedKalmanFilter()
// {
// }

// template <class M>
// M ExtendedKalmanFilter<M>::predictEkf(const M& X)
// {
//     double T = 6.0;               

//     M F(6, 6); // матрица процесса (state transition matrix). (F*x(k-1)) — это как раз модель эволюции процесса.
//     M G(6, 3); // матрица пересчета случайного воздействия в пространство траекторных параметров
//     M Q(3, 3);       // матрица ковариационных ошибок (шумовое влияние)


//     if(P.isZero())
//     {
//         double r_meas = sqrt(pow(X(0, 0), 2) + pow(X(2, 0), 2) + pow(X(4, 0), 2));
//         double az_meas =  atan2(X(2, 0), X(0, 0)) * (180/M_PI);
//         double um_meas = atan2(X(4, 0), sqrt(pow(X(0, 0), 2) + pow(X(2, 0), 2))) * (180/M_PI);
        
//         P = Utils<>::do_cart_P(Utils<>::sph2cartcov(R_sph_deg, r_meas, az_meas, um_meas));
//     }

//     F << 1.0, T, 0.0, 0.0, 0.0, 0.0,  // первая строка - уравнение движения (x)
//         0.0, 1.0, 0.0, 0.0, 0.0, 0.0, // вторая строка - изменение скорости (vx)
//         0.0, 0.0, 1.0, T, 0.0, 0.0,   // третья строка - уравнение движения (y)
//         0.0, 0.0, 0.0, 1.0, 0.0, 0.0, // четвертая строка- изменение скорости (vy)
//         0.0, 0.0, 0.0, 0.0, 1.0, T,
//         0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

//     G << (T * T) / 2.0, 0.0, 0.0, // матрица управлении для учета внешних воздействий (xyz)
//         T, 0.0, 0.0,
//         0.0, (T * T) / 2.0, 0.0,
//         0.0, T, 0.0,
//         0.0, 0.0, (T * T) / 2.0,
//         0.0, 0.0, T;

//     Q << sa , 0.0, 0.0, // дисперсия (то есть sigma^2)
//         0.0, sa, 0.0,
//         0.0, 0.0, sa;

//     Cv structCv;

//     structCv.Xe = F * X; // Здесь должно быть еще слагаемое матрицы управления (но мы не управляем, а следим)
//     structCv.Pe = (F * P) * F.transpose() + (G * Q) * G.transpose(); // последнее слагаемое показывает внешнее воздействие на систему

//     Hcv structHcv;
//     M H (3, 6);

//     double xm = structCv.Xe(0);
//     double ym = structCv.Xe(2);
//     double zm = structCv.Xe(4);

//     H << xm / sqrt(pow(xm,2) + pow(ym,2) + pow(zm,2)), 0.0, ym / sqrt(pow(xm,2) + pow(ym,2) + pow(zm,2)), 0.0, zm / sqrt(pow(xm,2) + pow(ym,2) + pow(zm,2)), 0.0,   //  матриц измерении, показывающая отношения измерений и состояния
//                   -ym , 0.0, xm , 0.0, 0.0, 0.0,
//                   -zm * ym/(sqrt(pow(xm,2) + pow(ym,2))), 0.0, -zm * xm/(sqrt(pow(xm,2) + pow(ym,2))), 0.0, sqrt(pow(xm,2) + pow(ym,2)), 0.0; 


//     structHcv.Ze.resize(3, 1);
//     structHcv.Ze << sqrt(pow(xm,2) + pow(ym,2)), atan2(ym,xm), atan2(zm, sqrt(pow(xm,2) + pow(ym,2))); // предсказанное измерение по заданной модели наблюдения.

//     M K;
//     K = (structCv.Pe * H.transpose()) * structHcv.Se.inverse(); // коэф. калмана. Большой, если экстраполяция грубая, а измерения точные.

//     predictStruct.Xe = structCv.Xe;
//     predictStruct.Pe = structCv.Pe;
//     predictStruct.Se = structHcv.Se;
//     predictStruct.Ze = structHcv.Ze;
//     predictStruct.K = K;

//     structHcv.Se = H * structCv.Pe * H.transpose() + R_sph_rad; //  R не понятно в каких ед. передавать
//     M X_p = predictStruct.Xe;
//     return X_p;
// }

// template <class M>
// M ExtendedKalmanFilter<M>::correctEkf(const M& Z)
// {
//     correctStruct.X = predictStruct.Xe + predictStruct.K * (Z - predictStruct.Ze);  
//     correctStruct.P = predictStruct.Pe - (predictStruct.K * predictStruct.Se) * predictStruct.K.transpose();
//     P = correctStruct.P;
//     M X_c = correctStruct.X;
//     return X_c;
// }
