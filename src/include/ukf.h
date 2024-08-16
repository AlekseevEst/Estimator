#pragma once

#include <array>
#include <cmath>
#include <vector>
#include <iostream>
#include "utils.h"
#include "models.h"
#include "ukfMath.h"
#include "ifilter.h"

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc,
          template <typename> class ControlFunc>

struct UnscentedKalmanfilter : public IFilter<M>
{
private:
    UnscentedKalmanFilterMath<M> UKfilterMath;
    M extrapolatedStateSigmaVectors;
    M Qp;
    StateFunc<M> stateFunc;
    MeasurementFunc<M> measFunc;
    ControlFunc<M> controlFunc;
    
public:

    using IFilter<M>::predictStruct;
    using IFilter<M>::correctStruct;

    // void initFilter();
    M predict(double dt) override;
    M correct(const M &Z) override;
    double computeDistance(const M&Z);

    std::type_index getTypeModelState() const override;

    UnscentedKalmanfilter(const M& X, const M& procNoise, const M& measNoise , ParamSigmaPoints paramSigmaPoints): 
                                                                                                            UKfilterMath(measNoise, procNoise, paramSigmaPoints)    
                                                                                                            {                                                                                                   
                                                                                                                correctStruct.X = X;
                                                                                                                correctStruct.P = UKfilterMath.make_P0_cart(X);
                                                                                                                std::cout << "Коструктор N Фильтра в Ukf";
                                                                                                                PRINTM(correctStruct.X);
                                                                                                                PRINTM(correctStruct.P);
                                                                                                                                                                                                                            
                                                                                                            }
};
// template <class M,
//           template <typename> class StateFunc,
//           template <typename> class MeasurementFunc,
//           template <typename> class ControlFunc>

// void UnscentedKalmanfilter<M, StateFunc, MeasurementFunc, ControlFunc>::init()
// {
// }


template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc,
          template <typename> class ControlFunc>

M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc, ControlFunc>::predict(double dt)
{
    // correctStruct.P = UKfilterMath.make_P_cart(correctStruct.P, correctStruct.X);
    M sigmaVectors = UKfilterMath.doSigmaVectors(correctStruct.X,correctStruct.P);
    std::cout<< "Проверка coorectStruct.X и Р В UKF predict после doSigmavec" << std::endl<< std::endl;
    PRINTM(correctStruct.X);
    PRINTM(correctStruct.P);
    extrapolatedStateSigmaVectors = stateFunc(sigmaVectors, dt);
    predictStruct.Xe = UKfilterMath.doExtrapolatedStateVector(extrapolatedStateSigmaVectors);
    M G = controlFunc(dt);
    predictStruct.Pe = UKfilterMath.doCovMatExtrapolatedStateVector(extrapolatedStateSigmaVectors, predictStruct.Xe, G);
    std::cout<< "Проверка predictStruct.Xe" << std::endl<< std::endl;
    PRINTM(predictStruct.Xe);
    std::cout<< "Проверка PredictStruct.Pe" << std::endl<< std::endl;
    PRINTM(predictStruct.Pe);

    return predictStruct.Xe;
}

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc,
          template <typename> class ControlFunc>
M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc, ControlFunc>::correct(const M& Z)
{   
    M extrapolatedMeasSigmaVectors = measFunc(extrapolatedStateSigmaVectors, Z); 
    predictStruct.Ze = UKfilterMath.doExtrapolatedMeasVector(extrapolatedMeasSigmaVectors);
    predictStruct.Pzz = UKfilterMath.doCovMatExtrapolatedMeasVector(extrapolatedMeasSigmaVectors, predictStruct.Ze);
    predictStruct.Se = UKfilterMath.doCovMatInnovation(predictStruct.Pzz);
    predictStruct.K = UKfilterMath.calcGainFilter(extrapolatedStateSigmaVectors, predictStruct.Xe, extrapolatedMeasSigmaVectors,predictStruct.Ze, predictStruct.Se);
    std::cout<< "Проверка coorectStruct.X и Р В UKF correct до correctState и correctCov" << std::endl<< std::endl;
    PRINTM(correctStruct.X);
    PRINTM(correctStruct.P);
    correctStruct.X = UKfilterMath.correctState(predictStruct.Xe,Z, predictStruct.Ze, predictStruct.K);
    correctStruct.P = UKfilterMath.correctCov(predictStruct.Pe, predictStruct.K, predictStruct.Se);
    std::cout<< "Проверка coorectStruct.X и Р В UKF correct После correctState и correctCov" << std::endl<< std::endl;
    PRINTM(correctStruct.X);
    PRINTM(correctStruct.P);

    return correctStruct.X;
}

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc,
          template <typename> class ControlFunc>
double UnscentedKalmanfilter<M, StateFunc, MeasurementFunc, ControlFunc>::computeDistance(const M& Z)
{   
    M v = Z - predictStruct.Ze; //невязка
    double distance = v.transpose() * predictStruct.Se.inverse() + predictStruct.Se.determinant();
    
    return distance;
}

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc,
          template <typename> class ControlFunc>
std::type_index UnscentedKalmanfilter<M, StateFunc, MeasurementFunc, ControlFunc>::getTypeModelState() const
{
    return typeid(StateFunc<M>);
}
