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
    M step(const M& Z, double dt);
    M step(double dt);
    double computeDistance(const M&Z);

    std::type_index getTypeModelState() const override;

    UnscentedKalmanfilter(const M& X, const M& procNoise, const M& measNoise , ParamSigmaPoints paramSigmaPoints): 
                                                                                                            UKfilterMath(measNoise, procNoise, paramSigmaPoints)    
                                                                                                            {                                                                                                   
                                                                                                                correctStruct.X = X;
                                                                                                                correctStruct.P = UKfilterMath.make_P0_cart(X);                                                                                 
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
    M sigmaVectors = UKfilterMath.doSigmaVectors(correctStruct.X,correctStruct.P);
    extrapolatedStateSigmaVectors = stateFunc(sigmaVectors, dt);
    predictStruct.Xe = UKfilterMath.doExtrapolatedStateVector(extrapolatedStateSigmaVectors);
    M G = controlFunc(dt);
    predictStruct.Pe = UKfilterMath.doCovMatExtrapolatedStateVector(extrapolatedStateSigmaVectors, predictStruct.Xe, G);

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

    correctStruct.X = UKfilterMath.correctState(predictStruct.Xe,Z, predictStruct.Ze, predictStruct.K);
    correctStruct.P = UKfilterMath.correctCov(predictStruct.Pe, predictStruct.K, predictStruct.Se);


    return correctStruct.X;
}

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc,
          template <typename> class ControlFunc>
M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc, ControlFunc>::step(const M &Z, double dt)
{
    predict(dt);
    return correct(Z);
}

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc,
          template <typename> class ControlFunc>
M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc, ControlFunc>::step(double dt)
{
    correctStruct.X = predict(dt);
    correctStruct.P = predictStruct.Pe;
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
