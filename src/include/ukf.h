#pragma once

#include <array>
#include <cmath>
#include <vector>
#include <iostream>
#include "utils.h"
#include "models.h"
#include "ukfMath.h"

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc,
          template <typename> class ControlFunc>

struct UnscentedKalmanfilter
{
private:
    UnscentedKalmanFilterMath<M> UKfilterMath;
    M extrapolatedStateSigmaVectors;
    M Qp;
    StateFunc<M> stateFunc;
    MeasurementFunc<M> measFunc;
    ControlFunc<M> controlFunc;
    

public:
    Predict<M> predictStruct;
    Correct<M> correctStruct;

    M predict(double dt);
    M correct(const M &Z);

    UnscentedKalmanfilter(const M& X, const M& procNoise, const M& measNoiseMatRadian , Points points): 
                                                                                                            UKfilterMath(measNoiseMatRadian, procNoise, points)    
                                                                                                            {                                                                                                   
                                                                                                                correctStruct.X = X;
                                                                                                                correctStruct.P.resize(X.rows(),X.rows());
                                                                                                                correctStruct.P.setZero();
                                                                                                            }
};

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc,
          template <typename> class ControlFunc>

M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc, ControlFunc>::predict(double dt)
{
    correctStruct.P = UKfilterMath.make_P_cart(correctStruct.P, correctStruct.X);
    M sigmaVectors = UKfilterMath.doSigmaVectors(correctStruct.X,correctStruct.P);
    extrapolatedStateSigmaVectors = stateFunc(sigmaVectors, dt);
    predictStruct.Xe = UKfilterMath.doExtrapolatedStateVector(extrapolatedStateSigmaVectors);
    M G = controlFunc(dt);
    predictStruct.Pe = UKfilterMath.doCovMatExtrapolatedStateVector(extrapolatedStateSigmaVectors, predictStruct.Xe, G, dt);

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
    predictStruct.Se = UKfilterMath.doCovMatExtrapolatedMeasVector(extrapolatedMeasSigmaVectors, predictStruct.Ze);
    predictStruct.K = UKfilterMath.calcGainFilter(extrapolatedStateSigmaVectors, predictStruct.Xe, extrapolatedMeasSigmaVectors,predictStruct.Ze, predictStruct.Se);
    correctStruct.X = UKfilterMath.correctState(predictStruct.Xe,Z, predictStruct.Ze, predictStruct.K);
    correctStruct.P = UKfilterMath.correctCov(predictStruct.Pe, predictStruct.K, predictStruct.Se);

    return correctStruct.X;
}
