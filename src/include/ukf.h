#pragma once

#include <array>
#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <boost/math/distributions/chi_squared.hpp>
#include "utils.h"
#include "models.h"
#include "ukfMath.h"

template <class M,
                template <typename> class StateFunc,
                template <typename> class MeasurementFunc>

struct UnscentKalmanfilter
{
    private:

    UnscentedKalmanFilterMath<M>  UKfilterMath;
    
    public:

    double T;
    StateFunc<M> stateFunc;
    MeasurementFunc<M> measFunc; 

    M state;
    M stateCovariance;
    M measurement;
    M predict(const M& X, const M& P);
    M correct(const M& Z);
    UnscentKalmanfilter(M X, M Z,double t, M procNoise, M measNoiseMatRadian ,double koef): state(X), measurement(Z),T(t), UKfilterMath(stateCovariance,
                                                                                             measNoiseMatRadian, procNoise , t, koef){}

};

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc>
M UnscentKalmanfilter<M, StateFunc, MeasurementFunc>::predict(const M& X, const M& P)
{
    M make_P_cart(const M &P, const M &X);
    M rootMat = UKfilterMath.sqrt_matrix_p(P);
    
    std::vector<M> sigmaVectors = UKfilterMath.doSigmaVectors(state,rootMat);
    std::vector<M> weightVectors = UKfilterMath.calculationVectorWeights();
    std::vector<M> ExtrapolatedSigmaVectors = stateFunc(sigmaVectors,T);
    UKfilterMath.doExtrapolatedStateVector(ExtrapolatedSigmaVectors, weightVectors);
    UKfilterMath.doCovMatExtrapolatedStateVector(ExtrapolatedSigmaVectors, weightVectors);
    std::vector<M> ExtrapolatedSigmaVectorsSph = UKfilterMath.doSigmaVectorsSph(ExtrapolatedSigmaVectors);
    UKfilterMath.doExtrapolatedStateVectorSph(ExtrapolatedSigmaVectorsSph, weightVectors);
    UKfilterMath.doCovMatExtrapolatedStateVectorSph(ExtrapolatedSigmaVectorsSph, weightVectors);
    UKfilterMath.calcGainFilter(ExtrapolatedSigmaVectors, ExtrapolatedSigmaVectorsSph, weightVectors);

    return UKfilterMath.predictStruct.Xe;
}

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc>
M UnscentKalmanfilter<M, StateFunc, MeasurementFunc>::correct(const M& Z)
{   
    M correctState = UKfilterMath.correctState(Z);
    M correctCov = UKfilterMath.correctCov();
    return correctState;
}
