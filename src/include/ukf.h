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

struct UnscentedKalmanfilter
{
    private:

    UnscentedKalmanFilterMath<M>  UKfilterMath;
    
    public:

    double T;
    StateFunc<M> stateFunc;
    MeasurementFunc<M> measFunc; 

    // M state;
    // M stateCovariance;
    M measurement;
    M predict();
    M correct(const M& Z);
    UnscentedKalmanfilter(const M& X, const M& Z, double t, const M& procNoise, const M& measNoiseMatRadian ,double koef): 
                                                                                            measurement(Z),T(t),
                                                                                            UKfilterMath(X,measNoiseMatRadian, procNoise, t, koef){}
};

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc>

M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc>::predict()
{
    UKfilterMath.make_P_cart();
    M rootMat = UKfilterMath.sqrt_matrix_p();
    std::vector<M> sigmaVectors = UKfilterMath.doSigmaVectors(rootMat);
    std::vector<double> weightVectors = UKfilterMath.calculationVectorWeights();
    std::vector<M> ExtrapolatedSigmaVectors = stateFunc(sigmaVectors,T);
    UKfilterMath.doExtrapolatedStateVector(ExtrapolatedSigmaVectors, weightVectors);
    UKfilterMath.doCovMatExtrapolatedStateVector(ExtrapolatedSigmaVectors, weightVectors);
    std::vector<M> ExtrapolatedSigmaVectorsSph = UKfilterMath.doSigmaVectorsSph(ExtrapolatedSigmaVectors, measurement);
    UKfilterMath.doExtrapolatedStateVectorSph(ExtrapolatedSigmaVectorsSph, weightVectors);
    UKfilterMath.doCovMatExtrapolatedStateVectorSph(ExtrapolatedSigmaVectorsSph, weightVectors);
    UKfilterMath.calcGainFilter(ExtrapolatedSigmaVectors, ExtrapolatedSigmaVectorsSph, weightVectors);

    return UKfilterMath.predictStruct.Xe;
}

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc>
M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc>::correct(const M& Z)
{   
    M correctState = UKfilterMath.correctState(Z);
    M correctCov = UKfilterMath.correctCov();
    return correctState;
}
