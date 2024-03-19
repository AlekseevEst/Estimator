
#pragma once
#include "utils.h"

template <class M, class TypeEstimator>
struct Track
{
    M Step(M &meas)
    {
        M xe = estimator.predict();
        M x = estimator.correct(meas);
        return x;
    }
    M Step(double dt)
    {
        M xe = estimator.predict();
        estimator.UKfilterMath.X = estimator.UKfilterMath.predictStruct.Xe;

        try
        {
            Utils<M>::СheckingСonditionsMat(estimator.UKfilterMath.predictStruct.Pe)
        }
        catch (const std::runtime_error &e)
        {
            std::cerr << e.what() << '\n';
        }

        estimator.UKfilterMath.P = estimator.UKfilterMath.predictStruct.Pe;
        return estimator.UKfilterMath.P;
    }

private:
    TypeEstimator estimator;
};