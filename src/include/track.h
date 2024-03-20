
#pragma once
#include "utils.h"

template <class M, class TypeEstimator>
struct Track
{

    Track(const M &X, double t, const M &procNoise, const M &measNoiseMatRadian, double koef) : estimator(X, t, procNoise, measNoiseMatRadian, koef) {}

    M Step(M &meas)
    {
        M xe = estimator.predict();
        M x = estimator.correct(meas);
        return x;
    }
    M Step(double dt)
    {
        auto FilterMath = estimator.getFilterMath();
        FilterMath.X = estimator.predict();

        try
        {
            Utils<M>::СheckingСonditionsMat(FilterMath.predictStruct.Pe);
        }
        catch (const std::runtime_error &e)
        {
            std::cerr << e.what() << '\n';
        }

        FilterMath.P = FilterMath.predictStruct.Pe;
        estimator.setFilterMath(FilterMath);
        return FilterMath.X;
    }

private:
    TypeEstimator estimator;
};