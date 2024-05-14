#pragma once
#include "utils.h"

template <class M, class TypeEstimator/*, class TypeDetection*/>
struct Track
{

    Track(const M &X, double t, const M &procNoise, const M &measNoiseMatRadian, Points points) : estimator(X, t, procNoise, measNoiseMatRadian, points) {}

    M Step(const M &meas)
    {
        try
        {
            M xe = estimator.predict();
            M x = estimator.correct(meas);

            return x;
        }
        catch (const std::runtime_error &e)
        {
            std::cerr << e.what() << '\n';
            return M();
        }
    }

    M Step()
    {

        try
        {
            estimator.correctStruct.X = estimator.predict();
            estimator.correctStruct.P = estimator.predictStruct.Pe;
            return estimator.correctStruct.X;
        }

        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return M(); 
        }
    }

private:
    TypeEstimator estimator;
};