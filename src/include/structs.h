#pragma once

template <class M>
struct Predict
{
    M Xe;
    M Pe;
    M Pzz;
    M Se;
    M v; 

    M Ze;
    M K;
};
template <class M>
struct Correct
{
    M X;
    M P;
};

struct Measurement
{
    double r_meas;
    double az_meas;
    double um_meas;
};
struct ParamSigmaPoints
{
    double alpha;
    double beta;
    double kappa;
};

template<class M>
struct Detection
{
    double timePoint;
    M point;
    //M SKO_measurement;
};