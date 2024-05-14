#pragma once

template <class M>
struct Predict
{
    M Xe;
    M Pe;
    M Se;

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
struct Points
{
    double alpha;
    double beta;
    double kappa;
};
