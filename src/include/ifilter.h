#pragma once
#include <typeindex>
#include <Eigen/Dense>
#include "structs.h"
template <class M>
struct IFilter
{
    virtual std::pair<M, M> predict(double dt) = 0;
    virtual std::pair<M, M> correct(const M &Z) = 0;


    virtual double likelihood(const M &Z) = 0;
    virtual double distance(const M &Z) = 0;

    virtual std::type_index getModelType() const = 0;

    virtual void setCorrectInfo(const M& X, const M& P) = 0;
    virtual Predict<M> getPredictInfo() = 0;
    virtual Correct<M> getCorrectInfo() = 0;

    virtual ~IFilter() {}
};