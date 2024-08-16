#pragma once
#include <typeindex>
#include <Eigen/Dense>
template<class M>
struct IFilter
{
public:
    Predict<M> predictStruct;
    Correct<M> correctStruct;

    virtual M predict(double dt) = 0;
    virtual M correct(const M &Z) = 0;

    virtual std::type_index getTypeModelState() const = 0;
    virtual ~IFilter() {}
};
