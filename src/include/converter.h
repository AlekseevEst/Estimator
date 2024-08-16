#pragma once
#include "models.h"
#include "keyMap.h"
#include <functional>
#include <vector>
// #include <Eigen/Sparse>
template <class M>
struct Converter
{
    FuncConstVel<M> modelCv;
    FuncConstAcceleration<M> modelCa;
    FuncConstTurn<M> modelCt;
    using SpMat = Eigen::SparseMatrix<double>;
    using T = Eigen::Triplet<double>;


    std::unordered_map<namespaceKeyMap::KeyMap, std::function<M(M)>, namespaceKeyMap::Hash_fn> m;

    Converter()

    {
        SpMat HVelTurn(6, 7);
        SpMat HVelAcc(6, 9);
        SpMat HTurnAcc(7, 9);
        std::vector<T> tripletList;
        tripletList.reserve(6);
        tripletList.push_back(T(0, 0, 1.0)); // использовать POS_X и т.д
        tripletList.push_back(T(1, 1, 1.0));
        tripletList.push_back(T(2, 2, 1.0));
        tripletList.push_back(T(3, 3, 1.0));
        tripletList.push_back(T(4, 4, 1.0));
        tripletList.push_back(T(5, 5, 1.0));
        HVelTurn.setFromTriplets(tripletList.begin(), tripletList.end());
        tripletList.clear();

        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 1, 1.0));
        tripletList.push_back(T(2, 3, 1.0));
        tripletList.push_back(T(3, 4, 1.0));
        tripletList.push_back(T(4, 6, 1.0));
        tripletList.push_back(T(5, 7, 1.0));
        HVelAcc.setFromTriplets(tripletList.begin(), tripletList.end());
        HTurnAcc.setFromTriplets(tripletList.begin(), tripletList.end());
        tripletList.clear();

        m[{typeid(modelCv), typeid(modelCt)}] = [HVelTurn](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {   
                res = HVelTurn.transpose() * matStateOrCov;
                return res;
            }
            res = HVelTurn.transpose() * matStateOrCov * HVelTurn;
            return res;
        };

        m[{typeid(modelCt), typeid(modelCv)}] = [HVelTurn](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HVelTurn * matStateOrCov;
                return res;
            }
            res = HVelTurn * matStateOrCov * HVelTurn.transpose();
            return res;
        };

        m[{typeid(modelCv), typeid(modelCa)}] = [HVelAcc](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HVelAcc.transpose() * matStateOrCov;
                return res;
            }
            res = HVelAcc.transpose() * matStateOrCov * HVelAcc;
            return res;
        };

        m[{typeid(modelCa), typeid(modelCv)}] = [HVelAcc](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HVelAcc * matStateOrCov;
                return res;
            }
            res = HVelAcc * matStateOrCov * HVelAcc.transpose();
            return res;
        };

        m[{typeid(modelCt), typeid(modelCa)}] = [HTurnAcc](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HTurnAcc.transpose() * matStateOrCov;
                return res;
            }
            res = HTurnAcc.transpose() * matStateOrCov * HTurnAcc;
            return res;
        };

        m[{typeid(modelCa), typeid(modelCt)}] = [HTurnAcc](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HTurnAcc * matStateOrCov;
                return res;
            }
            res = HTurnAcc * matStateOrCov * HTurnAcc.transpose();
            return res;
        };
    }
};
