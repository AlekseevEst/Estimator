#pragma once
#include "models.h"
#include "keyMap.h"
#include <functional>
#include <vector>
#include <Eigen/Sparse>

template <class M>
struct Converter
{
    FuncConstVel<M> modelCv;
    FuncConstAcceleration<M> modelCa;
    FuncConstTurnXZ<M> modelCtXz;
    FuncConstTurnXY<M> modelCtXy;
    FuncBalreentry<M> modelBal;
    using SpMat = Eigen::SparseMatrix<double>;
    using T = Eigen::Triplet<double>;

    std::unordered_map<namespaceKeyMap::KeyMap, std::function<M(M)>, namespaceKeyMap::Hash_fn> m;

    Converter()

    {
        SpMat HVelTurn(6, 7);
        SpMat HVelAcc(6, 9);
        SpMat HTurnAcc(7, 9);
        SpMat HBalAcc(7,9);
        SpMat HVelBal(6,7);
        SpMat HTurnBal(7,7);
       
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


        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 1, 1.0));
        tripletList.push_back(T(2, 3, 1.0));
        tripletList.push_back(T(3, 4, 1.0));
        tripletList.push_back(T(5, 6, 1.0));
        tripletList.push_back(T(6, 7, 1.0));
        HBalAcc.setFromTriplets(tripletList.begin(), tripletList.end());
        tripletList.clear();

        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 1, 1.0));
        tripletList.push_back(T(2, 2, 1.0));
        tripletList.push_back(T(3, 3, 1.0));
        tripletList.push_back(T(4, 5, 1.0));
        tripletList.push_back(T(5, 6, 1.0));
        HVelBal.setFromTriplets(tripletList.begin(), tripletList.end());
        HTurnBal.setFromTriplets(tripletList.begin(), tripletList.end());
        tripletList.clear();

        m[{typeid(modelCv), typeid(modelCv)}] = [](const M &matStateOrCov)
        {
            return matStateOrCov;
        };

        m[{typeid(modelCtXy), typeid(modelCtXy)}] = [](const M &matStateOrCov)
        {
            return matStateOrCov;
        };

        m[{typeid(modelCtXz), typeid(modelCtXz)}] = [](const M &matStateOrCov)
        {
            return matStateOrCov;
        };
        m[{typeid(modelCtXy), typeid(modelCtXz)}] = [](const M &matStateOrCov)
        {
            return matStateOrCov;
        };
        m[{typeid(modelCtXz), typeid(modelCtXy)}] = [](const M &matStateOrCov)
        {
            return matStateOrCov;
        };

        m[{typeid(modelCa), typeid(modelCa)}] = [](const M &matStateOrCov)
        {
            return matStateOrCov;
        };
        
        m[{typeid(modelBal), typeid(modelBal)}] = [](const M &matStateOrCov)
        {
            return matStateOrCov;
        };


        m[{typeid(modelCv), typeid(modelCtXy)}] = [HVelTurn](const M &matStateOrCov)
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

        m[{typeid(modelCv), typeid(modelCtXz)}] = [HVelTurn](const M &matStateOrCov)
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

        m[{typeid(modelCtXy), typeid(modelCv)}] = [HVelTurn](const M &matStateOrCov)
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

        m[{typeid(modelCtXz), typeid(modelCv)}] = [HVelTurn](const M &matStateOrCov)
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

        m[{typeid(modelCtXy), typeid(modelCa)}] = [HTurnAcc](const M &matStateOrCov)
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
        m[{typeid(modelCtXz), typeid(modelCa)}] = [HTurnAcc](const M &matStateOrCov)
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

        m[{typeid(modelCa), typeid(modelCtXy)}] = [HTurnAcc](const M &matStateOrCov)
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

        m[{typeid(modelCa), typeid(modelCtXz)}] = [HTurnAcc](const M &matStateOrCov)
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


        m[{typeid(modelCv), typeid(modelBal)}] = [HVelBal](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HVelBal.transpose() * matStateOrCov;
                return res;
            }
            res = HVelBal.transpose() * matStateOrCov * HVelBal;
            return res;

        };


        m[{typeid(modelBal), typeid(modelCv)}] = [HVelBal](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HVelBal * matStateOrCov;
                return res;
            }
            res = HVelBal * matStateOrCov * HVelBal.transpose();
            return res;
        };


        m[{typeid(modelBal), typeid(modelCa)}] = [HBalAcc](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HBalAcc.transpose() * matStateOrCov;
                return res;
            }
            res = HBalAcc.transpose() * matStateOrCov * HBalAcc;
            return res;
        };

        m[{typeid(modelCa), typeid(modelBal)}] = [HBalAcc](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HBalAcc * matStateOrCov;
                return res;
            }
            res = HBalAcc * matStateOrCov * HBalAcc.transpose();
            return res;
        };

        m[{typeid(modelCtXy), typeid(modelBal)}] = [HTurnBal](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HTurnBal.transpose() * matStateOrCov;
                return res;
            }
            res = HTurnBal.transpose() * matStateOrCov * HTurnBal;
            return res;
        };


        m[{typeid(modelBal), typeid(modelCtXy)}] = [HTurnBal](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HTurnBal * matStateOrCov;
                return res;
            }
            res = HTurnBal * matStateOrCov * HTurnBal.transpose();
            return res;
        };


        m[{typeid(modelCtXz), typeid(modelBal)}] = [HTurnBal](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HTurnBal.transpose() * matStateOrCov;
                return res;
            }
            res = HTurnBal.transpose() * matStateOrCov * HTurnBal;
            return res;
        };

        m[{typeid(modelBal), typeid(modelCtXz)}] = [HTurnBal](const M &matStateOrCov)
        {
            M res;
            if (matStateOrCov.cols() == 1)
            {
                res = HTurnBal * matStateOrCov;
                return res;
            }
            res = HTurnBal * matStateOrCov * HTurnBal.transpose();
            return res;
        };

    }
};
