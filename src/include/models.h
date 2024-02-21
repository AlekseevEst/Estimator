#pragma once

// #include "Eigen/Dense"

template <class M>
struct FuncConstVel
{

    std::vector<M> operator()(std::vector<M> Xu, double T)
    {

        M F(Xu[0].rows(), Xu[0].rows());

        F << 1.0, T, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, T, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, T,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

        //------------ЭКСТРАПОЛЯЦИЯ СИГМА-ВЕКТОРОВ-----------------
        std::vector<M> Xue(Xu.size());

        for (int i = 0; i < Xu.size(); i++)
        {
            Xue[i] = F * Xu[i]; 
           
        }
        
        return Xue;
    }
};

template <class M>
struct FuncMeasSph
{

    M operator()(const M &Z)
    {

    }
};
