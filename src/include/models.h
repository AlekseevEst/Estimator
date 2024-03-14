#pragma once
#include "utils.h"
#define ENUM_TO_INT(x) static_cast<int>(x)

enum class SizeMat {
    ROW1 = 1,
    ROW3 = 3,
    ROW6 = 6,
    COL1 = 1,
    COL6 = 6
};
template <class M>
struct FuncConstVel
{
        M operator()(M X, double T)
    {
        SizeMat rows = SizeMat::ROW6;
        SizeMat cols = SizeMat::COL6;
        M F(ENUM_TO_INT(rows),ENUM_TO_INT(cols));

        F << 1.0, T, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, T, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, T,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
        
        return F*X;
    }

    std::vector<M> operator()(std::vector<M> Xu, double T)
    {
        SizeMat rows = SizeMat::ROW6;
        SizeMat cols = SizeMat::COL6;
        M F(ENUM_TO_INT(rows),ENUM_TO_INT(cols));

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
        M operator()(const M X, M Z)
    {   
        SizeMat rows = SizeMat::ROW3;
        SizeMat cols = SizeMat::COL1;

        double range = sqrt(pow(X(0, 0), 2) + pow(X(2, 0), 2) + pow(X(4, 0), 2));
        double az = atan2(X(2, 0), X(0, 0));
        double el = atan2(X(4, 0), sqrt(pow(X(0, 0), 2) + pow(X(2, 0), 2)));
        
        az =  Z(1,0) + Utils<M>::ComputeAngleDifference(az, Z(1,0));
        M z (Z.rows(), Z.cols());
        z << range,az,el;
        return z;
    }


    std::vector<M> operator()(const std::vector<M> Xue, M Z)
    {   
        SizeMat rows = SizeMat::ROW3;
        SizeMat cols = SizeMat::COL1;
        std::vector<M> Zue(Xue.size());
        for (int i = 0; i < Xue.size(); i++)
        {
            double range = sqrt(pow(Xue[i](0, 0), 2) + pow(Xue[i](2, 0), 2) + pow(Xue[i](4, 0), 2));
            double az = atan2(Xue[i](2, 0), Xue[i](0, 0));
            double el = atan2(Xue[i](4, 0), sqrt(pow(Xue[i](0, 0), 2) + pow(Xue[i](2, 0), 2)));

            M zTmp(ENUM_TO_INT(rows), ENUM_TO_INT(cols));
            zTmp << range, az, el;
            Zue[i] = zTmp;

            Zue[i](1,0) =  Z(1,0) + Utils<M>::ComputeAngleDifference(Zue[i](1,0), Z(1,0));
        }
        return Zue;
    }
};


    