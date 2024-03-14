#pragma once
#include "utils.h"
#define ENUM_TO_INT(x) static_cast<int>(x)

enum class SizeMat
{
    ROW1 = 1,
    ROW3 = 3,
    ROW6 = 6,
    COL1 = 1,
    COL6 = 6
};
enum class MeasPositionMat
{
    AZ = 1,
    EL = 2
};
enum class CoordPositionMat
{
    X = 0,
    Y = 2,
    Z = 4
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
        MeasPositionMat azPos = MeasPositionMat::AZ;

        CoordPositionMat xPos = CoordPositionMat::X;
        CoordPositionMat yPos = CoordPositionMat::Y;
        CoordPositionMat zPos = CoordPositionMat::Z;

        double range = sqrt(pow(X(ENUM_TO_INT(xPos), 0), 2) + pow(X(ENUM_TO_INT(yPos), 0), 2) + pow(X(ENUM_TO_INT(zPos), 0), 2));
        double az = atan2(X(ENUM_TO_INT(yPos), 0), X(ENUM_TO_INT(xPos), 0));
        double el = atan2(X(ENUM_TO_INT(zPos), 0), sqrt(pow(X(ENUM_TO_INT(xPos), 0), 2) + pow(X(ENUM_TO_INT(yPos), 0), 2)));
        
        az =  Z(ENUM_TO_INT(azPos),0) + Utils<M>::ComputeAngleDifference(az, Z(ENUM_TO_INT(azPos),0));
        M z (Z.rows(), Z.cols());
        z << range,az,el;
        return z;
    }

    std::vector<M> operator()(const std::vector<M> Xue, M Z)
    {   
        MeasPositionMat azPos = MeasPositionMat::AZ;

        CoordPositionMat xPos = CoordPositionMat::X;
        CoordPositionMat yPos = CoordPositionMat::Y;
        CoordPositionMat zPos = CoordPositionMat::Z;

        std::vector<M> Zue(Xue.size());
        for (int i = 0; i < Xue.size(); i++)
        {
            double range = sqrt(pow(Xue[i](ENUM_TO_INT(xPos), 0), 2) + pow(Xue[i](ENUM_TO_INT(yPos), 0), 2) + pow(Xue[i](ENUM_TO_INT(zPos), 0), 2));
            double az = atan2(Xue[i](ENUM_TO_INT(yPos), 0), Xue[i](ENUM_TO_INT(xPos), 0));
            double el = atan2(Xue[i](ENUM_TO_INT(zPos), 0), sqrt(pow(Xue[i](ENUM_TO_INT(xPos), 0), 2) + pow(Xue[i](ENUM_TO_INT(yPos), 0), 2)));

            M zTmp(Z.rows(), Z.cols());
            zTmp << range, az, el;
            Zue[i] = zTmp;

            Zue[i](ENUM_TO_INT(azPos),0) =  Z(ENUM_TO_INT(azPos),0) + Utils<M>::ComputeAngleDifference(Zue[i](ENUM_TO_INT(azPos),0), Z(ENUM_TO_INT(azPos),0));
        }
        return Zue;
    }
};
