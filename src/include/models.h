#pragma once
#include "utils.h"

template <class M>
struct FuncConstVel
{
    M operator()(const M &Xu, double T)
    {

        M F(ENUM_TO_INT(SizeMat::ROW6),ENUM_TO_INT(SizeMat::COL6));

        F << 1.0, T, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, T, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, T,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

        M Xue(Xu.rows(),Xu.cols());

        for (int i = 0; i < Xu.cols(); i++)
        {
            Xue.col(i) = F * Xu.col(i); 
        }
        
        return Xue;
    }
};

template <class M>
struct FuncConstTurn
{

    M operator()(const M &Xu, double T)
    {
        M F(ENUM_TO_INT(SizeMat::ROW7),ENUM_TO_INT(SizeMat::COL7));

        M Xue(Xu.rows(),Xu.cols());
        double w = Xu.col(ENUM_TO_INT(SizeMat::COL0))(ENUM_TO_INT(CoordPositionMat::W));

        for (int i = 0; i < Xu.cols(); i++)
        {

            F <<1.0,  sin(w*T)/w,       0.0,   -(1-cos(w*T))/w,    0.0,    0.0,   0.0,
                0.0,  cos(w*T),         0.0,    -sin(w*T),         0.0,    0.0,   0.0,
                0.0,  (1-cos(w*T))/w,   1.0,    sin(w*T)/w,        0.0,    0.0,   0.0,
                0.0,  sin(w*T),         0.0,    cos(w*T),          0.0,    0.0,   0.0,
                0.0,  0.0,              0.0,    0.0,               1.0,      T,   0.0,
                0.0,  0.0,              0.0,    0.0,               0.0,    1.0,   0.0,
                0.0,  0.0,              0.0,    0.0,               0.0,    0.0,   1.0;

            
            Xue.col(i) = F * Xu.col(i); 
        }
        
        return Xue;
    }

};

template <class M>
struct FuncMeasSph
{ 

    M operator()(const M &Xue, const M &Z)
    {   
        MeasPositionMat azPos = MeasPositionMat::AZ;
        CoordPositionMat xPos = CoordPositionMat::X;
        CoordPositionMat yPos = CoordPositionMat::Y;
        CoordPositionMat zPos = CoordPositionMat::Z;
        CoordPositionMat vxPos = CoordPositionMat::VX;
        CoordPositionMat vyPos = CoordPositionMat::VY;
        CoordPositionMat vzPos = CoordPositionMat::VZ;

        M Zue(Z.rows(),Xue.cols());
        for (int i = 0; i < Xue.cols(); i++)
        {
            M zTmp(Z.rows(), 1);

                double range = sqrt(pow(Xue.col(i)(ENUM_TO_INT(xPos)), 2) + pow(Xue.col(i)(ENUM_TO_INT(yPos)), 2) + pow(Xue.col(i)(ENUM_TO_INT(zPos)), 2));
                double az = atan2(Xue.col(i)(ENUM_TO_INT(yPos)), Xue.col(i)(ENUM_TO_INT(xPos)));
                double el = atan2(Xue.col(i)(ENUM_TO_INT(zPos)), sqrt(pow(Xue.col(i)(ENUM_TO_INT(xPos)), 2) + pow(Xue.col(i)(ENUM_TO_INT(yPos)), 2)));

            if (Z.rows() == ENUM_TO_INT(SizeMat::ROW3))
            {
                zTmp << range, az, el;
            }    
            else
            {        
                double vr = (Xue.col(i)(ENUM_TO_INT(vxPos)) * Xue.col(i)(ENUM_TO_INT(xPos)) + \
                            Xue.col(i)(ENUM_TO_INT(vyPos)) * Xue.col(i)(ENUM_TO_INT(yPos)) + \
                            Xue.col(i)(ENUM_TO_INT(vzPos)) * Xue.col(i)(ENUM_TO_INT(zPos)))/range;

            
                zTmp << range, az, el, vr;
            }
            Zue.col(i) = zTmp; 

            Zue.col(i)(ENUM_TO_INT(azPos)) =  Z(ENUM_TO_INT(azPos)) + Utils<M>::ComputeAngleDifference(Zue.col(i)(ENUM_TO_INT(azPos)), Z(ENUM_TO_INT(azPos)));
        }
        return Zue;
    }
};


