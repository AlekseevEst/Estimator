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
        // PRINTM(Xue);
        return Xue;
    }
};

template <class M>
struct FuncConstTurn
{

    M operator()(M &Xu, double T)
    {
        M F(ENUM_TO_INT(SizeMat::ROW7),ENUM_TO_INT(SizeMat::COL7));

        M Xue(Xu.rows(),Xu.cols());
        
        for (int i = 0; i < Xu.cols(); i++)
        {

            double w = Xu.col(i)(ENUM_TO_INT(CoordPositionMat::W)) * (M_PI/180.0);
            if (w == 0)
                w = std::nextafter(0.0, 1.0);
            F <<1.0,  sin(w*T)/w,       0.0,   -(1-cos(w*T))/w,    0.0,    0.0,   0.0,
                0.0,  cos(w*T),         0.0,    -sin(w*T),         0.0,    0.0,   0.0,
                0.0,  (1-cos(w*T))/w,   1.0,    sin(w*T)/w,        0.0,    0.0,   0.0,
                0.0,  sin(w*T),         0.0,    cos(w*T),          0.0,    0.0,   0.0,
                0.0,  0.0,              0.0,    0.0,               1.0,      T,   0.0,
                0.0,  0.0,              0.0,    0.0,               0.0,    1.0,   0.0,
                0.0,  0.0,              0.0,    0.0,               0.0,    0.0,   1.0;

        
            Xue.col(i) = F * Xu.col(i); 
        }
        // PRINTM(Xue);
        return Xue;
    }

};

template <class M>
struct FuncConstAcceleration
{

    M operator()(M &Xu, double T)
    {
        double T2 = T*T/2.0;

        M F(ENUM_TO_INT(SizeMat::ROW9),ENUM_TO_INT(SizeMat::COL9));

        M Xue(Xu.rows(),Xu.cols());
        
        for (int i = 0; i < Xu.cols(); i++)
        {

            F <<1.0,   T,         T2,       0.0,    0.0,     0.0,      0.0,     0.0,       0.0,
                0.0,  1.0,         T,       0.0,    0.0,     0.0,      0.0,     0.0,       0.0,
                0.0,  0.0,       1.0,       0.0,    0.0,     0.0,      0.0,     0.0,       0.0,
                0.0,  0.0,       0.0,       1.0,     T,       T2,      0.0,     0.0,       0.0,
                0.0,  0.0,       0.0,       0.0,    1.0,       T,      0.0,     0.0,       0.0,
                0.0,  0.0,       0.0,       0.0,    0.0,     1.0,      0.0,     0.0,       0.0,
                0.0,  0.0,       0.0,       0.0,    0.0,     0.0,      1.0,       T,        T2,
                0.0,  0.0,       0.0,       0.0,    0.0,     0.0,      0.0,     1.0,         T,
                0.0,  0.0,       0.0,       0.0,    0.0,     0.0,      0.0,     0.0,       1.0;

        
            Xue.col(i) = F * Xu.col(i); 
        }
        // PRINTM(Xue);
        return Xue;
    }

};

template <class M>
struct FuncMeasSphCVCT
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
                double az = atan2(Xue.col(i)(ENUM_TO_INT(yPos)), Xue.col(i)(ENUM_TO_INT(xPos))) * (180.0/M_PI);
                double el = atan2(Xue.col(i)(ENUM_TO_INT(zPos)), sqrt(pow(Xue.col(i)(ENUM_TO_INT(xPos)), 2) + pow(Xue.col(i)(ENUM_TO_INT(yPos)), 2))) * (180.0/M_PI);

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
        
            Zue.col(i)(ENUM_TO_INT(azPos)) =  Z(ENUM_TO_INT(azPos)) + Utils<M>::ComputeAngleDifference(Zue.col(i)(ENUM_TO_INT(azPos)) * (M_PI/180.0), Z(ENUM_TO_INT(azPos)) * (M_PI/180.0)) * (180.0/M_PI);

        }
        // PRINTM(Zue); 
        return Zue;
    }
};

template <class M>
struct FuncMeasSphCA
{ 
    M operator()(const M &Xue, const M &Z)
    {   
        MeasPositionMat azPos = MeasPositionMat::AZ;
        CoordPositionMat xPos = CoordPositionMat::X_CA;
        CoordPositionMat yPos = CoordPositionMat::Y_CA;
        CoordPositionMat zPos = CoordPositionMat::Z_CA;
        CoordPositionMat vxPos = CoordPositionMat::VX_CA;
        CoordPositionMat vyPos = CoordPositionMat::VY_CA;
        CoordPositionMat vzPos = CoordPositionMat::VZ_CA;

        M Zue(Z.rows(),Xue.cols());
        for (int i = 0; i < Xue.cols(); i++)
        {
            M zTmp(Z.rows(), 1);

                double range = sqrt(pow(Xue.col(i)(ENUM_TO_INT(xPos)), 2) + pow(Xue.col(i)(ENUM_TO_INT(yPos)), 2) + pow(Xue.col(i)(ENUM_TO_INT(zPos)), 2));
                double az = atan2(Xue.col(i)(ENUM_TO_INT(yPos)), Xue.col(i)(ENUM_TO_INT(xPos))) * (180.0/M_PI);
                double el = atan2(Xue.col(i)(ENUM_TO_INT(zPos)), sqrt(pow(Xue.col(i)(ENUM_TO_INT(xPos)), 2) + pow(Xue.col(i)(ENUM_TO_INT(yPos)), 2))) * (180.0/M_PI);

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

            Zue.col(i)(ENUM_TO_INT(azPos)) =  Z(ENUM_TO_INT(azPos)) + Utils<M>::ComputeAngleDifference(Zue.col(i)(ENUM_TO_INT(azPos)) * (M_PI/180.0), Z(ENUM_TO_INT(azPos)) * (M_PI/180.0)) * (180.0/M_PI);
            
        }
        // PRINTM(Zue); 
        return Zue;
    }
};

template <class M>
struct FuncControlMatrix_XvXYvYZvZ
{
    M operator()(double T)
    {
        M G(ENUM_TO_INT(SizeMat::ROW6), ENUM_TO_INT(SizeMat::COL3));
        G << (T * T) / 2.0,          0.0,            0.0,
                    T,               0.0,            0.0,
                   0.0,         (T * T) / 2.0,       0.0,
                   0.0,               T,             0.0,
                   0.0,              0.0,       (T * T) / 2.0,
                   0.0,              0.0,             T;
        return G;
    }
};
template <class M>
struct FuncControlMatrix_XvXYvYZvZW
{
    M operator()(double T)
    {
        M G(ENUM_TO_INT(SizeMat::ROW7), ENUM_TO_INT(SizeMat::COL4));
        G << (T * T) / 2.0,    0.0,               0.0,          0.0,
                T,             0.0,               0.0,          0.0,
                0.0,      (T * T) / 2.0,          0.0,          0.0,
                0.0,            T,                0.0,          0.0,
                0.0,           0.0,          (T * T) / 2.0,     0.0,
                0.0,           0.0,                T,           0.0,
                0.0,           0.0,               0.0,          1.0;
            return G;
    }
};
template <class M>
struct FuncControlMatrix_XvXaXYvYaYZvZaZ
{
    M operator()(double T)
    {
        M G(ENUM_TO_INT(SizeMat::ROW9), ENUM_TO_INT(SizeMat::COL3));
        G << (T * T) / 2.0,    0.0,               0.0,
                 T,            0.0,               0.0,
                1.0,           0.0,               0.0,         
                0.0,      (T * T) / 2.0,          0.0,         
                0.0,            T,                0.0,
                0.0,           1.0,               0.0,         
                0.0,           0.0,          (T * T) / 2.0,    
                0.0,           0.0,                T,          
                0.0,           0.0,               1.0;
        return G;
    }
};