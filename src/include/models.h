#pragma once
#include "utils.h"

template <class M>
struct FuncConstVel
{
    enum class VelPos
    {
        POS_X = 0,
        POS_VX,
        POS_Y,
        POS_VY,
        POS_Z,
        POS_VZ,
        SIZE
    };

    M operator()(const M &Xu, double T)
    {

        if (Xu.rows() != ENUM_TO_INT(VelPos::SIZE) || Xu.cols() != 1)
        {
            throw std::invalid_argument("Xu.rows() != ENUM_TO_INT(VelPos::SIZE)  Xu.cols() != 1");
        }

        M F(Xu.rows(), Xu.rows());
        M Xue(Xu.rows(), Xu.cols());

        F << 1.0, T, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, T, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, T,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

        for (int i = 0; i < Xu.cols(); i++)
        {
            Xue.col(i) = F * Xu.col(i);
        }
        // PRINTM(Xue);
        return Xue;
    }
    int getSize()
    {
        return ENUM_TO_INT(VelPos::SIZE);
    }
};

template <class M>
struct FuncConstTurnXY
{
    enum class TurnPos
    {
        POS_X = 0,
        POS_VX,
        POS_Y,
        POS_VY,
        POS_Z,
        POS_VZ,
        POS_OMEGA,
        SIZE
    };

    M operator()(M &Xu, double T)
    {
        if (Xu.rows() != ENUM_TO_INT(TurnPos::SIZE) || Xu.cols() != 1)
        {
            throw std::invalid_argument("Xu.rows() != ENUM_TO_INT(TurnPos::SIZE)  Xu.cols() != 1");
        }
        M F(Xu.rows(), Xu.rows());
        M Xue(Xu.rows(), Xu.cols());

        for (int i = 0; i < Xu.cols(); i++)
        {

            double w = Xu.col(i)(ENUM_TO_INT(TurnPos::POS_OMEGA)) * (M_PI / 180.0);
            if (w == 0)
                w = std::nextafter(0.0, 1.0);
            F << 1.0,   sin(w * T) / w,        0.0, -(1 - cos(w * T)) / w,  0.0,    0.0, 0.0,
                 0.0,     cos(w * T),          0.0,      -sin(w * T),       0.0,    0.0, 0.0,
                 0.0,  (1 - cos(w * T)) / w,   1.0,       sin(w * T) / w,   0.0,    0.0, 0.0,
                 0.0,     sin(w * T),          0.0,       cos(w * T),       0.0,    0.0, 0.0,
                 0.0,           0.0,           0.0,          0.0,           1.0,     T,  0.0,
                 0.0,           0.0,           0.0,          0.0,           0.0,    1.0, 0.0,
                 0.0,           0.0,           0.0,          0.0,           0.0,    0.0, 1.0;

            Xue.col(i) = F * Xu.col(i);
        }
        // PRINTM(Xue);
        return Xue;
    }
    int getSize()
    {
        return ENUM_TO_INT(TurnPos::SIZE);
    }
};
template <class M>
struct FuncConstTurnXZ
{
    enum class TurnPos
    {
        POS_X = 0,
        POS_VX,
        POS_Y,
        POS_VY,
        POS_Z,
        POS_VZ,
        POS_OMEGA,
        SIZE
    };

    M operator()(M &Xu, double T)
    {

        if (Xu.rows() != ENUM_TO_INT(TurnPos::SIZE) || Xu.cols() != 1)
        {
            throw std::invalid_argument("Xu.rows() != ENUM_TO_INT(TurnPos::SIZE)  Xu.cols() != 1");
        }
        M F(Xu.rows(), Xu.rows());
        M Xue(Xu.rows(), Xu.cols());

        for (int i = 0; i < Xu.cols(); i++)
        {

            double w = Xu.col(i)(ENUM_TO_INT(TurnPos::POS_OMEGA)) * (M_PI / 180.0);
            if (w == 0)
                w = std::nextafter(0.0, 1.0);

            F << 1.0,   sin(w * T) / w,     0.0, 0.0, 0.0,  -(1 - cos(w * T)) / w,   0.0,
                 0.0,      cos(w * T),      0.0, 0.0, 0.0,      -sin(w * T),         0.0,
                 0.0,         0.0,          1.0,   T, 0.0,         0.0,              0.0,
                 0.0,         0.0,          0.0, 1.0, 0.0,         0.0,              0.0,
                 0.0, (1 - cos(w * T)) / w, 0.0, 0.0, 1.0,       sin(w * T) / w,     0.0,
                 0.0,      sin(w * T),      0.0, 0.0, 0.0,       cos(w * T),         0.0,
                 0.0,         0.0,          0.0, 0.0, 0.0,        0.0,               1.0;

            Xue.col(i) = F * Xu.col(i);
        }
        // PRINTM(Xue);
        return Xue;
    }
    int getSize()
    {
        return ENUM_TO_INT(TurnPos::SIZE);
    }
};
template <class M>
struct FuncConstAcceleration
{

    enum class AccPos
    {
        POS_X = 0,
        POS_VX,
        POS_AX,
        POS_Y,
        POS_VY,
        POS_AY,
        POS_Z,
        POS_VZ,
        POS_AZ,
        SIZE
    };

    M operator()(M &Xu, double T)
    {
        if (Xu.rows() != ENUM_TO_INT(AccPos::SIZE) || Xu.cols() != 1)
        {
            throw std::invalid_argument("Xu.rows() != ENUM_TO_INT(AccPos::SIZE)  Xu.cols() != 1");
        }
        M F(Xu.rows(), Xu.rows());
        M Xue(Xu.rows(), Xu.cols());

        double T2 = T * T / 2.0;

        for (int i = 0; i < Xu.cols(); i++)
        {

            F << 1.0, T, T2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 1.0, T, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 1.0, T, T2, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 1.0, T, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, T, T2,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, T,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

            Xue.col(i) = F * Xu.col(i);
        }
        // PRINTM(Xue);
        return Xue;
    }
    int getSize()
    {
        return ENUM_TO_INT(AccPos::SIZE);
    }
};

template <class M>
struct FuncMeasSph
{
    enum class VelPos
    {
        POS_X = 0,
        POS_VX,
        POS_Y,
        POS_VY,
        POS_Z,
        POS_VZ,
        SIZE
    };
    enum class TurnPos
    {
        POS_X = 0,
        POS_VX,
        POS_Y,
        POS_VY,
        POS_Z,
        POS_VZ,
        POS_OMEGA,
        SIZE
    };

    enum class AccPos
    {
        POS_X = 0,
        POS_VX,
        POS_AX,
        POS_Y,
        POS_VY,
        POS_AY,
        POS_Z,
        POS_VZ,
        POS_AZ,
        SIZE
    };

    enum class MeasPos
    {
        RANGE = 0,
        AZ,
        EL,
        SIZE
    };

    M operator()(const M &Xue, const M &Z)
    {

        if (Xue.rows() != ENUM_TO_INT(VelPos::SIZE) || Xue.rows() != ENUM_TO_INT(TurnPos::SIZE) || Xue.rows() != ENUM_TO_INT(AccPos::SIZE) || Z.rows() != ENUM_TO_INT(MeasPos::SIZE))
        {
            throw std::invalid_argument("Xue.rows() != VelPos::SIZE || Xue.rows() != TurnPos::SIZE|| Xue.rows() != AccPos::SIZE || Z.rows() != MeasPos::SIZE");
        }

        M Zue(Z.rows(), Xue.cols());
        M zTmp(Z.rows(), Z.cols());
        double range;
        double az;
        double el;

        if (Xue.rows() == ENUM_TO_INT(VelPos::SIZE) || Xue.rows() == ENUM_TO_INT(TurnPos::SIZE))
        {
            for (int i = 0; i < Xue.cols(); i++)
            {

                range = sqrt(pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_X)), 2) + pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Y)), 2) + pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Z)), 2));
                az = atan2(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Y)), Xue.col(i)(ENUM_TO_INT(VelPos::POS_X))) * (180.0 / M_PI);
                el = atan2(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Z)), sqrt(pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_X)), 2) + pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Y)), 2))) * (180.0 / M_PI);

                zTmp << range, az, el;

                Zue.col(i) = zTmp;
                Zue.col(i)(ENUM_TO_INT(MeasPos::AZ)) = Z(ENUM_TO_INT(MeasPos::AZ)) + Utils<M>::ComputeAngleDifference(Zue.col(i)(ENUM_TO_INT(MeasPos::AZ)) * (M_PI / 180.0), Z(ENUM_TO_INT(MeasPos::AZ)) * (M_PI / 180.0)) * (180.0 / M_PI);
            }
            return Zue;
        }

        for (int i = 0; i < Xue.cols(); i++)
        {

            double range = sqrt(pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_X)), 2) + pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Y)), 2) + pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Z)), 2));
            double az = atan2(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Y)), Xue.col(i)(ENUM_TO_INT(AccPos::POS_X))) * (180.0 / M_PI);
            double el = atan2(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Z)), sqrt(pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_X)), 2) + pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Y)), 2))) * (180.0 / M_PI);

            zTmp << range, az, el;

            Zue.col(i) = zTmp;
            Zue.col(i)(ENUM_TO_INT(MeasPos::AZ)) = Z(ENUM_TO_INT(MeasPos::AZ)) + Utils<M>::ComputeAngleDifference(Zue.col(i)(ENUM_TO_INT(MeasPos::AZ)) * (M_PI / 180.0), Z(ENUM_TO_INT(MeasPos::AZ)) * (M_PI / 180.0)) * (180.0 / M_PI);
        }

        // PRINTM(Zue);
        return Zue;
    }
    int getSize()
    {
        return ENUM_TO_INT(MeasPos::SIZE);
    }
};

template <class M>
struct FuncMeasSphVr
{
    enum class VelPos
    {
        POS_X = 0,
        POS_VX,
        POS_Y,
        POS_VY,
        POS_Z,
        POS_VZ,
        SIZE
    };
    enum class TurnPos
    {
        POS_X = 0,
        POS_VX,
        POS_Y,
        POS_VY,
        POS_Z,
        POS_VZ,
        POS_OMEGA,
        SIZE
    };

    enum class AccPos
    {
        POS_X = 0,
        POS_VX,
        POS_AX,
        POS_Y,
        POS_VY,
        POS_AY,
        POS_Z,
        POS_VZ,
        POS_AZ,
        SIZE
    };
    enum class MeasVrPos
    {
        RANGE = 0,
        AZ,
        EL,
        VR,
        SIZE
    };
    M operator()(const M &Xue, const M &Z)
    {

        if (Xue.rows() != ENUM_TO_INT(VelPos::SIZE) || Xue.rows() != ENUM_TO_INT(TurnPos::SIZE) || Xue.rows() != ENUM_TO_INT(AccPos::SIZE) || Z.rows() != ENUM_TO_INT(MeasVrPos::SIZE))
        {
            throw std::invalid_argument("Xue.rows() != VelPos::SIZE || Xue.rows() != TurnPos::SIZE|| Xue.rows() != AccPos::SIZE || Z.rows() != MeasVrPos::SIZE");
        }

        M Zue(Z.rows(), Xue.cols());
        M zTmp(Z.rows(), Z.cols());
        double range;
        double az;
        double el;
        double vr;
        if (Xue.rows() == ENUM_TO_INT(VelPos::SIZE) || Xue.rows() == ENUM_TO_INT(TurnPos::SIZE))
        {
            for (int i = 0; i < Xue.cols(); i++)
            {

                range = sqrt(pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_X)), 2) + pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Y)), 2) + pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Z)), 2));
                az = atan2(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Y)), Xue.col(i)(ENUM_TO_INT(VelPos::POS_X))) * (180.0 / M_PI);
                el = atan2(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Z)), sqrt(pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_X)), 2) + pow(Xue.col(i)(ENUM_TO_INT(VelPos::POS_Y)), 2))) * (180.0 / M_PI);

                vr = (Xue.col(i)(ENUM_TO_INT(VelPos::POS_VX)) * Xue.col(i)(ENUM_TO_INT(VelPos::POS_X)) +
                      Xue.col(i)(ENUM_TO_INT(VelPos::POS_VY)) * Xue.col(i)(ENUM_TO_INT(VelPos::POS_Y)) +
                      Xue.col(i)(ENUM_TO_INT(VelPos::POS_VZ)) * Xue.col(i)(ENUM_TO_INT(VelPos::POS_Z))) /
                     range;

                zTmp << range, az, el, vr;

                Zue.col(i) = zTmp;
                Zue.col(i)(ENUM_TO_INT(MeasVrPos::AZ)) = Z(ENUM_TO_INT(MeasVrPos::AZ)) + Utils<M>::ComputeAngleDifference(Zue.col(i)(ENUM_TO_INT(MeasVrPos::AZ)) * (M_PI / 180.0), Z(ENUM_TO_INT(MeasVrPos::AZ)) * (M_PI / 180.0)) * (180.0 / M_PI);
            }
            return Zue;
        }

        for (int i = 0; i < Xue.cols(); i++)
        {

            double range = sqrt(pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_X)), 2) + pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Y)), 2) + pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Z)), 2));
            double az = atan2(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Y)), Xue.col(i)(ENUM_TO_INT(AccPos::POS_X))) * (180.0 / M_PI);
            double el = atan2(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Z)), sqrt(pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_X)), 2) + pow(Xue.col(i)(ENUM_TO_INT(AccPos::POS_Y)), 2))) * (180.0 / M_PI);

            vr = (Xue.col(i)(ENUM_TO_INT(AccPos::POS_VX)) * Xue.col(i)(ENUM_TO_INT(AccPos::POS_X)) +
                  Xue.col(i)(ENUM_TO_INT(AccPos::POS_VY)) * Xue.col(i)(ENUM_TO_INT(AccPos::POS_Y)) +
                  Xue.col(i)(ENUM_TO_INT(AccPos::POS_VZ)) * Xue.col(i)(ENUM_TO_INT(AccPos::POS_Z))) /
                 range;

            zTmp << range, az, el, vr;

            Zue.col(i) = zTmp;
            Zue.col(i)(ENUM_TO_INT(MeasVrPos::AZ)) = Z(ENUM_TO_INT(MeasVrPos::AZ)) + Utils<M>::ComputeAngleDifference(Zue.col(i)(ENUM_TO_INT(MeasVrPos::AZ)) * (M_PI / 180.0), Z(ENUM_TO_INT(MeasVrPos::AZ)) * (M_PI / 180.0)) * (180.0 / M_PI);
        }

        // PRINTM(Zue);
        return Zue;
    }
    int getSize()
    {
        return ENUM_TO_INT(MeasVrPos::SIZE);
    }
};

template <class M>
struct FuncControlMatrix_XvXYvYZvZ
{
    enum class VelPos
    {
        POS_X = 0,
        POS_VX,
        POS_Y,
        POS_VY,
        POS_Z,
        POS_VZ,
        SIZE
    };
    enum class Pos
    {
        X = 0,
        Y,
        Z,
        SIZE
    };

    M operator()(double T)
    {
        M G(ENUM_TO_INT(VelPos::SIZE), ENUM_TO_INT(Pos::SIZE));
        G << (T * T) / 2.0, 0.0, 0.0,
            T, 0.0, 0.0,
            0.0, (T * T) / 2.0, 0.0,
            0.0, T, 0.0,
            0.0, 0.0, (T * T) / 2.0,
            0.0, 0.0, T;
        return G;
    }
    std::pair<int, int> getSize()
    {
        return std::make_pair(ENUM_TO_INT(VelPos::SIZE), ENUM_TO_INT(Pos::SIZE));
    }
};
template <class M>
struct FuncControlMatrix_XvXYvYZvZW
{
    enum class TurnPos
    {
        POS_X = 0,
        POS_VX,
        POS_Y,
        POS_VY,
        POS_Z,
        POS_VZ,
        POS_OMEGA,
        SIZE
    };
    enum class Pos
    {
        X = 0,
        Y,
        Z,
        OMEGA,
        SIZE
    };

    M operator()(double T)
    {
        M G(ENUM_TO_INT(TurnPos::SIZE), ENUM_TO_INT(Pos::SIZE));
        G << (T * T) / 2.0, 0.0, 0.0, 0.0,
            T, 0.0, 0.0, 0.0,
            0.0, (T * T) / 2.0, 0.0, 0.0,
            0.0, T, 0.0, 0.0,
            0.0, 0.0, (T * T) / 2.0, 0.0,
            0.0, 0.0, T, 0.0,
            0.0, 0.0, 0.0, 1.0;
        return G;
    }

    std::pair<int, int> getSize()
    {
        return std::make_pair(ENUM_TO_INT(TurnPos::SIZE), ENUM_TO_INT(Pos::SIZE));
    }
};
template <class M>
struct FuncControlMatrix_XvXaXYvYaYZvZaZ
{
    enum class AccPos
    {
        POS_X = 0,
        POS_VX,
        POS_AX,
        POS_Y,
        POS_VY,
        POS_AY,
        POS_Z,
        POS_VZ,
        POS_AZ,
        SIZE
    };
    enum class Pos
    {
        X = 0,
        Y,
        Z,
        SIZE
    };
    M operator()(double T)
    {
        M G(ENUM_TO_INT(AccPos::SIZE), ENUM_TO_INT(Pos::SIZE));
        G << (T * T) / 2.0, 0.0, 0.0,
            T, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, (T * T) / 2.0, 0.0,
            0.0, T, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, (T * T) / 2.0,
            0.0, 0.0, T,
            0.0, 0.0, 1.0;
        return G;
    }
    std::pair<int, int> getSize()
    {
        return std::make_pair(ENUM_TO_INT(AccPos::SIZE), ENUM_TO_INT(Pos::SIZE));
    }
};