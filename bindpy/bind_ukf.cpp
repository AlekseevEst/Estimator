#include "bind_ukf.h"

namespace py = pybind11;

class BindUkf
{
private:
    // UnscentedKalmanFilterMath<Eigen::MatrixXd> ukfMath;
    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel,FuncMeasSph> ukf;
    double T;
    Eigen::MatrixXd X;
    Eigen::MatrixXd Z;

public:
    
        BindUkf(Eigen::MatrixXd state,
                Eigen::MatrixXd covMat,
                Eigen::MatrixXd measurement,
                double t,
                Eigen::MatrixXd processNoise,
                Eigen::MatrixXd measureNoise,
                double k):              
                          ukf(state, covMat,measurement,t,processNoise,measureNoise,k){}
                              

    Eigen::MatrixXd predUkf(Eigen::MatrixXd &X, Eigen::MatrixXd &P)
    {
        return ukf.predict(X,P);
    }
    Eigen::MatrixXd corrUkf(Eigen::MatrixXd &Z)
    {
        return ukf.correct(Z);
    }
};

void bind_ukf(pybind11::module &m)
{
    py::class_<BindUkf>(m, "BindUkf")
        .def(py::init<const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, double, const Eigen::MatrixXd&, const Eigen::MatrixXd&, double>())
        .def("predictUkf",&BindUkf::predUkf)
            // py::arg("X"))
        .def("correctUkf",&BindUkf::corrUkf);
            // py::arg("Z"));
}