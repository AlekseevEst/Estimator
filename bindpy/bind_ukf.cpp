#include "bind_ukf.h"

namespace py = pybind11;

class BindUkf
{
private:

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel,FuncMeasSph> ukf;


public:
    
        BindUkf(Eigen::MatrixXd state,
                double t,
                Eigen::MatrixXd processNoise,
                Eigen::MatrixXd measureNoise,
                double k):              
                          ukf(state,t,processNoise,measureNoise,k){}
                              

    Eigen::MatrixXd predUkf(Eigen::MatrixXd &Z)
    {
        return ukf.predict(Z);
    }
    Eigen::MatrixXd corrUkf(Eigen::MatrixXd &Z)
    {
        return ukf.correct(Z);
    }
};

void bind_ukf(pybind11::module &m)
{
    py::class_<BindUkf>(m, "BindUkf")
        .def(py::init<const Eigen::MatrixXd&, double, const Eigen::MatrixXd&, const Eigen::MatrixXd&, double>())
        .def("predictUkf",&BindUkf::predUkf)
            // py::arg("X"))
        .def("correctUkf",&BindUkf::corrUkf);
            // py::arg("Z"));
}