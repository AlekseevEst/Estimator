#include "bind_ukf.h"

namespace py = pybind11;

class BindUkf
{
private:

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn,FuncMeasSph> ukf;
 
public:
    
        BindUkf(Eigen::MatrixXd state,
                double dt,
                Eigen::MatrixXd processNoise,
                Eigen::MatrixXd measureNoise,
                double k):              
                          ukf(state, dt,processNoise,measureNoise,k){}
                              

    Eigen::MatrixXd predUkf()
    {
        return ukf.predict();
    }
    Eigen::MatrixXd corrUkf(const Eigen::MatrixXd &Z)
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