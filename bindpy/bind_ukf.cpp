#include "bind_ukf.h"

namespace py = pybind11;

class BindUkf
{
private:

public:

    UnscentKalmanfilter<Eigen::MatrixXd,FuncConstVel,FuncMeasSph> ukf;
    
    // BindUkf()
    // {

    // }

    Eigen::MatrixXd predUkf(Eigen::MatrixXd &X, Eigen::MatrixXd &P)
    {
        // return ukf.predict(X,P);
        
    }
    Eigen::MatrixXd corrUkf(Eigen::MatrixXd &Z)
    {
    //    return ukf.correctUkf(Z);
    }
};

void bind_ukf(pybind11::module &m)
{
    // py::class_<BindUkf>(m, "BindUkf")
    //     .def(py::init<>())
    //     .def("predictUkf",&BindUkf::predUkf)
    //         // py::arg("X"))
    //     .def("correctUkf",&BindUkf::corrUkf);
    //         // py::arg("Z"));
}