#include "bind_ukf.h"
// #include "ukf.h"
// #include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
// #include <pybind11/eigen.h>
namespace py = pybind11;

class BindUkf
{
private:
    unscent_filter<> uf;
    
public:
    MatrixXd predUkf(MatrixXd &X)
    {
        return uf.predictUkf(X);
    }
    MatrixXd corrUkf(MatrixXd &Z)
    {
       return uf.correctUkf(Z);
    }
};

void bind_ukf(pybind11::module &m)
{
    py::class_<BindUkf>(m, "BindUkf")
        .def(py::init<>())
        .def("predictUkf",&BindUkf::predUkf)
            // py::arg("X"))
        .def("correctUkf",&BindUkf::corrUkf);
            // py::arg("Z"));
}