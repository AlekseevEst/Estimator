#include "bind_ukf.h"
namespace py = pybind11;

class BindUkf
{
private:

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn,FuncMeasSphCVCT, FuncControlMatrix_XvXaXYvYaYZvZaZ> ukf;
 
public:
    
        BindUkf(Eigen::MatrixXd state,
                Eigen::MatrixXd processNoise,
                Eigen::MatrixXd measureNoise,
                ParamSigmaPoints paramSigmaPoints):              
                          ukf(state, processNoise, measureNoise, paramSigmaPoints){}
                              

    Eigen::MatrixXd predUkf(double dt)
    {
        return ukf.predict(dt);
    }
    Eigen::MatrixXd corrUkf(const Eigen::MatrixXd &Z)
    {
        return ukf.correct(Z);
    }
};

void bind_ukf(pybind11::module &m)
{
    py::class_<BindUkf>(m, "BindUkf")
        .def(py::init<const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, ParamSigmaPoints>())
        .def("predictUkf",&BindUkf::predUkf)
            // py::arg("X"))
        .def("correctUkf",&BindUkf::corrUkf);
            // py::arg("Z"));
}