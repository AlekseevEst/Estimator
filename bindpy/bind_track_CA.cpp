#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf_CA
{
using TypeFilterCA = UnscentedKalmanFilter<Eigen::MatrixXd, InitUnscentedKalmanFilterCA,FuncConstAcceleration<Eigen::MatrixXd>, FuncMeasSph<Eigen::MatrixXd>, FuncControlMatrix_XvXaXYvYaYZvZaZ<Eigen::MatrixXd>>;
private:
    Track<Eigen::MatrixXd, TypeFilterCA>
        track;

public:

    void init(const Detection<Eigen::MatrixXd>& detection)
    {
        track.Initialization(detection);
    }

    Eigen::MatrixXd step(const Detection<Eigen::MatrixXd>& detection)
    {
        return track.step(detection);
    }

    Eigen::MatrixXd step(double dt)
    {
        return track.step(dt);
    }
};

void bind_track_CA(pybind11::module &m)
{

    py::class_<BindTrackUkf_CA>(m, "BindTrackUkf_CA")
        .def(py::init<>())
        .def("init", (void(BindTrackUkf_CA::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CA::init)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CA::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CA::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CA::*)(double)) & BindTrackUkf_CA::step);

    py::class_<Detection<Eigen::MatrixXd>>(m, "Detection")
        .def(py::init<>())
        .def_readwrite("time", &Detection<Eigen::MatrixXd>::time)
        .def_readwrite("measurement", &Detection<Eigen::MatrixXd>::measurement)
        .def_readwrite("measurementNoise", &Detection<Eigen::MatrixXd>::measurementNoise);

}
