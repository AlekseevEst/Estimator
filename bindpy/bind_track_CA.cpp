#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf_CA
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>,
          InitUKFStateModelCAMeasureModelSph<Eigen::MatrixXd, FuncConstAcceleration,
                                             FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>>
        track;

public:

    BindTrackUkf_CA(const Detection<Eigen::MatrixXd>& detection):track(detection){}

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
        .def(py::init<const Detection<Eigen::MatrixXd>&>())
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CA::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CA::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CA::*)(double)) & BindTrackUkf_CA::step);

    py::class_<Detection<Eigen::MatrixXd>>(m, "Detection")
        .def(py::init<>())
        .def_readwrite("point", &Detection<Eigen::MatrixXd>::point)
        .def_readwrite("timePoint", &Detection<Eigen::MatrixXd>::timePoint);

}
