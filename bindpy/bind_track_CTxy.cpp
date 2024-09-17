
#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf_CTxy
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanFilter<Eigen::MatrixXd, InitUnscentedKalmanFilterCT,FuncConstTurnXY<Eigen::MatrixXd>, FuncMeasSph<Eigen::MatrixXd>, FuncControlMatrix_XvXYvYZvZW<Eigen::MatrixXd>>>
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

void bind_track_CTxy(pybind11::module &m)
{
    py::class_<BindTrackUkf_CTxy>(m, "BindTrackUkf_CTxy")
        .def(py::init<>())
        .def("init", (void(BindTrackUkf_CTxy::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CTxy::init)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CTxy::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CTxy::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CTxy::*)(double)) & BindTrackUkf_CTxy::step);
}
