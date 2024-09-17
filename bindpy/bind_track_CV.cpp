
#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf_CV
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanFilter<Eigen::MatrixXd, InitUnscentedKalmanFilterCV,FuncConstVel<Eigen::MatrixXd>, FuncMeasSph<Eigen::MatrixXd>, FuncControlMatrix_XvXYvYZvZ<Eigen::MatrixXd>>>
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

void bind_track_CV(pybind11::module &m)
{
    py::class_<BindTrackUkf_CV>(m, "BindTrackUkf_CV")
        .def(py::init<>())
        .def("init", (void(BindTrackUkf_CV::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CV::init)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CV::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)(double)) & BindTrackUkf_CV::step);
}
