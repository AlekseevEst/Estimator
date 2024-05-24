#include "bind_track.h"
namespace py = pybind11;

// Points Points_CV;
class BindTrackUkf_CV
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSphCVCT>> track;

public:
    BindTrackUkf_CV(Eigen::MatrixXd state,
                    Eigen::MatrixXd processNoise,
                    Eigen::MatrixXd measureNoise,
                    Points points) : track(state, processNoise, measureNoise, points) {}

    Eigen::MatrixXd step(double dt, const Eigen::MatrixXd &meas)
    {
        return track.Step(dt, meas);
    }

    Eigen::MatrixXd step(double dt)
    {
        return track.Step(dt);
    }
};

void bind_track_CV(pybind11::module &m)
{
    py::class_<BindTrackUkf_CV>(m, "BindTrackUkf_CV")
        .def(py::init<const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, Points>())
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)(double, const Eigen::MatrixXd &)) & BindTrackUkf_CV::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)(double)) & BindTrackUkf_CV::step);

}

