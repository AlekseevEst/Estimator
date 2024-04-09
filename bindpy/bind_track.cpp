#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSph>> track;

public:
    BindTrackUkf(Eigen::MatrixXd state,
                 double t,
                 Eigen::MatrixXd processNoise,
                 Eigen::MatrixXd measureNoise,
                 double k) : track(state, t, processNoise, measureNoise, k) {}

    Eigen::MatrixXd step(const Eigen::MatrixXd &meas)
    {
        return track.Step(meas);
    }

    Eigen::MatrixXd step()
    {
        return track.Step();
    }
};

void bind_track(pybind11::module &m)
{
    py::class_<BindTrackUkf>(m, "BindTrackUkf")
        .def(py::init<const Eigen::MatrixXd&, double, const Eigen::MatrixXd&, const Eigen::MatrixXd&, double>())
        .def("step", (Eigen::MatrixXd(BindTrackUkf::*)(const Eigen::MatrixXd &)) & BindTrackUkf::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf::*)()) & BindTrackUkf::step);
}
