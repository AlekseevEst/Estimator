#include "bind_track.h"
namespace py = pybind11;

// Points Points_CV;
class BindTrackUkf_CV
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSph>> track;

public:
    BindTrackUkf_CV(Eigen::MatrixXd state,
                 double t,
                 Eigen::MatrixXd processNoise,
                 Eigen::MatrixXd measureNoise,
                 Points points) : track(state, t, processNoise, measureNoise, points) {}

    Eigen::MatrixXd step(const Eigen::MatrixXd &meas)
    {
        return track.Step(meas);
    }

    Eigen::MatrixXd step()
    {
        return track.Step();
    }
};

void bind_track_CV(pybind11::module &m)
{
    py::class_<BindTrackUkf_CV>(m, "BindTrackUkf_CV")
        .def(py::init<const Eigen::MatrixXd&, double, const Eigen::MatrixXd&, const Eigen::MatrixXd&, Points>())
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)(const Eigen::MatrixXd &)) & BindTrackUkf_CV::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)()) & BindTrackUkf_CV::step);

}

