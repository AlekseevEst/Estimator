#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf_CA
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstAcceleration, FuncMeasSph>> track;
    
public:
    
    BindTrackUkf_CA(Eigen::MatrixXd state,
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

void bind_track_CA(pybind11::module &m)
{
    py::class_<BindTrackUkf_CA>(m, "BindTrackUkf_CA")
        .def(py::init<const Eigen::MatrixXd&, double, const Eigen::MatrixXd&, const Eigen::MatrixXd&, Points>())
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CA::*)(const Eigen::MatrixXd &)) & BindTrackUkf_CA::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CA::*)()) & BindTrackUkf_CA::step);

}
