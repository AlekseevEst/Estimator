#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf_CA
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstAcceleration, FuncMeasSphCA>> track;
    
public:
    
    BindTrackUkf_CA(Eigen::MatrixXd state,
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

void bind_track_CA(pybind11::module &m)
{
    py::class_<BindTrackUkf_CA>(m, "BindTrackUkf_CA")
        .def(py::init<const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, Points>())
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CA::*)(double, const Eigen::MatrixXd &)) & BindTrackUkf_CA::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CA::*)(double)) & BindTrackUkf_CA::step);

}
