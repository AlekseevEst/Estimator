#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf_CT
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn, FuncMeasSph>> track;
    
public:
    
    BindTrackUkf_CT(Eigen::MatrixXd state,
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

void bind_track_CT(pybind11::module &m)
{
    py::class_<BindTrackUkf_CT>(m, "BindTrackUkf_CT")
        .def(py::init<const Eigen::MatrixXd&, double, const Eigen::MatrixXd&, const Eigen::MatrixXd&, Points>())
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CT::*)(const Eigen::MatrixXd &)) & BindTrackUkf_CT::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CT::*)()) & BindTrackUkf_CT::step);
    py::class_<Points>(m,"Points")
        .def(py::init<>())
        .def_readwrite("alpha", &Points::alpha) 
        .def_readwrite("beta", &Points::beta)
        .def_readwrite("kappa",&Points::kappa);   
}
