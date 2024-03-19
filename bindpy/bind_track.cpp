#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf
{
private:
    Track <Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSph>> track;
    

public :

    Eigen::MatrixXd step(Eigen::MatrixXd& meas)
    {
        return track.Step(meas);
    }

    Eigen::MatrixXd step(double dt)
    {
        return track.Step(dt);
    }

};

void bind_track(pybind11::module &m)
{
    py::class_<BindTrackUkf>(m, "BindTrack")
        .def(py::init<>())
        .def("step",(Eigen::MatrixXd(BindTrackUkf::*)(Eigen::MatrixXd&)) &BindTrackUkf::step)
        .def("step",(Eigen::MatrixXd(BindTrackUkf::*)(double)) &BindTrackUkf::step);
        
}
