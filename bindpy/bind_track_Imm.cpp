
#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkfImm
{
private:
    TrackImm<Eigen::MatrixXd, Imm<Eigen::MatrixXd, ConteinerCVCACT>,
          InitUkfImmCVCTCAMeasureModelSph<Eigen::MatrixXd, ConteinerCVCACT>>
        track;

public:

    BindTrackUkfImm(const Detection<Eigen::MatrixXd>& detection):track(detection){}

    Eigen::MatrixXd step(const Detection<Eigen::MatrixXd>& detection)
    {
        return track.step(detection);
    }

    Eigen::MatrixXd step(double dt)
    {
        // return track.step(dt);
    }
    
    Eigen::MatrixXd get_m_i()
    {
        return track.estimator->mu_i;
    }
};

void bind_track_Imm(pybind11::module &m)
{
    py::class_<BindTrackUkfImm>(m, "BindTrackUkfImm")
        .def(py::init<const Detection<Eigen::MatrixXd>&>())
        .def("step", (Eigen::MatrixXd(BindTrackUkfImm::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkfImm::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkfImm::*)(double)) & BindTrackUkfImm::step)
        .def("get_m_i", (Eigen::MatrixXd(BindTrackUkfImm::*)()) & BindTrackUkfImm::get_m_i);
}
