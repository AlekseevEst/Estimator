
#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkfImm_CTxy
{
private:
    TrackImm<Eigen::MatrixXd, Imm<Eigen::MatrixXd, ConteinerCVCACTxy>,
          InitUkfImmCVCTxyCAMeasureModelSph<Eigen::MatrixXd, ConteinerCVCACTxy>>
        track;

public:

    BindTrackUkfImm_CTxy(const Detection<Eigen::MatrixXd>& detection):track(detection){}

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

void bind_track_Imm_CTxy(pybind11::module &m)
{
    py::class_<BindTrackUkfImm_CTxy>(m, "BindTrackUkfImm_CTxy")
        .def(py::init<const Detection<Eigen::MatrixXd>&>())
        .def("step", (Eigen::MatrixXd(BindTrackUkfImm_CTxy::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkfImm_CTxy::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkfImm_CTxy::*)(double)) & BindTrackUkfImm_CTxy::step)
        .def("get_m_i", (Eigen::MatrixXd(BindTrackUkfImm_CTxy::*)()) & BindTrackUkfImm_CTxy::get_m_i);
}
