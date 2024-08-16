// #include "bind_imm.h"
// namespace py = pybind11;

// class BindImm
// {
// private:

//     // ConteinerCVCACT<Eigen::MatrixXd> conteiner;
//     Imm <Eigen::MatrixXd, ConteinerCVCACT> imm;
  

// public:
    
//     BindImm(    const  Eigen::MatrixXd& modeProbability,
//                 const  Eigen::MatrixXd& transmitProbability,
//                 const  Eigen::MatrixXd& state,
//                 const  Eigen::MatrixXd& covariance,
//                 const  Eigen::MatrixXd& processNoise,
//                 const  Eigen::MatrixXd& measureNoise,
//                 ParamSigmaPoints paramSigmaPoints):imm(modeProbability, transmitProbability, state, covariance, processNoise, measureNoise, paramSigmaPoints)
//                 {

//                 }

// // conteiner(state, covariance, processNoise, measureNoise, paramSigmaPoints)


//     Eigen::MatrixXd step(const Eigen::MatrixXd& Z, double dt)
//     {
//         return imm.step(Z, dt);
//     }

// };

// void bind_imm(pybind11::module &m)
// {
//     py::class_<BindImm>(m, "BindImm")
//         .def(py::init<const  Eigen::MatrixXd&, const  Eigen::MatrixXd&, const  Eigen::MatrixXd&, const  Eigen::MatrixXd&, const  Eigen::MatrixXd&, const  Eigen::MatrixXd&, ParamSigmaPoints>())
//         .def("step", (Eigen::MatrixXd(BindImm::*)(const  Eigen::MatrixXd&, double)) & BindImm::step);


// }
