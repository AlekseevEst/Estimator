#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include "track.h"
#include "initFilters.h"
#include "imm.h"
#include "models.h"
#include "structs.h"

void bind_track_CV(pybind11::module &m);
void bind_track_CTxy(pybind11::module &m);
// void bind_track_CTxz(pybind11::module &m);
void bind_track_CA(pybind11::module &m);
void bind_track_Imm_ConteinerCVCTxyCA(pybind11::module &m);
// void bind_track_Imm_ConteinerCVCACTxz(pybind11::module &m);
