#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include "track.h"
#include "models.h"
#include "ukf.h"

void bind_track(pybind11::module &m);