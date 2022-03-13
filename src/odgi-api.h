// odgi-api.h ODGI C-API for FFIs
//
// Copyright (c) 2022 Erik Garrison and Pjotr Prins
//
// ODGI is published under the MIT license

#pragma once

#include <cstddef>
// #include <cstdint>

extern "C" {

  // These functions are exposed as a C API
  const char *odgi_version();
  size_t odgi_graph_nodes();

}
