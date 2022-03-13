// odgi-api exports the C-API for FFIs
//
// Copyright (c) 2022 Erik Garrison and Pjotr Prins
//
// ODGI is published under the MIT license

#include "odgi-api.h"
#include "odgi.hpp"
#include "version.hpp"

namespace odgi {

size_t odgi_graph_nodes() {
  return 100; // odgi::graph_t::get_node_count();
}

const char *odgi_version() {
  return Version::get_version().c_str();
}

}
