#pragma once
#include <iostream>
#include <chrono>
#include <vector>
#include <omp.h>
#include <random>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include <sstream>
#include <iomanip>

#include "odgi.hpp"
#include "XoshiroCpp.hpp"
#include "dirty_zipfian_int_distribution.h"

namespace cuda {


struct __align__(8) node_t {
    float coords[4];
    int32_t seq_length;
};
struct node_data_t {
    uint32_t node_count;
    node_t *nodes;
};


struct __align__(8) path_element_t {
    uint32_t pidx;
    uint32_t node_id;
    int64_t pos;    // if position negative: reverse orientation
};

struct path_t {
    uint32_t step_count;
    uint64_t first_step_in_path;  // precomputed position in path
    path_element_t *elements;
};

struct path_data_t {
    uint32_t path_count;
    uint64_t total_path_steps;
    path_t *paths;
    path_element_t *element_array;
};


// #define SM_COUNT 84
// #define SM_COUNT 108
#define BLOCK_SIZE 1024
#define WARP_SIZE 32
struct curandStateXORWOWCoalesced_t {
    unsigned int d[BLOCK_SIZE];
    unsigned int w0[BLOCK_SIZE];
    unsigned int w1[BLOCK_SIZE];
    unsigned int w2[BLOCK_SIZE];
    unsigned int w3[BLOCK_SIZE];
    unsigned int w4[BLOCK_SIZE];
};
typedef struct curandStateXORWOWCoalesced_t curandStateCoalesced_t;


struct layout_config_t {
    uint64_t iter_max;
    uint64_t min_term_updates;
    double eta_max;
    double eps;
    int32_t iter_with_max_learning_rate;
    uint32_t first_cooling_iteration;
    double theta;
    uint32_t space;
    uint32_t space_max;
    uint32_t space_quantization_step;
    int nthreads;
};


void gpu_layout(layout_config_t config, const odgi::graph_t &graph, std::vector<std::atomic<double>> &X, std::vector<std::atomic<double>> &Y);

}