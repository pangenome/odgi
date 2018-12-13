#ifndef btype_h
#define btypes_h
#include <cstdint>
#include <vector>
#include "spp.h"

namespace betagraph{
    struct b_node_t {
        bool orientation;
        std::int64_t id;
        char* seq;
    };
    struct b_edge_t{
        bool forward;
        std::int64_t id;
        std::int64_t from_id;
        std::int64_t to_id;
    };
    struct b_path_t{
        char* name;
        std::vector<int64_t> ids; 
    };
    struct b_graph_t{
        spp::sparse_hash_map<int64_t, b_node_t> id_to_node;
        spp::sparse_hash_map<int64_t, vector<b_edge_t>> node_to_edges;
        spp::sparse_hash_map<char*, b_path_t> name_to_path;
    };

}

#endif