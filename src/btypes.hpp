#ifndef btype_h
#define btypes_h
#include <cstdint>
#include <vector>
#include "spp.h"

namespace betagraph{

    struct bholder_t{
        bool orientation;
        uint32_t length;
        vector<dankgraph::edge_t> in_edges;
        vector<dankgraph::edge_t> out_edges;
        string sequence;
    };
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
    struct bgraph_t{
        // spp::sparse_hash_map<int64_t, b_node_t> id_to_node;
        // spp::sparse_hash_map<int64_t, vector<b_edge_t>> node_to_edges;
        // spp::sparse_hash_map<char*, b_path_t> name_to_path;
        spp::sparse_hash_map<id_t, bholder_t> backer;
        id_t min_node_id = 0;
        id_t max_node_id = UINT64_MAX;
    };
    struct bpath_occurrence_t{
        int32_t rank = 0;
        id_t id;
    };
    struct bpathstore_t{
        spp::sparse_hash_map<string, bpath_occurrence_t> paths;

    };

}

#endif