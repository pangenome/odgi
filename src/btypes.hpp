#ifndef btype_h
#define btypes_h
#include <cstdint>
#include <vector>
#include <string>
#include "spp.h"

namespace betagraph{

    struct bholder_t{
        bool orientation;
        uint32_t length;
        std::vector<dankgraph::edge_t> in_edges;
        std::vector<dankgraph::edge_t> out_edges;
        std::string sequence;
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
        spp::sparse_hash_map<dankgraph::id_t, bholder_t> backer;
        id_t min_node_id = 0;
        id_t max_node_id = std::numeric_limits<id_t>::max();
    };
    struct bpath_occurrence_t{
        int32_t rank = 0;
        id_t id;
    };
    struct bpath_key_t{
        dankgraph::id_t id;
        std::string name;
        bool operator==(const bpath_key_t &o) const 
        { return name == o.name && name == o.name; }
    };

    struct bpath_val_t{
        dankgraph::id_t id;
        std::vector<bpath_occurrence_t> occurrences;
    };



    struct bpathstore_t{
        spp::sparse_hash_map<std::string, bpath_val_t> paths;
        spp::sparse_hash_map<dankgraph::id_t, std::string> id_to_name;
        bool has_path(const std::string& s) const{
            return (paths.count(s) != 0);
        };
        dankgraph::id_t get_id(std::string s) const{
            return paths.at(s).id;
        };
        std::string get_name(dankgraph::path_handle_t p) const{
            return id_to_name.at(as_integer(p));
        };
        size_t get_path_occurrence_count(dankgraph::path_handle_t p) const{
            return paths.at(id_to_name.at(as_integer(p))).occurrences.size();
        };

    };

}

    namespace std
{
    // inject specialization of std::hash for Person into namespace std
    // ----------------------------------------------------------------
    template<> 
    struct hash<betagraph::bpath_key_t>
    {
        std::size_t operator()(betagraph::bpath_key_t const &p) const
        {
            std::size_t seed = 0;
            spp::hash_combine(seed, p.name);
            spp::hash_combine(seed, p.name);
            return seed;
        }
    };
};

#endif
