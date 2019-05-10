#ifndef btype_h
#define btypes_h
#include <cstdint>
#include <vector>
#include <string>
#include "spp.h"
#include <handlegraph/types.hpp>

namespace betagraph{

    struct bholder_t{
        bool orientation;
        uint32_t length;
        std::vector<odgi::edge_t> in_edges;
        std::vector<odgi::edge_t> out_edges;
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
        // spp::sparse_hash_map<int64_t, b_node_t> nid_to_node;
        // spp::sparse_hash_map<int64_t, vector<b_edge_t>> node_to_edges;
        // spp::sparse_hash_map<char*, b_path_t> name_to_path;
        spp::sparse_hash_map<odgi::nid_t, bholder_t> backer;
        handlegraph::nid_t min_node_id = 0;
        handlegraph::nid_t max_node_id = std::numeric_limits<handlegraph::nid_t>::max();
    };
    struct bpath_step_t{
        int32_t rank = 0;
        handlegraph::nid_t id;
    };
    struct bpath_key_t{
        handlegraph::nid_t id;
        std::string name;
        bool operator==(const bpath_key_t &o) const 
        { return name == o.name && name == o.name; }
    };

    struct bpath_val_t{
        handlegraph::nid_t id;
        std::vector<bpath_step_t> steps;
    };



    struct bpathstore_t{
        spp::sparse_hash_map<std::string, bpath_val_t> paths;
        spp::sparse_hash_map<handlegraph::nid_t, std::string> nid_to_name;
        bool has_path(const std::string& s) const{
            return (paths.count(s) != 0);
        };
        handlegraph::nid_t get_id(std::string s) const{
            return paths.at(s).id;
        };
        std::string get_name(odgi::path_handle_t p) const{
            return nid_to_name.at(as_integer(p));
        };
        size_t get_path_step_count(odgi::path_handle_t p) const{
            return paths.at(nid_to_name.at(as_integer(p))).steps.size();
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
