#ifndef btype_h
#define btypes_h
#include <cstdint>
#include <vector>

namespace betagraph{
    struct node_t {
        bool orientation;
        std::int64_t id;
        char* seq;
    };
    struct edge_t{
        std::int64_t id;
        std::int64_t from_id;
        std::int64_t to_id;
    };
    struct path_t{
        char* name;
        std::vector<int64_t> ids; 
    };
    struct graph_t{

    };

}

#endif