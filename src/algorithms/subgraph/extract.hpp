#ifndef ODGI_EXTRACT_H
#define ODGI_EXTRACT_H

#include <regex>
#include <unordered_map>
#include <stack>
#include <map>
#include "odgi.hpp"
#include "handlegraph/path_position_handle_graph.hpp"
#include "weakly_connected_components.hpp"
#include "progress.hpp"
#include "utils.hpp"
#include "position.hpp"
#include "src/algorithms/subgraph/region.hpp"

namespace odgi {
    namespace algorithms {
        /* Several functions were inspired by https://github.com/vgteam/vg */

        void add_full_paths_to_component(const graph_t &source, graph_t &component, const uint64_t num_threads);

        std::string make_path_name(const string &path_name, size_t offset, size_t end_offset);

        path_handle_t create_subpath(graph_t &subgraph, const string &subpath_name, const bool is_circular);

        /// add subpaths to the subgraph, providing a concatenation of subpaths that are discontiguous over the subgraph
        void add_subpaths_to_subgraph(const graph_t &source, const std::vector<path_handle_t> source_paths,
                                      graph_t &subgraph, uint64_t num_threads,
                                      const std::string &progress_message = "");

        void extract_path_range(const graph_t &source, path_handle_t path_handle, int64_t start, int64_t end,
                                graph_t &subgraph);

        void for_handle_in_path_range(const graph_t &source, path_handle_t path_handle, int64_t start, int64_t end,
                                      const std::function<void(const handle_t&)>& lambda);

        void add_connecting_edges_to_subgraph(const graph_t &source, graph_t &subgraph,
                                              const std::string &progress_message = "");

        void expand_subgraph_by_length(const graph_t &source, graph_t &subgraph, const uint64_t &length, bool forward_only);

        void expand_subgraph_by_steps(const graph_t &source, graph_t &subgraph, const uint64_t &steps, bool forward_only);

        void extract_id_range(const graph_t &source, const nid_t &id1, const nid_t &id2, graph_t &subgraph,
                              const std::string &progress_message = "");

        void embed_lace_paths(graph_t &source, graph_t &subgraph,
                              const std::vector<path_handle_t>& lace_paths);

        void windows_in_out(
                const PathHandleGraph& graph,
                const std::vector<path_handle_t>& paths,
                const std::function<bool(handle_t)>& in_bounds,
                const uint64_t& length,
                const std::function<void(const std::vector<path_range_t>&)>& output,
                const uint64_t& num_threads);

        bool check_and_get_windows_in_out_parameter(
                const std::string& parameter,
                uint64_t &windows_len, uint64_t &windows_min, uint64_t &windows_max, bool &only_tips);

    }
}
#endif //ODGI_EXTRACT_H
