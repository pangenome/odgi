#include "gfa_to_handle.hpp"

namespace odgi {

    void gfa_to_handle(const string& gfa_filename, handlegraph::MutablePathMutableHandleGraph* graph,
                       bool show_progress) {
        
        char* filename = (char*) gfa_filename.c_str();
        //std::cerr << "filename is " << filename << std::endl;
        gfak::GFAKluge gg;
        //double version = gg.detect_version_from_file(filename);
        //std::cerr << version << " be version" << std::endl;
        //assert(version == 1.0);
        /*
         uint64_t num_nodes = 0;
         gg.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
         ++num_nodes;
         });
         graph_t graph(num_nodes+1); // include delimiter
         */
        uint64_t i = 0;
        uint64_t min_id = std::numeric_limits<uint64_t>::max();
        uint64_t max_id = std::numeric_limits<uint64_t>::min();
        gg.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
                uint64_t id = stol(s.name);
                min_id = std::min(min_id, id);
                max_id = std::max(max_id, id);
            });
        graph->set_id_increment(min_id);
        // set the min id as an offset
        gg.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
            uint64_t id = stol(s.name);
            graph->create_handle(s.sequence, id);
            if (show_progress) {
                if (i % 1000 == 0) std::cerr << "node " << i << "\r";
                ++i;
            }
        });
        if (show_progress) {
            i = 0; std::cerr << std::endl;
        }
        gg.for_each_edge_line_in_file(filename, [&](gfak::edge_elem e) {
            if (e.source_name.empty()) return;
            handlegraph::handle_t a = graph->get_handle(stol(e.source_name), !e.source_orientation_forward);
            handlegraph::handle_t b = graph->get_handle(stol(e.sink_name), !e.sink_orientation_forward);
            graph->create_edge(a, b);
            if (show_progress) {
                if (i % 1000 == 0) std::cerr << "edge " << i << "\r";
                ++i;
            }
        });
        if (show_progress) {
            i = 0; std::cerr << std::endl;
        }
        gg.for_each_path_element_in_file(filename, [&](const std::string& path_name_raw, const std::string& node_id, bool is_rev, const std::string& cigar) {
            handlegraph::path_handle_t path;
            std::string path_name = path_name_raw;
            path_name.erase(std::remove_if(path_name.begin(), path_name.end(), [](char c) { return std::isspace(c); }), path_name.end());
            if (!graph->has_path(path_name)) {
                if (show_progress) {
                    std::cerr << "path " << ++i << "\r";
                }
                path = graph->create_path_handle(path_name);
            } else {
                path = graph->get_path_handle(path_name);
            }
            handlegraph::handle_t occ = graph->get_handle(stol(node_id), is_rev);
            graph->append_step(path, occ);
            // ignores overlaps
        });
    }
}
