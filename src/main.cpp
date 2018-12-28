#include "graph.hpp"
#include "gfakluge.hpp"
//#include "io_helper.hpp"

using namespace dankgraph;

void make_graph(void) {
    graph_t graph;
    handle_t a = graph.create_handle("A");
    handle_t t = graph.create_handle("T");
    handle_t g = graph.create_handle("G");
    handle_t c = graph.create_handle("C");
    graph.create_edge(a, t);
    graph.create_edge(a, c);
    graph.create_edge(t, g);
    graph.create_edge(c, g);
}

int main(int argc, char** argv) {
    //make_graph();
    assert(argc > 0);
    char* filename = argv[1];
    //std::cerr << "filename is " << filename << std::endl;
    gfak::GFAKluge gg;
    //double version = gg.detect_version_from_file(filename);
    //std::cerr << version << " be version" << std::endl;
    //assert(version == 1.0);
    graph_t graph;
    gg.for_each_sequence_line_in_file(filename, [&](gfak::sequence_elem s) {
            uint64_t id = stol(s.name);
            graph.create_handle(s.sequence, id);
        });
    gg.for_each_edge_line_in_file(filename, [&](gfak::edge_elem e) {
            if (e.source_name.empty()) return;
            handle_t a = handle_helper::pack(stol(e.source_name), !e.source_orientation_forward);
            handle_t b = handle_helper::pack(stol(e.sink_name), !e.sink_orientation_forward);
            graph.create_edge(a, b);
        });
    gg.for_each_path_line_in_file(filename, [&](gfak::path_elem p) {
            path_handle_t path = graph.create_path_handle(p.name);
            for (uint64_t i = 0; i < p.segment_names.size(); ++i) {
                handle_t occ = graph.get_handle(stol(p.segment_names[i]),
                                                !p.orientations[i]);
                graph.append_occurrence(path, occ);
                // ignores overlaps
            }
        });
    // here we should measure memory usage etc.
    graph.display();
    return 0;
}
