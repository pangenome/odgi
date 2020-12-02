#include "gfa_to_handle.hpp"

namespace odgi {

void gfa_to_handle(const string& gfa_filename,
                   handlegraph::MutablePathMutableHandleGraph* graph,
                   uint64_t n_threads,
                   bool show_progress) {
        
    char* filename = (char*) gfa_filename.c_str();
    //std::cerr << "filename is " << filename << std::endl;
    gfak::GFAKluge gg;
    uint64_t i = 0;
    uint64_t min_id = std::numeric_limits<uint64_t>::max();
    uint64_t max_id = std::numeric_limits<uint64_t>::min();
    gg.for_each_sequence_line_in_file(
        filename,
        [&](gfak::sequence_elem s) {
            uint64_t id = stol(s.name);
            min_id = std::min(min_id, id);
            max_id = std::max(max_id, id);
        });
    uint64_t id_increment = min_id - 1;
    graph->set_id_increment(id_increment);
    // set the min id as an offset
    gg.for_each_sequence_line_in_file(
        filename,
        [&](gfak::sequence_elem s) {
            uint64_t id = stol(s.name);
            graph->create_handle(s.sequence, id - id_increment);
            if (show_progress) {
                if (i % 1000 == 0) std::cerr << "node " << i << "\r";
                ++i;
            }
        });
    if (show_progress) {
        i = 0; std::cerr << std::endl;
    }
    gg.for_each_edge_line_in_file(
        filename,
        [&](gfak::edge_elem e) {
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

    gfa_path_queue_t path_queue;
    std::atomic<bool> work_todo;
    auto worker =
        [&](uint64_t tid) {
            while (work_todo.load()) {
                gfak::path_elem* p;
                if (path_queue.try_pop(p)) {
                    handlegraph::path_handle_t path = graph->create_path_handle(p->name);
                    uint64_t i = 0;
                    //std::cerr << "adding path " << p->name << std::endl;
                    for (auto& s : p->segment_names) {
                        graph->append_step(path,
                                           graph->get_handle(std::stoi(s),
                                                             p->orientations[i++]));
                    }
                    delete p;
                } else {
                    std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                }
            }
        };
    std::vector<std::thread> workers;
    workers.reserve(n_threads);
    work_todo.store(true);
    for (uint64_t t = 0; t < n_threads; ++t) {
        workers.emplace_back(worker, t);
    }
    
    gg.for_each_path_line_in_file(
        filename,
        [&](const gfak::path_elem& path) {
            gfak::path_elem* p = new gfak::path_elem(path);
            path_queue.push(p);
        });

    while (!path_queue.was_empty()) {
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
    }
    work_todo.store(false);
    for (uint64_t t = 0; t < n_threads; ++t) {
        workers[t].join();
    }

}
}
