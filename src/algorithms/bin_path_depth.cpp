#include "bin_path_depth.hpp"

// #define  debug_bin_path_depth

namespace odgi {
    namespace algorithms {

        void bin_path_depth(const PathHandleGraph &graph,
                           const bool progress,
                           const uint64_t min_paths,
                           const uint64_t min_depth) {
            std::unique_ptr<progress_meter::ProgressMeter> progress_meter;
            if (progress) {
                progress_meter = std::make_unique<progress_meter::ProgressMeter>(
                        graph.get_node_count(), "[odgi::bin] bin_path_depth:");
            }
            const uint64_t path_count = graph.get_path_count();
            // write the header
            std::cout << "bin_id";
            graph.for_each_path_handle([&](const path_handle_t &path) {
                std::cout << "\t" << graph.get_path_name(path);
            });
            std::cout << std::endl;
            // the graph must be compacted for this to work
            graph.for_each_handle([&](const handle_t &h) {
                vector<uint64_t> p_i = vector<uint64_t>(path_count);
                graph.for_each_step_on_handle(h, [&](const step_handle_t &occ) {
                    const path_handle_t p_h = graph.get_path_handle_of_step(occ);
                    p_i[as_integer(p_h) - 1] += 1;
                });

                if (progress) {
                    progress_meter->increment(1);
                }
                std::string row = std::to_string(number_bool_packing::unpack_number(h) + 1);
                uint64_t paths_with_cov = 0;
                uint64_t idx = 0;
                for (const auto& path_cov : p_i) {
                    row += "\t";
                    if (path_cov < min_depth) {
                        row += std::to_string(0);
                    } else {
                        paths_with_cov++;
                        row += std::to_string(path_cov);
                    }
                    idx++;
                }
                row += "\n";
                if ((!(paths_with_cov == path_count)) && (paths_with_cov >= min_paths)) {
                    std::cout << row;
                }
            });
            if (progress) {
                progress_meter->finish();
            }
        }
    }
}
