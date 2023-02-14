#include "inject.hpp"

namespace odgi {

namespace algorithms {

using namespace handlegraph;

void inject_ranges(MutablePathDeletableHandleGraph& graph,
                   const ska::flat_hash_map<path_handle_t, std::vector<std::pair<interval_t, std::string>>>& path_intervals,
                   const std::vector<std::string>& ordered_intervals, const bool show_progress) {

    // we collect cut points based on where our intervals start and end
    ska::flat_hash_map<handle_t, std::vector<size_t>> cut_points;

    // parallel over paths requires collecting path handles in a vector
    std::vector<path_handle_t> paths;
    graph.for_each_path_handle([&](const path_handle_t& path) {
        paths.push_back(path);
    });

    std::atomic<uint64_t> num_intervals;
    num_intervals.store(0);

    std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress;
    if (show_progress) {
        // Count intervals to process
#pragma omp parallel for
        for (auto& path : paths) {
            auto x = path_intervals.find(path);
            if (x != path_intervals.end()) {
                num_intervals += x->second.size();
            }
        }

        progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                num_intervals.load(), "[odgi::inject] collecting cut points");
    }

#pragma omp parallel for
    for (auto& path : paths) {
        if (path_intervals.find(path) != path_intervals.end()) {
            auto& intervals = path_intervals.find(path)->second;
            auto ival = intervals.begin();
            std::set<uint64_t> interval_ends;
            uint64_t pos = 0;
            handle_t last_h;
            graph.for_each_step_in_path(
                path,
                [&](const step_handle_t& step) {
                    // remove intervals that ended before this node
                    while (interval_ends.size()
                           && *interval_ends.begin() <= pos) {
                        auto last_length = graph.get_length(last_h);
                        auto last_offset =
                            graph.get_is_reverse(last_h) ?
                            (pos - *interval_ends.begin())
                            : last_length - (pos - *interval_ends.begin());
                        if (last_offset > 0 && last_offset < last_length) {
                            // store cut points on the forward strand to avoid dups
                            auto last_h_fwd = graph.get_is_reverse(last_h) ?
                                graph.flip(last_h) : last_h;
#pragma omp critical (cut_points)
                            cut_points[last_h_fwd].push_back(last_offset);
                        }
                        interval_ends.erase(interval_ends.begin());
                    }
                    auto h = graph.get_handle_of_step(step);
                    auto len = graph.get_length(h);
                    // add intervals that start on this node
                    while (ival != intervals.end()
                           && ival->first.first >= pos
                           && ival->first.first < pos + len) {
                        interval_ends.insert(ival->first.second);
                        // mark a cut point
                        auto offset = graph.get_is_reverse(h) ?
                            len - (ival->first.first - pos)
                            : ival->first.first - pos;
                        if (offset > 0 && offset < len) {
                            // store cut points on the forward strand to avoid dups
                            auto h_fwd = graph.get_is_reverse(h) ?
                                graph.flip(h) : h;
#pragma omp critical (cut_points)
                            cut_points[h_fwd].push_back(offset);
                        }
                        ++ival;
                        if (show_progress) {
                            progress->increment(1);
                        }
                    }
                    pos += len;
                    last_h = h;
                });

            // Manage last offsets
            while (interval_ends.size() && *interval_ends.begin() <= pos) {
                auto last_length = graph.get_length(last_h);
                auto last_offset =
                        graph.get_is_reverse(last_h) ?
                        (pos - *interval_ends.begin())
                                                     : last_length - (pos - *interval_ends.begin());
                if (last_offset > 0 && last_offset < last_length) {
                    // store cut points on the forward strand to avoid dups
                    auto last_h_fwd = graph.get_is_reverse(last_h) ?
                                      graph.flip(last_h) : last_h;
#pragma omp critical (cut_points)
                    cut_points[last_h_fwd].push_back(last_offset);
                }
                interval_ends.erase(interval_ends.begin());
            }
        }
    }

    if (show_progress) {
        progress->finish();
    }

    for (auto& c : cut_points) {
        auto& v = c.second;
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    /*
    for (auto& c : cut_points) {
        std::cerr << "cutting at " << graph.get_id(c.first) << " -> ";
        for (auto& p : c.second) {
            std::cerr << p << " ";
        }
        std::cerr << std::endl;
    }
    */

    // then we cut the nodes in the graph at the interval starts and ends
    chop_at(graph, cut_points);

    if (show_progress) {
        progress = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                num_intervals.load(), "[odgi::inject] injecting intervals");
    }

    ska::flat_hash_map<std::string, path_handle_t> injected_paths;
    for (auto& n : ordered_intervals) {
        auto f = injected_paths.find(n);
        if (f != injected_paths.end()) {
            std::cerr << "[odgi::algorithms::inject_ranges] "
                      << "duplicate annotation path, unable to inject "
                      << n << std::endl;
            exit(1);
        } else {
            injected_paths[n] = graph.create_path_handle(n);
        }
    }

    // then we iterate back through the sorted path intervals and add paths at the appropriate points
#pragma omp parallel for
    for (auto& path : paths) {
        if (path_intervals.find(path) != path_intervals.end()) {
            auto& intervals = path_intervals.find(path)->second;
            auto ival = intervals.begin();
            // IMPORTANT: in the key of the map, the end position (uint64_t) has
            // to come first to keep the open intervals sorted by end coordinates.
            // The annotation name (std::string) is necessary to avoid losing
            // annotations having the same end coordinates.
            std::map<std::pair<uint64_t, std::string>, step_handle_t> open_intervals_by_end;
            uint64_t pos = 0;
            handle_t last_h;
            graph.for_each_step_in_path(
                path,
                [&](const step_handle_t& step) {
                    // remove intervals that ended before this node
                    while (open_intervals_by_end.size()
                           && open_intervals_by_end.begin()->first.first <= pos) {
                        // get reference to name
                        auto& name = open_intervals_by_end.begin()->first.second;
                        if (open_intervals_by_end.begin()->first.first != pos) {
                            std::cerr << "[odgi::algorithms::inject_ranges] "
                                      << "injection end point for interval " << name
                                      << " is not at a node boundary: "
                                      << "off by " << open_intervals_by_end.begin()->first.first
                                      << " vs " << pos
                                      << " on a node "
                                      << graph.get_length(graph.get_handle_of_step(step))
                                      << "bp long"
                                      << std::endl;
                            exit(1);
                        }
                        // add the path
                        path_handle_t p;
#pragma omp critical (get_path)
                        {
                            auto f = injected_paths.find(name);
                            assert(f != injected_paths.end());
                            p = f->second;
                        }
                        auto& c = open_intervals_by_end.begin()->second;
                        auto end = step;
                        do {
                            graph.append_step(p, graph.get_handle_of_step(c));
                            c = graph.get_next_step(c);
                        } while (c != end);
                        // clean up
                        open_intervals_by_end.erase(open_intervals_by_end.begin());
                    }
                    auto h = graph.get_handle_of_step(step);
                    auto len = graph.get_length(h);
                    // add intervals that start on this node
                    while (ival != intervals.end()
                           && ival->first.first >= pos
                           && ival->first.first < pos + len) {
                        // the intervals must start at the node start
                        open_intervals_by_end[std::make_pair(ival->first.second, ival->second)] = step;
                        ++ival;
                        if (show_progress) {
                            progress->increment(1);
                        }
                    }
                    pos += len;
                    last_h = h;
                });

            // Add last paths
            while (open_intervals_by_end.size()
                   && open_intervals_by_end.begin()->first.first <= pos) {
                // get reference to name
                auto& name = open_intervals_by_end.begin()->first.second;
                // add the path
                path_handle_t p;
#pragma omp critical (get_path)
                {
                    auto f = injected_paths.find(name);
                    assert(f != injected_paths.end());
                    p = f->second;
                }
                auto& c = open_intervals_by_end.begin()->second;
                auto end = graph.path_end(path);
                do {
                    graph.append_step(p, graph.get_handle_of_step(c));
                    c = graph.get_next_step(c);
                } while (c != end);
                // clean up
                open_intervals_by_end.erase(open_intervals_by_end.begin());
            }
        }
    }

    if (show_progress) {
        progress->finish();
    }
}

void chop_at(MutablePathDeletableHandleGraph &graph,
             const ska::flat_hash_map<handle_t, std::vector<size_t>>& cut_points) {

    std::vector<std::tuple<uint64_t, uint64_t, handle_t>> originalRank_inChoppedNodeRank_handle;
    std::vector<std::pair<uint64_t, handle_t>> originalRank_handleToChop;

    uint64_t rank = 0;
    graph.for_each_handle([&](const handle_t &handle) {
        if (cut_points.find(handle) != cut_points.end()) {
            originalRank_handleToChop.push_back(std::make_pair(rank, handle));
        } else {
            originalRank_inChoppedNodeRank_handle.push_back(std::make_tuple(rank, 0, handle));
        }
        rank++;
    });

    //std::cerr << "[odgi::inject::chop_at] " << originalRank_handleToChop.size() << " node(s) to chop." << std::endl;

    for (auto rank_handle : originalRank_handleToChop) {
        // get divide points
        auto& handle = rank_handle.second;
        uint64_t length = graph.get_length(handle);
        auto f = cut_points.find(handle);
        assert(f != cut_points.end());
        const std::vector<size_t>& cut_offsets = f->second;
        rank = 0;
        for (auto chopped_handle : graph.divide_handle(rank_handle.second, cut_offsets)) {
            originalRank_inChoppedNodeRank_handle.push_back(
                std::make_tuple(rank_handle.first, rank, chopped_handle));
            rank++;
        }
    }

    std::sort(originalRank_inChoppedNodeRank_handle.begin(), originalRank_inChoppedNodeRank_handle.end());

    std::vector<handle_t> new_handles;
    for (auto x_y_z : originalRank_inChoppedNodeRank_handle) {
        new_handles.push_back(std::get<2>(x_y_z));
    }

    graph.apply_ordering(new_handles, true);
}


}

}
