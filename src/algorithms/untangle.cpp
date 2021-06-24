#include "untangle.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

std::vector<step_handle_t> untangle_cuts(
    const PathHandleGraph& graph,
    const step_handle_t& start,
    const step_handle_t& end,
    const ska::flat_hash_map<step_handle_t, uint64_t>& step_pos) {
    //const uint64_t& min_length) {
    auto path = graph.get_path_handle_of_step(start);
    auto path_name = graph.get_path_name(path);
    // this assumes that the end is not inclusive
    uint64_t start_pos = step_pos.find(start)->second;
    uint64_t end_pos = step_pos.find(end)->second;
    std::cerr << "untangle_cuts(" << path_name << ", "
              << start_pos << ", "
              << end_pos << ")" << std::endl;
    std::vector<step_handle_t> cut_points;
    cut_points.push_back(start);
    // we go forward until we see a loop, where the other step has position < end_pos and > start_pos
    for (step_handle_t step = start; step != end; step = graph.get_next_step(step)) {
        //  we take the first and shortest loop we find
        // TODO change this, it can be computed based on the node length
        const auto& curr_pos = step_pos.find(step)->second;
        handle_t handle = graph.get_handle_of_step(step);
        bool found_loop = false;
        step_handle_t other;
        graph.for_each_step_on_handle(
            handle,
            [&](const step_handle_t& s) {
                if (step != s // not the step we're on
                    && graph.get_path_handle_of_step(s) == path) {
                    const auto& other_pos = step_pos.find(s)->second;
                    if (other_pos > start_pos
                        && other_pos < end_pos
                        && (!found_loop && other_pos > curr_pos
                            || (found_loop && (step_pos.find(other)->second > other_pos)))) {
                        found_loop = true;
                        other = s;
                    }
                }
            });
        if (found_loop) {
            //  recurse this function into it, taking start as our current handle other side of the loop as our end
            //  to cut_points we add the start position, the result from recursion, and our end position
            /*
            if (step_pos.find(other)->second < step_pos.find(step)->second) {
                std::cerr << "impossible" << std::endl;
                abort();
            }
            */
            auto internal_cuts = untangle_cuts(graph, step, other, step_pos);
            cut_points.insert(cut_points.end(),
                              internal_cuts.begin(),
                              internal_cuts.end());
            //  we then step forward to the loop end and continue iterating
            step = other;
        }
    }
    // TODO this block is the same as the previous algorithm, but in reverse
    // the differences in how positions are managed are subtle, making it hard to fold the
    // forward and reverse version together
    // now we reverse it
    step_handle_t path_begin = graph.path_begin(path);
    if (end == path_begin || !graph.has_previous_step(end)) {
        return cut_points;
    }
    //std::cerr << "reversing" << std::endl;
    for (step_handle_t step = end;
         step_pos.find(step)->second > start_pos;
         step = graph.get_previous_step(step)) {
        //  we take the first and shortest loop we find
        // TODO change this, it can be computed based on the node length
        const auto& curr_pos = step_pos.find(step)->second;
        handle_t handle = graph.get_handle_of_step(step);
        //std::cerr << "rev on step " << graph.get_id(handle) << " " << curr_pos << std::endl;
        bool found_loop = false;
        step_handle_t other;
        graph.for_each_step_on_handle(
            handle,
            [&](const step_handle_t& s) {
                if (step != s // not the step we're on
                    && graph.get_path_handle_of_step(s) == path) {
                    const auto& other_pos = step_pos.find(s)->second;
                    if (other_pos > start_pos
                        && other_pos < end_pos
                        && other_pos < curr_pos
                        && (!found_loop
                            || (found_loop && (step_pos.find(other)->second < other_pos)))) {
                        found_loop = true;
                        other = s;
                    }
                }
            });
        if (found_loop) {
            //  recurse this function into it, taking start as our current handle other side of the loop as our end
            //  to cut_points we add the start position, the result from recursion, and our end position
            /*
            if (step_pos.find(other)->second > step_pos.find(step)->second) {
                std::cerr << "impossible" << std::endl;
                abort();
            }
            */
            auto internal_cuts = untangle_cuts(graph, other, step, step_pos);
            cut_points.insert(cut_points.end(),
                              internal_cuts.begin(),
                              internal_cuts.end());
            //  we then step forward to the loop end and continue iterating
            step = other;
        }
    }
    cut_points.push_back(end);
    // and sort
    std::sort(cut_points.begin(),
              cut_points.end(),
              [&](const step_handle_t& a,
                  const step_handle_t& b) {
                  return step_pos.find(a)->second < step_pos.find(b)->second;
              });
    return cut_points;
}

void write_cuts(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const std::vector<step_handle_t>& cuts,
    const ska::flat_hash_map<step_handle_t, uint64_t>& step_pos) {
    auto path_name = graph.get_path_name(path);
    std::cout << "name\tcut" << std::endl;
    for (auto& step : cuts) {
        std::cout << path_name << "\t" << step_pos.find(step)->second << std::endl;
    }
}

std::vector<step_handle_t> merge_cuts(
    const std::vector<step_handle_t>& cuts,
    const uint64_t& dist,
    const ska::flat_hash_map<step_handle_t, uint64_t>& step_pos) {
    std::vector<step_handle_t> merged;
    uint64_t last = 0;
    std::cerr << "dist is " << dist << std::endl;
    for (auto& step : cuts) {
        auto& pos = step_pos.find(step)->second;
        if (pos == 0 || pos > (last + dist)) {
            merged.push_back(step);
            last = pos;
        }
    }
    return merged;
}

void self_dotplot(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const ska::flat_hash_map<step_handle_t, uint64_t>& step_pos) {
    auto path_name = graph.get_path_name(path);
    std::cout << "name\tfrom\tto" << std::endl;
    graph.for_each_step_in_path(
        path,
        [&](const step_handle_t& step) {
            // TODO this is given by the walk, no need for a hash table lookup
            const auto& curr_pos = step_pos.find(step)->second;
            handle_t handle = graph.get_handle_of_step(step);
            graph.for_each_step_on_handle(
                handle,
                [&](const step_handle_t& s) {
                    if (graph.get_path_handle_of_step(s) == path) {
                        const auto& other_pos = step_pos.find(s)->second;
                        std::cout << path_name << "\t"
                                  << curr_pos << "\t"
                                  << other_pos << std::endl;
                        /*
                        if (other_pos > curr_pos
                            && other_pos > start_pos
                            && other_pos < end_pos
                            && other_pos - start_pos > min_length) {
                        }
                        */
                    }
                });
        });
}

ska::flat_hash_map<step_handle_t, uint64_t> make_step_index(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& paths) {
    ska::flat_hash_map<step_handle_t, uint64_t> step_pos;
    for (auto& path : paths) {
        uint64_t pos = 0;
        graph.for_each_step_in_path(
            path,
            [&](const step_handle_t& step) {
                step_pos[step] = pos;
                handle_t handle = graph.get_handle_of_step(step);
                pos += graph.get_length(handle);
            });
        step_pos[graph.path_end(path)] = pos; // record the end position
    }
    return step_pos;
}

void show_steps(
    const PathHandleGraph& graph,
    const ska::flat_hash_map<step_handle_t, uint64_t>& steps) {
    for (auto& pos : steps) {
        auto h = graph.get_handle_of_step(pos.first);
        auto p = graph.get_path_handle_of_step(pos.first);
        auto name = graph.get_path_name(p);
        auto id = graph.get_id(h);
        auto rev = graph.get_is_reverse(h);
        std::cerr << name << " " << id << (rev?"-":"+") << " " << pos.second << std::endl;
    }
}

// pair-BED projection of the graph
// that describe nonlinear query : target relationships
// but can be sorted according to their query or target axes
void untangle(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& queries,
    const std::vector<path_handle_t>& targets,
    const uint64_t& merge_dist,
    const size_t& num_threads) {
    std::vector<path_handle_t> paths;
    paths.insert(paths.end(), queries.begin(), queries.end());
    paths.insert(paths.end(), targets.begin(), targets.end());
    auto step_pos = make_step_index(graph, paths);

    // compute the reference segmentations
    // and map them onto the graph using a static multiset index structure based on two arrays
    // we'll build up a big vector of node -> path segment pairings
    // then we'll sort them and build an index

    std::vector<std::vector<step_handle_t>> ref_cuts;
    for (auto& target : targets) {
        ref_cuts.push_back(
            merge_cuts(
                untangle_cuts(graph,
                              graph.path_begin(target),
                              graph.path_back(target),
                              step_pos),
                merge_dist,
                step_pos));
        auto& cuts = ref_cuts.back();
        
        //std::vector<>
        //write_cuts(graph, query, step_pos);
        // 
    }
    
    //
    // 
    // 
    //show_steps(graph, step_pos);
    //std::cout << "path\tfrom\tto" << std::endl;

    for (auto& query : queries) {
        std::vector<step_handle_t> cuts
            = merge_cuts(
                untangle_cuts(graph,
                              graph.path_begin(query),
                              graph.path_back(query),
                              step_pos),
                merge_dist,
                step_pos);
        //std::vector<>
        
        write_cuts(graph, query, cuts, step_pos);
        // 
    }
    //self_dotplot(graph, query, step_pos);
}

}
}
