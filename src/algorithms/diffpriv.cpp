#include "diffpriv.hpp"

namespace odgi {
namespace algorithms {

void diff_priv(
    const PathHandleGraph& graph,
    PathHandleGraph& priv,
    double epsilon,
    double target_coverage,
    uint64_t bp_limit) {

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<uint64_t> dist(1, graph.get_node_count());
    std::uniform_int_distribution<uint64_t> coin(0, 1);
    std::uniform_real_distribution<double> unif(0,1);

    typedef std::vector<std::pair<step_handle_t, step_handle_t>> step_ranges_t;

    bool todo = true;
    double step_count = 0;

    // algorithm
    while (todo) {
        // we randomly sample a starting node (todo: step)
        uint64_t rand_id = dist(mt);
        // we collect all steps on the node, picking a random orientation
        handle_t h = graph.get_handle(rand_id, (bool)coin(mt));
        // we collect all potential forward extensions
        step_ranges_t ranges;
        graph.for_each_step_on_handle(
            h, [&](const step_handle_t& s) {
                ranges.push_back(std::make_pair(s, s));
            });
        double initial_count = ranges.size();
        double walk_length = 0;
        // sampling loop
        while (!ranges.empty()) {
            // next handles
            std::map<handle_t, step_ranges_t> nexts;
            for (auto& range : ranges) {
                auto& s = range.second;
                if (graph.has_next_step(s)) {
                    step_handle_t q = graph.get_next_step(s);
                    handle_t n = graph.get_handle_of_step(q);
                    nexts[n].push_back(std::make_pair(range.first, q));
                }
            }
            // compute weights:
            // calculate the utility of each potential path range group extension
            // calculate delta utility
            std::vector<std::pair<double, handle_t>> weights;
            double sum_weights = 0;
            for (auto& n : nexts) {
                double u = (double)n.second.size() / initial_count;
                double d_u = (double)(n.second.size()-1)
                    / std::max(1.0,initial_count-1);
                double w = (epsilon * u) / (2 * d_u);
                weights.push_back(std::make_pair(w, n.first));
                sum_weights += w;
            }
            // normalize weights
            for (auto& w : weights) {
                w.first /= sum_weights;
            }
            std::sort(weights.begin(), weights.end());
            // apply the exponential mechanism using weighted sampling
            double d = unif(mt);
            handle_t opt;
            for (auto& w : weights) {
                if (w.first >= d) {
                    opt = w.second;
                    break;
                }
            }
            // set our ranges to the selected group
            ranges = nexts[opt];
            // check stopping conditions
            // something weighted by utility
            // if we draw >
            walk_length += graph.get_length(opt);
            if (walk_length > bp_limit) {
                // write the walk
                std::cerr << "got walk : ";
                for (step_handle_t s = ranges.front().first;
                     ;
                     s = graph.get_next_step(s)) {
                    handle_t h = graph.get_handle_of_step(s);
                    std::cerr << (graph.get_is_reverse(h) ? "<" : ">")
                              << graph.get_id(h);
                    if (s == ranges.front().second) break;
                }
                std::cerr << std::endl;
                break;
            }
            //ranges.size() / initial_count;
            /*
            if (ranges.size() == 1) {
                break;
            }
            */
        }
        todo = false; // one iteration for testing
    }

    // and look for potential extensions of open path intervals

}

}
}
