#include "diffpriv.hpp"

namespace odgi {
namespace algorithms {

void diff_priv_worker(const uint64_t tid,
                      const PathHandleGraph& graph,
                      const std::function<void(step_handle_t, step_handle_t)> callback,
                      const std::atomic<uint64_t>& step_count,
                      const double target_steps,
                      const double epsilon,
                      const double min_haplotype_freq,
                      const uint64_t bp_limit) {

    std::random_device rd;
    std::mt19937 mt(rd()); // fully random seed
    std::uniform_int_distribution<uint64_t> dist(1, graph.get_node_count());
    std::uniform_int_distribution<uint64_t> coin(0, 1);
    std::uniform_real_distribution<double> unif(0,1);
    random_selector<> selector{};

    typedef std::vector<std::pair<step_handle_t, step_handle_t>> step_ranges_t;

    // algorithm
    while (step_count < target_steps) {
        //std::cerr << "steps vs " << step_count << " < " << target_steps << std::endl;
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
        //std::cerr << "doing a thing " << initial_count << std::endl;
        double walk_length = 0;
        // sampling loop
        while (!ranges.empty()) {
            // next handles
            //std::cerr << "step! with " << ranges.size() << " ranges" << std::endl;
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
                double u = (double)n.second.size();
                // infs to remove
                double d_u = u - (double)(n.second.size()-1);
                double w = exp((epsilon * u) / (2 * d_u));
                weights.push_back(std::make_pair(w, n.first));
                sum_weights += w;
                //std::cerr << "u=" << u << " d_u=" << d_u << " w=" << w << std::endl;
            }
            //std::cerr << "sum weights " << sum_weights << std::endl;
            //std::sort(weights.begin(), weights.end());
            /*
            for (auto& w : weights) {
                std::cerr << "option " << w.first << " to " << graph.get_id(w.second) << std::endl;
            }
            */
            // apply the exponential mechanism using weighted sampling
            // first we sample within the range of the sum of weights
            double d = unif(mt) * sum_weights;
            handle_t opt;
            // respect ranges
            double x = 0;
            for (auto& w : weights) {
                if (x + w.first >= d) {
                    //std::cerr << "taking " << w.first << " " << graph.get_id(w.second) << std::endl;
                    opt = w.second;
                    break;
                }
                // they areas, areas
                x += w.first;
            }
            // set our ranges to the selected group
            ranges = nexts[opt];
            // check stopping conditions
            // 1) depth < min_haplotype_freq (2 by default)
            // 2) length > threshold
            // something weighted by utility
            walk_length += graph.get_length(opt);
            if (ranges.size() < min_haplotype_freq) {
                break; // do nothing
            }
            if (ranges.size() > min_haplotype_freq
                && walk_length > bp_limit) {
                // get a random range
                auto& r = selector(ranges);
                callback(r.first, r.second);
                break;
            }
        }
    }
}

void diff_priv(
    const PathHandleGraph& graph,
    MutablePathDeletableHandleGraph& priv,
    const double epsilon,
    const double target_coverage,
    const double min_haplotype_freq,
    const uint64_t bp_limit,
    const uint64_t nthreads,
    const bool write_samples) {

    // copy the sequence space of the graph into priv
    graph.for_each_handle([&](const handle_t& h) {
        priv.create_handle(graph.get_sequence(h), graph.get_id(h));
    });

    std::atomic<uint64_t> step_count(0);
    //std::cerr << "target coverage " << target_coverage << std::endl;
    uint64_t target_steps = graph.get_node_count() * target_coverage;
    //std::cerr << "target steps = " << target_steps << std::endl;
    std::atomic<size_t> written_paths(0);

    //
    auto writer_callback = [&](step_handle_t a, step_handle_t b) {
        // write the walk
        uint64_t range_step_count = 0;
        std::stringstream ss;
        ss << "hap" << ++written_paths;
        std::string name = ss.str();
        path_handle_t p;
#pragma omp critical (priv_path_create)
        {
            p = priv.create_path_handle(name);
        }
        for (step_handle_t s = a;
             ;
             s = graph.get_next_step(s)) {
            handle_t h = graph.get_handle_of_step(s);
            handle_t j = priv.get_handle(graph.get_id(h), graph.get_is_reverse(h));
            priv.append_step(p, j);
            ++range_step_count;
            if (s == b) break;
        }
        step_count.fetch_add(range_step_count);
        if (write_samples) {
            std::stringstream ss;
            ss << name << "\t";
            for (step_handle_t s = a;
                 ;
                 s = graph.get_next_step(s)) {
                handle_t h = graph.get_handle_of_step(s);
                ss << (graph.get_is_reverse(h) ? "<" : ">")
                   << graph.get_id(h);
                if (s == b) break;
            }
#pragma omp critical (cout)
            std::cout << ss.str() << std::endl;
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(nthreads);
    for (uint64_t t = 0; t < nthreads; ++t) {
        workers.emplace_back(&diff_priv_worker,
                             t,
                             std::cref(graph),
                             writer_callback,
                             std::cref(step_count),
                             target_steps,
                             epsilon,
                             min_haplotype_freq,
                             bp_limit);
    }

    // stuff happens

    for (uint64_t t = 0; t < nthreads; ++t) {
        workers[t].join();
    }

    // embed edges
    priv.for_each_path_handle([&](const path_handle_t& p) {
        priv.for_each_step_in_path(
            p,
            [&](const step_handle_t& s) {
                if (priv.has_next_step(s)) {
                    priv.create_edge(
                        priv.get_handle_of_step(s),
                        priv.get_handle_of_step(priv.get_next_step(s)));
                }
            });
    });

    // the emitted graph is not "differentially private" as it may leak information
    // due to its topology
    // further filtering and normalization are indicated to achieve a differentially private result
    // or, the path sequences can be extracted in FASTA and the graph rebuilt

}

}
}
