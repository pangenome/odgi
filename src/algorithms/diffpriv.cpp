#include "diffpriv.hpp"
#include <sdsl/bit_vectors.hpp>

namespace odgi {
namespace algorithms {

void diff_priv_worker(const uint64_t tid,
                      const PathHandleGraph& graph,
                      const std::function<void(step_handle_t, step_handle_t)> callback,
                      const std::function<handle_t(std::mt19937&)> sample_handle,
                      std::atomic<uint64_t>& sampled_length,
                      const uint64_t target_length,
                      const double epsilon,
                      const double min_haplotype_freq,
                      const uint64_t bp_limit) {

    std::random_device rd;
    std::mt19937 mt(rd()); // fully random seed
    std::uniform_real_distribution<double> unif(0,1);
    random_selector<> selector{};

    typedef std::vector<std::pair<step_handle_t, step_handle_t>> step_ranges_t;

    // algorithm
    while (sampled_length < target_length) {
        //std::cerr << "length vs target " << sampled_length << " < " << target_length << std::endl;
        // we randomly sample a starting node and orientation, weighted by node length
        handle_t h = sample_handle(mt);
        // we collect all potential forward extensions
        step_ranges_t ranges;
        graph.for_each_step_on_handle(
            h, [&](const step_handle_t& s) {
                ranges.push_back(std::make_pair(s, s));
            });
        uint64_t walk_length = graph.get_length(h);
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
                double u = std::log1p((double)n.second.size()); // utility == log(count)
                double d_u = u - std::log1p((double)n.second.size()-1); // sensitivity
                double w = exp((epsilon * u) / (2 * d_u)); // our weight
                weights.push_back(std::make_pair(w, n.first));
                sum_weights += w;
                //std::cerr << "f=" << n.second.size() << " u=" << u << " d_u=" << d_u << " w=" << w << std::endl;
            }
            // apply the exponential mechanism using weighted sampling
            // first we sample within the range of the sum of weights
            double d = unif(mt) * sum_weights;
            handle_t opt;
            // respect ranges
            double x = 0;
            for (auto& w : weights) {
                if (x + w.first >= d) {
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
            walk_length += graph.get_length(opt);
            if (ranges.size() < min_haplotype_freq) {
                break; // do nothing
            }
            if (ranges.size() >= min_haplotype_freq
                && walk_length >= bp_limit) {
                // get a random range to avoid orientation bias
                auto& r = selector(ranges);
                sampled_length.fetch_add(walk_length);
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
    const bool progress_reporting,
    const bool write_samples) {

    // copy the sequence space of the graph into priv
    graph.for_each_handle([&](const handle_t& h) {
        priv.create_handle(graph.get_sequence(h), graph.get_id(h));
    });

    uint64_t graph_bp = 0;
    // make a rank select dictionary over our sequence space
    // to get a node randomly distributed in the total length of the graph
    graph.for_each_handle(
        [&](const handle_t& h) {
            graph_bp += graph.get_length(h);
        });

    std::atomic<uint64_t> sampled_length(0);
    uint64_t target_length = graph_bp * target_coverage;

    sdsl::bit_vector graph_bv(graph_bp);
    graph_bp = 0;
    graph.for_each_handle(
        [&](const handle_t& h) {
            graph_bv[graph_bp] = 1;
            graph_bp += graph.get_length(h);
        });
    sdsl::bit_vector::rank_1_type graph_bv_rank;
    sdsl::util::assign(graph_bv_rank, sdsl::bit_vector::rank_1_type(&graph_bv));
    std::uniform_int_distribution<uint64_t> dis_graph_pos = std::uniform_int_distribution<uint64_t>(0, graph_bp-1);
    std::uniform_int_distribution<uint64_t> flip(0, 1);
    auto sample_handle = [&](std::mt19937& gen) {
        uint64_t pos = dis_graph_pos(gen)+1;
        uint64_t id = graph_bv_rank(pos);
        //std::cerr << "pos=" << pos << " id=" << id << std::endl;
        handle_t h = graph.get_handle(id, flip(gen));
        return h;
    };

    std::unique_ptr<progress_meter::ProgressMeter> sampling_progress;
    if (progress_reporting) {
        std::string banner = "[odgi::priv] exponential mechanism sampling subpaths:";
        sampling_progress = std::make_unique<progress_meter::ProgressMeter>(target_length, banner);
    }
    std::atomic<size_t> written_paths(0);
    // handles our sampled ranges
    auto writer_callback = [&](step_handle_t a, step_handle_t b) {
        // write the walk
        std::stringstream ss;
        ss << "hap" << ++written_paths;
        std::string name = ss.str();
        path_handle_t p;
#pragma omp critical (priv_path_create)
        {
            p = priv.create_path_handle(name);
        }
        uint64_t sample_length = 0;
        for (step_handle_t s = a;
             ;
             s = graph.get_next_step(s)) {
            handle_t h = graph.get_handle_of_step(s);
            sample_length += graph.get_length(h);
            handle_t j = priv.get_handle(graph.get_id(h), graph.get_is_reverse(h));
            priv.append_step(p, j);
            if (s == b) break;
        }
        if (progress_reporting) {
            sampling_progress->increment(sample_length);
        }
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
                             sample_handle,
                             std::ref(sampled_length),
                             target_length,
                             epsilon,
                             min_haplotype_freq,
                             bp_limit);
    }

    // stuff happens

    for (uint64_t t = 0; t < nthreads; ++t) {
        workers[t].join();
    }

    if (progress_reporting) {
        sampling_progress->finish();
    }

    // embed edges
    std::vector<path_handle_t> paths;
    priv.for_each_path_handle([&](const path_handle_t& p) {
        paths.push_back(p);
    });
#pragma omp parallel for
    for (auto& p : paths) {
        priv.for_each_step_in_path(
            p,
            [&](const step_handle_t& s) {
                if (priv.has_next_step(s)) {
                    priv.create_edge(
                        priv.get_handle_of_step(s),
                        priv.get_handle_of_step(priv.get_next_step(s)));
                }
            });
    }

    // the emitted graph is not "differentially private" as it may leak information due to its topology
    // further filtering and normalization are indicated to achieve a differentially private result
    // or, the path sequences can be extracted in FASTA and the graph rebuilt

}

}
}
