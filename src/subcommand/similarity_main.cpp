#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "split.hpp"
#include <omp.h>
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

// Functions used for the distance matrix generation
inline uint64_t encode_pair(uint32_t v, uint32_t h) {
    return ((uint64_t)v << 32) | (uint64_t)h;
}
inline void decode_pair(uint64_t pair, uint32_t *v, uint32_t *h) {
    *v = pair >> 32;
    *h = pair & 0x00000000FFFFFFFF;
}

int main_similarity(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi similarity";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Provides a sparse similarity matrix for paths or groups of paths. Each line prints in a tab-delimited format to stdout.");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group path_investigation_opts(parser, "[ Path Investigation Options ]");
    args::ValueFlag<std::string> path_delim(path_investigation_opts, "CHAR", "The part of each path name before this delimiter CHAR is a group identifier.",
                                                    {'D', "delim"});
    args::ValueFlag<std::uint16_t> path_delim_pos(path_investigation_opts, "N", "Consider the N-th occurrence of the delimiter specified with **-D, --delim**"
                                                    " to obtain the group identifier. Specify 1 for the 1st occurrence (default).",
                                                        {'p', "delim-pos"});   
    args::Flag distances(path_investigation_opts, "distances", "Provide distances (dissimilarities) instead of similarities. "
                                                             "Outputs additional columns with the Euclidean and Manhattan distances." , {'d', "distances"});
args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> threads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi similarity.", {'h', "help"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (!dg_in_file) {
        std::cerr << "[odgi::similarity] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (path_delim_pos && args::get(path_delim_pos) < 1) {
		std::cerr << "[odgi::similarity] error: -p,--delim-pos has to specify a value greater than 0." << std::endl;
		return 1;
	}

	const uint64_t num_threads = args::get(threads) ? args::get(threads) : 1;
    omp_set_num_threads(num_threads);

	graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(dg_in_file);
        if (infile.size()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                utils::handle_gfa_odgi_input(infile, "similarity", args::get(progress), num_threads, graph);
            }
        }
    }

    const uint16_t delim_pos = path_delim_pos ? args::get(path_delim_pos) - 1 : 0;

    const bool emit_distances = args::get(distances);

    auto group_identified_pos = [](const std::string& path_name, char delim, uint16_t delim_pos) -> std::pair<int32_t, int32_t> {
        int32_t pos = -1;
        int32_t cnt = -1;

        while (cnt != delim_pos) {
            ++pos;
            const int32_t current_pos = path_name.find(delim, pos);
            if (current_pos == std::string::npos) {
                return std::make_pair(cnt, pos - 1);
            }
            pos = current_pos;
            ++cnt;
        }

        return std::make_pair(cnt, pos);
    };

    // We support up to 4 billion paths (there are uint32_t variables in the implementation)

    bool using_delim = !args::get(path_delim).empty();
    char delim = '\0';
    ska::flat_hash_map<std::string, uint32_t> path_group_ids;
    ska::flat_hash_map<path_handle_t, uint32_t> path_handle_group_ids;
    std::vector<std::string> path_groups;
    if (using_delim) {
        delim = args::get(path_delim).at(0);
        uint32_t i = 0;
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                const std::pair<int32_t, int32_t> cnt_pos = delim ? group_identified_pos(graph.get_path_name(p), delim, delim_pos) : std::make_pair(0, 0);
                if (cnt_pos.first < 0) {
                    std::cerr << "[odgi::similarity] error: path name '" << graph.get_path_name(p) << "' has not occurrences of '" << delim << "'." << std::endl;
                    exit(-1);
                } else if (cnt_pos.first != delim_pos) {
                    std::cerr << "[odgi::similarity] warning: path name '" << graph.get_path_name(p) << "' has too few occurrences of '" << delim << "'. "
                            << "The " << cnt_pos.first + 1 << "-th occurrence is used." << std::endl;
                }
                std::string group_name = graph.get_path_name(p).substr(0, cnt_pos.second);
                auto f = path_group_ids.find(group_name);
                if (f == path_group_ids.end()) {
                    path_group_ids[group_name] = i++;
                    path_groups.push_back(group_name);
                }
                path_handle_group_ids[p] = path_group_ids[group_name];
            });
    }

    auto get_path_name
        = (using_delim ?
            (std::function<std::string(const uint32_t&)>)
            [&](const uint32_t& id) { return path_groups[id]; }
            :
            (std::function<std::string(const uint32_t&)>)
            [&](const uint32_t& id) { return graph.get_path_name(as_path_handle(id)); });

    auto get_path_id
        = (using_delim ?
            (std::function<uint32_t(const path_handle_t&)>)
            [&](const path_handle_t& p) {
                return path_handle_group_ids[p];
            }
            :
            (std::function<uint32_t(const path_handle_t&)>)
            [&](const path_handle_t& p) {
                return (uint32_t)as_integer(p);
            });

    std::vector<uint64_t> bp_count;
    if (using_delim) {
        bp_count.resize(path_groups.size());
    } else {
        bp_count.resize(graph.get_path_count() + 1);
    }

    uint32_t path_max = 0;
    graph.for_each_path_handle(
        [&](const path_handle_t& p) {
            path_max = std::max(path_max, (uint32_t)as_integer(p));
        });

    // Initialize first to avoid possible bugs later
    for (uint32_t i = 0; i < path_max; ++i) {
        path_handle_t p = as_path_handle(i + 1);
        bp_count[get_path_id(p)] = 0;
    }

#pragma omp parallel for
    for (uint32_t i = 0; i < path_max; ++i) {
        path_handle_t p = as_path_handle(i + 1);
        uint64_t path_length = 0;
        graph.for_each_step_in_path(
            p,
            [&](const step_handle_t& s) {
                path_length += graph.get_length(graph.get_handle_of_step(s));
            });
#pragma omp critical (bp_count)
        bp_count[get_path_id(p)] += path_length;
    }

    const bool show_progress = args::get(progress);
    std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
    if (show_progress) {
        progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                graph.get_node_count(), "[odgi::similarity] collecting path intersection lengths");
    }

    // ska::flat_hash_map<std::pair<uint64_t, uint64_t>, uint64_t> leads to huge memory usage with deep graphs
    ska::flat_hash_map<uint64_t, uint64_t> path_intersection_length;
    graph.for_each_handle(
        [&](const handle_t& h) {
            ska::flat_hash_map<uint32_t, uint64_t> local_path_lengths;
            size_t l = graph.get_length(h);
            graph.for_each_step_on_handle(
                h,
                [&](const step_handle_t& s) {
                    local_path_lengths[get_path_id(graph.get_path_handle_of_step(s))] += l;
                });

#pragma omp critical (path_intersection_length)
            for (auto& p : local_path_lengths) {
                for (auto& q : local_path_lengths) {
                    path_intersection_length[encode_pair(p.first, q.first)] += std::min(p.second, q.second);
                }
            }

            if (show_progress) {
                progress_meter->increment(1);
            }
        }, true);

    if (show_progress) {
        progress_meter->finish();
    }

    /*if (using_delim) {
        std::cout << "group.a" << "\t"
                    << "group.b" << "\t"
                    << "group.a.length" << "\t"
                    << "group.b.length" << "\t";
    } else {
        std::cout << "path.a" << "\t"
                    << "path.b" << "\t"
                    << "path.a.length" << "\t"
                    << "path.b.length" << "\t";
    }*/
    // Avoid changing column names (we use the more generic ones)
    std::cout << "group.a" << "\t"
            << "group.b" << "\t"
            << "group.a.length" << "\t"
            << "group.b.length" << "\t"
            << "intersection" << "\t";
    
    if (emit_distances) {
        std::cout << "jaccard.distance" << "\t"
                  << "cosine.distance" << "\t"
                  << "dice.distance" << "\t"
                  << "estimated.difference.rate" << "\t"
                  << "euclidean.distance" << "\t"
                  << "manhattan.distance";
    } else {
        std::cout << "jaccard.similarity" << "\t"
                  << "cosine.similarity" << "\t"
                  << "dice.similarity" << "\t"
                  << "estimated.identity";
    }

    std::cout << std::endl;
    for (auto& p : path_intersection_length) {
        uint32_t id_a, id_b;
        decode_pair(p.first, &id_a, &id_b);

        auto& intersection = p.second;

        // From https://stats.stackexchange.com/questions/58706/distance-metrics-for-binary-vectors
        const double jaccard = (double)intersection / (double)(bp_count[id_a] + bp_count[id_b] - intersection);
        const double cosine = (double)intersection / std::sqrt((double)(bp_count[id_a] * bp_count[id_b]));
        const double dice = 2.0 * ((double) intersection / (double)(bp_count[id_a] + bp_count[id_b]));
        const double estimated_identity = 2.0 * jaccard / (1.0 + jaccard);

        std::cout << get_path_name(id_a) << "\t"
                  << get_path_name(id_b) << "\t"
                  << bp_count[id_a] << "\t"
                  << bp_count[id_b] << "\t"
                  << intersection << "\t";

        if (emit_distances) {
            const double euclidian_distance = std::sqrt((double)((bp_count[id_a] + bp_count[id_b] - intersection) - intersection));
            const uint64_t manhattan_distance = (bp_count[id_a] + bp_count[id_b] - intersection) - intersection;
            std::cout << (1.0 - jaccard) << "\t"
                      << (1.0 - cosine) << "\t"
                      << (1.0 - dice) << "\t"
                      << (1.0 - estimated_identity) << "\t"
                      << euclidian_distance << "\t"
                      << manhattan_distance << std::endl;
        } else {
            std::cout << jaccard << "\t"
                      << cosine << "\t"
                      << dice << "\t"
                      << estimated_identity << std::endl;
        }
    }

    return 0;
}

static Subcommand odgi_similarity("similarity", "Provides a sparse similarity matrix for paths or groups of paths.",
                              PIPELINE, 3, main_similarity);


}
