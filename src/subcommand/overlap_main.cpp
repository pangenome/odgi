#include "subcommand.hpp"
#include "odgi.hpp"
#include "position.hpp"
#include "args.hxx"
#include "subgraph/region.hpp"
#include <omp.h>
#include <mutex>
#include "utils.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    int main_overlap(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi overlap";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("Find the paths touched by given input paths.");
        args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
        args::ValueFlag<std::string> og_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "input"});
        args::Group overlap_opts(parser, "[ Overlap Options ]");
        args::ValueFlag<std::string> subset_paths(overlap_opts, "FILE",
                                                  "Perform the search considering only the paths specified in the FILE. "
                                                  "The file must contain one path name per line and a subset of all paths can be specified."
                                                  "When searching the overlaps, only these paths will be considered.",
                                                  {'s', "subset-paths"});

        args::ValueFlag<std::string> path_name(overlap_opts, "PATH_NAME", "Perform the search of the given path STRING in the graph.",
                                               {'r', "path"});
        args::ValueFlag<std::string> path_file(overlap_opts, "FILE", "Report the search results only for the paths listed in FILE.",
                                               {'R', "paths"});

        args::ValueFlag<std::string> bed_input(overlap_opts, "FILE", "A BED FILE of ranges in paths in the graph.",
                                               {'b', "bed-input"});
        args::Group threading_opts(parser, "[ Threading ]");
        args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
		args::Group processing_info_opts(parser, "[ Processing Information ]");
		args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
        args::Group program_info_opts(parser, "[ Program Information ]");
        args::HelpFlag help(program_info_opts, "help", "display this help summary", {'h', "help"});

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
        if (argc == 1) {
            std::cout << parser;
            return 1;
        }

        if (!og_file || args::get(og_file).empty()) {
            std::cerr << "[odgi::overlap] error: please specify an input graph via -i=[FILE], --idx=[FILE]." << std::endl;
            return 1;
        }

        if ((!path_name || args::get(path_name).empty()) && (!path_file || args::get(path_file).empty()) && (!bed_input || args::get(bed_input).empty())) {
            std::cerr << "[odgi::overlap] error: please specify an input path (-r/--path), a list of paths (with -R/--paths), or a list of path ranges (-b/--bed-input)." << std::endl;
            return 1;
        }

		const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

        odgi::graph_t graph;
        assert(argc > 0);
        std::string infile = args::get(og_file);
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "overlap", args::get(progress), num_threads, graph);
        }

        omp_set_num_threads(num_threads);

        std::vector<path_handle_t> paths_to_consider;
        if (subset_paths && !args::get(subset_paths).empty()) {
            std::ifstream refs(args::get(subset_paths));
            std::string line;
            while (std::getline(refs, line)) {
                if (!line.empty()) {
                    if (!graph.has_path(line)) {
                        std::cerr << "[odgi::overlap] error: path " << line << " not found in graph" << std::endl;
                        exit(1);
                    }

                    paths_to_consider.push_back(graph.get_path_handle(line));
                }
            }
        } else {
            // All paths
            paths_to_consider.reserve(graph.get_path_count());
            graph.for_each_path_handle([&](const path_handle_t path) {
                paths_to_consider.push_back(path);
            });
        }


        std::vector<odgi::path_range_t> path_ranges;

        if (path_name) {
            add_bed_range(path_ranges, graph, args::get(path_name));
        } else if (path_file) {
            // for thing in things
            std::ifstream refs(args::get(path_file));
            std::string line;
            while (std::getline(refs, line)) {
                add_bed_range(path_ranges, graph, line);
            }
        } else {// if (bed_input) {
            std::ifstream bed_in(args::get(bed_input));
            std::string line;
            while (std::getline(bed_in, line)) {
                add_bed_range(path_ranges, graph, line);
            }
        }

        if (!path_ranges.empty()) {
            std::cout << "#path\tstart\tend\tpath.touched" << std::endl;

//#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
            for (auto &path_range : path_ranges) {
                // Collect handles crossed by the path in the specified range
                uint64_t start = path_range.begin.offset;
                uint64_t end = path_range.end.offset;
                path_handle_t path_handle = path_range.begin.path;

                std::unordered_set<handle_t> handles;

                uint64_t walked = 0;
                auto path_end = graph.path_end(path_handle);
                for (step_handle_t cur_step = graph.path_begin(path_handle);
                     cur_step != path_end && walked < end; cur_step = graph.get_next_step(cur_step)) {
                    handle_t cur_handle = graph.get_handle_of_step(cur_step);
                    uint64_t cur_length = graph.get_length(cur_handle);
                    walked += cur_length;
                    if (walked >= start) {
                        handles.insert(cur_handle);
                    }
                }

                // Collect paths that cross the collected handles
                std::vector<path_handle_t> touched_path_handles;

                std::mutex touched_path_handles_mutex;
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
                for (path_handle_t p_h : paths_to_consider) {
                    if (p_h != path_handle) {
                        bool stop = false;
                        graph.for_each_step_in_path(p_h, [&](const step_handle_t &step) {
                            if (!stop && handles.count(graph.get_handle_of_step(step)) > 0) {
                                {
                                    std::lock_guard<std::mutex> guard(touched_path_handles_mutex);
                                    touched_path_handles.push_back(p_h);
                                }

                                stop = true;
                            }
                        });
                    }
                }

                std::sort(touched_path_handles.begin(), touched_path_handles.end(),
                        [](const path_handle_t& a, const path_handle_t& b) {
                            return as_integer(a) < as_integer(b);
                        });

//#pragma omp critical (cout)
                for (auto touched_path_handle : touched_path_handles) {
                    std::cout << (graph.get_path_name(path_handle)) << "\t" << start << "\t" << end << "\t"
                              << graph.get_path_name(touched_path_handle) << std::endl;
                }
            }
        }

        return 0;
    }

    static Subcommand odgi_overlap("overlap", "Find the paths touched by given input paths.",
                                 PIPELINE, 3, main_overlap);

}
