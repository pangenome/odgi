#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "utils.hpp"

namespace odgi {

	using namespace odgi::subcommand;

	int main_tips(int argc, char **argv) {

		// trick argumentparser to do the right thing with the subcommand
		for (uint64_t i = 1; i < argc - 1; ++i) {
			argv[i] = argv[i + 1];
		}
		std::string prog_name = "odgi tips";
		argv[0] = (char *) prog_name.c_str();
		--argc;

		args::ArgumentParser parser(
				"Identifying break points positions relative to given (reference) path(s) of all the tips in the graph or of given path(s).");
		args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
		args::ValueFlag<std::string> og_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "input"});
		args::Group tips_opts(parser, "[ Tips Options ]");
		args::ValueFlag<std::string> path_in(tips_opts, "STRING", "Specify the (reference) path to which to walk from all path tips in the graph.", {'p', "path"});
		args::Group threading(parser, "[ Threading ]");
		args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
		args::Group processing_info_opts(parser, "[ Processing Information ]");
		args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
		args::Group program_information(parser, "[ Program Information ]");
		args::HelpFlag help(program_information, "help", "Print a help message for odgi tips.", {'h', "help"});

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

		if (!og_file) {
			std::cerr << "[odgi::tips] error: please specify a graph to walk the tips in via -i=[FILE], --idx=[FILE]."
					  << std::endl;
			return 1;
		}

		const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

		odgi::graph_t graph;
		assert(argc > 0);
		std::string infile = args::get(og_file);
		if (!infile.empty()) {
			if (infile == "-") {
				graph.deserialize(std::cin);
			} else {
				utils::handle_gfa_odgi_input(infile, "tips", args::get(progress), num_threads, graph);
			}
		}

		omp_set_num_threads(num_threads);

		std::string ref_path = args::get(path_in);
		// does the path exist in the graph?
		if (!graph.has_path(ref_path)) {
			std::cerr << "[odgi::tips] error: the given (reference) path '" << ref_path << "' is not in the graph! Please specify a valid path via -p=[FILE], --path=[FILE]."
					  << std::endl;
			return 1;
		}
		path_handle_t ref_path_t = graph.get_path_handle(ref_path);
		// make bit vector across nodes to tell us if we have a hit
		std::vector<bool> ref_handles;
		ref_handles.reserve(graph.get_node_count());
		graph.for_each_step_in_path(ref_path_t, [&](const step_handle_t &step) {
			handle_t h = graph.get_handle_of_step(step);
			ref_handles[as_integer(h)] = true;
		});

		std::vector<path_handle_t> paths;
		paths.reserve(graph.get_path_count());
		graph.for_each_path_handle([&](const path_handle_t &path) {
			paths.push_back(path);
		});

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
		for (auto path : paths) {
			// prevent self tips
			if (path == ref_path_t) {
				continue;
			}
			bool tip_reached_reference = false;
			// collecting tips from the front
			step_handle_t cur_step = graph.path_begin(path);
			handle_t cur_h = graph.get_handle_of_step(cur_step);
			while (!tip_reached_reference) {
				// did we already hit the given reference path?
				if (ref_handles[as_integer(cur_h)]) {
					graph.for_each_step_on_handle(
							cur_h,
							[&](const step_handle_t& s) {
								if (graph.get_path_handle_of_step(s) == ref_path_t) {
#pragma omp critical (cout)
									std::cerr << "[odgi::tips] error: We reached a tip from the front!" << std::endl;

									// TODO we need to collect all steps and get the lowest one?!
									// TODO use the XP index to get the path position quickly
									tip_reached_reference = true;
								}
							});
				}
				if (graph.has_next_step(cur_step)) {
					cur_step = graph.get_next_step(cur_step);
					cur_h = graph.get_handle_of_step(cur_step);
				} else {
					// did we iterate over all steps and we did not hit the reference?
					tip_reached_reference = true;
					// TODO emit a warning here
				}
			}

			// TODO come from the back
		}

		exit(0);
	}

	static Subcommand odgi_tips("tips",
									"Identifying break points positions relative to given (reference) path(s) of all the tips in the graph or of given path(s).",
									PIPELINE, 3, main_tips);

}
