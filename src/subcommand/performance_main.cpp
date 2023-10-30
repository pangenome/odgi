#include "subcommand.hpp"

#include "odgi.hpp"
#include "handlegraph/path_position_handle_graph.hpp"
#include "progress.hpp"

#include "args.hxx"
#include "utils.hpp"

namespace odgi {

	using namespace odgi::subcommand;

	int main_performance(int argc, char **argv) {

		// trick argumentparser to do the right thing with the subcommand
		for (uint64_t i = 1; i < argc - 1; ++i) {
			argv[i] = argv[i + 1];
		}
		std::string prog_name = "odgi performance";
		argv[0] = (char *) prog_name.c_str();
		--argc;

		args::ArgumentParser parser(
				"Runs through all paths in parallel and counts the steps.");
		args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
		args::ValueFlag<std::string> input_graph(mandatory_opts, "FILE",
												 "Input file containing the list of graphs to squeeze into the same\n"
												 "  file. The file must contain one graph per line. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!",
												 {'i', "input-graph"});
		args::Group threading_opts(parser, "[ Threading ]");
		args::ValueFlag<uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.",
										   {'t', "threads"});
		args::Group processing_info_opts(parser, "[ Processing Information ]");
		args::Flag progress(parser, "progress", "Print information about the progress to stderr.",
							{'P', "progress"});
		args::Group program_info_opts(parser, "[ Program Information ]");
		args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi squeeze.", {'h', "help"});

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

		const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
		omp_set_num_threads(num_threads);

		odgi::graph_t graph;
		assert(argc > 0);
		{
			const std::string infile = args::get(input_graph);
			if (!infile.empty()) {
				if (infile == "-") {
					graph.deserialize(std::cin);
				} else {
					utils::handle_gfa_odgi_input(infile, "position", args::get(progress), num_threads, graph);
				}
			}
		}

//		std::cerr << graph.get_node_count() << std::endl;

		std::vector<path_handle_t> paths;
		paths.reserve(graph.get_path_count());
		graph.for_each_path_handle([&](const path_handle_t &path) {
			paths.push_back(path);
		});

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
		for (auto path : paths) {
			uint64_t num_steps = 0;
			graph.for_each_step_in_path(path, [&](const step_handle_t &step) {
				num_steps++;
			});
#pragma omp critical (cout)
			std::cerr << graph.get_path_name(path) << ": " << num_steps << std::endl;
		}

		return 0;
	}

	static Subcommand odgi_squeeze("performance", "LJlsjföl multiple graphs in ODGI format into the same file in ODGI format.",
								   PIPELINE, 3, main_performance);


}
