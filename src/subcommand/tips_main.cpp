#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "utils.hpp"
#include "algorithms/stepindex.hpp"
#include "algorithms/tips.hpp"
#include "algorithms/tips_bed_writer_thread.hpp"

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
				"Identifying break point positions relative to given query (reference) path(s) of all the tips in the graph or of tips of given path(s).");
		args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
		args::ValueFlag<std::string> og_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "input"});
		args::Group tips_opts(parser, "[ Tips Options ]");
		args::ValueFlag<std::string> path_in(tips_opts, "STRING", "Specify the query (reference) path to which to walk from all path tips in the graph.", {'p', "path"});
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

		std::string query_path = args::get(path_in);
		// does the path exist in the graph?
		if (!graph.has_path(query_path)) {
			std::cerr << "[odgi::tips] error: the given query (reference) path '" << query_path << "' is not in the graph! Please specify a valid path via -p=[FILE], --path=[FILE]."
					  << std::endl;
			return 1;
		}
		path_handle_t query_path_t = graph.get_path_handle(query_path);
		// make bit vector across nodes to tell us if we have a hit
		// this is a speed up compared to iterating through all steps of a potential node for each walked step
		std::vector<bool> query_handles;
		query_handles.reserve(graph.get_node_count());
		graph.for_each_step_in_path(query_path_t, [&](const step_handle_t &step) {
			handle_t h = graph.get_handle_of_step(step);
			query_handles[number_bool_packing::unpack_number(h)] = true;
		});
		std::vector<path_handle_t> query_paths;
		query_paths.reserve(1);
		query_paths.push_back(query_path_t);

		std::vector<path_handle_t> paths;
		paths.reserve(graph.get_path_count());
		graph.for_each_path_handle([&](const path_handle_t &path) {
			paths.push_back(path);
		});
		algorithms::step_index_t step_index(graph, paths, num_threads, progress);

		// open the bed writer thread
		algorithms::tips_bed_writer bed_writer_thread;
		bed_writer_thread.open_writer();

		auto get_path_begin = [&](const path_handle_t& path) {
			return graph.path_begin(path);
		};
		auto get_next_step = [&](const step_handle_t& step) {
			return graph.get_next_step(step);
		};
		/// walk from the front
		algorithms::walk_tips(graph, paths, query_path_t, query_handles, step_index, num_threads, get_path_begin,
						get_next_step, bed_writer_thread, progress, true);

		auto get_path_back = [&](const path_handle_t& path) {
			return graph.path_back(path);
		};
		auto get_prev_step = [&](const step_handle_t& step) {
			return graph.get_previous_step(step);
		};
		/// walk from the back
		algorithms::walk_tips(graph, paths, query_path_t, query_handles, step_index, num_threads, get_path_back,
						get_prev_step, bed_writer_thread, progress, false);
		bed_writer_thread.close_writer();
		exit(0);
	}

	static Subcommand odgi_tips("tips",
									"Identifying break point positions relative to given query (reference) path(s) of all the tips in the graph or of tips of given path(s).",
									PIPELINE, 3, main_tips);

}
