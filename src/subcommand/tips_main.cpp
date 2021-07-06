#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "utils.hpp"
#include "algorithms/stepindex.hpp"

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
				"Identifying break points positions relative to given query (reference) path(s) of all the tips in the graph or of given path(s).");
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
		algorithms::step_index_t step_index(graph, paths, num_threads);

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
		for (auto path : paths) {
			// prevent self tips
			if (path == query_path_t) {
				continue;
			}
			bool tip_reached_query = false;
			// collecting tips from the front
			step_handle_t cur_step = graph.path_begin(path);
			handle_t cur_h = graph.get_handle_of_step(cur_step);
			while (!tip_reached_query) {
				// did we already hit the given reference path?
				if (query_handles[number_bool_packing::unpack_number(cur_h)]) {
					std::vector<uint64_t> query_path_positions;
					graph.for_each_step_on_handle(
							cur_h,
							[&](const step_handle_t& s) {
								if (graph.get_path_handle_of_step(s) == query_path_t) {
//#pragma omp critical (cout)
//									std::cerr << "[odgi::tips] error: We reached a tip from the front!" << std::endl;

									// we collect all positions
									// later we will get the smallest, the highest and the median from these
									query_path_positions.push_back(step_index.get_position(s));
								}
							});
					std::sort(query_path_positions.begin(),
							  query_path_positions.end(),
							  [&](const uint64_t & pos_a,
								  const uint64_t & pos_b) {
								  return pos_a < pos_b;
							  });
					double query_pos_median = utils::median_of_sorted_vec(query_path_positions);
					uint64_t query_min_pos = query_path_positions[0]; // 0-based starting position in BED
					uint64_t query_max_pos = query_path_positions[query_path_positions.size() - 1] + 1; // 1-based ending position in BED
#pragma omp critical (cout)
					std::cout << query_path << "\t" << query_min_pos << "\t" << query_max_pos << "\t"
						<< query_pos_median << "\t" << graph.get_path_name(path) << "\t" << step_index.get_position(cur_step) << std::endl;
					// TODO add BED record to queue of BED_writer_thread
					tip_reached_query = true;
				}
				if (graph.has_next_step(cur_step)) {
					cur_step = graph.get_next_step(cur_step);
					cur_h = graph.get_handle_of_step(cur_step);
				} else {
					// did we iterate over all steps and we did not hit the reference?
					tip_reached_query = true;
					// TODO emit a warning here on stderr
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
