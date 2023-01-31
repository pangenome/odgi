#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/stepindex.hpp"
#include "algorithms/progress.hpp"

namespace odgi {

	using namespace odgi::subcommand;

	int main_stepindex(int argc, char **argv) {

		// trick argumentparser to do the right thing with the subcommand
		for (uint64_t i = 1; i < argc - 1; ++i) {
			argv[i] = argv[i + 1];
		}
		std::string prog_name = "odgi stepindex";
		argv[0] = (char *) prog_name.c_str();
		--argc;

		args::ArgumentParser parser(
				"Generate a step index from a given graph. If no output file is provided via *-o, --out*, the index will be directly written to *INPUT_GRAPH.stpidx*.");
		args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
		args::ValueFlag<std::string> og_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "input"});
		args::Group step_index_opts(parser, "[ Step Index Options ]");
		args::ValueFlag<uint64_t> _step_index_sample_rate(step_index_opts, "N", "The sample rate when building the step index. We index a node only if mod(node_id, step-index-sample-rate) == 0! Number must be dividable by 2 or 0 to disable sampling. (default: 8).",
														  {'a', "step-index-sample-rate"});
		args::ValueFlag<std::string> dg_out_file(mandatory_opts, "FILE", "Write the created step index to the specified file. A file"
																		 " ending with *.stpidx* is recommended. (default: *INPUT_GRAPH.stpidx*).", {'o', "out"});
		args::Group threading(parser, "[ Threading ]");
		args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
		args::Group processing_info_opts(parser, "[ Processing Information ]");
		args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
		args::Group program_information(parser, "[ Program Information ]");
		args::HelpFlag help(program_information, "help", "Print a help message for odgi validate.", {'h', "help"});

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
			std::cerr << "[odgi::stepindex] error: please specify a graph to index via -i=[FILE], --idx=[FILE]."
					  << std::endl;
			return 1;
		}

		uint64_t step_index_sample_rate = args::get(_step_index_sample_rate);
		if (step_index_sample_rate > 0 && step_index_sample_rate % 2 != 0) {
			std::cerr << "[odgi::stepindex] error: The given sample rate of " << step_index_sample_rate << " is not dividable by 2. Please provide a different sample rate." << std::endl;
			exit(1);
		}

		std::string step_index_out_file = dg_out_file ? args::get(dg_out_file) : (args::get(og_file) + ".stpidx");

		const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

		odgi::graph_t graph;
		assert(argc > 0);
		std::string infile = args::get(og_file);
		if (!infile.empty()) {
			if (infile == "-") {
				graph.deserialize(std::cin);
			} else {
				utils::handle_gfa_odgi_input(infile, "stepindex", args::get(progress), num_threads, graph);
			}
		}

		omp_set_num_threads(num_threads);


		std::vector<path_handle_t> paths;
		paths.reserve(graph.get_path_count());
		graph.for_each_path_handle([&](const path_handle_t &path) {
			paths.push_back(path);
		});

		algorithms::step_index_t step_index(graph, paths, num_threads, progress, step_index_sample_rate);
		step_index.save(step_index_out_file);

		return 0;
	}

	static Subcommand odgi_stepindex("stepindex",
									"Generate a step index and access the position of each step of each path once.",
									PIPELINE, 3, main_stepindex);

}