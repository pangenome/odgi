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
				"Generate a step index and access the position of each step of each path once. Prints the run time of each the index generation and the position fetching to stdout.");
		args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
		args::ValueFlag<std::string> og_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "input"});
		args::Group step_index_opts(parser, "[ Step Index Options ]");
		args::ValueFlag<uint64_t> _step_index_sample_rate(step_index_opts, "N", "The sample rate when building the step index. We index a node only if mod(node_id, step-index-sample-rate) == 0! Number must be dividable by 2. (default: 8).",
														  {'a', "step-index-sample-rate"});
		args::Flag naked_run(step_index_opts, "naked-run", "Only read in the given graph and then exit gracefully", {'n', "naked-run"});
		args::ValueFlag<uint64_t> _iterations(step_index_opts, "N", "The number of position fetching rounds. For each path step a position is fetched N times. (default: 10).",
														  {'b', "iterations"});
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

		uint64_t  step_index_sample_rate = args::get(_step_index_sample_rate) ? args::get(_step_index_sample_rate) : 8;
		if (step_index_sample_rate % 2 != 0) {
			std::cerr << "[odgi::stepindex] error: The given sample rate of " << step_index_sample_rate << " is not dividable by 2. Please provide a different sample rate." << std::endl;
			exit(1);
		}

		uint64_t iterations = args::get(_iterations) ? args::get(_iterations) : 10;

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

		if (naked_run) {
			exit(0);
		}

		omp_set_num_threads(num_threads);

		std::vector<path_handle_t> paths;
		paths.reserve(graph.get_path_count());
		graph.for_each_path_handle([&](const path_handle_t &path) {
			paths.push_back(path);
		});

		auto distribute_seconds = [&](int& days, int& hours, int& minutes, int& seconds, const double& input_seconds) {
			const int cseconds_in_day = 86400;
			const int cseconds_in_hour = 3600;
			const int cseconds_in_minute = 60;
			const int cseconds = 1;
			days = std::floor(input_seconds / cseconds_in_day);
			hours = std::floor(((int)input_seconds % cseconds_in_day) / cseconds_in_hour);
			minutes = std::floor((((int)input_seconds % cseconds_in_day) % cseconds_in_hour) / cseconds_in_minute);
			seconds = ((((int)input_seconds % cseconds_in_day) % cseconds_in_hour) % cseconds_in_minute) / cseconds;
		};

		auto print_time = [&](const double& _seconds) {
			int days = 0, hours = 0, minutes = 0, seconds = 0;
			distribute_seconds(days, hours, minutes, seconds, _seconds);
			std::stringstream buffer;
			buffer << std::setfill('0') << std::setw(2) << days << ":"
				   << std::setfill('0') << std::setw(2) << hours << ":"
				   << std::setfill('0') << std::setw(2) << minutes << ":"
				   << std::setfill('0') << std::setw(2) << seconds;
			return buffer.str();
		};

		std::chrono::time_point<std::chrono::steady_clock> start_time_index_building = std::chrono::steady_clock::now();
		algorithms::step_index_t step_index(graph, paths, num_threads, progress, step_index_sample_rate);
		auto cur_time_index_building = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds_index_building = cur_time_index_building-start_time_index_building;
		std::cout << "[odgi::stepindex] info: Building step index elapsed time: " << print_time(elapsed_seconds_index_building.count()) << std::endl;

		std::chrono::time_point<std::chrono::steady_clock> start_time_position_fetching = std::chrono::steady_clock::now();

		std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
		if (progress) {
			progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(paths.size()*iterations, "[odgi::stepindex::position_fetching] Progress:");
		}

		for (int i = 0; i < iterations; i++) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
			for (auto path: paths) {
				graph.for_each_step_in_path(path, [&](const step_handle_t &step) {
					std::cout << step_index.get_position(step, graph) << std::endl;
				});
				if (progress) {
					progress_meter->increment(1);
				}
			}
		}
		if (progress) {
			progress_meter->finish();
		}
		auto cur_time_position_fetching = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds_position_fetching = cur_time_position_fetching-start_time_position_fetching;
		std::cout << "[odgi::stepindex] info: Position fetching elapsed time: " << print_time(elapsed_seconds_position_fetching.count()) << std::endl;

		return 0;
	}

	static Subcommand odgi_stepindex("stepindex",
									"Generate a step index and access the position of each step of each path once. Prints the run time of each the index generation and the position fetching to stdout.",
									PIPELINE, 3, main_stepindex);

}