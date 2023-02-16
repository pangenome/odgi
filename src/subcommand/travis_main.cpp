#include "subcommand.hpp"
#include "odgi.hpp"
#include "gfa_to_handle.hpp"
#include "args.hxx"
#include <cstdio>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <fstream>
#include "utils.hpp"

namespace odgi {

	using namespace odgi::subcommand;

	int main_travis(int argc, char** argv) {

		// trick argumentparser to do the right thing with the subcommand
		for (uint64_t i = 1; i < argc-1; ++i) {
			argv[i] = argv[i+1];
		}
		std::string prog_name = "odgi travis";
		argv[0] = (char*)prog_name.c_str();
		--argc;

		args::ArgumentParser parser("Prepare files for the building of a graph BWT for Travis.");
		args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
		args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
		args::ValueFlag<std::string> seq_out_file(mandatory_opts, "FILE", "Write path name and the sequence to *FILE*. A file ending with *.tsv* is recommended.", {'s', "seq"});
		args::ValueFlag<std::string> node_out_file(mandatory_opts, "FILE", "Write path name and the node identifier for each base to *FILE*. A file ending with *.tsv* is recommended.", {'n', "node"});
		args::Group threading(parser, "[ Threading ]");
		args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
		args::Group processing_information(parser, "[ Processing Information ]");
		args::Flag progress(processing_information, "progress", "Write the current progress to stderr.", {'P', "progress"});
		args::Group program_information(parser, "[ Program Information ]");
		args::HelpFlag help(program_information, "help", "Print a help message for odgi build.", {'h', "help"});
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

		const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;
		graph_t graph;
		assert(argc > 0);
		{
			const std::string infile = args::get(dg_in_file);
			if (!infile.empty()) {
				if (infile == "-") {
					graph.deserialize(std::cin);
				} else {
					utils::handle_gfa_odgi_input(infile, "travis", args::get(progress), num_threads, graph);
				}
			}
		}

		graph.set_number_of_threads(num_threads);

		ofstream seq_out_stream(args::get(seq_out_file));
		ofstream node_out_stream(args::get(node_out_file));

		std::vector<path_handle_t> paths;
		paths.reserve(graph.get_path_count());
		graph.for_each_path_handle([&](const path_handle_t &path) {
			paths.push_back(path);
		});

		std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
		if (progress) {
			progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
					paths.size(), "[odgi::travis::main] Path Progress:");
		}

		for (auto path : paths) {
			graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
				handle_t h = graph.get_handle_of_step(step);
				uint64_t h_id = graph.get_id(h);
				uint64_t h_len = graph.get_length(h);
				std::string seq = graph.get_sequence(h);
				for (uint64_t i = 0; i < h_len; i++) {
					if (i != 0) {
						node_out_stream << " ";
					}
					node_out_stream << h_id;
					seq_out_stream << seq[i];
				}
				if (graph.has_next_step(step)) {
					node_out_stream << "\t";
				}
			});
			node_out_stream << std::endl;
			if (progress) {
				progress_meter->increment(1);
			}
		}

		seq_out_stream.close();
		node_out_stream.close();

		if (progress) {
			progress_meter->finish();
		}

		return 0;
	}

	static Subcommand odgi_travis("travis", "Prepare files for the building of a graph BWT for Travis.",
								 PIPELINE, 3, main_travis);


}
