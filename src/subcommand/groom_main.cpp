#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/groom.hpp"
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_groom(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    const std::string prog_name = "odgi groom";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Harmonize node orientations.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> og_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::ValueFlag<std::string> og_out_file(mandatory_opts, "FILE", "Write the groomed succinct variation graph in ODGI format to *FILE*. A file ending with *.og* is recommended.", {'o', "out"});
    args::Group grooming_opts(parser, "[ Grooming Options ]");
    args::Flag use_dfs(grooming_opts, "use-dfs", "Use depth-first search for grooming.", {'d', "use-dfs"});
	args::ValueFlag<std::string> _target_paths(grooming_opts, "FILE", "Read the paths that should be considered as target paths (references) from this *FILE*. odgi groom tries to force a forward orientation of all steps for the given paths. A path's rank determines it's weight for decision making and is given by its position in the given *FILE*.", {'R', "target-paths"});
	args::Group threading(parser, "[ Threading ]");
	args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi groom.", {'h', "help"});

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

    if (!og_in_file) {
        std::cerr << "[odgi::groom] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!og_out_file) {
        std::cerr << "[odgi::groom] error: please specify an output file to where to store the groomped graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

	const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

	graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(og_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "groom", args::get(progress), num_threads, graph);
            }
        }
    }

	std::vector<path_handle_t> target_paths;
	if (_target_paths) {
		target_paths = utils::load_paths(args::get(_target_paths), graph, "groom");
	}

    graph.apply_ordering(algorithms::groom(graph, progress, target_paths, !args::get(use_dfs)));

    {
        const std::string outfile = args::get(og_out_file);
        if (!outfile.empty()) {
            if (outfile == "-") {
                graph.serialize(std::cout);
            } else {
                ofstream f(outfile.c_str());
                graph.serialize(f);
                f.close();
            }
        }
    }

    return 0;
}

static Subcommand odgi_groom("groom", "Harmonize node orientations.",
                              PIPELINE, 3, main_groom);


}
