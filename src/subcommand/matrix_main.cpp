#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/matrix_writer.hpp"
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_matrix(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    const std::string prog_name = "odgi matrix";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Write the graph topology in sparse matrix format.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group matrix_opts(parser, "[ Matrix Options ]");
    args::Flag weight_by_edge_depth(matrix_opts, "edge-depth-weight", "Weigh edges by their path depth.", {'e', "edge-depth-weight"});
    args::Flag weight_by_edge_delta(matrix_opts, "delta-weight", "Weigh edges by the inverse id delta.", {'d', "delta-weight"});
	args::Group threading(parser, "[ Threading ]");
	args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi matrix.", {'h', "help"});

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
        std::cerr << "[odgi::matrix] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
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
				utils::handle_gfa_odgi_input(infile, "matrix", args::get(progress), num_threads, graph);
            }
        }
    }

    algorithms::write_as_sparse_matrix(std::cout, graph, args::get(weight_by_edge_depth), args::get(weight_by_edge_delta));

    return 0;
}

static Subcommand odgi_matrix("matrix", "Write the graph topology in sparse matrix format.",
                              PIPELINE, 3, main_matrix);


}
