#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_view(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi view";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Project a graph into other formats.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group out_opts(parser, "[ Output Options ]");
    args::Flag to_gfa(out_opts, "to_gfa", "Write the graph in GFAv1 format to standard output.", {'g', "to-gfa"});
    args::Flag emit_node_annotation(out_opts, "node_annotation", "Emit node annotations for the graph in GFAv1 format.", {'a', "node-annotation"});
    args::Flag display(out_opts, "display", "Show the internal structures of a graph. Print to stderr the maximum"
                                          " node identifier, the minimum node identifier, the nodes vector, the"
                                          " delete nodes bit vector and the path metadata, each in a separate"
                                          " line.", {'d', "display"});
	args::Group threading(parser, "[ Threading ]");
	args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_information(parser, "[ Program Information ]");
    args::HelpFlag help(program_information, "help", "Print a help message for odgi view.", {'h', "help"});

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
        std::cerr << "[odgi::view] error: Please specify an input file to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
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
                utils::handle_gfa_odgi_input(infile, "view", args::get(progress), num_threads, graph);
            }
        }
    }

    if (args::get(display)) {
        graph.display();
    }
    if (args::get(to_gfa)) {
        graph.to_gfa(std::cout, args::get(emit_node_annotation));
    }

    return 0;
}

static Subcommand odgi_view("view", "Project a graph into other formats.",
                             PIPELINE, 3, main_view);


}
