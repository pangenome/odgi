#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/xp.hpp"

namespace odgi {

    using namespace odgi::subcommand;
    using namespace xp;

    int main_path_index(int argc, char** argv) {

        for (uint64_t i = 1; i < argc-1; ++i) {
            argv[i] = argv[i+1];
        }
        const std::string prog_name = "odgi path index";
        argv[0] = (char*)prog_name.c_str();
        --argc;

        args::ArgumentParser parser("create a path index for a given graph");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> idx_out_file(parser, "FILE", "store the index in this file", {'o', "out"});

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

        // read in the graph
        graph_t graph;
        assert(argc > 0);
        std::string infile = args::get(dg_in_file);
        if (infile.size()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                ifstream f(infile.c_str());
                graph.deserialize(f);
                f.close();
            }
        }
        // TODO Get graph_handle. @ekg Is the graph_t already a graph_handle?
        // TODO Build the index.
        XP path_index;
        path_index.from_handle_graph(graph);
        // TODO Open stream to idx_out_file.
        // index.serialize_members(out);
        std::cout << "The current graph has " << graph.get_node_count() << " number of nodes." << std::endl;
        return 0;
    }

    static Subcommand odgi_path_index("path_index", "create a path index for a given graph",
                               PIPELINE, 3, main_path_index);


}
