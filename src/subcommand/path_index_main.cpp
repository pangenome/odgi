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
        std::cout << "The current graph has " << graph.get_node_count() << " number of nodes." << std::endl;

        XP path_index;
        path_index.from_handle_graph(graph);
        std::cout << "Indexed " << path_index.path_count << " paths." << std::endl;

        // writ out the index
        std::ofstream out;
        out.open(args::get(idx_out_file));
        std::cout << "Writing index to disk..." << std::endl;
        path_index.serialize_members(out);
        out.close();

        std::cout << "Reading index from disk " << args::get(idx_out_file) << std::endl;
        XP path_index_1;
        std::ifstream in;
        in.open(args::get(idx_out_file));
        path_index_1.load(in);
        in.close();
        std::cout << "Loaded index has " << path_index_1.path_count << " paths." << std::endl;

        size_t bin_id = path_index.get_bin_id("2196", 2, 2);
        std::cout << "Bin id for input \"2196\":2:2 and constructed index is: " << bin_id << std::endl;

        size_t bin_id_1 = path_index_1.get_bin_id("2196", 2, 2);
        std::cout << "Bin id for input \"2196\":2:2 and loaded index is: " << bin_id_1 << std::endl;

        return 0;
    }

    static Subcommand odgi_path_index("path_index", "create a path index for a given graph",
                               PIPELINE, 3, main_path_index);


}
