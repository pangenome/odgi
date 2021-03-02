#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"

namespace odgi {

    using namespace odgi::subcommand;

    int main_explode(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi explode";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("breaks a graph into connected components in their own files in the given directory");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> dg_out_file(parser, "FILE","store the graph with the generated paths in this file", {'o', "out"});
        args::ValueFlag<uint64_t> nthreads(parser, "N", "number of threads to use", {'t', "threads"});
        args::Flag debug(parser, "debug", "print information about the components and the progress to stderr",{'d', "debug"});

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

        if (!dg_in_file) {
            std::cerr
                    << "[odgi::explode] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                    << std::endl;
            return 1;
        }

        if (!dg_out_file) {
            std::cerr
                    << "[odgi::explode] error: please specify an output file to where to store the graph via -o=[FILE], --out=[FILE]."
                    << std::endl;
            return 1;
        }

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

        uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

        // TO DO

        std::string outfile = args::get(dg_out_file);
        if (outfile.size()) {
            if (outfile == "-") {
                graph.serialize(std::cout);
            } else {
                ofstream f(outfile.c_str());
                graph.serialize(f);
                f.close();
            }
        }
        return 0;
    }

    static Subcommand odgi_explode("explode", "breaks a graph into connected components",
                                 PIPELINE, 3, main_explode);


}
