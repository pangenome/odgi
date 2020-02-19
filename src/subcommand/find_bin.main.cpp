#include "subcommand.hpp"
#include "args.hxx"

namespace odgi {

    using namespace odgi::subcommand;

    int main_find_bin(int argc, char** argv) {

        for (uint64_t i = 1; i < argc-1; ++i) {
            argv[i] = argv[i+1];
        }
        const std::string prog_name = "odgi find_bin";
        argv[0] = (char*)prog_name.c_str();
        --argc;

        args::ArgumentParser parser("retrieving a bin id by given path:position");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> index_in_file(parser, "FILE", "load the bin index from this file", {'b', "bin-index"});
        args::ValueFlag<std::string> path_position(parser, "STRING", "find the bin for path:position", {'p', "path-position"});
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

        const std::string pos = args::get(path_position);

        std::cout << "\"path_position entered\": " << pos << std::endl;

        // TODO throw an error if not both arguments were entered.

        // TODO Verify that the path:position argument is correct.

        // TODO @ekg I need to verify that the path and its position are in the index.

        // TODO @ekg Read in the index file and return the desired bin id.

        return 0;
    }

    static Subcommand odgi_find_bin("find_bin", "find the bin of a given path:position",
                               PIPELINE, 3, main_find_bin);
}
