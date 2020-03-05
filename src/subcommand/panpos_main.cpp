#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/xp.hpp"

namespace odgi {

    using namespace odgi::subcommand;
    using namespace xp;

    int main_panpos(int argc, char** argv) {

        for (uint64_t i = 1; i < argc-1; ++i) {
            argv[i] = argv[i+1];
        }
        const std::string prog_name = "odgi panpos";
        argv[0] = (char*)prog_name.c_str();
        --argc;

        args::ArgumentParser parser("get the pangenome position of a given path and nucleotide position");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the index from this file", {'i', "idx"});
        args::ValueFlag<std::string> path_name(parser, "PATH_NAME", "get the pangenome positon of this path", {'p', "path"});
        args::ValueFlag<uint64_t> nuc_pos(parser, "NUC_POS", "get the pangenome position of this nucleotide position", {'n', "nuc_pos"});

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

        // TODO KEEP IN MIND THAT THE WHOLE THING IS BASED ON A 0 POSITIONING!
        /*

        // writ out the index
        std::ofstream out;
        out.open(args::get(idx_out_file));
        std::cout << "Writing index to " << args::get(idx_out_file) << std::endl;
        path_index.serialize_members(out);
        out.close();

        std::cout << "Reading index from " << args::get(idx_out_file) << std::endl;
        XP path_index_1;
        std::ifstream in;
        in.open(args::get(idx_out_file));
        path_index_1.load(in);
        in.close();

        size_t bin_id_1 = path_index_1.get_pangenome_pos("5", 5);
        std::cout << "Pangenome position \"5\":5 in loaded index is: " << bin_id_1 << std::endl;
        bin_id_1 = path_index_1.get_pangenome_pos("5-", 5);
        std::cout << "Pangenome position \"5-\":5 in loaded index is: " << bin_id_1 << std::endl;
        bin_id_1 = path_index_1.get_pangenome_pos("5", 12);
        std::cout << "Pangenome position \"5-\":12 in loaded index is: " << bin_id_1 << std::endl;

        std::cout << path_index_1.has_path("5") << std::endl;
        std::cout << path_index_1.has_path("5-") << std::endl;
        std::cout << path_index.has_path("34adf") << std::endl;
        std::cout << path_index.has_position("5", 3) << std::endl;
        std::cout << path_index.has_position("5", 45) << std::endl;
        std::cout << path_index.has_position("543", 3) << std::endl;
         */

        return 0;
    }

    static Subcommand odgi_panpos("panpos",
            "get the pangenome position for a given path and nucleotide position",
                                     PIPELINE, 3, main_panpos);

}
