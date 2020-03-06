#include "subcommand.hpp"
#include "args.hxx"
#include "algorithms/xp.hpp"
#include <httplib.h>

namespace odgi {

    using namespace odgi::subcommand;
    using namespace xp;
    using namespace httplib;

    int main_server(int argc, char** argv) {

        for (uint64_t i = 1; i < argc-1; ++i) {
            argv[i] = argv[i+1];
        }
        const std::string prog_name = "odgi server";
        argv[0] = (char*)prog_name.c_str();
        --argc;

        args::ArgumentParser parser("start a HTTP server on localhost://3000 with a given index file to query a pangenome position");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the index from this file", {'i', "idx"});

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
            std::cerr << "Please enter a file to read the index from." << std::endl;
            exit(1);
        }

        XP path_index;
        std::ifstream in;
        in.open(args::get(dg_in_file));
        path_index.load(in);
        in.close();

        Server svr;

        svr.Get("/hi", [](const Request& req, Response& res) {
            res.set_content("Hello World!", "text/plain");
        });

        svr.Get(R"(/(\w+)/(\d+))", [&](const Request& req, Response& res) {
            /*
            for (size_t i = 0; i < req.matches.size(); i++) {
                std::cout << req.matches[i] << std::endl;
            }
             */
            auto path_name = req.matches[1];
            auto nuc_pos_1 = req.matches[2];
            std::cout << "GOT REQUEST : path name: " << path_name << "; 1-based nucleotide position: " << nuc_pos_1 << std::endl;
            size_t nuc_pos_0 = std::stoi(nuc_pos_1) - 1;
            size_t pan_pos = 0;
            if (path_index.has_path(path_name)) {
                if (path_index.has_position(path_name, nuc_pos_0)) {
                    pan_pos = path_index.get_pangenome_pos(path_name, nuc_pos_0) + 1;
                }
            }
            std::cout << "SEND RESPONSE: pangenome position: " << pan_pos << std::endl;
            res.set_content(std::to_string(pan_pos), "text/plain");
        });

        svr.Get("/stop", [&](const Request& req, Response& res) {
            svr.stop();
        });

        std::cout << "http server listening on http://localhost:3000" << std::endl;
        svr.listen("localhost", 3000);

        /*
        // we have a 0-based positioning
        uint64_t nucleotide_pos = args::get(nuc_pos) - 1;
        std::string p_name = args::get(path_name);

        if (!path_index.has_path(p_name)) {
            std::cerr << "The given path name " << p_name << " is not in the index." << std::endl;
            exit(1);
        }

        if (!path_index.has_position(p_name, nucleotide_pos)) {
            std::cerr << "The given path " << p_name << " with nucleotide position " << nuc_pos << " is not in the index." << std::endl;
            exit(1);
        }

        size_t pangenome_pos = path_index.get_pangenome_pos(p_name, nucleotide_pos) + 1;
        cout << pangenome_pos << endl;

        */
        return 0;
    }

    static Subcommand odgi_server("server",
                                  "start a HTTP server on localhost://3000 with a given index file to query a pangenome position",
                                  PIPELINE, 3, main_server);

}
