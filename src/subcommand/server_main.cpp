#include "subcommand.hpp"
#include "args.hxx"
#include "algorithms/xp.hpp"
#include <httplib.h>
#include <filesystem>

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

        args::ArgumentParser parser("Start a basic HTTP server with a given path index file to go from *path:position* to *pangenome:position* very efficiently.");
        args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
        args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph index from this *FILE*. The file name usually ends with *.xp*.", {'i', "idx"});
        args::ValueFlag<std::string> port(mandatory_opts, "N", "Run the server under this port.", {'p', "port"});
        args::Group http_opts(parser, "[ HTTP Options ]");
        args::ValueFlag<std::string> ip_address(http_opts, "IP", "Run the server under this IP address. If not specified, *IP* will be *localhost*.", {'a', "ip"});
        args::Group program_information(parser, "[ Program Information ]");
        args::HelpFlag help(program_information, "help", "Print a help message for odgi server.", {'h', "help"});

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
            std::cerr << "[odgi::server]: please enter a file to read the index from via -i=[FILE], --idx=[FILE]." << std::endl;
            exit(1);
        }

        if (!port) {
            std::cerr << "[odgi::server]: please enter a port for the server via -p=[N], --port=[N]." << std::endl;
            exit(1);
        }

        XP path_index;
		if (!std::filesystem::exists(args::get(dg_in_file))) {
			std::cerr << "[odgi::" << "panpos" << "] error: the given file \"" << args::get(dg_in_file) << "\" does not exist. Please specify an existing input file in xp format via -i=[FILE], --idx=[FILE]." << std::endl;
			return 1;
		}
        std::ifstream in;
        in.open(args::get(dg_in_file));
        path_index.load(in);
        in.close();

        /*
        const char* pattern = R"(/(\d+)/(\w+))";
        std::regex regexi = std::regex(pattern);
        std::cmatch cm;    // same as std::match_results<const char*> cm;
        std::regex_match ("/3/test",cm,regexi);
        std::cout << "the matches were: " << std::endl;
        for (unsigned i=0; i<cm.size(); ++i) {
            std::cout << "[" << cm[i] << "] " << std::endl;
        }
        std::cout << std::endl;

        // const char* pattern1 = R"(/(\w+)/(\d+))";
        // const char* pattern1 = R"(/([a-zA-Z]*[0-9]*)/(\d+))";
        const char* pattern1 = R"(/(\w*.*)/(\d+))";
        std::regex regexi1 = std::regex(pattern1);
        std::cmatch cm1;    // same as std::match_results<const char*> cm;
        std::regex_match ("/5-/3",cm1,regexi1); // /1741.hr2/3
        std::cout << "the matches were: " << std::endl;
        for (unsigned i=0; i<cm1.size(); ++i) {
            std::cout << "[" << cm1[i] << "] " << std::endl;
        }
        */

        Server svr;

        svr.Get("/hi", [](const Request& req, Response& res) {
            res.set_header("Access-Control-Allow-Origin", "*");
            res.set_header("Access-Control-Expose-Headers", "text/plain");
            res.set_header("Access-Control-Allow-Methods", "GET, POST, DELETE, PUT");
            res.set_content("Hello World!", "text/plain");
            std::cout << "GOT REQUEST : HELLO WORLD!" << std::endl;
        });

        svr.Get(R"(/(\w*.*)/(\d+))", [&](const Request& req, Response& res) {
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
            res.set_header("Access-Control-Allow-Origin", "*");
            res.set_header("Access-Control-Expose-Headers", "text/plain");
            res.set_header("Access-Control-Allow-Methods", "GET, POST, DELETE, PUT");
            res.set_content(std::to_string(pan_pos), "text/plain");
        });

        svr.Get("/stop", [&](const Request& req, Response& res) {
            svr.stop();
        });

        const int p = std::stoi(args::get(port));
        std::string ip;
        if (!ip_address) {
            ip = "localhost";
        } else {
            ip = args::get(ip_address);
        }

        std::cout << "http server listening on http://" << ip << ":" << args::get(port) << std::endl;
        svr.listen(ip.c_str(), p);

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
                                  "Start a basic HTTP server to lift coordinates between path and pangenomic positions.",
                                  PIPELINE, 3, main_server);

}
