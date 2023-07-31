#include "subcommand.hpp"
#include "args.hxx"
#include "../version.hpp"
#include <cstdint>

namespace odgi {

    using namespace odgi::subcommand;

    int main_version(int argc, char** argv) {

        for (uint64_t i = 1; i < argc-1; ++i) {
            argv[i] = argv[i+1];
        }
        const std::string prog_name = "odgi version";
        argv[0] = (char*)prog_name.c_str();
        --argc;

        args::ArgumentParser parser("Print the version of ODGI to stdout.");
        args::Group version_opts(parser, "[ Version Options ]");
        args::Flag version(parser, "version", "Print only the version (like *v0.4.0-44-g89d022b*).", {'v', "version"});
        args::Flag codename(parser, "codename", "Print only the codename (like *back to old ABI*).", {'c', "codename"});
        args::Flag release(parser, "release", "Print only the release (like *v0.4.0*)", {'r', "release"});
        args::Group program_information(parser, "[ Program Information ]");
        args::HelpFlag help(program_information, "help", "Print a help message for odgi version.", {'h', "help"});

        if (argc == 1) {
            std::cout << Version::get_short() << endl;
            return 0;
        }

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

        if (version) {
            std::cout << Version::get_version() << endl;
        } else if (codename) {
            std::cout << Version::get_codename() << endl;
        } else if (release) {
            std::cout << Version::get_release() << endl;
        } else {
            std::cout << Version::get_short() << endl;
        }

        return 0;
    }

    static Subcommand odgi_version("version", "Print the version of ODGI to stdout.",
                                   PIPELINE, 3, main_version);

}
