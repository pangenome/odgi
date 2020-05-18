#include "subcommand.hpp"
#include "args.hxx"
#include "../version.hpp"

namespace odgi {

    using namespace odgi::subcommand;

    int main_version(int argc, char** argv) {

        for (uint64_t i = 1; i < argc-1; ++i) {
            argv[i] = argv[i+1];
        }
        const std::string prog_name = "odgi version";
        argv[0] = (char*)prog_name.c_str();
        --argc;

        args::ArgumentParser parser("get the git version of odgi");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::Flag version(parser, "version", "print only the version (like v1.7.0-68-g224e7625)", {'v', "version"});
        args::Flag codename(parser, "codename", "print only the codename (like edgy)", {'c', "codename"});
        args::Flag release(parser, "release", "print only the release (like v1.7.0)", {'r', "release"});

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
        if (argc==1) {
            std::cout << parser;
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

    static Subcommand odgi_version("version",
            "get the git version of odgi",
                                     PIPELINE, 3, main_version);

}
