#include <iostream>
#include <fstream>

// New subcommand system provides all the subcommands that used to live here
#include "subcommand/subcommand.hpp"
#include "version.hpp"

using namespace std;
using namespace odgi;

void odgi_help(char** argv) {
    cerr << "odgi: optimized dynamic genome/graph implementation, version " << Version::get_short() << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << odgi::subcommand::PIPELINE << ":" << endl;

    odgi::subcommand::Subcommand::for_each(odgi::subcommand::PIPELINE, [](const odgi::subcommand::Subcommand& command) {
        // Announce every subcommand we have

        // Pad all the names so the descriptions line up
        string name = command.get_name();
        name.resize(14, ' ');
        cerr << "  -- " << name << command.get_description() << endl;
     });

     cerr << endl;// << "For more commands, type `odgi help`." << endl;
 }

// We make sure to compile main for the lowest common denominator architecture.
// This works on GCC and Clang. But we have to declare main and then define it.
int main(int argc, char *argv[]) __attribute__((__target__("arch=x86-64")));

int main(int argc, char *argv[]) {
    // set a higher value for tcmalloc warnings
    setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD", "1000000000000000", 1);

    if (argc == 1) {
        odgi_help(argv);
        return 0;
    }

    const auto* subcommand = odgi::subcommand::Subcommand::get(argc, argv);
    if (subcommand != nullptr) {
        // We found a matching subcommand, so run it
        return (*subcommand)(argc, argv);
    } else {
        // No subcommand found
        cerr << "[odgi] error: command '" << argv[1] << "' not found.\n\nType `odgi` to list the available commands.\n" << endl;
        return 1;
    }
}
