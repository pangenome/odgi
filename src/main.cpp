// Needed for crash.hpp to work because it uses newer types
#define _POSIX_C_SOURCE 200809L

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <csignal>
#include <getopt.h>
#include <sys/stat.h>

// New subcommand system provides all the subcommands that used to live here
#include "subcommand/subcommand.hpp"

#include "crash.hpp"

using namespace std;
using namespace dsgvg;


void dsgvg_help(char** argv) {
    cerr << "dsgvg: dynamic succinct variation graph tool" << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << dsgvg::subcommand::PIPELINE << ":" << endl;
         
    dsgvg::subcommand::Subcommand::for_each(dsgvg::subcommand::PIPELINE, [](const dsgvg::subcommand::Subcommand& command) {
        // Announce every subcommand we have
        
        // Pad all the names so the descriptions line up
        string name = command.get_name();
        name.resize(14, ' ');
        cerr << "  -- " << name << command.get_description() << endl;
     });
     
     cerr << endl << "For more commands, type `dg help`." << endl;
 }

// We make sure to compile main for the lowest common denominator architecture.
// This works on GCC and Clang. But we have to decalre main and then define it.
int main(int argc, char *argv[]) __attribute__((__target__("arch=x86-64")));

int main(int argc, char *argv[]) {

    // Make sure the system meets system requirements (i.e. has all the instructions we need)
    //preflight_check();
    dsgvg::enable_crash_handling();

    // Set up stack trace support from crash.hpp
    //enable_crash_handling();
    
    // set a higher value for tcmalloc warnings
    setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD", "1000000000000000", 1);

    if (argc == 1) {
        dsgvg_help(argv);
        return 1;
    }
    
    auto* subcommand = dsgvg::subcommand::Subcommand::get(argc, argv);
    if (subcommand != nullptr) {
        // We found a matching subcommand, so run it
        return (*subcommand)(argc, argv);
    } else {
        // No subcommand found
        string command = argv[1];
        cerr << "error:[dg] command " << command << " not found" << endl;
        dsgvg_help(argv);
        return 1;
    }

}
