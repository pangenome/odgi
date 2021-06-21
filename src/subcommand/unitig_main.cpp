#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/simple_components.hpp"
#include <random>
#include <deque>
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_unitig(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    const std::string prog_name = "odgi unitig";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Output unitigs of the graph.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group fastq_opts(parser, "[ FASTQ Options ]");
    args::Flag fake_fastq(fastq_opts, "fake", "Write the unitigs in FASTQ format to stdout with a fixed quality value of *I*.", {'f', "fake-fastq"});
    args::Group unitig_opts(parser, "[ Unitig Options ]");
    args::ValueFlag<uint64_t> unitig_to(unitig_opts, "N", "Continue unitigs with a random walk in the graph so that they have at least the given *N* length.", {'t', "sample-to"});
    args::ValueFlag<uint64_t> unitig_plus(unitig_opts, "N", "Continue unitigs with a random walk in the graph by *N* past their natural end.", {'p', "sample-plus"});
    args::ValueFlag<uint64_t> min_begin_node_length(unitig_opts, "N", "Only begin unitigs collection from nodes which have at least length *N*.", {'l', "min-begin-node-length"});
	args::Group threading(parser, "[ Threading ]");
	args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_information(parser, "[ Program Information ]");
    args::HelpFlag help(program_information, "help", "Print a help message for odgi unitig.", {'h', "help"});

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
        std::cerr << "[odgi::viz] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

	const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

	graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(dg_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "unitig", args::get(progress), num_threads, graph);
            }
        }
    }

    uint64_t max_id = 0;
    graph.for_each_handle([&](const handle_t& handle) {
        max_id = std::max((uint64_t)graph.get_id(handle), max_id);
    });

    std::vector<bool> seen_handles(max_id+1);
    if (args::get(min_begin_node_length)) {
        const uint64_t min_length = args::get(min_begin_node_length);
        graph.for_each_handle([&](const handle_t& handle) {
            if (graph.get_length(handle) < min_length) {
                seen_handles[graph.get_id(handle)] = true;
            }
        });
    }

    std::random_device rseed;
    std::mt19937 rgen(rseed()); // mersenne_twister

    uint64_t unitig_num = 0;
    graph.for_each_handle([&](const handle_t& handle) {
        if (!seen_handles.at(graph.get_id(handle))) {
            seen_handles[graph.get_id(handle)] = true;
            std::unordered_set<handle_t> seen_in_unitig;
            // extend the unitig
            std::deque<handle_t> unitig;
            unitig.push_back(handle);
            handle_t curr = handle;
            seen_in_unitig.insert(curr);
            while (graph.get_degree(curr, false) == 1) {
                graph.follow_edges(curr, false, [&](const handle_t& n) {
                    curr = n;
                });

                if (seen_in_unitig.count(curr)) {
                    break;
                }

                unitig.push_back(curr);
                seen_handles[graph.get_id(curr)] = true;
                seen_in_unitig.insert(curr);
            }
            curr = handle;
            while (graph.get_degree(curr, true) == 1) {
                graph.follow_edges(curr, true, [&](const handle_t& n) {
                    curr = n;
                });

                if (seen_in_unitig.count(curr)) {
                    break;
                }

                unitig.push_front(curr);
                seen_handles[graph.get_id(curr)] = true;
                seen_in_unitig.insert(curr);
            }
            // if we should extend further, do it
            uint64_t unitig_length = 0;
            for (auto& h : unitig) {
                unitig_length += graph.get_length(h);
            }
            uint64_t to_add = 0;
            if (args::get(unitig_plus)) {
                to_add = args::get(unitig_plus) * 2; // bi-ended extension
            }
            if (args::get(unitig_to) > unitig_length) {
                to_add = args::get(unitig_to) - unitig_length;
            }
            uint64_t added_fwd = 0;
            curr = unitig.back();
            uint64_t i = 0;
            while (added_fwd < to_add/2 && (i = graph.get_degree(curr, false)) > 0) {
                std::uniform_int_distribution<uint64_t> idist(0,i);
                uint64_t j = idist(rgen);
                graph.follow_edges(curr, false, [&](const handle_t& h) {
                    if (j == 0) {
                        unitig.push_back(h);
                        added_fwd += graph.get_length(h);
                        curr = h;
                        return false;
                    } else {
                        --j;
                        return true;
                    }
                });
            }
            curr = unitig.front();
            uint64_t added_rev = 0;
            i = 0;
            while (added_rev < to_add/2 && (i = graph.get_degree(curr, true)) > 0) {
                std::uniform_int_distribution<uint64_t> idist(0,i);
                uint64_t j = idist(rgen);
                graph.follow_edges(curr, true, [&](const handle_t& h) {
                    if (j == 0) {
                        unitig.push_front(h);
                        added_rev += graph.get_length(h);
                        curr = h;
                        return false;
                    } else {
                        --j;
                        return true;
                    }
                });
            }
            unitig_length += added_fwd + added_rev;
            if (args::get(fake_fastq)) {
                std::cout << "@";
            } else {
                std::cout << ">";
            }
            std::cout << "unitig" << ++unitig_num;
            std::cout << " length=" << unitig_length;
            std::cout << " path=";
            for (uint64_t i = 0; i < unitig.size(); ++i) {
                auto& h = unitig.at(i);
                std::cout << graph.get_id(h) << (graph.get_is_reverse(h) ? "-" : "+") << (i+1 < unitig.size() ? "," : "");
            }
            std::cout << std::endl;
            for (auto& h : unitig) {
                std::cout << graph.get_sequence(h);
            }
            std::cout << std::endl;
            if (args::get(fake_fastq)) {
                std::cout << "+" << std::endl;
                for (auto& h : unitig) {
                    for (auto& c : graph.get_sequence(h)) {
                        std::cout << "I";
                    }
                }
                std::cout << std::endl;
            }
        }
    });

    return 0;
}

static Subcommand odgi_unitig("unitig", "Output unitigs of the graph.",
                              PIPELINE, 3, main_unitig);


}
