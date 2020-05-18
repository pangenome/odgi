#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/simple_components.hpp"
#include <random>
#include <deque>

namespace odgi {

using namespace odgi::subcommand;

int main_unitig(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi unitig";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("output unitigs of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    //args::ValueFlag<std::string> dg_out_file(parser, "FILE", "store the graph self index in this file", {'o', "out"});
    args::Flag fake_fastq(parser, "fake", "write unitigs as FASTQ with fixed quality", {'f', "fake-fastq"});
    args::ValueFlag<uint64_t> unitig_to(parser, "N", "continue unitigs with a random walk in the graph so that they are at least this length", {'t', "sample-to"});
    args::ValueFlag<uint64_t> unitig_plus(parser, "N", "continue unitigs with a random walk in the graph this far past their natural end", {'p', "sample-plus"});
    args::ValueFlag<uint64_t> min_begin_node_length(parser, "N", "only begin unitigs from nodes of this length", {'l', "min-begin-node-length"});

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
        std::cerr << "Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
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

    uint64_t max_id = 0;
    graph.for_each_handle([&](const handle_t& handle) {
            max_id = std::max((uint64_t)graph.get_id(handle), max_id);
        });
    std::vector<bool> seen_handles(max_id+1);

    std::random_device rseed;
    std::mt19937 rgen(rseed()); // mersenne_twister
    uint64_t unitig_num = 0;

    if (args::get(min_begin_node_length)) {
        uint64_t min_length = args::get(min_begin_node_length);
        graph.for_each_handle([&](const handle_t& handle) {
                if (graph.get_length(handle) < min_length) {
                    seen_handles[graph.get_id(handle)] = true;
                }
            });
    }

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
                while (added_fwd < to_add/2
                       && (i = graph.get_degree(curr, false)) > 0) {
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
                while (added_rev < to_add/2
                       && (i = graph.get_degree(curr, true)) > 0) {
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

static Subcommand odgi_unitig("unitig", "emit the unitigs of the graph",
                              PIPELINE, 3, main_unitig);


}
