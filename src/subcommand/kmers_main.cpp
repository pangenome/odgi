#include "subcommand.hpp"
#include "odgi.hpp"
#include "algorithms/kmer.hpp"
#include "args.hxx"
#include <omp.h>
#include "algorithms/hash.hpp"
#include "phf.hpp"
#include "algorithms/prune.hpp"
#include "algorithms/remove_high_degree.hpp"
#include <chrono>
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_kmers(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    const std::string prog_name = "odgi kmers";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Display and characterize the kmer space of a graph.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::ValueFlag<uint64_t> kmer_length(mandatory_opts, "K", "The kmer length to generate kmers from.", {'k', "kmer-length"});
    args::Group kmer_opts(parser, "[ Kmer Options ]");
    args::ValueFlag<uint64_t> max_furcations(kmer_opts, "N", "Break at edges that would be induce this many furcations in a kmer.", {'e', "max-furcations"});
    args::ValueFlag<uint64_t> max_degree(kmer_opts, "N", "Don't take nodes into account that have a degree greater than N.", {'D', "max-degree"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<int> threads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Flag kmers_stdout(kmer_opts, "", "Write the kmers to stdout. Kmers are line-separated.", {'c', "stdout"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi kmers.", {'h', "help"});

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

    if (!kmer_length) {
        std::cerr << "Please specify a kmer length via -k=[N], --kmer-lenght=[N]." << std::endl;
        return 1;
    }
    assert(args::get(kmer_length));

	const uint64_t num_threads = args::get(threads) ? args::get(threads) : 1;

	graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(dg_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "kmers", args::get(progress), num_threads, graph);
            }
        }
    }

    omp_set_num_threads(num_threads);

    if (args::get(max_degree)) {
        algorithms::remove_high_degree_nodes(graph, args::get(max_degree));
    }

    /*
    if (args::get(max_furcations)) {
        std::vector<edge_t> to_prune = algorithms::find_edges_to_prune(graph, args::get(kmer_length), args::get(max_furcations));
        std::cerr << "edges to prune: " << to_prune.size() << std::endl;
        for (auto& edge : to_prune) {
            graph.destroy_edge(edge);
        }
        std::cerr << "done prune" << std::endl;
    }
    */

    if (args::get(kmers_stdout)) {
        std::vector<std::vector<kmer_t>> buffers(num_threads);

        algorithms::for_each_kmer(graph, args::get(kmer_length), args::get(max_furcations), [&](const kmer_t& kmer) {
                const int tid = omp_get_thread_num();
                auto& buffer = buffers.at(tid);
                buffer.push_back(kmer);
                if (buffer.size() > 1e5) {
#pragma omp critical (cout)
                    {
                        for (auto& kmer : buffer) {
                            std::cout << kmer << "\n";
                        }
                        buffer.clear();
                    }
                }
            });

        // last kmers in the buffer
        for (auto& buffer : buffers) {
            for (auto& kmer : buffer) {
                std::cout << kmer << "\n";
            }
            buffer.clear();
        }

        std::cout.flush();
    } else {
        //ska::flat_hash_map<uint32_t, uint32_t> kmer_table;
        vector<uint64_t> kmers;
        uint64_t seen_kmers = 0;
        algorithms::for_each_kmer(graph, args::get(kmer_length), args::get(max_furcations), [&](const kmer_t& kmer) {
                //int tid = omp_get_thread_num();
//#pragma omp atomic
                uint64_t hash = djb2_hash64(kmer.seq.c_str());
                //if (hash % 31 == 0) {
#pragma omp critical (kmers)
                    kmers.push_back(hash);
                    //}
                /*
#pragma omp critical (cerr)
                if (seen_kmers % 100000 == 0) {
                    std::cerr << seen_kmers << " " << kmer_table.size() << "\r";
                } else {
                    ++seen_kmers;
                }
                */
            });

        std::cerr << std::endl;
        std::sort(kmers.begin(), kmers.end());
        kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
        boophf_t* bphf = new boomphf::mphf<uint64_t,hasher_t>(kmers.size(),kmers,num_threads);
        //kmers.clear();
        std::cerr << "querying kmers" << std::endl;
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        for (auto& hash : kmers) {
            bphf->lookup(hash);
        }
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        auto used_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        std::cerr << "done with " << kmers.size() << " @ " << (double)used_time/(double)kmers.size() << "ns/kmer" << std::endl;
        /*
        algorithms::for_each_kmer(graph, args::get(kmer_length), [&](const kmer_t& kmer) {
                uint64_t hash = djb2_hash64(kmer.seq.c_str());
                bphf->lookup(hash);
#pragma omp atomic
                ++seen_kmers;
            });
        std::cerr << "done with " << seen_kmers << " kmers" << std::endl;
        */
        delete bphf;
    }
    return 0;
}

static Subcommand odgi_kmers("kmers", "Display and characterize the kmer space of a graph.",
                              PIPELINE, 3, main_kmers);


}
