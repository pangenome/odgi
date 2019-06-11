#include "subcommand.hpp"
#include "odgi.hpp"
#include "algorithms/kmer.hpp"
#include "args.hxx"
#include "threads.hpp"
#include "algorithms/hash.hpp"
#include "phf.hpp"
#include "algorithms/prune.hpp"
#include "algorithms/remove_high_degree.hpp"
#include <chrono>
#include "mmmultimap.hpp"
#include "mmmultiset.hpp"
#include "kmer.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_kmers(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi kmers";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("show and characterize the kmer space of the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<uint64_t> kmer_length(parser, "K", "the length of the kmers to generate", {'k', "kmer-length"});
    args::ValueFlag<uint64_t> max_furcations(parser, "N", "break at edges that would be induce this many furcations in a kmer", {'e', "max-furcations"});
    args::ValueFlag<uint64_t> max_degree(parser, "N", "remove nodes that have degree greater that this level", {'D', "max-degree"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});
    args::ValueFlag<std::string> temp_base(parser, "FILE", "use this temporary file base during index generation", {'b', "tmp-base"});
    args::Flag kmers_stdout(parser, "", "write the kmers to stdout", {'c', "stdout"});

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

    graph_t graph;
    assert(argc > 0);
    assert(args::get(kmer_length));
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.load(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.load(f);
            f.close();
        }
    }
    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    }
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
        std::vector<std::vector<kmer_t>> buffers(get_thread_count());
        algorithms::for_each_kmer(graph, args::get(kmer_length), args::get(max_furcations), [&](const kmer_t& kmer) {
                int tid = omp_get_thread_num();
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
        for (auto& buffer : buffers) {
            for (auto& kmer : buffer) {
                std::cout << kmer << "\n";
            }
            buffer.clear();
        }
        std::cout.flush();
    } else {
        //ska::flat_hash_map<uint32_t, uint32_t> kmer_table;
        //vector<uint64_t> kmers;
        std::vector<uint64_t> unique_kmers;
        {
            mmmulti::set<uint64_t> kmers(args::get(temp_base));
            uint64_t seen_kmers = 0;
            uint64_t k = args::get(kmer_length);
            algorithms::for_each_kmer(graph, k, args::get(max_furcations), [&](const kmer_t& kmer) {
                    //uint64_t hash = djb2_hash64(kmer.seq.c_str());
                    bool is_dna = true;
                    uint64_t kint = kmer::seq2bit(kmer.seq.c_str(), k, is_dna);
                    if (is_dna) kmers.append(kint);
                });
            kmers.index();
            kmers.for_each_value_count([&](const uint64_t& v, const uint64_t& count) {
                    unique_kmers.push_back(v);
                });
        }
        std::remove(args::get(temp_base).c_str());
        //std::cerr << std::endl;
        //std::sort(kmers.begin(), kmers.end());
        //kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
        boophf_t* bphf = new boomphf::mphf<uint64_t,hasher_t>(unique_kmers.size(),unique_kmers,get_thread_count());
        //kmers.clear();
        std::cerr << "querying kmers" << std::endl;
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        for (auto& hash : unique_kmers) {
            bphf->lookup(hash);
        }
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        auto used_time = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        std::cerr << "done with " << unique_kmers.size() << " @ " << (double)used_time/(double)unique_kmers.size() << "ns/kmer" << std::endl;
        delete bphf;
    }
    return 0;
}

static Subcommand odgi_kmers("kmers", "process and dump the kmers of the graph",
                              PIPELINE, 3, main_kmers);


}
