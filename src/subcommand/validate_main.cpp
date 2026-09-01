#include "subcommand.hpp"
#include "odgi.hpp"
#include "position.hpp"
#include "args.hxx"
#include "split.hpp"
#include "algorithms/bfs.hpp"
#include <omp.h>
#include "utils.hpp"
#include <fstream>
#include <unordered_set>
#include <cctype>

namespace odgi {

    using namespace odgi::subcommand;

    // Reject compressed FASTA: odgi does not link zlib.
    static void reject_compressed_fasta(const std::string &filename) {
        std::ifstream probe(filename, std::ios::binary);
        if (!probe) {
            std::cerr << "[odgi::validate] error: cannot open FASTA file " << filename << std::endl;
            exit(1);
        }
        unsigned char magic[2] = {0, 0};
        probe.read((char *) magic, 2);
        if (probe.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b) {
            std::cerr << "[odgi::validate] error: " << filename
                      << " is gzip-compressed; odgi cannot read it. Decompress it first "
                      << "(for example: zcat in.fa.gz > in.fa)." << std::endl;
            exit(1);
        }
    }

    // Walk a path and concatenate its oriented segment sequences.
    static std::string path_sequence(const odgi::graph_t &graph, const path_handle_t &path) {
        std::string seq;
        graph.for_each_step_in_path(path, [&](const step_handle_t &step) {
            seq.append(graph.get_sequence(graph.get_handle_of_step(step)));
        });
        return seq;
    }

    // Report the first symbol that differs; returns true when identical.
    static bool first_divergence(const std::string &source, const std::string &walked, uint64_t &pos) {
        const uint64_t shared = std::min(source.size(), walked.size());
        for (uint64_t i = 0; i < shared; ++i) {
            if (std::toupper((unsigned char) source[i]) != std::toupper((unsigned char) walked[i])) {
                pos = i;
                return false;
            }
        }
        if (source.size() != walked.size()) {
            pos = shared;
            return false;
        }
        return true;
    }

    int main_validate(int argc, char **argv) {

        // trick argumentparser to do the right thing with the subcommand
        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi validate";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser(
                "Validate a graph checking if the paths are consistent with the graph topology, and optionally with their source sequences.");
		args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
		args::ValueFlag<std::string> og_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1 or GFAz (compressed GFA), but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "input"});
		args::Group seq_opts(parser, "[ Sequence Validation ]");
		args::ValueFlag<std::string> fasta_file(seq_opts, "FILE", "Also check that every path spells out its source sequence in this uncompressed FASTA *FILE*. Sequences are paired to paths by name and compared symbol-exactly (case-insensitive), so N does not match A.", {'r', "fasta"});
        args::Group threading(parser, "[ Threading ]");
        args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
		args::Group processing_info_opts(parser, "[ Processing Information ]");
		args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
        args::Group program_information(parser, "[ Program Information ]");
        args::HelpFlag help(program_information, "help", "Print a help message for odgi validate.", {'h', "help"});

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
        if (argc == 1) {
            std::cout << parser;
            return 1;
        }

        if (!og_file) {
            std::cerr << "[odgi::validate] error: please specify a graph to validate via -i=[FILE], --idx=[FILE]."
                      << std::endl;
            return 1;
        }

		const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

		odgi::graph_t graph;
        assert(argc > 0);
        std::string infile = args::get(og_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
				utils::handle_gfa_odgi_input(infile, "validate", args::get(progress), num_threads, graph);
            }
        }

        omp_set_num_threads(num_threads);

        bool valid_graph = true;

        std::vector<path_handle_t> paths;
        paths.reserve(graph.get_path_count());
        graph.for_each_path_handle([&](const path_handle_t &path) {
            paths.push_back(path);
        });

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (auto path : paths) {
            graph.for_each_step_in_path(path, [&](const step_handle_t &step) {
                if (graph.has_next_step(step)) {
                    step_handle_t next_step = graph.get_next_step(step);
                    handle_t h = graph.get_handle_of_step(step);
                    handle_t next_h = graph.get_handle_of_step(next_step);

                    if (!graph.has_edge(h, next_h)) {
#pragma omp critical (cout)
                        std::cerr << "[odgi::validate] error: the path " << graph.get_path_name(path) << " does not "
                                  << "respect the graph topology: the link "
                                  << graph.get_id(h) << (graph.get_is_reverse(h) ? "-" : "+")
                                  << ","
                                  << graph.get_id(next_h) << (graph.get_is_reverse(next_h) ? "-" : "+")
                                  << " is missing." << std::endl;

                        valid_graph = false;
                    }
                }
            });
        }

        if (fasta_file) {
            const std::string fasta_name = args::get(fasta_file);
            reject_compressed_fasta(fasta_name);

            std::ifstream in(fasta_name);
            if (!in) {
                std::cerr << "[odgi::validate] error: cannot open FASTA file " << fasta_name << std::endl;
                return 1;
            }

            uint64_t identical = 0, divergent = 0, missing_path = 0;
            std::unordered_set<std::string> paired;
            std::string line, name, source;

            // Compare one record as soon as it is complete, then discard it.
            auto check_record = [&]() {
                if (name.empty()) {
                    return;
                }
                if (!graph.has_path(name)) {
                    std::cerr << "[odgi::validate] error: the source " << name
                              << " has no path in the graph." << std::endl;
                    ++missing_path;
                    valid_graph = false;
                    return;
                }
                if (!paired.insert(name).second) {
                    std::cerr << "[odgi::validate] error: the source " << name
                              << " appears more than once in " << fasta_name << "." << std::endl;
                    valid_graph = false;
                    return;
                }
                const path_handle_t path = graph.get_path_handle(name);
                const std::string walked = path_sequence(graph, path);
                uint64_t pos = 0;
                if (first_divergence(source, walked, pos)) {
                    ++identical;
                } else {
                    std::cerr << "[odgi::validate] error: the path " << name
                              << " does not spell out its source sequence: first difference at "
                              << "source position " << (pos + 1)
                              << " (source length " << source.size()
                              << ", path length " << walked.size() << ")";
                    if (pos < source.size() && pos < walked.size()) {
                        std::cerr << ", source symbol " << source[pos]
                                  << ", path symbol " << walked[pos];
                    }
                    std::cerr << "." << std::endl;
                    ++divergent;
                    valid_graph = false;
                }
            };

            while (std::getline(in, line)) {
                if (!line.empty() && line.back() == '\r') {
                    line.pop_back();
                }
                if (!line.empty() && line[0] == '>') {
                    check_record();
                    const size_t end = line.find_first_of(" \t");
                    name = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
                    source.clear();
                } else {
                    source.append(line);
                }
            }
            check_record();

            // Paths the FASTA never named are reported too.
            uint64_t missing_source = 0;
            graph.for_each_path_handle([&](const path_handle_t &path) {
                const std::string path_name = graph.get_path_name(path);
                if (!paired.count(path_name)) {
                    std::cerr << "[odgi::validate] error: the path " << path_name
                              << " has no source in " << fasta_name << "." << std::endl;
                    ++missing_source;
                    valid_graph = false;
                }
            });

            std::cerr << "[odgi::validate] sequence check against " << fasta_name << ": "
                      << identical << " identical, "
                      << divergent << " divergent, "
                      << missing_path << " missing path, "
                      << missing_source << " missing source." << std::endl;
        }

        return (valid_graph ? 0 : 1);
    }

    static Subcommand odgi_validate("validate",
                                    "Validate a graph checking if the paths are consistent with the graph topology, and optionally with their source sequences.",
                                    PIPELINE, 3, main_validate);

}
