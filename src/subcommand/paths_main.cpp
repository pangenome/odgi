#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "split.hpp"
#include "position.hpp"
#include <omp.h>
#include "utils.hpp"
#include "algorithms/path_keep.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_paths(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi paths";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Interrogate the embedded paths of a graph. Does not print anything to stdout by default!");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group path_investigation_opts(parser, "[ Path Investigation Options ]");
    args::ValueFlag<std::string> overlaps_file(path_investigation_opts, "FILE", "Read in the path grouping *FILE* to generate the overlap statistics"
                                                               " from. The file must be tab-delimited. The first column lists a"
                                                               " grouping and the second the path itself. Each line has one path entry."
                                                               " For each group the pairwise overlap statistics for each pairing will"
                                                               " be calculated and printed to stdout.", {'O', "overlaps"});
    args::Flag list_names(path_investigation_opts, "list-names", "Print the paths in the graph to stdout. Each path is printed in its"
                                                " own line.", {'L', "list-paths"});
	args::Flag list_path_start_end(path_investigation_opts, "list-path-start-end", "If -L,--list-paths was specified, this additionally prints the start and end positions of each path in additional, tab-delimited coloumns."
						   , {'l', "list-path-start-end"});
    args::Flag write_fasta(path_investigation_opts, "fasta", "Print paths in FASTA format to stdout. One line for the FASTA header, another line for the whole sequence.", {'f', "fasta"});
    args::Flag haplo_matrix(path_investigation_opts, "haplo", "Print to stdout the paths in a path coverage haplotype matrix"
                                                              " based on the graph’s sort order. The output is tab-delimited:"
                                                              " *path.name*, *path.length*, *path.step.count*, *node.1*,"
                                                              " *node.2*, *node.n*. Each path entry is printed in its own line.", {'H', "haplotypes"});
    args::Flag scale_by_node_length(path_investigation_opts, "haplo", "Scale the haplotype matrix cells by node length.", {'N', "scale-by-node-len"});
    args::ValueFlag<std::string> path_delim(path_investigation_opts, "CHAR", "The part of each path name before this delimiter CHAR is a group"
                                                    " identifier. For use with -H, --haplotypes**: it prints an additional, first column **group.name** to stdout.",
                                                    {'D', "delim"});
    args::ValueFlag<std::uint16_t> path_delim_pos(path_investigation_opts, "N", "Consider the N-th occurrence of the delimiter specified with **-D, --delim**"
                                                    " to obtain the group identifier. Specify 1 for the 1st occurrence (default).",
                                                    {'p', "delim-pos"});
    args::ValueFlag<std::string> non_reference_paths(path_investigation_opts, "FILE", "Print to stdout (in BED format) path ranges that are not in the paths listed (by line) in *FILE*.", {'n', "non-reference-paths"});
    args::Group path_modification_opts(parser, "[ Path Modification Options ]");
    args::ValueFlag<std::string> keep_paths_file(path_modification_opts, "FILE", "Keep paths listed (by line) in *FILE*.", {'K', "keep-paths"});
    args::ValueFlag<std::string> drop_paths_file(path_modification_opts, "FILE", "Drop paths listed (by line) in *FILE*.", {'X', "drop-paths"});
    args::ValueFlag<std::string> dg_out_file(path_modification_opts, "FILE", "Write the dynamic succinct variation graph to this file (e.g. *.og*).", {'o', "out"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> threads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi paths.", {'h', "help"});

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
        std::cerr << "[odgi::paths] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!dg_out_file && (keep_paths_file || drop_paths_file)) {
        std::cerr << "[odgi::paths] error: when keeping or dropping paths please specify an output file to store the graph via -o=[FILE], --out=[FILE]." << std::endl;
        return 1;
    }

	if (list_path_start_end && !list_names) {
		std::cerr << "[odgi::paths] error: please specify also -L,--list-path with the -l,--list-path-start-end option!" << std::endl;
		return 1;
	}

    if (path_delim_pos && args::get(path_delim_pos) < 1) {
		std::cerr << "[odgi::paths] error: -p,--delim-pos has to specify a value greater than 0." << std::endl;
		return 1;
	}

	const uint64_t num_threads = args::get(threads) ? args::get(threads) : 1;
    omp_set_num_threads(num_threads);

	graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "paths", args::get(progress), num_threads, graph);
        }
    }

    if (list_path_start_end && list_names) {
    	std::vector<path_handle_t> paths;
		graph.for_each_path_handle([&](const path_handle_t& p) {
			paths.push_back(p);
		});
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
		for (auto path : paths) {
			uint64_t path_len = 0;
			graph.for_each_step_in_path(path, [&](const step_handle_t& s) {
				handle_t h = graph.get_handle_of_step(s);
				path_len += graph.get_length(h);
			});
#pragma omp critical (cout)
			std::cout << graph.get_path_name(path) << "\t" << 1 << "\t" << path_len << std::endl;
		}
	} else if (args::get(list_names)) {
        graph.for_each_path_handle([&](const path_handle_t& p) {
                std::cout << graph.get_path_name(p) << std::endl;
            });
    }

    if (args::get(write_fasta)) {
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                std::cout << ">" << graph.get_path_name(p) << std::endl;
                graph.for_each_step_in_path(
                    p, [&](const step_handle_t& s) {
                           std::cout << graph.get_sequence(graph.get_handle_of_step(s));
                       });
                std::cout << std::endl;
            });
    }

    const uint16_t delim_pos = path_delim_pos ? args::get(path_delim_pos) - 1 : 0;

    auto group_identified_pos = [](const std::string& path_name, char delim, uint16_t delim_pos) -> std::pair<int32_t, int32_t> {
        int32_t pos = -1;
        int32_t cnt = -1;

        while (cnt != delim_pos) {
            ++pos;
            const int32_t current_pos = path_name.find(delim, pos);
            if (current_pos == std::string::npos) {
                return std::make_pair(cnt, pos - 1);
            }
            pos = current_pos;
            ++cnt;
        }

        return std::make_pair(cnt, pos);
    };

    if (args::get(haplo_matrix)) {
        char delim = '\0';
        if (!args::get(path_delim).empty()) {
            delim = args::get(path_delim).at(0);
        }
        
        { // write the header
            stringstream header;
            if (delim) {
                header << "group.name" << "\t";
            }
            header << "path.name" << "\t"
                   << "path.length" << "\t"
                   << "path.step.count";
            graph.for_each_handle(
                [&](const handle_t& handle) {
                    header << "\t" << "node." << graph.get_id(handle);
                });
            std::cout << header.str() << std::endl;
        }
        bool node_length_scale = args::get(scale_by_node_length);
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                std::string full_path_name = graph.get_path_name(p);
                const std::pair<int32_t, int32_t> cnt_pos = delim ? group_identified_pos(full_path_name, delim, delim_pos) : std::make_pair(0, 0);
                if (cnt_pos.first < 0) {
                    std::cerr << "[odgi::paths] error: path name '" << full_path_name << "' has not occurrences of '" << delim << "'." << std::endl;
                    exit(-1);
                } else if (cnt_pos.first != delim_pos) {
                    std::cerr << "[odgi::paths] warning: path name '" << full_path_name << "' has too few occurrences of '" << delim << "'. "
                              << "The " << cnt_pos.first + 1 << "-th occurrence is used." << std::endl;
                }
                std::string group_name = (delim ? full_path_name.substr(0, cnt_pos.second) : "");
                std::string path_name = (delim ? full_path_name.substr(cnt_pos.second+1) : full_path_name);
                uint64_t path_length = 0;
                uint64_t path_step_count = 0;
                std::vector<uint64_t> row(graph.get_node_count());
                graph.for_each_step_in_path(
                    p,
                    [&](const step_handle_t& s) {
                        const handle_t& h = graph.get_handle_of_step(s);
                        path_length += graph.get_length(h);
                        ++path_step_count;
                        row[graph.get_id(h)-1]++;
                    });
                if (delim) {
                    std::cout << group_name << "\t";
                }
                std::cout << path_name << "\t"
                          << path_length << "\t"
                          << path_step_count;
                if (node_length_scale) {
                    for (uint64_t i = 0; i < row.size(); ++i) {
                        std::cout << "\t" << row[i] * graph.get_length(graph.get_handle(i+1));
                    }
                } else {
                    for (uint64_t i = 0; i < row.size(); ++i) {
                        std::cout << "\t" << row[i];
                    }
                }
                std::cout << std::endl;
            });
    }

    if (!args::get(overlaps_file).empty()) {
        std::string line;
        ska::flat_hash_map<std::string, std::vector<std::string> > path_sets;
        ska::flat_hash_map<std::string, std::vector<pos_t> > path_decomposition;
        auto& x = args::get(overlaps_file);
        std::ifstream overlaps_in(x);
        while (std::getline(overlaps_in, line)) {
            // This file should contain tab-delimited lists of path names, one per line
            // the first field is an identifier
            std::vector<string> fields = split(line, '\t');
            assert(fields.size()==2);
            auto& name = fields[0];
            auto& pathset = path_sets[name];
            pathset.push_back(fields[1]);
            path_decomposition[fields[1]] = {};
        }

        std::vector<std::string> path_names;
        for (auto& p : path_decomposition) {
            path_decomposition[p.first] = {};
            path_names.push_back(p.first);
        }

        // the header
        std::cout << "group.name" << "\t"
                  << "query" << "\t"
                  << "target" << "\t"
                  << "overlap" << "\t"
                  << "overlap.frac" << std::endl;

#pragma omp parallel for
        for (uint32_t k = 0; k < path_names.size(); ++k) {
            auto& path_name = path_names.at(k);
            auto& decomposition = path_decomposition[path_name];
            // walk the path, adding each position to the decomposition
            path_handle_t path = graph.get_path_handle(path_name);
            uint64_t pos = 0;
            graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    nid_t id = graph.get_id(h);
                    uint64_t len = graph.get_length(h);
                    for (uint64_t i = 0; i < len; ++i) {
                        decomposition.push_back(make_pos_t(id, i, graph.get_is_reverse(h)));
                    }
                });
        }

        for (auto& s : path_sets) {
            auto& group_name = s.first;
            auto& paths = s.second;
            for (uint32_t i = 0; i < paths.size(); ++i) {
                for (uint32_t j = i+1; j < paths.size(); ++j) {
                    auto& v1 = path_decomposition[paths[i]];
                    auto& v2 = path_decomposition[paths[j]];
                    std::vector<pos_t> v3;
                    std::sort(v1.begin(), v1.end());
                    std::sort(v2.begin(), v2.end());
                    std::set_intersection(v1.begin(),v1.end(),
                                          v2.begin(),v2.end(),
                                          back_inserter(v3));
                    //ska::flat_hash_map<std::string, std::vector<pos_t> > path_decomposition;
                    std::cout << group_name << "\t" << paths[i] << "\t" << paths[j]
                              << "\t" << v3.size()
                              << "\t" << (float)v3.size()/((float)(v1.size()+v2.size())/2) << std::endl;
                }
            }
        }
    }

    if (keep_paths_file || drop_paths_file) {
        //to_keep
        ska::flat_hash_set<path_handle_t> to_keep;
        if (keep_paths_file) {
            std::string line;
            auto& x = args::get(keep_paths_file);
            std::ifstream infile(x);
            while (std::getline(infile, line)) {
                // This file should contain path names, one per line
                auto& name = line;
                if (graph.has_path(name)) {
                    to_keep.insert(graph.get_path_handle(name));
                } else {
                    std::cerr << "[odgi::paths] error: path'" << name
                              << "' does not exist in graph." << std::endl;
                    return 1;
                }
            }
        } else {
            graph.for_each_path_handle([&to_keep](const path_handle_t& p) {
                to_keep.insert(p);
            });
        }
        if (drop_paths_file) {
            std::string line;
            auto& x = args::get(drop_paths_file);
            std::ifstream infile(x);
            while (std::getline(infile, line)) {
                // This file should contain path names, one per line
                auto& name = line;
                if (graph.has_path(name)) {
                    auto p = graph.get_path_handle(name);
                    assert(to_keep.count(p) == 1);
                    to_keep.erase(p);
                } else {
                    std::cerr << "[odgi::paths] error: path'" << name
                              << "' does not exist in graph." << std::endl;
                    return 1;
                }
            }
        }
        graph_t into;
        algorithms::keep_paths(graph, into, to_keep);
        // write the graph
        if (dg_out_file) {
            const std::string outfile = args::get(dg_out_file);
            if (outfile == "-") {
                into.serialize(std::cout);
            } else {
                ofstream f(outfile.c_str());
                into.serialize(f);
                f.close();
            }
        }
    }

    if (non_reference_paths && !args::get(non_reference_paths).empty()) {
        // Check if the node IDs are compacted
        const uint64_t shift = graph.min_node_id();
        if (graph.max_node_id() - shift >= graph.get_node_count()){
            std::cerr << "[odgi::paths] error: the node IDs are not compacted. Please run 'odgi sort' using -O, --optimize to optimize the graph." << std::endl;
            exit(1);
        }

        // Read paths to use as reference paths
        std::vector<path_handle_t> reference_paths;
        std::string line;
        auto& x = args::get(non_reference_paths);
        std::ifstream infile(x);
        while (std::getline(infile, line)) {
            // This file should contain path names, one per line
            auto& name = line;
            if (graph.has_path(name)) {
                reference_paths.push_back(graph.get_path_handle(name));
            } else {
                std::cerr << "[odgi::paths] error: path'" << name
                            << "' does not exist in graph." << std::endl;
                return 1;
            }
        }

        // Set the reference nodes
        atomicbitvector::atomic_bv_t reference_nodes(graph.get_node_count()+1);
#pragma omp parallel for schedule(dynamic,1)
        for (auto &path : reference_paths) {
			graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
                const handle_t handle = graph.get_handle_of_step(step);
                reference_nodes.set(graph.get_id(handle) - shift);
            });
        }

        std::vector<path_handle_t> non_reference_paths;
		graph.for_each_path_handle([&non_reference_paths](const path_handle_t& path) {
			non_reference_paths.push_back(path);
        });

        // Prepare non reference path handles for parallel processing
        std::sort(non_reference_paths.begin(), non_reference_paths.end());
        std::sort(reference_paths.begin(), reference_paths.end());

        non_reference_paths.erase(
            std::remove_if(non_reference_paths.begin(), non_reference_paths.end(), 
            [&reference_paths](const auto &x) { 
                return std::binary_search(reference_paths.begin(), reference_paths.end(), x); 
            }), non_reference_paths.end());

		// Traverse non reference paths to emit non reference ranges
#pragma omp parallel for schedule(dynamic, 1)
		for (auto& path : non_reference_paths) {
            uint64_t start = 0, end = 0;
			graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
                const handle_t handle = graph.get_handle_of_step(step);
                const uint64_t index = graph.get_id(handle) - shift;
                if (reference_nodes.test(index)) {
                    // Emit the previous non reference range, if any
                    if (end > start) {
                        #pragma omp critical (cout)
                            std::cout << graph.get_path_name(path) << "\t" << start << "\t" << end << std::endl;
                    }
                    end += graph.get_length(handle);
                    start = end;
                } else {
                    end += graph.get_length(handle);
                }
            });
        }
    }

    return 0;
}

static Subcommand odgi_paths("paths", "Interrogate the embedded paths of a graph.",
                              PIPELINE, 3, main_paths);


}
