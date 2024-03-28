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


bool isNumberAndLessThanOrEqualToX(const std::string& s, double x) {
    char* end = nullptr;
    double val = strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0' && val <= x;
}

bool isValid(const std::string& levels, double max_val) {
    const std::vector<string> levels_vector = split(levels, ',');
    for (auto l : levels_vector) {
        if (!isNumberAndLessThanOrEqualToX(l, max_val)) {
            return false;
        }
    }
    return true;
}

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
   
    args::Group non_ref_opts(parser, "[ Non-ref. Sequence Options ]");
    args::ValueFlag<std::string> non_reference_nodes(non_ref_opts, "FILE", "Print to stdout IDs of nodes that are not in the paths listed (by line) in *FILE*.", {"non-reference-nodes"});
    args::ValueFlag<std::string> non_reference_ranges(non_ref_opts, "FILE", "Print to stdout (in BED format) path ranges that are not in the paths listed (by line) in *FILE*.", {"non-reference-ranges"});

    args::Group seq_type_opts(parser, "[ Sequence Type Options ]");
    args::ValueFlag<std::string> coverage_levels(seq_type_opts, "c1,c2,...,cN", "List of coverage thresholds (number of paths that pass through the node). Print to stdout a TAB-separated table with node ID, node length, and class.", {"coverage-levels"});
    args::ValueFlag<std::string> fraction_levels(seq_type_opts, "f1,f2,...,fN", "List of fraction thresholds (fraction of paths that pass through the node). Print to stdout a TAB-separated table with node ID , node length, and class.')]", {"fraction-levels"});
    args::Flag path_range_class(seq_type_opts, "N", "Instead of node information, print to stdout a TAB-separated table with path range and class.", {"path-range-class"});
    
    args::Group common_opts(parser, "[ Common Options ]");
    args::ValueFlag<uint64_t> min_size(common_opts, "N", "Minimum size (in bps) of nodes (with --non-reference-nodes or --coverage-levels/fraction-levels) or ranges (with --non-reference-ranges or --coverage-levels/fraction-levels together with --path-range-class).", {"min-size"});
    args::Flag show_step_ranges(common_opts, "N", "Show steps (that is, node IDs and strands) (with --non-reference-ranges or --coverage-levels/fraction-levels together with --path-range-class).", {"show-step-ranges"});
    args::ValueFlag<std::string> path_delim(common_opts, "CHAR", "The part of each path name before this delimiter CHAR is a group"
                                                    " identifier. For use with -H/--haplotypes, --non-reference-ranges or --coverage-levels/fraction-levels. "
                                                    "With -H/--haplotypes it prints an additional, first column **group.name** to stdout.",
                                                    {'D', "delim"});
    args::ValueFlag<std::uint16_t> path_delim_pos(common_opts, "N", "Consider the N-th occurrence of the delimiter specified with **-D, --delim**"
                                                    " to obtain the group identifier. Specify 1 for the 1st occurrence (default).",
                                                    {'p', "delim-pos"});

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

    if (non_reference_nodes && non_reference_ranges) {
		std::cerr << "[odgi::paths] error: specify --non-reference-nodes or --non-reference-ranges, not both." << std::endl;
		return 1;
    }

    if (coverage_levels && fraction_levels) {
        std::cerr << "[odgi::paths] error: specify --coverage-levels or --fraction-levels, not both." << std::endl;
        return 1;
    }

    if ((non_reference_nodes || non_reference_ranges) && (coverage_levels || fraction_levels)) {
        std::cerr << "[odgi::paths] error: specify only one of --non-reference-nodes, --non-reference-range, --coverage-levels or --fraction-levels." << std::endl;
        return 1;
    }

    if (coverage_levels && !isValid(args::get(coverage_levels), std::numeric_limits<double>::max())) {
        std::cerr << "[odgi::paths] error: invalid coverage levels. Only comma-separated numbers are accepted." << std::endl;
        return 1;
    }

    if (fraction_levels && !isValid(args::get(fraction_levels), 1.0)) {
        std::cerr << "[odgi::paths] error: invalid fraction levels. Only comma-separated numbers <= 1.0 are accepted." << std::endl;
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

    const uint64_t shift = graph.min_node_id();
    if (
        args::get(haplo_matrix) || 
        (non_reference_nodes && !args::get(non_reference_nodes).empty()) ||
        (non_reference_ranges && !args::get(non_reference_ranges).empty())
        ) {
        // Check if the node IDs are compacted
        if (graph.max_node_id() - shift >= graph.get_node_count()){
            std::cerr << "[odgi::paths] error: the node IDs are not compacted. Please run 'odgi sort' using -O, --optimize to optimize the graph." << std::endl;
            exit(1);
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
                // Initialize first to avoid possible bugs later
                for (uint32_t i = 0; i < graph.get_node_count(); ++i) {
                    row[i] = 0;
                }

                graph.for_each_step_in_path(
                    p,
                    [&](const step_handle_t& s) {
                        const handle_t& h = graph.get_handle_of_step(s);
                        path_length += graph.get_length(h);
                        ++path_step_count;
                        row[graph.get_id(h)-shift]++;
                    });
                if (delim) {
                    std::cout << group_name << "\t";
                }
                std::cout << path_name << "\t"
                          << path_length << "\t"
                          << path_step_count;
                if (node_length_scale) {
                    for (uint64_t i = 0; i < row.size(); ++i) {
                        std::cout << "\t" << row[i] * graph.get_length(graph.get_handle(i+shift));
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

    if (
        (non_reference_nodes && !args::get(non_reference_nodes).empty()) ||
        (non_reference_ranges && !args::get(non_reference_ranges).empty())
        ) {
        const uint64_t min_size_in_bp = min_size ? args::get(min_size) : 0;

        // Read paths to use as reference paths
        std::vector<path_handle_t> reference_paths;
        std::string line;
        auto& x = non_reference_nodes && !args::get(non_reference_nodes).empty() ? args::get(non_reference_nodes) : args::get(non_reference_ranges);
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

        if (non_reference_nodes && !args::get(non_reference_nodes).empty()){
            // Emit non-reference nodes

            // Set non-reference nodes
            atomicbitvector::atomic_bv_t non_reference_nodes(graph.get_node_count());
            for(uint64_t x = 0; x < non_reference_nodes.size(); ++x) {
                if (min_size_in_bp == 0 || graph.get_length(graph.get_handle(x + shift)) >= min_size_in_bp) {
                    non_reference_nodes.set(x);
                }
            }
#pragma omp parallel for schedule(dynamic,1)
            for (auto &path : reference_paths) {
                graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
                    const handle_t handle = graph.get_handle_of_step(step);
                    non_reference_nodes.reset(graph.get_id(handle) - shift);
                });
            }

            // Emit non-reference nodes
            std::cout << "#node.id\tnode.len\tnum.uncalled.bases\tpaths" << std::endl;
            for (auto x : non_reference_nodes) {
                const handle_t handle = graph.get_handle(x + shift);
            
                // Check paths that go through this node, if any
                std::unordered_set<path_handle_t> unique_path_handles;
                graph.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
                    unique_path_handles.insert(graph.get_path_handle_of_step(step));
                });
                std::string result;
                for (const auto& path : unique_path_handles) {
                    if (!result.empty()) {
                        result += ",";
                    }
                    result += graph.get_path_name(path);
                }

                uint64_t n_count = 0;
                for (auto c : graph.get_sequence(handle)) {
                    if (c == 'N' || c == 'n') { // Increment n_count if character is 'N' or 'n'
                        ++n_count;
                    }
                }

                std::cout << graph.get_id(handle) << "\t" << graph.get_length(handle) << "\t" << n_count << "\t" << result << std::endl;
            }
        } else {
            // Emit non-reference ranges

            const bool _show_step_ranges = args::get(show_step_ranges);

            // Set the reference nodes
            atomicbitvector::atomic_bv_t reference_nodes(graph.get_node_count());
    #pragma omp parallel for schedule(dynamic,1)
            for (auto &path : reference_paths) {
                graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
                    const handle_t handle = graph.get_handle_of_step(step);
                    reference_nodes.set(graph.get_id(handle) - shift);
                });
            }

            // Prepare non reference path handles for parallel processing
            std::vector<path_handle_t> non_reference_paths;
            graph.for_each_path_handle([&non_reference_paths](const path_handle_t& path) {
                non_reference_paths.push_back(path);
            });
            std::sort(non_reference_paths.begin(), non_reference_paths.end());
            std::sort(reference_paths.begin(), reference_paths.end());
            non_reference_paths.erase(
                std::remove_if(non_reference_paths.begin(), non_reference_paths.end(), 
                [&reference_paths](const auto &x) { 
                    return std::binary_search(reference_paths.begin(), reference_paths.end(), x); 
                }), non_reference_paths.end());

            // Traverse non reference paths to emit non-reference ranges
            if (_show_step_ranges) {
                std::cout << "#path.name\tstart\tend\tsteps" << std::endl;
            } else {
                std::cout << "#path.name\tstart\tend" << std::endl;
            }
    #pragma omp parallel for schedule(dynamic, 1)
            for (auto& path : non_reference_paths) {
                uint64_t start = 0, end = 0;
                std::vector<step_handle_t> step_range;
                graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
                    const handle_t handle = graph.get_handle_of_step(step);
                    const uint64_t index = graph.get_id(handle) - shift;
                    if (reference_nodes.test(index)) {
                        // Emit the previous non reference range, if any
                        if (end > start && (end - start) >= min_size_in_bp) {
                            if (_show_step_ranges) {
                                std::string step_range_str = "";
                                for (auto& step : step_range) {
                                    const handle_t handle = graph.get_handle_of_step(step);
                                    step_range_str += std::to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+") + ",";
                                }
                                #pragma omp critical (cout)
                                    std::cout << graph.get_path_name(path) << "\t" << start << "\t" << end << "\t" << step_range_str.substr(0, step_range_str.size() - 1) << std::endl; // trim the trailing comma from step_range
                            } else {
                                #pragma omp critical (cout)
                                    std::cout << graph.get_path_name(path) << "\t" << start << "\t" << end << std::endl;
                            }
                        }
                        end += graph.get_length(handle);
                        start = end;
                        if (_show_step_ranges) {
                            step_range.clear();
                        }
                    } else {
                        end += graph.get_length(handle);
                    }
                    if (_show_step_ranges) {
                        step_range.push_back(step);
                    }
                });

                // Emit last non reference range, if any
                if (end > start && (end - start) >= min_size_in_bp) {
                    if (_show_step_ranges) {
                        std::string step_range_str = "";
                        for (auto& step : step_range) {
                            const handle_t handle = graph.get_handle_of_step(step);
                            step_range_str += std::to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+") + ",";
                        }
                        #pragma omp critical (cout)
                            std::cout << graph.get_path_name(path) << "\t" << start << "\t" << end << "\t" << step_range_str.substr(0, step_range_str.size() - 1) << std::endl; // trim the trailing comma from step_range
                    } else {
                        #pragma omp critical (cout)
                            std::cout << graph.get_path_name(path) << "\t" << start << "\t" << end << std::endl;
                    }
                }
            }
        }
    }

    if (coverage_levels || fraction_levels) {
        char delim = '\0';
        if (!args::get(path_delim).empty()) {
            delim = args::get(path_delim).at(0);
        }

        const uint64_t min_size_in_bp = min_size ? args::get(min_size) : 0;

        // Read levels and sort them
        std::vector<double> sorted_levels;
        {
            const std::string levels = coverage_levels ? args::get(coverage_levels) : args::get(fraction_levels);
            const std::vector<string> levels_vector = split(levels, ',');
            for (auto l : levels_vector) {
                sorted_levels.push_back(std::stod(l));
            }
            std::sort(sorted_levels.begin(), sorted_levels.end());
            // Duplicate the first value by inserting it at the beginning (this to manage values < min_value)
            sorted_levels.insert(sorted_levels.begin(), sorted_levels.front());
        }

        std::vector<ska::flat_hash_set<handle_t>> set_of_handles_for_level(sorted_levels.size());

        ska::flat_hash_map<path_handle_t, std::string> phandle_2_name;
		ska::flat_hash_set<std::string> sample_names;
		graph.for_each_path_handle([&](const path_handle_t &path) {
			const std::string path_name = graph.get_path_name(path);

            std::string sample_name;
            if (delim) {
                const std::pair<int32_t, int32_t> cnt_pos = group_identified_pos(path_name, delim, delim_pos);
                if (cnt_pos.first < 0) {
                    std::cerr << "[odgi::paths] error: path name '" << path_name << "' has not occurrences of '" << delim << "'." << std::endl;
                    exit(-1);
                } else if (cnt_pos.first != delim_pos) {
                    std::cerr << "[odgi::paths] warning: path name '" << path_name << "' has too few occurrences of '" << delim << "'. "
                              << "The " << cnt_pos.first + 1 << "-th occurrence is used." << std::endl;
                }           
                sample_name = path_name.substr(0, cnt_pos.second);
            } else {
                sample_name = path_name;
            }

			phandle_2_name[path] = sample_name;
			sample_names.insert(sample_name);
		});

        // Classify nodes in parallel
		graph.for_each_handle([&](const handle_t &h) {
			const uint64_t hl = graph.get_length(h);
			ska::flat_hash_set<std::string> samples;
			graph.for_each_step_on_handle(h, [&](const step_handle_t &occ) {
				const path_handle_t p_h = graph.get_path_handle_of_step(occ);
				samples.insert(phandle_2_name[p_h]);
			});
            // Number of samples or fraction of samples
            const double value = coverage_levels ? samples.size() : (double) samples.size() / (double) sample_names.size();

            for (int64_t i = sorted_levels.size() - 1; i >= 0; --i) {
                // if i == 0, then value < min_level
                if (i == 0 || value >= sorted_levels[i]) {
                    #pragma omp critical (set_of_handles_for_level)
                        set_of_handles_for_level[i].insert(h); // Save the forward handle
                    // #pragma omp critical (cout)
                    //     std::cout << graph.get_id(h) << "\t" << hl << "\t" << value << "\t" << ">= " << sorted_levels[i] << std::endl;
                    break;
                }
            }
		}, true);

        const std::string symbol = coverage_levels ? "c" : "f";

        if (args::get(path_range_class)) {
            const bool _show_step_ranges = args::get(show_step_ranges);
            
            // Traverse all paths to emit classified path ranges
            if (_show_step_ranges) {
                std::cout << "#path.name\tstart\tend\tclass\tsteps" << std::endl;
            } else {
                std::cout << "#path.name\tstart\tend\tclass" << std::endl;
            }
            std::vector<path_handle_t> all_paths;
            graph.for_each_path_handle([&all_paths](const path_handle_t& path) {
                all_paths.push_back(path);
            });

    #pragma omp parallel for schedule(dynamic, 1)
            for (auto& path : all_paths) {
                uint64_t start = 0, end = 0;
                int64_t last_class = -1;
                std::vector<step_handle_t> step_range;
                graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
                    handle_t handle = graph.get_handle_of_step(step);
                    if (graph.get_is_reverse(handle)) {
                        handle = graph.flip(handle); // We saved the forward handle
                    }

                    // Check the class of the node
                    int64_t current_class = -1;
                    for (uint64_t i = 0; i < set_of_handles_for_level.size(); ++i) {
                        if (set_of_handles_for_level[i].count(handle)) {
                            current_class = i;
                            break;
                        }
                    }

                    if (last_class != -1 && (last_class != current_class)) {
                        // Emit the previous range, if any
                        if (end > start && (end - start) >= min_size_in_bp) {
                            std::string seq_class;
                            if (last_class == 0){
                                seq_class = symbol + "<" + utils::to_string_custom(sorted_levels[last_class]);
                            } else if (last_class == set_of_handles_for_level.size() - 1){
                                seq_class = symbol + ">=" + utils::to_string_custom(sorted_levels[last_class]);
                            } else {
                                seq_class = utils::to_string_custom(sorted_levels[last_class]) + "<=" + symbol + "<" + utils::to_string_custom(sorted_levels[last_class + 1]);
                            }
                            if (_show_step_ranges) {
                                std::string step_range_str = "";
                                for (auto& step : step_range) {
                                    const handle_t handle = graph.get_handle_of_step(step);
                                    step_range_str += std::to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+") + ",";
                                }
                                
                                #pragma omp critical (cout)
                                    std::cout << graph.get_path_name(path) << "\t" << start << "\t" << end << "\t" << seq_class << "\t" << step_range_str.substr(0, step_range_str.size() - 1) << std::endl; // trim the trailing comma from step_range
                            } else {
                                #pragma omp critical (cout)
                                    std::cout << graph.get_path_name(path) << "\t" << start << "\t" << end << "\t" << seq_class << std::endl;
                            }
                        }
                        start = end;
                        end += graph.get_length(handle);
                        if (_show_step_ranges) {
                            step_range.clear();
                        }
                    } else {
                        end += graph.get_length(handle);
                    }
                    if (_show_step_ranges) {
                        step_range.push_back(step);
                    }
                    last_class = current_class;
                });

                // Emit last range, if any
                if (end > start && (end - start) >= min_size_in_bp) {
                    std::string seq_class;
                    if (last_class == 0){
                        seq_class = symbol + "<" + utils::to_string_custom(sorted_levels[last_class]);
                    } else if (last_class == set_of_handles_for_level.size() - 1){
                        seq_class = symbol + ">=" + utils::to_string_custom(sorted_levels[last_class]);
                    } else {
                        seq_class = utils::to_string_custom(sorted_levels[last_class]) + "<=" + symbol + "<" + utils::to_string_custom(sorted_levels[last_class + 1]);
                    }
                    if (_show_step_ranges) {
                        std::string step_range_str = "";
                        for (auto& step : step_range) {
                            const handle_t handle = graph.get_handle_of_step(step);
                            step_range_str += std::to_string(graph.get_id(handle)) + (graph.get_is_reverse(handle) ? "-" : "+") + ",";
                        }
                        
                        #pragma omp critical (cout)
                            std::cout << graph.get_path_name(path) << "\t" << start << "\t" << end << "\t" << seq_class << "\t" << step_range_str.substr(0, step_range_str.size() - 1) << std::endl; // trim the trailing comma from step_range
                    } else {
                        #pragma omp critical (cout)
                            std::cout << graph.get_path_name(path) << "\t" << start << "\t" << end << "\t" << seq_class << std::endl;
                    }
                }
            }
        } else {
            std::cout << "#node.id\tnode.len\tclass" << std::endl;
            for (uint64_t i = 0; i < set_of_handles_for_level.size(); ++i) {
                for (auto h : set_of_handles_for_level[i]) {
                    const uint64_t hl = graph.get_length(h);

                    if (hl >= min_size_in_bp) {
                        std::string seq_class;
                        if (i == 0){
                            seq_class = symbol + "<" + utils::to_string_custom(sorted_levels[i]);
                        } else if (i == set_of_handles_for_level.size() - 1){
                            seq_class = symbol + ">=" + utils::to_string_custom(sorted_levels[i]);
                        } else {
                            seq_class = utils::to_string_custom(sorted_levels[i]) + "<=" + symbol + "<" + utils::to_string_custom(sorted_levels[i + 1]);
                        }
                        std::cout << graph.get_id(h) << "\t" << hl << "\t" << seq_class << std::endl;
                    }
                }
            }
        }
    }

    return 0;
}

static Subcommand odgi_paths("paths", "Interrogate the embedded paths of a graph.",
                              PIPELINE, 3, main_paths);


}
