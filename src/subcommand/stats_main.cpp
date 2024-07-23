#include "subcommand.hpp"
#include "odgi.hpp"
#include "IntervalTree.h"
#include "args.hxx"
#include "split.hpp"
#include <omp.h>
#include "algorithms/layout.hpp"
#include "algorithms/weakly_connected_components.hpp"
#include "cover.hpp"
#include "utils.hpp"
#include <filesystem>
#include "split.hpp"

//#define debug_odgi_stats

namespace odgi {

using namespace odgi::subcommand;

int main_stats(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi stats";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Metrics describing a variation graph and its path relationship.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group summary_opts(parser, "[ Summary Options ]");
    args::Flag _summarize(summary_opts, "summarize", "Summarize the graph properties and dimensions. Print to stdout the #nucleotides, #nodes, #edges, #paths, #steps in a tab-delimited format.", {'S', "summarize"});

    args::Flag _weakly_connected_components(summary_opts, "show-components", "Shows the properties of the weakly connected components.", {'W', "weak-connected-components"});

    args::Flag _num_self_loops(summary_opts, "show-self-loops", "Number of nodes with a self-loop.", {'L', "self-loops"});
    args::Flag _show_nondeterministic_edges(summary_opts, "show-edges", "Show nondeterministic edges (those that extend to the same next base).", {'N', "nondeterministic-edges"});
    args::Flag base_content(summary_opts, "base-content", "Describe the base content of the graph. Print to stdout the #A, #C, #G\n"
                                                          "  and #T in a tab-delimited format.", {'b', "base-content"});
    //args::Flag path_coverage(parser, "coverage", "provide a histogram of path coverage over bases in the graph", {'C', "coverage"});
    //args::Flag path_setcov(parser, "setcov", "provide a histogram of coverage over unique sets of paths", {'V', "set-coverage"});
    //args::Flag path_setcov_count(parser, "setcountcov", "provide a histogram of coverage over counts of unique paths", {'Q', "set-count-coverage"});
    //args::Flag path_multicov(parser, "multicov", "provide a histogram of coverage over unique multisets of paths", {'M', "multi-coverage"});
    //args::Flag path_multicov_count(parser, "multicountcov", "provide a histogram of coverage over counts of paths", {'L', "multi-count-coverage"});
    //args::ValueFlag<std::string> path_bedmulticov(parser, "BED", "for each BED entry, provide a table of path coverage over unique multisets of paths in the graph. Each unique multiset of paths overlapping a given BED interval is described in terms of its length relative to the total interval, the number of path traversals, and unique paths involved in these traversals.", {'B', "bed-multicov"});
    args::ValueFlag<std::string> path_delim(summary_opts, "STRING", "The part of each path name before this delimiter is a group identifier, which when specified will ensure that odgi stats collects the summary information per group and not per path.", {'D', "delim"});
	args::Flag _file_size(summary_opts, "file-size", "Show the file size in bytes.", {'f', "file-size"});
	args::ValueFlag<std::string> _pangenome_sequence_class_counts(summary_opts, "DELIM,POS", "Show counted pangenome sequence class counts of all samples. "
                                                                                             "Classes are Private (only one sample visiting the node), Core (all samples visiting the node), and Shell (not Core or Private). "
                                                                                             "The given String determines how to find the sample name in the path names: DELIM,POS. "
                                                                                             "Split the whole path name by DELIM and access the actual sample name at POS of the split result. "
                                                                                             "If the full path name is the sample name, select a DELIM that is not in the path names and set POS to 0. "
                                                                                             "If -m,--multiqc was set, this OPTION has to be set implicitly.", {'a', "pangenome-sequence-class-counts"});
	args::Group sorting_goodness_evaluation_opts(parser, "[ Sorting Goodness Eval Options ]");
	args::ValueFlag<std::string> layout_in_file(sorting_goodness_evaluation_opts, "FILE", "Load the 2D layout coordinates in binary layout format from this *FILE*. The file name usually ends with *.lay*. The sorting goodness evaluation will then be performed for this *FILE*. When the layout coordinates are provided, the mean links length and the sum path nodes distances statistics are evaluated in 2D, else in 1D. Such a file can be generated with *odgi layout*.", {'c', "coords-in"});
	args::Flag mean_links_length(sorting_goodness_evaluation_opts, "mean_links_length", "Calculate the mean links length. This metric is path-guided and"
                                                              " computable in 1D and 2D.", {'l', "mean-links-length"});
	args::Flag dont_penalize_gap_links(sorting_goodness_evaluation_opts, "dont-penalize-gap-links", "Don’t penalize gap links in the mean links length. A gap link is a"
                                                                                                    " link which connects two nodes that are consecutive in the linear"
                                                                                                    " pangenomic order. This option is specifiable only to compute the mean"
                                                                                                    " links length in 1D.", {'g', "no-gap-links"});
	args::Flag sum_of_path_node_distances(sorting_goodness_evaluation_opts, "sum_of_path_node_distances", "Calculate the sum of path nodes distances. This metric is path-guided"
                                                                                                          " and computable in 1D and 2D. For each path, it iterates from node to"
                                                                                                          " node, summing their distances, and normalizing by the path length. In"
                                                                                                          " 1D, if a link goes back in the linearized viewpoint of the graph, this"
                                                                                                          " is penalized (adding 3 times its length in the sum).", {'s', "sum-path-nodes-distances"});
	args::Flag penalize_diff_orientation(sorting_goodness_evaluation_opts, "penalize_diff_orientation", "If a link connects two nodes which have different orientations, this"
                                                                                                        " is penalized (adding 2 times its length in the sum).", {'d', "penalize-different-orientation"});
	args::Flag path_statistics(sorting_goodness_evaluation_opts, "path_statistics", "Display the statistics (mean links length or sum path nodes distances) for each path.", {'p', "path-statistics"});
	
    args::Flag weighted_feedback_arc(sorting_goodness_evaluation_opts, "weighted_feedback_arc", "Compute the sum of weights of all feedback arcs, i.e. backward pointing edges the statistics (the weight is the number of times the edge is traversed by paths).", {'w', "weighted-feedback-arc"});
	args::Flag weighted_reversing_join(sorting_goodness_evaluation_opts, "weighted_reversing_join", "Compute the sum of weights of all reversing joins, i.e. edges joining two in- or two out-sides (the weight is the number of times the edge is traversed by paths).", {'j', "weighted-reversing-join"});
	args::Flag links_length_per_nuc(sorting_goodness_evaluation_opts, "links_length_per_nuc", "Compute the links length per nucleotide, i.e. sum up the links lengths of all paths and divide this value by the nucleotide lengths of all paths. This metric can be used to compare the linearity of different graphs. By default we don't count gap links.", {'q', "links_length_per_nuc"});

    args::Group io_format_opts(parser, "[ IO Format Options ]");
    args::Flag _multiqc(io_format_opts, "multiqc", "Setting this option prints all! statistics in YAML format instead of pseudo TSV to stdout. This includes *-S,--summarize*, *-W,--weak-connected-components*, *-L,--self-loops*, *-b,--base-content*, *-l,--mean-links-length*, *-g,--no-gap-links*, *-s,--sum-path-nodes-distances*, *-f,--file-size*, and *-d,--penalize-different-orientation*. *-p,path-statistics* is still optional. Not applicable to *-N,--nondeterministic-edges*. Overwrites all other given OPTIONs! The output is perfectly curated for the ODGI MultiQC module.", {'m', "multiqc"});
	args::Flag _yaml(io_format_opts, "yaml", "Setting this option prints all selected statistics in YAML format instead of pseudo TSV to stdout.", {'y', "yaml"});
    args::Group processing_information(parser, "[ Processing Information ]");
    args::ValueFlag<uint64_t> threads(processing_information, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_information(parser, "[ Program Information ]");
    args::HelpFlag help(program_information, "help", "Print a help message for odgi stats.", {'h', "help"});

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
        std::cerr << "[odgi::stats] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (!args::get(mean_links_length) && !args::get(sum_of_path_node_distances) && !args::get(weighted_feedback_arc) && !args::get(weighted_reversing_join)){
        if (args::get(path_statistics)){
            std::cerr
                    << "[odgi::stats] error: please specify the -l/--mean-links-length and/or the -s/--sum-path-nodes-distances and/or "
                    "-w/--weighted-feedback-arc and/or -j/--weighted-reversing-join options to use the -p/--path-statistics option."
                    << std::endl;
            return 1;
        }

        if (layout_in_file){
            std::cerr
                    << "[odgi::stats] error: please specify the -l/--mean-links-length and/or the -s/--sum-path-nodes-distances options to specify the layout coordinates (-c/--coords-in option)."
                    << std::endl;
            return 1;
        }
    }

    if (args::get(penalize_diff_orientation) && !args::get(sum_of_path_node_distances)){
        std::cerr
                << "[odgi::stats] error: please specify the -s/--sum-path-nodes-distances option to use the -d/--penalize-different-orientation option."
                << std::endl;
        return 1;
    }
    if (args::get(dont_penalize_gap_links)){
        if (!args::get(mean_links_length)) {
            std::cerr
                    << "[odgi::stats] error: please specify the -l/--mean-links-length option to use the -g/--no-gap-links option."
                    << std::endl;
            return 1;
        } else if (layout_in_file){
            std::cerr
                    << "[odgi::stats] error: The -g/--no-gap-links option can be specified together with the layout coordinates (-c/--coords-in option)."
                    << std::endl;
            return 1;
        }
    }

	const std::string pangenome_sequence_class_counts = args::get(_pangenome_sequence_class_counts);
	std::vector<string> delim_pos = split(pangenome_sequence_class_counts, ',');
    if (_pangenome_sequence_class_counts) {
		if (delim_pos.size() != 2) {
			std::cerr << "[odgi::stats] error: Argument for -a,--pangenome-sequence-classes malformed. Please follow DEL,POS. "
						 "DEL is the delimiter in the paths and POS the position of the sample name after the split with the delimiter." << std::endl;
			exit(1);
		}
    }

	const uint64_t num_threads = args::get(threads) ? args::get(threads) : 1;
	omp_set_num_threads(num_threads);

    graph_t graph;
    assert(argc > 0);
    const std::string infile = args::get(dg_in_file);
    if (!infile.empty()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "stats", args::get(progress), num_threads, graph);
        }
    }

	const uint64_t shift = number_bool_packing::unpack_number(graph.get_handle(graph.min_node_id()));

    if (args::get(mean_links_length) || args::get(sum_of_path_node_distances)) {
		if (number_bool_packing::unpack_number(graph.get_handle(graph.max_node_id())) - shift >= graph.get_node_count()){
			std::cerr << "[odgi::stats] error: the node IDs are not compacted. Please run 'odgi sort' using -O, --optimize to optimize the graph." << std::endl;
			exit(1);
		}
    }
	if (_multiqc || _yaml) {
    	std::cout << "---" << std::endl;
    }

	// per default, we want option -S to be shown
	const bool no_args = !(args::get(_weakly_connected_components) ||
			args::get(_num_self_loops) ||
			args::get(_show_nondeterministic_edges) ||
			args::get(base_content) ||
			path_delim ||
			args::get(_file_size) ||
			_pangenome_sequence_class_counts ||
			args::get(mean_links_length) ||
			args::get(dont_penalize_gap_links) ||
			args::get(sum_of_path_node_distances) ||
			args::get(penalize_diff_orientation) ||
			args::get(path_statistics) ||
			args::get(weighted_feedback_arc) ||
			args::get(weighted_reversing_join) ||
			args::get(links_length_per_nuc) ||
			args::get(_multiqc) ||
			args::get(_yaml));

    if (args::get(_summarize) || _multiqc || no_args) {
        uint64_t length_in_bp = 0, node_count = 0;
        graph.for_each_handle([&](const handle_t& h) {
                length_in_bp += graph.get_length(h);
                ++node_count;
            });

        uint64_t edge_count = graph.get_edge_count();
        uint64_t path_count = graph.get_path_count();

		uint64_t step_count = 0;
		graph.for_each_path_handle([&](const path_handle_t& p) {
			step_count += graph.get_step_count(p);
		});

        if (_multiqc || _yaml) {
        	std::cout << "length: " << length_in_bp << std::endl;
        	std::cout << "nodes: " << node_count << std::endl;
        	std::cout << "edges: " << edge_count << std::endl;
        	std::cout << "paths: " << path_count << std::endl;
			std::cout << "steps: " << step_count << std::endl;
        } else {
			std::cout << "#length\tnodes\tedges\tpaths\tsteps" << std::endl;
			std::cout << length_in_bp << "\t" << node_count << "\t" << edge_count << "\t" << path_count << "\t" << step_count << std::endl;
		}
    }

    if (args::get(_weakly_connected_components) || _multiqc) {
        std::vector<ska::flat_hash_set<handlegraph::nid_t>> weak_components = algorithms::weakly_connected_components(&graph);
		if (_multiqc || _yaml) {
			std::cout << "num_weakly_connected_components: " << weak_components.size() << std::endl;
			std::cout << "weakly_connected_components: " << std::endl;
		} else {
			std::cout << "##num_weakly_connected_components: " << weak_components.size() << std::endl;
			std::cout << "#component\tnodes\tis_acyclic" << std::endl;
		}
        for(uint64_t i = 0; i < weak_components.size(); ++i) {
            auto& weak_component = weak_components[i];

            ska::flat_hash_set<handlegraph::nid_t> head_nodes = algorithms::is_nice_and_acyclic(graph, weak_components[i]);
            bool acyclic = !(head_nodes.empty());
			if (_multiqc || _yaml) {
				std::cout << "  - component:" << std::endl;
				std::cout << "      id: " << i << std::endl;
				std::cout << "      nodes: " << weak_components[i].size() << std::endl;
				std::cout << "      is_acyclic: " << (acyclic ? "'yes'" : "'no'") << std::endl;
			} else {
				std::cout << i << "\t" << weak_components[i].size() << "\t" << (acyclic ? "yes" : "no") << std::endl;
			}
        }
    }

    if (_num_self_loops || _multiqc) {
        uint64_t total_self_loops = 0;
        std::unordered_set<nid_t> loops;
        graph.for_each_edge([&](const edge_t& e) {
            if (graph.get_id(e.first) == graph.get_id(e.second)) {
                ++total_self_loops;
                loops.insert(graph.get_id(e.first));
            }
        });

        // Should be these always equal?
        if (_multiqc || _yaml) {
        	std::cout << "num_nodes_self_loops:" << std::endl;
        	std::cout << "  total: " << total_self_loops << std::endl;
        	std::cout << "  unique: " << loops.size() << std::endl;
        } else {
			cout << "#type\tnum" << endl;
			cout << "total" << "\t" << total_self_loops << endl;
			cout << "unique" << "\t" << loops.size() << endl;
		}
    }

	/// we don't do this when `-y, --_multiqc` was specified
    if (_show_nondeterministic_edges) {
        // This edges could be compressed in principle

        std::cout << "#from_node\tto_node" << std::endl;
        graph.for_each_handle([&](const handle_t& handle) {
            nid_t id = graph.get_id(handle);
            for (bool is_reverse : { false, true }) {
                ska::flat_hash_map<char, std::vector<handle_t>> edges;
                graph.follow_edges(graph.get_handle(id, is_reverse), false, [&](const handle_t& to) {
                    edges[graph.get_base(to, 0)].push_back(to);
                });
                for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
                    if (iter->second.size() > 1) {
                        for (const handle_t& to : iter->second) {
                            std::cout << id << (is_reverse ? "-" : "+") << "\t" << graph.get_id(to) << (graph.get_is_reverse(to) ? "-" : "+") << std::endl;
                        }
                    }
                }
            }
        });
    }

    if (args::get(base_content) || _multiqc) {
        std::vector<uint64_t> chars(256);
        graph.for_each_handle([&](const handle_t& h) {
                std::string seq = graph.get_sequence(h);
                for (auto c : seq) {
                    ++chars[c];
                }
            });
        for (uint64_t i = 0; i < 256; ++i) {
            if (chars[i]) {
            	if (_multiqc || _yaml) {
					std::cout << (char)i << ": " << chars[i] << std::endl;
            	} else {
					std::cout << (char)i << "\t" << chars[i] << std::endl;
				}
            }
        }
    }

	if (_file_size || _multiqc) {
		// 1. get the file size with error handling
		const filesystem::path path_infile = infile;
		std::error_code err_code;
		std::uintmax_t file_size = filesystem::file_size(path_infile, err_code);
		if (err_code) {
			std::cerr << "[odgi::stats] error: " << infile << " : " << err_code.message() << std::endl;
			exit(1);
		} else {
			if (_multiqc || _yaml) {
				std::cout << "file_size_in_bytes: " << file_size << std::endl;
			} else {
				std::cout << file_size << std::endl;
			}
		}
	}

	if (_pangenome_sequence_class_counts) {
		const char delim = delim_pos[0][0];
		uint64_t sample_pos = std::stoul(delim_pos[1]);
		struct pangenome_seq_class_counts {
			uint64_t core = 0;
			uint64_t priv = 0;
			uint64_t shell = 0;
		};
		ska::flat_hash_map<std::string, pangenome_seq_class_counts> sample_counts;
		/// collect the path names in a map which goes *path_handle* -> *sample_name*
		ska::flat_hash_map<path_handle_t, std::string> p_h_s_n_set;
		ska::flat_hash_set<std::string> sample_names;
		graph.for_each_path_handle([&](const path_handle_t &path) {
			const std::string path_name = graph.get_path_name(path);
			std::vector<string> path_name_split = split(path_name, delim);
			const std::string sample_name = path_name_split[sample_pos];
			p_h_s_n_set[path] = sample_name;
			sample_names.insert(sample_name);
		});
		graph.for_each_handle([&](const handle_t &h) {
			uint64_t hl = graph.get_length(h);
			ska::flat_hash_set<std::string> samples;
			graph.for_each_step_on_handle(h, [&](const step_handle_t &occ) {
				const path_handle_t p_h = graph.get_path_handle_of_step(occ);
				samples.insert(p_h_s_n_set[p_h]);
			});
			/// Private
			if (samples.size() == 1) {
				for (auto& sample_name : samples) {
					sample_counts[sample_name].priv = sample_counts[sample_name].priv + hl;
				}
				/// Shell
			} else if (samples.size() < sample_names.size()) {
				for (auto& sample_name : samples) {
					sample_counts[sample_name].shell = sample_counts[sample_name].shell + hl;
				}
				/// Core
			} else {
				for (auto& sample_name : samples) {
					sample_counts[sample_name].core = sample_counts[sample_name].core + hl;
				}
			}
		});
		if (_yaml || _multiqc) {
			std::cout << "pangenome_sequence_class_counts:" << std::endl;
			for (auto& sample_count : sample_counts) {
				std::cout << "  - sample: " << std::endl;
				std::cout << "      name: " << sample_count.first << std::endl;
				std::cout << "      core: " << sample_count.second.core << std::endl;
				std::cout << "      private: " << sample_count.second.priv << std::endl;
				std::cout << "      shell: " << sample_count.second.shell << std::endl;
			}
		} else {
			for (auto& sample_count : sample_counts) {
				std::cout << sample_count.first << "\t" << sample_count.second.core << "\t" << sample_count.second.priv
				<< "\t" << sample_count.second.shell << std::endl;
			}
		}
		// TODO clear all sets?
	}

    if (args::get(mean_links_length) || args::get(sum_of_path_node_distances)) {
        // This vector is needed for computing the metrics in 1D and for detecting gap-links
        std::vector<uint64_t> position_map(graph.get_node_count() + 1);

        // These vectors are needed for computing the metrics in 2D
        std::vector<double> X, Y;

        if (layout_in_file) {
            const auto& infile = args::get(layout_in_file);
            if (!infile.empty()) {
                algorithms::layout::Layout layout;

                if (infile == "-") {
                    layout.load(std::cin);
                } else {
                    ifstream f(infile.c_str());
                    layout.load(f);
                    f.close();
                }

                X = layout.get_X();
                Y = layout.get_Y();
            }
        }

        {
            uint64_t len = 0;
            graph.for_each_handle([&](const handle_t &h) {
                position_map[number_bool_packing::unpack_number(h) - shift] = len;
                uint64_t hl = graph.get_length(h);
                len += hl;
#ifdef debug_odgi_viz
                std::cerr << "SEGMENT ID: " << graph.get_id(h) << " - " << as_integer(h) << " - index_in_position_map (" << number_bool_packing::unpack_number(h) << ") = " << len << std::endl;
#endif
            });
            position_map[position_map.size() - 1] = len;
        }

        if (args::get(mean_links_length)){
            bool _dont_penalize_gap_links = args::get(dont_penalize_gap_links);

            uint64_t sum_all_node_space = 0;
            uint64_t sum_all_nt_space = 0;
            double sum_all_2D_space = 0.0;
            uint64_t num_all_links = 0;
            uint64_t num_all_gap_links = 0;

            if (_multiqc || _yaml) {
            	std::cout << "mean_links_length:" << std::endl;
            } else {
				std::cout << "#mean_links_length" << std::endl;
				if (layout_in_file) {
					std::cout << "path\tin_2D_space\tnum_links_considered" << std::endl;
				}else{
					std::cout << "path\tin_node_space\tin_nucleotide_space\tnum_links_considered";

					if (dont_penalize_gap_links){
						std::cout << "\tnum_gap_links_not_penalized" << std::endl;
					}else{
						std::cout << std::endl;
					}
				}
			}

            // TODO Could we run this in parallel?
			graph.for_each_path_handle([&](const path_handle_t &path) {
				std::set<uint64_t> ordered_unpacked_numbers_in_path; // ordered set to retain the positions order

                uint64_t sum_node_space = 0;
                uint64_t sum_nt_space = 0;
                double sum_2D_space = 0.0;
                uint64_t num_links = 0;
                uint64_t num_gap_links = 0;

                if (_dont_penalize_gap_links){
                    graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                        ordered_unpacked_numbers_in_path.insert(number_bool_packing::unpack_number(graph.get_handle_of_step(occ)));
                    });
                }

                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);

                    if (graph.has_next_step(occ)){
                        handle_t i = graph.get_handle_of_step(graph.get_next_step(occ));

                        uint64_t unpacked_h = number_bool_packing::unpack_number(h);
                        uint64_t unpacked_i = number_bool_packing::unpack_number(i);

                        // The position map includes the start and end of the node in successive entries.
                        // Edges leave from the end (or start) of one node (depending on whether they are on the
                        // forward or reverse strand) and go to the start (or end) of the other side of the link.
                        uint64_t _info_a = unpacked_h + !number_bool_packing::unpack_bit(h);
                        uint64_t _info_b = unpacked_i + number_bool_packing::unpack_bit(i);

                        if (!dont_penalize_gap_links || next(ordered_unpacked_numbers_in_path.find(unpacked_h)) != ordered_unpacked_numbers_in_path.find(unpacked_i)){
                            if (_info_b < _info_a){
                                _info_a = _info_b;
                                _info_b = unpacked_h + !number_bool_packing::unpack_bit(h);
                            }

                            if (layout_in_file) {
                                // 2D metric
                                double dx = X[2 * (unpacked_h - shift) + number_bool_packing::unpack_bit(h)] - X[2 * (unpacked_i - shift) + number_bool_packing::unpack_bit(i)];
                                double dy = Y[2 * (unpacked_h - shift) + number_bool_packing::unpack_bit(h)] - Y[2 * (unpacked_i - shift) + number_bool_packing::unpack_bit(i)];

                                sum_2D_space += sqrt(dx * dx + dy * dy);
                            }else{
                                // 1D metric (in node space and in nucleotide space)
                                sum_node_space += _info_b - _info_a;
                                sum_nt_space += position_map[_info_b - shift] - position_map[_info_a - shift];

#ifdef debug_odgi_stats
                                std::cerr << _info_b << " - " << _info_a << ": " << position_map[_info_b - shift] - position_map[_info_a - shift] << std::endl;
#endif
                            }
                        }else if (dont_penalize_gap_links){
                            num_gap_links++;
                        }

                        num_links++;
                    }
                });

                /// this could land in the YAML, but we don't force it, because we don't need it for the MultiQC module
                if (args::get(path_statistics)) {
                    double ratio_node_space = 0;
                    double ratio_nt_space = 0;
                    double ratio_2D_space = 0;
                    if (num_links > 0){
                        if (layout_in_file) {
                            ratio_2D_space = sum_2D_space / (double)num_links;
                        } else{
                            ratio_node_space = (double)sum_node_space / (double)num_links;
                            ratio_nt_space = (double)sum_nt_space / (double)num_links;
                        }
                    }
					if (_multiqc || _yaml) {
						std::cout << "  - length:" << std::endl;
						std::cout << "      path: " << graph.get_path_name(path) << std::endl;
						if (layout_in_file) {
							std::cout << "      in_2D_space: " << ratio_2D_space << std::endl;
						} else {
							std::cout << "      in_node_space: " << ratio_node_space << std::endl;
							std::cout << "      in_nucleotide_space: " << ratio_nt_space << std::endl;
						}
						std::cout << "      num_links_considered: " << num_links << std::endl;
						if (dont_penalize_gap_links) {
							std::cout << "      num_gap_links_not_penalized: " << num_gap_links << std::endl;
						}
					} else {
						if (layout_in_file) {
							std::cout << graph.get_path_name(path) << "\t" << ratio_2D_space << "\t" << num_links << std::endl;
						}else{
							std::cout << graph.get_path_name(path) << "\t" << ratio_node_space << "\t" << ratio_nt_space << "\t" << num_links;

							if (dont_penalize_gap_links){
								std::cout << "\t" << num_gap_links << std::endl;
							}else{
								std::cout << std::endl;
							}
						}
					}
                }

                sum_all_node_space += sum_node_space;
                sum_all_nt_space += sum_nt_space;
                sum_all_2D_space += sum_2D_space;
                num_all_links += num_links;
                num_all_gap_links += num_gap_links;
            });

            double ratio_node_space = 0;
            double ratio_nt_space = 0;
            double ratio_2D_space = 0;
            if (num_all_links > 0) {
                if (layout_in_file) {
                    ratio_2D_space = sum_all_2D_space / (double)num_all_links;
                }else{
                    ratio_node_space = (double)sum_all_node_space / (double)num_all_links;
                    ratio_nt_space = (double)sum_all_nt_space / (double)num_all_links;
                }
            }
			if (_multiqc || _yaml) {
				std::cout << "  - length:" << std::endl;
				std::cout << "      path: " << "all_paths" << std::endl;
				if (layout_in_file) {
					std::cout << "      in_2D_space: " << ratio_2D_space << std::endl;
				} else {
					std::cout << "      in_node_space: " << ratio_node_space << std::endl;
					std::cout << "      in_nucleotide_space: " << ratio_nt_space << std::endl;
				}
				std::cout << "      num_links_considered: " << num_all_links << std::endl;
				if (dont_penalize_gap_links || _multiqc) {
					std::cout << "      num_gap_links_not_penalized: " << num_all_gap_links << std::endl;
				}
			} else {
				if (layout_in_file) {
					std::cout << "all_paths\t" << ratio_2D_space << "\t" << num_all_links << std::endl;
				}else{
					std::cout << "all_paths\t" << ratio_node_space << "\t" << ratio_nt_space << "\t" << num_all_links;

					if (_dont_penalize_gap_links){
						std::cout << "\t" << num_all_gap_links << std::endl;
					}else{
						std::cout << std::endl;
					}
				}
			}
        }

        if (args::get(sum_of_path_node_distances)){
            bool _penalize_diff_orientation = args::get(penalize_diff_orientation);

            uint64_t sum_all_path_node_dist_node_space = 0;
            uint64_t sum_all_path_node_dist_nt_space = 0;
            double sum_all_path_node_dist_2D_space = 0.0;
            uint64_t len_all_path_node_space = 0;
            uint64_t len_all_path_nt_space = 0;
            uint64_t num_all_penalties = 0;
            uint64_t num_all_penalties_diff_orientation = 0;

            if (_multiqc || _yaml) {
            	std::cout << "sum_of_path_node_distances:" << std::endl;
            } else {
				std::cout << "#sum_of_path_node_distances" << std::endl;

				if (layout_in_file) {
					std::cout << "path\tin_2D_space_by_nodes\tin_2D_space_by_nucleotides\tnodes\tnucleotides";
				}else{
					std::cout << "path\tin_node_space\tin_nucleotide_space\tnodes\tnucleotides\tnum_penalties";
				}

				if (_penalize_diff_orientation){
					std::cout << "\tnum_penalties_different_orientation" << std::endl;
				}else{
					std::cout << std::endl;
				}
			}

            // TODO Could we run this in parallel?
			handle_t last_handle;
            graph.for_each_path_handle([&](const path_handle_t &path) {
#ifdef debug_odgi_stats
                std::cerr << "path_name: " << graph.get_path_name(path) << std::endl;
#endif
                uint64_t sum_path_node_dist_node_space = 0;
                uint64_t sum_path_node_dist_nt_space = 0;
                double sum_path_node_dist_2D_space = 0.0;
                uint64_t len_path_node_space = 0;
                uint64_t len_path_nt_space = 0;
                uint64_t num_penalties = 0;
                uint64_t num_penalties_diff_orientation = 0;

                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
					last_handle = h;

                    if (graph.has_next_step(occ)){
                        handle_t i = graph.get_handle_of_step(graph.get_next_step(occ));
						last_handle = i;

                        uint64_t unpacked_a = number_bool_packing::unpack_number(h);
                        uint64_t unpacked_b = number_bool_packing::unpack_number(i);

                        double euclidean_distance_2D;

                        if (layout_in_file) {
                            // 2D metric
                            double dx = X[2 * (unpacked_a - shift) + number_bool_packing::unpack_bit(h)] - X[2 * (unpacked_b - shift) + number_bool_packing::unpack_bit(i)];
                            double dy = Y[2 * (unpacked_a - shift) + number_bool_packing::unpack_bit(h)] - Y[2 * (unpacked_b - shift) + number_bool_packing::unpack_bit(i)];

                            euclidean_distance_2D = sqrt(dx * dx + dy * dy);
                            sum_path_node_dist_2D_space += euclidean_distance_2D;
                        }else{
                            uint8_t weight = 1;
                            if (unpacked_b < unpacked_a){
                                unpacked_a = unpacked_b;
                                unpacked_b = number_bool_packing::unpack_number(h);

                                // When a path goes back in terms of pangenomic order, this is punished
                                weight = 3;
                                num_penalties++;
                            }

#ifdef debug_odgi_stats
                            std::cerr << unpacked_b << " - " << unpacked_a << ": " << position_map[unpacked_b - shift] - position_map[unpacked_a - shift] << " * " << weight << std::endl;
#endif

                            sum_path_node_dist_node_space += weight * (unpacked_b - unpacked_a);
                            sum_path_node_dist_nt_space += weight * (position_map[unpacked_b - shift] - position_map[unpacked_a - shift]);
                        }

                        if (_penalize_diff_orientation && (number_bool_packing::unpack_bit(h) != number_bool_packing::unpack_bit(i))){
                            if (layout_in_file) {
                                sum_path_node_dist_2D_space += 2 * euclidean_distance_2D;
                            }else{
                                sum_path_node_dist_node_space += 2 * (unpacked_b - unpacked_a);
                                sum_path_node_dist_nt_space += 2 * (position_map[unpacked_b - shift] - position_map[unpacked_a - shift]);
                            }

                            num_penalties_diff_orientation++;
                        }
                    }

                    len_path_node_space++;
                    len_path_nt_space += graph.get_length(h);
                });

				// add end of path so the best metric equals 1
				sum_path_node_dist_node_space++;
				sum_path_node_dist_nt_space += graph.get_length(last_handle);

				/// this could land in the YAML, but we don't force it, because we don't need it for the MultiQC module
                if (args::get(path_statistics)) {
                	if (_multiqc || _yaml) {
						std::cout << "  - distance:" << std::endl;
						std::cout << "      path: " << graph.get_path_name(path) << std::endl;
						if (layout_in_file) {
							std::cout << "      in_2D_space_by_nodes: " << (double)sum_path_node_dist_2D_space / (double)len_path_node_space << std::endl;
							std::cout << "      in_2D_space_by_nucleotides: " << (double)sum_path_node_dist_2D_space / (double)len_path_nt_space << std::endl;
							std::cout << "      nodes: " << len_path_node_space << std::endl;
							std::cout << "      nucleotides: " << len_path_nt_space << std::endl;
						} else {
							std::cout << "      in_node_space: " << (double)sum_path_node_dist_node_space / (double)len_path_node_space << std::endl;
							std::cout << "      in_nucleotide_space: " << (double)sum_path_node_dist_nt_space / (double)len_path_nt_space << std::endl;
							std::cout << "      nodes: " << len_path_node_space << std::endl;
							std::cout << "      nucleotides: " << len_path_nt_space << std::endl;
							std::cout << "      num_penalties: " << num_penalties << std::endl;
						}
						if (_penalize_diff_orientation || _multiqc) {
							std::cout << "      num_penalties_different_orientation: " << num_penalties_diff_orientation << std::endl;
						}
                	} else {
						if (layout_in_file) {
							std::cout << graph.get_path_name(path) << "\t" << (double)sum_path_node_dist_2D_space / (double)len_path_node_space << "\t" << (double)sum_path_node_dist_2D_space / (double)len_path_nt_space << "\t" << len_path_node_space << "\t" << len_path_nt_space;
						}else{
							std::cout << graph.get_path_name(path) << "\t" << (double)sum_path_node_dist_node_space / (double)len_path_node_space << "\t" << (double)sum_path_node_dist_nt_space / (double)len_path_nt_space << "\t" << len_path_node_space << "\t" << len_path_nt_space  << "\t" << num_penalties;
						}

						if (_penalize_diff_orientation){
							std::cout << "\t" << num_penalties_diff_orientation << std::endl;
						}else{
							std::cout << std::endl;
						}
					}
                }

                sum_all_path_node_dist_node_space += sum_path_node_dist_node_space;
                sum_all_path_node_dist_nt_space += sum_path_node_dist_nt_space;
                sum_all_path_node_dist_2D_space += sum_path_node_dist_2D_space;
                len_all_path_node_space += len_path_node_space;
                len_all_path_nt_space += len_path_nt_space;
                num_all_penalties += num_penalties;
                num_all_penalties_diff_orientation += num_penalties_diff_orientation;
            });

            if (_multiqc || _yaml) {
				std::cout << "  - distance:" << std::endl;
				std::cout << "      path: " << "all_paths" << std::endl;
				if (layout_in_file) {
					std::cout << "      in_2D_space_by_nodes: " << (double)sum_all_path_node_dist_2D_space / (double)len_all_path_node_space << std::endl;
					std::cout << "      in_2D_space_by_nucleotides: " << (double)sum_all_path_node_dist_2D_space / (double)len_all_path_nt_space << std::endl;
					std::cout << "      nodes: " << len_all_path_node_space << std::endl;
					std::cout << "      nucleotides: " << len_all_path_nt_space << std::endl;
				} else {
					std::cout << "      in_node_space: " << (double)sum_all_path_node_dist_node_space / (double)len_all_path_node_space << std::endl;
					std::cout << "      in_nucleotide_space: " << (double)sum_all_path_node_dist_nt_space / (double)len_all_path_nt_space << std::endl;
					std::cout << "      nodes: " << len_all_path_node_space << std::endl;
					std::cout << "      nucleotides: " << len_all_path_nt_space << std::endl;
					std::cout << "      num_penalties: " << num_all_penalties << std::endl;
				}
				if (_penalize_diff_orientation || _multiqc) {
					std::cout << "      num_penalties_different_orientation: " << num_all_penalties_diff_orientation << std::endl;
				}
            } else {
				if (layout_in_file) {
					std::cout << "all_paths\t" << (double)sum_all_path_node_dist_2D_space / (double)len_all_path_node_space << "\t" << (double)sum_all_path_node_dist_2D_space / (double)len_all_path_nt_space << "\t" << len_all_path_node_space << "\t" << len_all_path_nt_space;
				}else{
					std::cout << "all_paths\t" << (double)sum_all_path_node_dist_node_space / (double)len_all_path_node_space << "\t" << (double)sum_all_path_node_dist_nt_space / (double)len_all_path_nt_space << "\t" << len_all_path_node_space << "\t" << len_all_path_nt_space << "\t" << num_all_penalties;
				}

				if (_penalize_diff_orientation){
					std::cout << "\t" << num_all_penalties_diff_orientation << std::endl;
				}else{
					std::cout << std::endl;
				}
			}
        }
    }


    // Put path handles in a vector to work on them in parallel
    std::vector<path_handle_t> paths;
    paths.reserve(graph.get_path_count());

    // Check if we are going to use the vector with the path handles
    if (args::get(weighted_feedback_arc) || args::get(weighted_reversing_join) || args::get(links_length_per_nuc)) {
        graph.for_each_path_handle([&](const path_handle_t path) {
            paths.push_back(path);
        });
    }

    if (args::get(weighted_feedback_arc)) {
		if (_multiqc || _yaml) {
			std::cout << "weighted_feedback_arc: ";
		} else {
			std::cout << "path\tweighted_feedback_arc" << std::endl;
		}

        uint64_t wfa_all_paths = 0;

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (auto& path : paths) {
            uint64_t wfa_current_path = 0;

            graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                handle_t h = graph.get_handle_of_step(occ);

                if (graph.has_next_step(occ)){
                    handle_t i = graph.get_handle_of_step(graph.get_next_step(occ));

                    uint64_t unpacked_a = number_bool_packing::unpack_number(h);
                    uint64_t unpacked_b = number_bool_packing::unpack_number(i);

                    // Check if it is a feedback arc (edge joining out-sides with in-sides such that the outside node does not precede the inside node)
                    if (
                        (!graph.get_is_reverse(h) && !graph.get_is_reverse(i) && unpacked_a >= unpacked_b) ||
                        ( graph.get_is_reverse(h) &&  graph.get_is_reverse(i) && unpacked_a <= unpacked_b)
                    ) {
                        wfa_current_path += 1;
                    }
                }
            });

            if (args::get(path_statistics)) {
#pragma omp critical (cout)
                std::cout << graph.get_path_name(path) << "\t" << wfa_current_path << std::endl;
            }

#pragma omp critical (wfa_all_paths)
            wfa_all_paths += wfa_current_path;
        }
		if (_multiqc || _yaml) {
			std::cout << wfa_all_paths << std::endl;
		} else {
			std::cout << "all_paths" << "\t" << wfa_all_paths << std::endl;
		}
    }

    if (args::get(weighted_reversing_join)) {
		if (_multiqc || _yaml) {
			std::cout << "weighted_reversing_join: ";
		} else {
			std::cout << "path\tweighted_reversing_join" << std::endl;
		}

        uint64_t wrj_all_paths = 0;

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
        for (auto& path : paths) {
            uint64_t wrj_current_path = 0;

            graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                handle_t h = graph.get_handle_of_step(occ);

                if (graph.has_next_step(occ)){
                    handle_t i = graph.get_handle_of_step(graph.get_next_step(occ));

                    uint64_t unpacked_a = number_bool_packing::unpack_number(h);
                    uint64_t unpacked_b = number_bool_packing::unpack_number(i);

                    // Check if it is a reversing arc (edges joining two in- or two out-sides)
                    if (graph.get_is_reverse(h) != graph.get_is_reverse(i)) {
                        wrj_current_path += 1;
                    }
                }
            });

            if (args::get(path_statistics)) {
#pragma omp critical (cout)
                std::cout << graph.get_path_name(path) << "\t" << wrj_current_path << std::endl;
            }

#pragma omp critical (cout)
            wrj_all_paths += wrj_current_path;
        }
		if (_multiqc || _yaml) {
			std::cout << wrj_all_paths << std::endl;
		} else {
			std::cout << "all_paths" << "\t" << wrj_all_paths << std::endl;
		}
    }

	std::atomic<uint64_t> total_links_length = 0;
	std::atomic<uint64_t> total_num_nuc = 0;

	if (args::get(links_length_per_nuc)) {
		if (_multiqc || _yaml) {
			std::cout << "links_length_per_nuc: ";
		} else {
			std::cout << "path\tlinks_length_per_nuc" << std::endl;
		}
		// TODO maybe combine this with other position_maps?
		std::vector<uint64_t> pan_pos(graph.get_node_count());
		uint64_t pos = 0;
		graph.for_each_handle([&](const handle_t &h) {
			pan_pos[number_bool_packing::unpack_number(h) - shift] = pos;
			uint64_t hl = graph.get_length(h);
			pos += hl;
#ifdef debug_odgi_viz
			std::cerr << "SEGMENT ID: " << graph.get_id(h) << " - " << as_integer(h) << " - index_in_position_map (" << number_bool_packing::unpack_number(h) << ") = " << len << std::endl;
#endif
		});

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
		for (auto& path : paths) {
			uint64_t links_len_curr_path = 0;
			uint64_t nuc_curr_path = 0;
			std::set<uint64_t> ordered_unpacked_numbers_in_path; // ordered set to retain the positions order

			graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
				ordered_unpacked_numbers_in_path.insert(number_bool_packing::unpack_number(graph.get_handle_of_step(occ)));
			});

			graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
				handle_t h = graph.get_handle_of_step(occ);
				nuc_curr_path += graph.get_length(h);

				if (graph.has_next_step(occ)) {
					handle_t i = graph.get_handle_of_step(graph.get_next_step(occ));
					uint64_t unpacked_h = number_bool_packing::unpack_number(h);
					uint64_t unpacked_i = number_bool_packing::unpack_number(i);
					uint64_t pos_h = pan_pos[unpacked_h - shift];
					uint64_t pos_i = pan_pos[unpacked_i - shift];
					uint64_t len_h = graph.get_length(h);
					uint64_t len_i = graph.get_length(i);
					uint64_t nid_h = graph.get_id(h);
					uint64_t nid_i = graph.get_id(i);
					bool is_rev_h = graph.get_is_reverse(h);
					bool is_rev_i = graph.get_is_reverse(i);

					/// f means forward oriented step, r means reverse oriented step

					/// handle_h: f, handle_i: f
					if (!is_rev_h && !is_rev_i) {
						/// nid_h <= nid_i
						if (nid_h <= nid_i) {
							/// If we have a gap link, we don't count up
							if (next(ordered_unpacked_numbers_in_path.find(unpacked_h)) != ordered_unpacked_numbers_in_path.find(unpacked_i)) {
								links_len_curr_path += pos_i - (pos_h + len_h);
							}
							/// nid_h > nid_i
						} else {
							links_len_curr_path += pos_h - pos_i + len_h;
						}
						/// handle_h: f, handle_i: r
					} else if (!is_rev_h && is_rev_i) {
						if (nid_h <= nid_i) {
							links_len_curr_path += pos_i + len_i - (pos_h + len_h);
							/// nid_h > nid_i
						} else {
							links_len_curr_path += pos_h - pos_i - len_i + len_h;
						}
						/// handle_h: r, handle_i: f
					} else if (is_rev_h && !is_rev_i) {
						if (nid_h <= nid_i) {
							links_len_curr_path += pos_i - pos_h;
							/// nid_h > nid_i
						} else {
							links_len_curr_path += pos_h - pos_i + len_h + len_i;
						}
						/// handle_h: r, handle_i: r
					} else {
						if (nid_h <= nid_i) {
							links_len_curr_path += pos_i - pos_h + len_h + len_i;
							/// nid_h > nid_i
						} else {
							links_len_curr_path += pos_h - (pos_i + len_i);
						}
					}
				}
			});

			total_links_length += links_len_curr_path;
			total_num_nuc += nuc_curr_path;
		}

		double result_links_length_per_nuc = (double)total_links_length / (double) total_num_nuc;
		if (_multiqc || _yaml) {
			std::cout << result_links_length_per_nuc << std::endl;
		} else {
			std::cout << "all_paths\t" << result_links_length_per_nuc << std::endl;
		}
	}

    bool using_delim = !args::get(path_delim).empty();
    char delim = '\0';
    ska::flat_hash_map<std::string, uint64_t> path_group_ids;
    std::vector<std::string> path_groups;
    if (using_delim) {
        delim = args::get(path_delim).at(0);
        uint64_t i = 0;
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                std::string group_name = split(graph.get_path_name(p), delim)[0];
                auto f = path_group_ids.find(group_name);
                if (f == path_group_ids.end()) {
                    path_group_ids[group_name] = i++;
                    path_groups.push_back(group_name);
                }
            });
    }

    auto get_path_name
        = (using_delim ?
           (std::function<std::string(const uint64_t&)>)
           [&](const uint64_t& id) { return path_groups[id]; }
           :
           (std::function<std::string(const uint64_t&)>)
           [&](const uint64_t& id) { return graph.get_path_name(as_path_handle(id)); });

    auto get_path_id
        = (using_delim ?
           (std::function<uint64_t(const path_handle_t&)>)
           [&](const path_handle_t& p) {
               std::string group_name = split(graph.get_path_name(p), delim)[0];
               return path_group_ids[group_name];
           }
           :
           (std::function<uint64_t(const path_handle_t&)>)
           [&](const path_handle_t& p) {
               return as_integer(p);
           });

    /*
if (args::get(path_setcov)) {
    uint64_t total_length = 0;
    std::map<std::set<uint64_t>, uint64_t> setcov;
    graph.for_each_handle(
        [&](const handle_t& h) {
            std::set<uint64_t> paths_here;
            graph.for_each_step_on_handle(
                h,
                [&](const step_handle_t& occ) {
                    paths_here.insert(get_path_id(graph.get_path(occ)));
                });
            size_t l = graph.get_length(h);
#pragma omp critical (setcov)
            {
                setcov[paths_here] += l;
                total_length += l;
            }
        }, true);
    std::cout << "length\tgraph.frac\tn.paths\tpath.set" << std::endl;
    for (auto& p : setcov) {
        std::cout << p.second << "\t"
                  << (double)p.second/(double)total_length << "\t"
                  << p.first.size() << "\t";
        for (auto& i : p.first) {
            std::cout << get_path_name(i) << ",";
        }
        std::cout << std::endl;
    }
}


if (args::get(path_setcov_count)) {
uint64_t total_length = 0;
std::map<uint64_t, uint64_t> setcov_count;
graph.for_each_handle(
    [&](const handle_t& h) {
        std::set<uint64_t> paths_here;
        graph.for_each_step_on_handle(
            h,
            [&](const step_handle_t& occ) {
                paths_here.insert(get_path_id(graph.get_path(occ)));
            });
        size_t l = graph.get_length(h);
#pragma omp critical (setcov)
        {
            setcov_count[paths_here.size()] += l;
            total_length += l;
        }
    }, true);
std::cout << "length\tgraph.frac\tunique.path.count" << std::endl;
for (auto& p : setcov_count) {
    std::cout << p.second << "\t"
              << (double)p.second/(double)total_length
              << "\t" << p.first << std::endl;
}
}

if (args::get(path_multicov)) {
uint64_t total_length = 0;
std::map<std::vector<uint64_t>, uint64_t> multisetcov;
graph.for_each_handle(
    [&](const handle_t& h) {
        std::vector<uint64_t> paths_here;
        graph.for_each_step_on_handle(
            h,
            [&](const step_handle_t& occ) {
                paths_here.push_back(get_path_id(graph.get_path(occ)));
            });
        std::sort(paths_here.begin(), paths_here.end());
        size_t l = graph.get_length(h);
#pragma omp critical (multisetcov)
        {
            multisetcov[paths_here] += l;
            total_length += l;
        }
    }, true);
std::cout << "length\tgraph.frac\tn.paths\tpath.multiset" << std::endl;
for (auto& p : multisetcov) {
    std::cout << p.second << "\t"
              << (double)p.second/(double)total_length << "\t"
              << p.first.size() << "\t";
    bool first = true;
    for (auto& i : p.first) {
        std::cout << (first ? (first=false, "") :",") << get_path_name(i);
    }
    std::cout << std::endl;
}
}

if (args::get(path_multicov_count)) {
uint64_t total_length = 0;
std::map<uint64_t, uint64_t> multisetcov_count;
graph.for_each_handle(
    [&](const handle_t& h) {
        uint64_t paths_here = 0;
        graph.for_each_step_on_handle(
            h,
            [&](const step_handle_t& occ) {
                ++paths_here;
            });
        size_t l = graph.get_length(h);
#pragma omp critical (multisetcov)
        {
            multisetcov_count[paths_here] += l;
            total_length += l;
        }
    }, true);
std::cout << "length\tgraph.frac\tpath.step.count" << std::endl;
for (auto& p : multisetcov_count) {
    std::cout << p.second << "\t"
              << (double)p.second/(double)total_length << "\t"
              << p.first << std::endl;
}
}


    if (!args::get(path_bedmulticov).empty()) {
        std::string line;
        typedef std::map<std::vector<uint64_t>, uint64_t> setcov_t;
        typedef IntervalTree<uint64_t, std::pair<std::string, setcov_t*> > itree_t;
        map<std::string, itree_t::interval_vector> intervals;
        auto& x = args::get(path_bedmulticov);
        std::ifstream bed_in(x);
        while (std::getline(bed_in, line)) {
            // BED is base-numbered, 0-origin, half-open.  This parse turns that
            // into base-numbered, 0-origin, fully-closed for internal use.  All
            // coordinates used internally should be in the latter, and coordinates
            // from the user in the former should be converted immediately to the
            // internal format.
            std::vector<string> fields = split(line, '\t');
            intervals[fields[0]].push_back(
                itree_t::interval(
                    std::stoul(fields[1]),
                    std::stoul(fields[2]),
                    make_pair(fields[3], new setcov_t())));
        }

        std::vector<std::string> path_names;
        path_names.reserve(intervals.size());
        for(auto const& i : intervals) {
            path_names.push_back(i.first);
        }

        // the header
        std::cout << "path.name" << "\t"
                  << "bed.name" << "\t"
                  << "bed.start" << "\t"
                  << "bed.stop" << "\t"
                  << "bed.len" << "\t"
                  << "path.set.state.len" << "\t"
                  << "path.set.frac" << "\t"
                  << "path.traversals" << "\t"
                  << "uniq.paths.in.state" << "\t"
                  << "path.multiset" << std::endl;

        #pragma omp parallel for
        for (uint64_t k = 0; k < path_names.size(); ++k) {
            auto& path_name = path_names.at(k);
            auto& path_ivals = intervals[path_name];
            // build the intervals for each path we'll query
            itree_t itree(std::move(path_ivals)); //, 16, 1);
            uint64_t pos = 0;
            graph.for_each_step_in_path(graph.get_path_handle(path_name), [&](const step_handle_t& occ) {
                    std::vector<uint64_t> paths_here;
                    handle_t h = graph.get_handle_of_step(occ);
                    graph.for_each_step_on_handle(
                        h,
                        [&](const step_handle_t& occ) {
                            paths_here.push_back(get_path_id(graph.get_path(occ)));
                        });
                    uint64_t len = graph.get_length(h);
                    // check each position in the node
                    auto hits = itree.findOverlapping(pos, pos+len);
                    if (hits.size()) {
                        for (auto& h : hits) {
                            auto& q = *h.value.second;
                            // adjust length for overlap length
                            uint64_t ovlp = len - (h.start > pos ? h.start - pos : 0) - (h.stop < pos+len ? pos+len - h.stop : 0);
                            q[paths_here] += ovlp;
                        }
                    }
                    pos += len;
                });
            itree.visit_all([&](const itree_t::interval& ival) {
                    auto& name = ival.value.first;
                    auto& setcov = *ival.value.second;
                    for (auto& p : setcov) {
                        std::set<uint64_t> u(p.first.begin(), p.first.end());
        #pragma omp critical (cout)
                        {
                            std::cout << path_name << "\t" << name << "\t" << ival.start << "\t" << ival.stop << "\t"
                                      << ival.stop - ival.start << "\t"
                                      << p.second << "\t"
                                      << (float)p.second/(ival.stop-ival.start) << "\t"
                                      << p.first.size() << "\t"
                                      << u.size() << "\t";
                            bool first = true;
                            for (auto& i : p.first) {
                                std::cout << (first ? (first=false, "") :",") << get_path_name(i);
                            }
                            std::cout << std::endl;
                        }
                    }
            });
            for (auto& ival : path_ivals) {
                delete ival.value.second;
            }
        }

        //todo there is still a memory leak somewhere
    }
*/
    return 0;
}

static Subcommand odgi_stats("stats", "Metrics describing a variation graph and its path relationship.",
                              PIPELINE, 3, main_stats);


}
