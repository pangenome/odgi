#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "split.hpp"
#include "subgraph/region.hpp"
#include <omp.h>
#include "algorithms/degree.hpp"

#include "src/algorithms/subgraph/extract.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_degree(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi degree";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Describe the graph in terms of node degree.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> og_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "input"});
    args::Group degree_opts(parser, "[ Degree Options ]");
	args::ValueFlag<std::string> _subset_paths(degree_opts, "FILE",
											   "Compute the degree considering only the paths specified in the FILE. "
											   "The file must contain one path name per line and a subset of all paths can be specified; "
											   "If a step is of a path of the given list, it is taken into account when calculating a node's degree. Else not.",
											   {'s', "subset-paths"});

	args::ValueFlag<std::string> path_name(degree_opts, "PATH_NAME", "Compute the degree of the given path PATH_NAME in the graph.",
										   {'r', "path"});
	args::ValueFlag<std::string> path_file(degree_opts, "FILE", "Report the degree only for the paths listed in FILE.",
										   {'R', "paths"});
	args::ValueFlag<std::string> graph_pos(degree_opts, "[node_id][,offset[,(+|-)]*]*",
										   "Compute the degree at the given node, e.g. 7 or 3,4 or 42,10,+ or 302,0,-.",
										   {'g', "graph-pos"});
	args::ValueFlag<std::string> graph_pos_file(degree_opts, "FILE", "A file with one graph position per line.",
												{'G', "graph-pos-file"});
	args::ValueFlag<std::string> path_pos(degree_opts, "[path_name][,offset[,(+|-)]*]*",
										  "Return degree at the given path position e.g. chrQ or chr3,42 or chr8,1337,+ or chrZ,3929,-.",
										  {'p', "path-pos"});
	args::ValueFlag<std::string> path_pos_file(degree_opts, "FILE", "A file with one path position per line.",
											   {'F', "path-pos-file"});
	args::ValueFlag<std::string> bed_input(degree_opts, "FILE", "A BED file of ranges in paths in the graph.",
										   {'b', "bed-input"});

	args::Flag graph_degree_table(degree_opts, "graph-degree-table",
								  "Compute the degree on each node in the graph, writing a table by node: node.id, degree.",
								  {'d', "graph-degree-table"});

	args::Flag graph_degree_vec(degree_opts, "graph-degree-vec",
								"Compute the degree on each node in the graph, writing a vector by base in one line.",
								{'v', "graph-degree-vec"});

	args::Flag path_degree(degree_opts, "path-degree",
						   "Compute a vector of degree on each base of each path. Each line consists of a path name and subsequently the space-separated degree of each base.",
						   {'D', "path-degree"});

	args::Flag self_degree(degree_opts, "self-degree",
						   "Compute the degree of the path versus itself on each base in each path. Each line consists of a path name and subsequently the space-separated degree of each base.",
						   {'a', "self-degree"});
    args::Flag summarize_degree(degree_opts, "summarize", "Summarize the graph properties and dimensions. Print in tab-delimited format to stdout the node.count, edge.count, avg.degree, min.degree, max.degree.", {'S', "summarize-graph-degree"});
    args::ValueFlag<std::string> _windows_in(degree_opts, "LEN:MIN:MAX",
                                             "Print to stdout a BED file of path intervals where the degree is between MIN and MAX, "
                                             "merging regions not separated by more than LEN bp.",
											 {'w', "windows-in"});
    args::ValueFlag<std::string> _windows_out(degree_opts, "LEN:MIN:MAX",
                                              "Print to stdout a BED file of path intervals where the degree is outside of MIN and MAX, "
                                              "merging regions not separated by more than LEN bp.",
											  {'W', "windows-out"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> _num_threads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi degree.", {'h', "help"});
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

    if (!og_file) {
        std::cerr << "[odgi::degree] error: please specify a target graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    if (_windows_in && _windows_out) {
        std::cerr << "[odgi::degree] error: please specify -w/--windows-in or -W/--windows-out, not both." << std::endl;
        return 1;
    }

    if (summarize_degree && (_windows_in || _windows_out)) {
        std::cerr << "[odgi::degree] error: please specify -S/--summarize without specifying windows-in or -W/--windows-out." << std::endl;
        return 1;
    }

    uint64_t windows_in_len = 0, windows_in_min = 0, windows_in_max = 0;
    bool windows_in_only_tips = false;
    if (_windows_in) {
        if (!algorithms::check_and_get_windows_in_out_parameter(args::get(_windows_in)+":0", windows_in_len, windows_in_min, windows_in_max, windows_in_only_tips)) {
            std::cerr << "[odgi::degree] error: please specify a valid string (LEN:MIN:MAX) for the -w/--windows-in option." << std::endl;
            return 1;
        }
    }

    uint64_t windows_out_len = 0, windows_out_min = 0, windows_out_max = 0;
    bool windows_out_only_tips = false;
    if (_windows_out) {
        if (!algorithms::check_and_get_windows_in_out_parameter(args::get(_windows_out)+":0", windows_out_len, windows_out_min, windows_out_max, windows_out_only_tips)) {
            std::cerr << "[odgi::degree] error: please specify a valid string (LEN:MIN:MAX) for the -W/--windows-out option." << std::endl;
            return 1;
        }
    }

	const uint64_t num_threads = args::get(_num_threads) ? args::get(_num_threads) : 1;

	odgi::graph_t graph;
    assert(argc > 0);
    if (!args::get(og_file).empty()) {
        std::string infile = args::get(og_file);
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "degree", args::get(progress), num_threads, graph);
        }
    }

    omp_set_num_threads((int) num_threads);
	const uint64_t shift = graph.min_node_id();
	if (_windows_in || _windows_out) {
		if (graph.max_node_id() - shift >= graph.get_node_count()){
			std::cerr << "[odgi::degree] error: the node IDs are not compacted. Please run 'odgi sort' using -O, --optimize to optimize the graph." << std::endl;
			exit(1);
		}
	}

	std::vector<bool> paths_to_consider;
	if (_subset_paths) {
		paths_to_consider.resize(graph.get_path_count() + 1, false);

		std::ifstream refs(args::get(_subset_paths).c_str());
		std::string line;
		while (std::getline(refs, line)) {
			if (!line.empty()) {
				if (!graph.has_path(line)) {
					std::cerr << "[odgi::degree] error: path " << line << " not found in graph" << std::endl;
					exit(1);
				}

				paths_to_consider[as_integer(graph.get_path_handle(line))] = true;
			}
		}
	} else {
		paths_to_consider.resize(graph.get_path_count() + 1, true);
	}

	// these options are exclusive (probably we should say with a warning)
	std::vector<odgi::pos_t> graph_positions;
	std::vector<odgi::path_pos_t> path_positions;
	std::vector<odgi::path_range_t> path_ranges;

	// TODO refactor this into a file, we are using this in odgi depth, odgi position, ....
	auto add_graph_pos = [&graph_positions](const odgi::graph_t &graph,
											const std::string &buffer) {
		auto vals = split(buffer, ',');
		uint64_t id = std::stoi(vals[0]);
		if (!graph.has_node(id)) {
			std::cerr << "[odgi::degree] error: no node " << id << " in graph" << std::endl;
			exit(1);
		}
		uint64_t offset = 0;
		if (vals.size() >= 2) {
			offset = std::stoi(vals[1]);
			handle_t h = graph.get_handle(id);
			if (graph.get_length(h) < offset) {
				std::cerr << "[odgi::degree] error: offset of " << offset << " lies beyond the end of node " << id
						  << std::endl;
				exit(1);
			}
		}
		bool is_rev = false;
		if (vals.size() == 3) {
			is_rev = vals[2] == "-";
		}
		graph_positions.push_back(make_pos_t(id, is_rev, offset));
	};

	// TODO refactor this into a file, we are using this in odgi depth, odgi position, ....
	auto add_path_pos = [&path_positions](const odgi::graph_t &graph,
										  const std::string &buffer) {
		if (!buffer.empty()) {
			auto vals = split(buffer, ',');
			auto &path_name = vals[0];
			if (!graph.has_path(path_name)) {
				std::cerr << "[odgi::degree] error: path " << path_name << " not found in graph" << std::endl;
				exit(1);
			} else {
				path_positions.push_back({
												 graph.get_path_handle(path_name),
												 (vals.size() > 1 ? (uint64_t) std::stoi(vals[1]) : 0),
												 (vals.size() == 3 ? vals[2] == "-" : false)
										 });
			}
		}
	};



	if (summarize_degree) {
		// we do nothing here, we iterate over the handles in the graph later, spitting out the information by default
	} else if (graph_degree_table) {
		graph.for_each_handle([&](const handle_t &h) {
			add_graph_pos(graph, std::to_string(graph.get_id(h)));
		});
	} else if (graph_degree_vec) {
		std::cout << (og_file ? args::get(og_file) : "graph") << "_vec";
		graph.for_each_handle(
				[&](const handle_t &h) {
					uint64_t degree = 0;
					bool consider = false;
					graph.for_each_step_on_handle(
							h,
							[&](const step_handle_t &occ) {
								if (paths_to_consider[
										as_integer(
												graph.get_path_handle_of_step(occ))
								]) {
									consider = true;
									return;
								}
							});
					if (consider) {
						degree += graph.get_degree(h, false) + graph.get_degree(h, true);
					}
					auto length = graph.get_length(h);
					for (uint64_t i = 0; i < length; ++i) {
						std::cout << " " << degree;
					}
				});
		std::cout << std::endl;
	} else if (path_degree) {
		std::vector<path_handle_t> paths;
		graph.for_each_path_handle(
				[&paths,&paths_to_consider](const path_handle_t& path) {
					if (paths_to_consider[as_integer(path)]) {
						paths.push_back(path);
					}
				});
		// for each path handle
#pragma omp parallel for schedule(dynamic, 1)
		for (auto& path : paths) {
			std::stringstream ss;
			ss << graph.get_path_name(path);
			// for each step
			uint64_t pos = 0;
			graph.for_each_step_in_path(
					path,
					[&](const step_handle_t& step) {
						handle_t handle = graph.get_handle_of_step(step);
						auto degree = graph.get_degree(handle, false) + graph.get_degree(handle, true);
						auto next_pos = pos + graph.get_length(handle);
						while (pos++ < next_pos) {
							ss << " " << degree;
						}
					});
#pragma omp critical (cout)
			std::cout << ss.str() << std::endl;
		}
	} else if (self_degree) {
		std::vector<path_handle_t> paths;
		graph.for_each_path_handle(
				[&paths,&paths_to_consider](const path_handle_t& path) {
					if (paths_to_consider[as_integer(path)]) {
						paths.push_back(path);
					}
				});
		// for each path handle
#pragma omp parallel for schedule(dynamic, 1)
		for (auto& path : paths) {
			std::stringstream ss;
			ss << graph.get_path_name(path);
			// for each step
			uint64_t pos = 0;
			graph.for_each_step_in_path(
					path,
					[&](const step_handle_t& step) {
						handle_t handle = graph.get_handle_of_step(step);
						uint64_t degree = 0;
						graph.for_each_step_on_handle(
								handle,
								[&](const step_handle_t& other) {
									if (path == graph.get_path_handle_of_step(other)) {
										degree += graph.get_degree(handle, false) + graph.get_degree(handle, true);
									}
								});
						auto next_pos = pos + graph.get_length(handle);
						while (pos++ < next_pos) {
							ss << " " << degree;
						}
					});
#pragma omp critical (cout)
			std::cout << ss.str() << std::endl;
		}
	} else if (graph_pos) {
		add_graph_pos(graph, args::get(graph_pos));
	} else if (graph_pos_file) {
		std::ifstream gpos(args::get(graph_pos_file));
		std::string buffer;
		while (std::getline(gpos, buffer)) {
			add_graph_pos(graph, buffer);
		}
	} else if (path_pos) {
		// if given a path pos, we convert it into a path pos in our reference set
		add_path_pos(graph, args::get(path_pos));
	} else if (path_pos_file) {
		// if we're given a file of path degrees, we'll convert them all
		std::ifstream refs(args::get(path_pos_file));
		std::string buffer;
		while (std::getline(refs, buffer)) {
			add_path_pos(graph, buffer);
		}
	} else if (bed_input) {
		std::ifstream bed_in(args::get(bed_input));
		std::string buffer;
		while (std::getline(bed_in, buffer)) {
		    add_bed_range(path_ranges, graph, buffer);
		}
	} else if (path_name) {
	    add_bed_range(path_ranges, graph, args::get(path_name));
	} else if (path_file) {
		// for thing in things
		std::ifstream refs(args::get(path_file));
		std::string line;
		while (std::getline(refs, line)) {
		    add_bed_range(path_ranges, graph, line);
		}
	} else if (!_windows_in && !_windows_out) {
		// using all the paths in the graph
		graph.for_each_path_handle(
		        [&](const path_handle_t &path) { add_bed_range(path_ranges, graph, graph.get_path_name(path)); });
	}

	auto get_graph_pos = [](const odgi::graph_t &graph,
							const path_pos_t &pos) {
		const auto path_end = graph.path_end(pos.path);
		uint64_t walked = 0;
		for (step_handle_t s = graph.path_begin(pos.path);
			 s != path_end; s = graph.get_next_step(s)) {
			handle_t h = graph.get_handle_of_step(s);
			uint64_t node_length = graph.get_length(h);
			if (walked + node_length > pos.offset) {
				return make_pos_t(graph.get_id(h), graph.get_is_reverse(h), pos.offset - walked);
			}
			walked += node_length;
		}

#pragma omp critical (cout)
		std::cerr << "[odgi::degree] warning: position " << graph.get_path_name(pos.path) << ":" << pos.offset
				  << " outside of path" << std::endl;
		return make_pos_t(0, false, 0);
	};

	auto get_offset_in_path = [](const odgi::graph_t &graph,
								 const path_handle_t &path, const step_handle_t &target) {
		const auto path_end = graph.path_end(path);
		uint64_t walked = 0;
		step_handle_t s = graph.path_begin(path);
		for (; s != target; s = graph.get_next_step(s)) {
			handle_t h = graph.get_handle_of_step(s);
			walked += graph.get_length(h);
		}
		assert(s != path_end);
		return walked;
	};

	auto get_graph_node_degree = [](const odgi::graph_t &graph, const nid_t node_id,
									const std::vector<bool>& paths_to_consider) {

		uint64_t node_degree = 0;
		std::set<uint64_t> unique_paths;

		const handle_t h = graph.get_handle(node_id);
		bool consider = false;

		graph.for_each_step_on_handle(
				h,
				[&](const step_handle_t &occ) {
					if (paths_to_consider[
							as_integer(graph.get_path_handle_of_step(occ))]) {
						consider = true;
						unique_paths.insert(as_integer(graph.get_path(occ)));
						return;
					}
				});
		if (consider) {
			node_degree +=graph.get_degree(h, false) + graph.get_degree(h, true);
		}

		return make_pair(node_degree, unique_paths.size());
	};

	if (_windows_in || _windows_out) {
		std::vector<path_handle_t> paths;
		paths.reserve(graph.get_path_count());
		graph.for_each_path_handle([&](const path_handle_t path) {
			paths.push_back(path);
		});

		// precompute degrees for all handles in parallel
		std::vector<uint64_t> degrees(graph.get_node_count() + 1);
		graph.for_each_handle(
				[&](const handle_t& h) {
					auto id = graph.get_id(h);
					degrees[id - shift] = graph.get_degree(h, false) + graph.get_degree(h, true);
				}, true);

		auto in_bounds =
				[&](const handle_t &handle) {
					uint64_t degree = degrees[graph.get_id(handle) - shift];
					return _windows_in ? (degree >= windows_in_min && degree <= windows_in_max) : (degree < windows_out_min || degree > windows_out_max);
				};

		std::cout << "#path\tstart\tend" << std::endl;

		algorithms::windows_in_out(graph, paths, in_bounds, _windows_in ? windows_in_len : windows_out_len,
								   [&](const std::vector<path_range_t>& path_ranges) {
#pragma omp critical (cout)
									   for (auto path_range : path_ranges) {
										   std::cout << graph.get_path_name(path_range.begin.path) << "\t"
													 << path_range.begin.offset << "\t"
													 << path_range.end.offset << std::endl;
									   }
								   }, num_threads);
	}

// TODO modify from here the results
// TODO summarize_degree needs a major overhaul!
    if (summarize_degree) {
        uint64_t total_edges = 0;
        uint64_t min_degree = std::numeric_limits<uint64_t>::max();
        uint64_t max_degree = std::numeric_limits<uint64_t>::min();
        graph.for_each_handle(
                [&](const handle_t& handle) {
                    uint64_t degree = graph.get_degree(handle, false) + graph.get_degree(handle, true);
                    total_edges += degree;
                    min_degree = std::min(min_degree, degree);
                    max_degree = std::max(max_degree, degree);
                });
        std::cout << "#node.count\tedge.count\tavg.degree\tmin.degree\tmax.degree" << std::endl
                  << graph.get_node_count() << "\t"
                  << total_edges / 2 << "\t" // we double-count edges
                  << (double) total_edges / (double)graph.get_node_count() << "\t"
                  << min_degree << "\t"
                  << max_degree
                  << std::endl;
    }

	if (!graph_positions.empty()) {
		std::cout << "#node.id\tnode.degree" << std::endl;
		for (auto& graph_pos : graph_positions) {
#pragma omp critical (cout)
			std::cout << id(graph_pos) << "\t"
					  << graph.get_degree(graph.get_handle(id(graph_pos)), false) +
					  graph.get_degree(graph.get_handle(id(graph_pos)), true)
					  << std::endl;
		}
	}

	if (!path_positions.empty()) {
		std::cout << "#path.position\tdegree\tdegree.uniq" << std::endl;
#pragma omp parallel for schedule(dynamic, 1)
		for (auto &path_pos : path_positions) {
			const pos_t pos = get_graph_pos(graph, path_pos);

			const nid_t node_id = id(pos);
			const auto degree = get_graph_node_degree(graph, node_id, paths_to_consider);

#pragma omp critical (cout)
			std::cout << (graph.get_path_name(path_pos.path)) << "," << path_pos.offset << ","
					  << (path_pos.is_rev ? "-" : "+") << "\t"
					  << degree.first << "\t" << degree.second << std::endl;
		}
	}

	if (!path_ranges.empty()) {
		std::cout << "#path\tstart\tend\tmean.degree" << std::endl;
		algorithms::for_each_path_range_degree(
				graph,
				path_ranges,
				paths_to_consider,
				[&](const path_range_t& range,
					const double& degree) {
#pragma omp critical (cout)
					std::cout << (graph.get_path_name(range.begin.path)) << "\t"
							  << range.begin.offset << "\t"
							  << range.end.offset << "\t"
							  << degree << std::endl;
				});
	}

    return 0;
}

static Subcommand odgi_degree("degree", "Describe the graph in terms of node degree.",
                              PIPELINE, 3, main_degree);


}
