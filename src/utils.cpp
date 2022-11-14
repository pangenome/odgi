#include <string>
#include <algorithm>
#include "utils.hpp"

namespace utils {
    bool is_number(const std::string &s) {
        return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
    }

    void graph_deep_copy(const odgi::graph_t& source,
                         odgi::graph_t* target) {
        // copy the now-compacted graph to our output_graph
        source.for_each_handle(
                [&](const handle_t& old_handle) {
                    target->create_handle(
                            source.get_sequence(old_handle),
                            source.get_id(old_handle));
                });

        source.for_each_handle(
                [&](const handle_t& curr) {
                    source.follow_edges(
                            curr, false,
                            [&](const handle_t& next) {
                                target->create_edge(
                                        target->get_handle(source.get_id(curr),
                                                           source.get_is_reverse(curr)),
                                        target->get_handle(source.get_id(next),
                                                           source.get_is_reverse(next)));
                            });
                    source.follow_edges(
                            curr, true,
                            [&](const handle_t& prev) {
                                target->create_edge(
                                        target->get_handle(source.get_id(prev),
                                                           source.get_is_reverse(prev)),
                                        target->get_handle(source.get_id(curr),
                                                           source.get_is_reverse(curr)));
                            });
                });

        source.for_each_path_handle(
                [&](const path_handle_t& old_path) {
                    path_handle_t new_path = target->create_path_handle(source.get_path_name(old_path));
                    source.for_each_step_in_path(old_path, [&](const step_handle_t& step) {
                        handle_t old_handle = source.get_handle_of_step(step);
                        handle_t new_handle = target->get_handle(
                                source.get_id(old_handle),
                                source.get_is_reverse(old_handle));
                        target->append_step(new_path, new_handle);
                    });
                });
    }

	bool ends_with(const std::string &fullString, const std::string &ending) {
		if (fullString.length() >= ending.length()) {
			return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
		} else {
			return false;
		}
	}

	int handle_gfa_odgi_input(const std::string infile, const std::string subcommmand_name, const bool progress,
							const uint64_t num_threads, odgi::graph_t &graph) {
		if (!std::filesystem::exists(infile)) {
			std::cerr << "[odgi::" << subcommmand_name << "] error: the given file \"" << infile << "\" does not exist. Please specify an existing input file in ODGI format via -i=[FILE], --idx=[FILE]." << std::endl;
			exit(1);
		}
		if (utils::ends_with(infile, "gfa")) {
			if (progress) {
				std::cerr << "[odgi::" << subcommmand_name << "] warning: the given file \"" << infile << "\" is not in ODGI format. "
																				   "To save time in the future, please use odgi build -i=[FILE], --idx=[FILE] -o=[FILE], --out=[FILE] "
																				   "to generate a graph in ODGI format. Such a graph can be supplied to all ODGI subcommands. Building graph in ODGI format from given GFA." << std::endl;
			}
			gfa_to_handle(infile, &graph, false, num_threads, progress);
			graph.set_number_of_threads(num_threads);
		} else {
			ifstream f(infile.c_str());
			graph.deserialize(f);
			f.close();
		}
		return 0;
    }

	uint64_t modulo(const uint64_t n, const uint64_t d) {
		return (n & (d - 1));
	}

	std::vector<path_handle_t> load_paths(const std::string& path_names_file,
										  odgi::graph_t &graph,
										  const std::string subcommmand_name) {
		std::ifstream path_names_in(path_names_file);
		uint64_t num_of_paths_in_file = 0;
		std::vector<bool> path_already_seen;
		path_already_seen.resize(graph.get_path_count(), false);
		std::string line;
		std::vector<path_handle_t> paths;
		while (std::getline(path_names_in, line)) {
			if (!line.empty()) {
				if (graph.has_path(line)) {
					const path_handle_t path = graph.get_path_handle(line);
					const uint64_t path_rank = as_integer(path) - 1;
					if (!path_already_seen[path_rank]) {
						path_already_seen[path_rank] = true;
						paths.push_back(path);
					} else {
						std::cerr << "[odgi::" << subcommmand_name << "] error: in the path list there are duplicated path names."
								  << std::endl;
						exit(1);
					}
				}
				++num_of_paths_in_file;
			}
		}
		path_names_in.close();
		std::cerr << "[odgi::" << subcommmand_name << "] found " << paths.size() << "/" << num_of_paths_in_file
				  << " paths to consider." << std::endl;
		if (paths.empty()) {
			std::cerr << "[odgi::" << subcommmand_name << "] error: no path to consider." << std::endl;
			exit(1);
		}
		return paths;
	}
}
