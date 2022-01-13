/**
 * \file
 * unittest/stepindex.cpp: test cases for the implementations of the step index.
 */

#include "catch.hpp"

#include <handlegraph/util.hpp>
#include "odgi.hpp"
#include "stepindex.hpp"

#include "algorithms/xp.hpp"

namespace odgi {
	namespace unittest {

		using namespace std;
		using namespace handlegraph;
		using namespace algorithms;

		TEST_CASE("step index construction, and position retrieval.", "[stepindex]") {

			graph_t graph;
			handle_t n1 = graph.create_handle("CAA");
			handle_t n2 = graph.create_handle("A");
			handle_t n3 = graph.create_handle("G");
			handle_t n4 = graph.create_handle("T");
			handle_t n5 = graph.create_handle("C");
			handle_t n6 = graph.create_handle("TTG");
			handle_t n7 = graph.create_handle("A");
			handle_t n8 = graph.create_handle("GAT");
			handle_t n9 = graph.create_handle("GTC");
			handle_t n10 = graph.create_handle("CATG");
			graph.create_edge(n1, n2);
			graph.create_edge(n1, n3);
			graph.create_edge(n2, n4);
			graph.create_edge(n2, n5);
			graph.create_edge(n3, n5);
			graph.create_edge(n4, n6);
			graph.create_edge(n5, n6);
			graph.create_edge(n6, n6);
			graph.create_edge(n6, n7);
			graph.create_edge(n6, n8);
			graph.create_edge(n8, n9);

			graph.create_path_handle("target", false);
			path_handle_t target = graph.get_path_handle("target");
			graph.append_step(target, n3);
			graph.append_step(target, n5);
			graph.append_step(target, n6);
			graph.append_step(target, n6);
			graph.append_step(target, n8);
			graph.append_step(target, n9);

			graph.create_path_handle("query1", false);
			path_handle_t query1 = graph.get_path_handle("query1");
			graph.append_step(query1, n2);
			graph.append_step(query1, n4);
			graph.append_step(query1, n6);
			graph.append_step(query1, n7);

			graph.create_path_handle("query2", false);
			path_handle_t query2 = graph.get_path_handle("query2");
			graph.append_step(query2, n10);

			graph.create_path_handle("query3", false);
			path_handle_t query3 = graph.get_path_handle("query3");
			graph.append_step(query3, n1);
			graph.append_step(query3, n2);
			graph.append_step(query3, n5);
			graph.append_step(query3, n6);
			graph.append_step(query3, n8);
			graph.append_step(query3, n9);

			vector<path_handle_t> paths;
			graph.for_each_path_handle([&](const path_handle_t path) {
				paths.push_back(path);
			});

			SECTION("The index delivers the correct positions for a given step. Sample rate: 1.") {
				step_index_t step_index_1(graph, paths, 1, false, 1);
				graph.for_each_path_handle([&](const path_handle_t path) {
					std::string cur_path = graph.get_path_name(path);
					if (cur_path == "target") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_1.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_1.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_1.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_1.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_1.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_1.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
						step_handle_t path_end = graph.path_end(path);
						step_handle_t path_back_from_end = graph.get_previous_step(path_end);
						REQUIRE(step_index_1.get_position(path_end, graph) == step_index_1.get_path_len(path));
						REQUIRE(step_index_1.get_position(path_back_from_end, graph) == step_index_1.get_position(graph.path_back(path), graph));
					}

					if (cur_path == "query1") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_1.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_1.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_1.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_1.get_position(occ, graph) == 5);
									break;
							}
							cur_step_rank++;
						});
						step_handle_t path_end = graph.path_end(path);
						step_handle_t path_back_from_end = graph.get_previous_step(path_end);
						REQUIRE(step_index_1.get_position(path_end, graph) == step_index_1.get_path_len(path));
						REQUIRE(step_index_1.get_position(path_back_from_end, graph) == step_index_1.get_position(graph.path_back(path), graph));
					}

					if (cur_path == "query2") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_1.get_position(occ, graph) == 0);
									break;
							}
							cur_step_rank++;
						});
						step_handle_t path_end = graph.path_end(path);
						step_handle_t path_back_from_end = graph.get_previous_step(path_end);
						REQUIRE(step_index_1.get_position(path_end, graph) == step_index_1.get_path_len(path));
						REQUIRE(step_index_1.get_position(path_back_from_end, graph) == step_index_1.get_position(graph.path_back(path), graph));
					}

					if (cur_path == "query3") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_1.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_1.get_position(occ, graph) == 3);
									break;
								case 2:
									REQUIRE(step_index_1.get_position(occ, graph) == 4);
									break;
								case 3:
									REQUIRE(step_index_1.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_1.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_1.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
						step_handle_t path_end = graph.path_end(path);
						step_handle_t path_back_from_end = graph.get_previous_step(path_end);
						REQUIRE(step_index_1.get_position(path_end, graph) == step_index_1.get_path_len(path));
						REQUIRE(step_index_1.get_position(path_back_from_end, graph) == step_index_1.get_position(graph.path_back(path), graph));
					}
				});
			}

			SECTION("The index delivers the correct positions for a given step. Sample rate: 2.") {
				step_index_t step_index_2(graph, paths, 1, false, 2);
				graph.for_each_path_handle([&](const path_handle_t path) {
					std::string cur_path = graph.get_path_name(path);
					if (cur_path == "target") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_2.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_2.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_2.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_2.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_2.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_2.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query1") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_2.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_2.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_2.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_2.get_position(occ, graph) == 5);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query2") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_2.get_position(occ, graph) == 0);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query3") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_2.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_2.get_position(occ, graph) == 3);
									break;
								case 2:
									REQUIRE(step_index_2.get_position(occ, graph) == 4);
									break;
								case 3:
									REQUIRE(step_index_2.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_2.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_2.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}
				});
			}

			SECTION("The index delivers the correct positions for a given step. Sample rate: 4.") {
				step_index_t step_index_4(graph, paths, 1, false, 4);
				graph.for_each_path_handle([&](const path_handle_t path) {
					std::string cur_path = graph.get_path_name(path);
					if (cur_path == "target") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_4.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_4.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_4.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_4.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_4.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_4.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query1") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_4.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_4.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_4.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_4.get_position(occ, graph) == 5);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query2") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_4.get_position(occ, graph) == 0);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query3") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_4.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_4.get_position(occ, graph) == 3);
									break;
								case 2:
									REQUIRE(step_index_4.get_position(occ, graph) == 4);
									break;
								case 3:
									REQUIRE(step_index_4.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_4.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_4.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}
				});
			}

			SECTION("The index delivers the correct positions for a given step. Sample rate: 4.") {
				step_index_t step_index_4(graph, paths, 1, false, 4);
				graph.for_each_path_handle([&](const path_handle_t path) {
					std::string cur_path = graph.get_path_name(path);
					if (cur_path == "target") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_4.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_4.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_4.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_4.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_4.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_4.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query1") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_4.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_4.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_4.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_4.get_position(occ, graph) == 5);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query2") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_4.get_position(occ, graph) == 0);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query3") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_4.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_4.get_position(occ, graph) == 3);
									break;
								case 2:
									REQUIRE(step_index_4.get_position(occ, graph) == 4);
									break;
								case 3:
									REQUIRE(step_index_4.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_4.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_4.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}
				});
			}

			SECTION("The index delivers the correct positions for a given step. Sample rate: 8.") {
				step_index_t step_index_8(graph, paths, 1, false, 8);
				graph.for_each_path_handle([&](const path_handle_t path) {
					std::string cur_path = graph.get_path_name(path);
					if (cur_path == "target") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_8.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_8.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_8.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_8.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_8.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_8.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query1") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_8.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_8.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_8.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_8.get_position(occ, graph) == 5);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query2") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_8.get_position(occ, graph) == 0);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query3") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_8.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_8.get_position(occ, graph) == 3);
									break;
								case 2:
									REQUIRE(step_index_8.get_position(occ, graph) == 4);
									break;
								case 3:
									REQUIRE(step_index_8.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_8.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_8.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}
				});
			}

			SECTION("The index delivers the correct positions for a given step. Sample rate: 16.") {
				step_index_t step_index_16(graph, paths, 1, false, 16);
				graph.for_each_path_handle([&](const path_handle_t path) {
					std::string cur_path = graph.get_path_name(path);
					if (cur_path == "target") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_16.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_16.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_16.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_16.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_16.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_16.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query1") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_16.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_16.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_16.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_16.get_position(occ, graph) == 5);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query2") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_16.get_position(occ, graph) == 0);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query3") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_16.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_16.get_position(occ, graph) == 3);
									break;
								case 2:
									REQUIRE(step_index_16.get_position(occ, graph) == 4);
									break;
								case 3:
									REQUIRE(step_index_16.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_16.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_16.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}
				});
			}

			SECTION("The index can be saved and loaded. After loading, we can retrieve all the correct positions again.") {
				step_index_t step_index_to_save(graph, paths, 1, false, 8);
				// Write index to temporary file in preparation for the next tests.
				std::string basename = xp::temp_file::create();
				step_index_to_save.save(basename + "unittest.stpidx");

				step_index_t step_index_loaded;
				step_index_loaded.load(basename + "unittest.stpidx");

				graph.for_each_path_handle([&](const path_handle_t path) {
					std::string cur_path = graph.get_path_name(path);
					if (cur_path == "target") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query1") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 1);
									break;
								case 2:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 2);
									break;
								case 3:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 5);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query2") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 0);
									break;
							}
							cur_step_rank++;
						});
					}

					if (cur_path == "query3") {
						uint64_t cur_step_rank = 0;
						graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
							switch(cur_step_rank) {
								case 0:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 0);
									break;
								case 1:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 3);
									break;
								case 2:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 4);
									break;
								case 3:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 5);
									break;
								case 4:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 8);
									break;
								case 5:
									REQUIRE(step_index_loaded.get_position(occ, graph) == 11);
									break;
							}
							cur_step_rank++;
						});
					}
				});
			}

		}
	}
}
