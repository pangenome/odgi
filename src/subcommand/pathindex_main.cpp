#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/xp.hpp"
#include "utils.hpp"

namespace odgi {

    using namespace odgi::subcommand;
    using namespace xp;

    int main_pathindex(int argc, char **argv) {

        for (uint64_t i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        const std::string prog_name = "odgi pathindex";
        argv[0] = (char *) prog_name.c_str();
        --argc;

        args::ArgumentParser parser("Create a path index for a given graph.");
        args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
        args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*.", {'i', "idx"});
        args::ValueFlag<std::string> idx_out_file(mandatory_opts, "FILE", "Write the succinct variation graph index to this FILE. A file ending with *.xp* is recommended.", {'o', "out"});
        args::Group threading_opts(parser, "[ Threading ]");
        args::ValueFlag<std::uint64_t> nthreads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
		args::Group processing_info_opts(parser, "[ Processing Information ]");
		args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
        args::Group program_info_opts(parser, "[ Program Information ]");
        args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi pathindex.", {'h', "help"});

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

        if (!dg_in_file) {
            std::cerr << "[odgi::pathindex] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]."
                      << std::endl;
            return 1;
        }

        if (!idx_out_file) {
            std::cerr << "[odgi::pathindex] error: please specify an output file to where to store the path index via -o=[FILE], --out=[FILE]."
                      << std::endl;
            return 1;
        }

		const uint64_t num_threads = nthreads ? args::get(nthreads) : 1;

		// read in the graph
        graph_t graph;
        assert(argc > 0);
        {
            const std::string infile = args::get(dg_in_file);
            if (!infile.empty()) {
                if (infile == "-") {
                    graph.deserialize(std::cin);
                } else {
					utils::handle_gfa_odgi_input(infile, "pathindex", args::get(progress), num_threads, graph);

                }
            }
        }

        XP path_index;
        path_index.from_handle_graph(graph, num_threads);
		if (progress) {
			std::cout << "Indexed " << path_index.path_count << " path(s)." << std::endl;
		}

#ifdef debug_pathindex
        size_t pangenome_pos = path_index.get_pangenome_pos("5", 1);
        std::cerr << "Pangenome position for input \"5\":1 in constructed index is: " << pangenome_pos << std::endl;
        pangenome_pos = path_index.get_pangenome_pos("5", 2);
        std::cerr << "Pangenome position for input \"5\":2 in constructed index is: " << pangenome_pos << std::endl;
        pangenome_pos = path_index.get_pangenome_pos("5", 13);
        std::cerr << "Pangenome position for input \"5\":13 in constructed index is: " << pangenome_pos << std::endl;
        pangenome_pos = path_index.get_pangenome_pos("5", 5);
        std::cerr << "Pangenome position for input \"5\":5 in constructed index is: " << pangenome_pos << std::endl;
        pangenome_pos = path_index.get_pangenome_pos("5", 12);
        std::cerr << "Pangenome position for input \"5\":12 in constructed index is: " << pangenome_pos << std::endl;
#endif

        // writ out the index
        std::ofstream out;
        out.open(args::get(idx_out_file));
        if (progress) {
			std::cout << "Writing index to " << args::get(idx_out_file) << "." << std::endl;
		}
        path_index.serialize_members(out);
        out.close();

        return 0;
    }

    static Subcommand odgi_pathindex("pathindex", "Create a path index for a given graph.",
                                     PIPELINE, 3, main_pathindex);

}
