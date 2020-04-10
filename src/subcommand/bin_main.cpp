#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/bin_path_info.hpp"

#include <arrow/io/file.h>
#include <parquet/exception.h>
#include <parquet/stream_writer.h>

#include <filesystem>
#include <regex>

namespace odgi {

using namespace odgi::subcommand;

struct TsvSerializer : public algorithms::BinSerializer {
    TsvSerializer(const std::string& path_delim, bool aggregate_delim) :
        algorithms::BinSerializer(path_delim, aggregate_delim)
    {}

    void write_header(const uint64_t pangenome_length, const uint64_t bin_width) override {
        std::cout << "path.name" << "\t"
                  << "path.prefix" << "\t"
                  << "path.suffix" << "\t"
                  << "bin" << "\t"
                  << "mean.cov" << "\t"
                  << "mean.inv" << "\t"
                  << "mean.pos" << "\t"
                  << "first.nucl" << "\t"
                  << "last.nucl" << std::endl;
    }

    void write_seq(const uint64_t& bin_id, const std::string& seq) override {}

    void write_path(const std::string& path_name, const link_vec_t& links, const bin_map_t &bins) override {
        std::string name_prefix = this->get_path_prefix(path_name);
        std::string name_suffix = this->get_path_suffix(path_name);
        for (auto& entry : bins) {
            auto& bin_id = entry.first;
            auto& info = entry.second;
            if (info.mean_cov) {
                std::cout << path_name << "\t"
                          << name_prefix << "\t"
                          << name_suffix << "\t"
                          << bin_id << "\t"
                          << info.mean_cov << "\t"
                          << info.mean_inv << "\t"
                          << info.mean_pos << "\t"
                          << info.first_nucleotide << "\t"
                          << info.last_nucleotide << std::endl;
            }
        }
    }
};

struct JsonSerializer : public algorithms::BinSerializer {
    static const uint64_t ODGI_JSON_VERSION = 10;

    bool write_seqs;

    JsonSerializer(const std::string& path_delim, bool aggregate_delim, bool write_seqs) :
        algorithms::BinSerializer(path_delim, aggregate_delim),
        write_seqs(write_seqs)
    {}

    void write_header(const uint64_t pangenome_length, const uint64_t bin_width) override {
        std::cout << "{\"odgi_version\": " << ODGI_JSON_VERSION << ",";
        std::cout << "\"bin_width\": " << bin_width << ",";
        std::cout << "\"pangenome_length\": " << pangenome_length << "}" << std::endl;
    };

    void write_seq(const uint64_t& bin_id, const std::string& seq) override {
        if (!this->write_seqs) {
            std::cout << "{\"bin_id\":" << bin_id << "}" << std::endl;
        } else {
            std::cout << "{\"bin_id\":" << bin_id << ","
                      << "\"sequence\":\"" << seq << "\"}" << std::endl;
        }
    }

    void write_path(const std::string& path_name, const link_vec_t& links, const bin_map_t& bins) override {
        std::string name_prefix = this->get_path_prefix(path_name);
        std::string name_suffix = this->get_path_suffix(path_name);
        std::cout << "{\"path_name\":\"" << path_name << "\",";
        if (!this->path_delim.empty()) {
            std::cout << "\"path_name_prefix\":\"" << name_prefix << "\","
                      << "\"path_name_suffix\":\"" << name_suffix << "\",";
        }
        std::cout << "\"bins\":[";
        auto entry_it = bins.begin();
        for (uint64_t i = 0; i < bins.size(); ++i) {
            auto& bin_id = entry_it->first;
            auto& info = entry_it->second;
            std::cout << "[" << bin_id << ","
                      << info.mean_cov << ","
                      << info.mean_inv << ","
                      << info.mean_pos << ","
                      << info.first_nucleotide << ","
                      << info.last_nucleotide << "]";
            if (i+1 != bins.size()) {
                std::cout << ",";
            }
            ++entry_it;
        }
        std::cout << "],";
        std::cout << "\"links\":[";
        for (uint64_t i = 0; i < links.size(); ++i) {
            auto& link = links[i];
            std::cout << "[" << link.first << "," << link.second << "]";
            if (i+1 < links.size()) std::cout << ",";
        }
        std::cout << "]}" << std::endl;
    }
};

parquet::schema::NodePtr UInt64(const std::string& name) {
    using namespace parquet;
    return schema::PrimitiveNode::Make(name, Repetition::REQUIRED, Type::INT64, ConvertedType::UINT_64);
}

parquet::schema::NodePtr String(const std::string& name) {
    using namespace parquet;
    return schema::PrimitiveNode::Make(name, Repetition::REQUIRED, Type::BYTE_ARRAY, ConvertedType::UTF8);
}

struct Link {
    uint64_t path_id, upstream, downstream;

    static std::shared_ptr<parquet::schema::GroupNode> GetSchema() {
        using namespace parquet;
        schema::NodeVector fields{UInt64("path_id"), UInt64("upstream"), UInt64("downstream")};
        auto group_node = schema::GroupNode::Make("schema", Repetition::REQUIRED, fields);
        return std::static_pointer_cast<schema::GroupNode>(group_node);
    }
};

parquet::StreamWriter& operator<<(parquet::StreamWriter& os, const Link& link) {
    os << link.path_id << link.upstream << link.downstream << parquet::EndRow;
    return os;
}

struct PathBin {
    uint64_t path_id;
    uint64_t bin_id;
    algorithms::path_info_t path_info;

    std::vector<std::pair<uint64_t, uint64_t>> ranges() {
        return {{path_info.first_nucleotide, path_info.last_nucleotide}};
    }

    static std::shared_ptr<parquet::schema::GroupNode> GetSchema() {
        using namespace parquet;
        schema::NodeVector fields;
        fields.push_back(UInt64("path_id"));
        fields.push_back(UInt64("bin_id"));
        fields.push_back(schema::Double("mean_cov"));
        fields.push_back(schema::Double("mean_inv"));
        fields.push_back(schema::Double("mean_pos"));
        auto group_node = schema::GroupNode::Make("schema", Repetition::REQUIRED, fields);
        return std::static_pointer_cast<schema::GroupNode>(group_node);
    }
};

parquet::StreamWriter& operator<<(parquet::StreamWriter& os, const PathBin& bin) {
    os << bin.path_id
       << bin.bin_id
       << bin.path_info.mean_cov
       << bin.path_info.mean_inv
       << bin.path_info.mean_pos
       << parquet::EndRow;
    return os;
}

struct Path {
    uint64_t path_id;
    std::string name, prefix, suffix;

    static std::shared_ptr<parquet::schema::GroupNode> GetSchema() {
        using namespace parquet;
        schema::NodeVector fields;
        fields.push_back(UInt64("path_id"));
        fields.push_back(String("name"));
        fields.push_back(String("prefix"));
        fields.push_back(String("suffix"));
        auto group_node = schema::GroupNode::Make("schema", Repetition::REQUIRED, fields);
        return std::static_pointer_cast<schema::GroupNode>(group_node);
    }
};

parquet::StreamWriter& operator<<(parquet::StreamWriter& os, const Path& path) {
    os << path.path_id << path.name << path.prefix << path.suffix << parquet::EndRow;
    return os;
}

struct PathBinRange {
    uint64_t path_id, bin_id, start, end;

    static std::shared_ptr<parquet::schema::GroupNode> GetSchema() {
        using namespace parquet;
        schema::NodeVector fields;
        fields.push_back(UInt64("path_id"));
        fields.push_back(UInt64("bin_id"));
        fields.push_back(UInt64("start"));
        fields.push_back(UInt64("end"));
        auto group_node = schema::GroupNode::Make("schema", Repetition::REQUIRED, fields);
        return std::static_pointer_cast<schema::GroupNode>(group_node);
    }
};

parquet::StreamWriter& operator<<(parquet::StreamWriter& os, const PathBinRange& pbr) {
    os << pbr.path_id << pbr.bin_id << pbr.start << pbr.end << parquet::EndRow;
    return os;
}

class ParquetSerializer : public algorithms::BinSerializer {
    parquet::StreamWriter paths_;
    parquet::StreamWriter path_bins_;
    parquet::StreamWriter path_bin_ranges_;
    parquet::StreamWriter links_;

    uint64_t path_id;
    std::filesystem::path output_dir_;

    static parquet::StreamWriter create_writer(const std::string& filename,
                                               const std::shared_ptr<parquet::schema::GroupNode>& schema) {
        std::shared_ptr<arrow::io::FileOutputStream> outfile;
        PARQUET_ASSIGN_OR_THROW(outfile, arrow::io::FileOutputStream::Open(filename));

        parquet::StreamWriter w{parquet::ParquetFileWriter::Open(outfile, schema)};
        w.SetMaxRowGroupSize(1000000);
        return w;
    }

    static const int ODGI_PARQUET_VERSION = 1;

public:
    ParquetSerializer(const std::string& path_delim, bool aggregate_delim, const std::filesystem::path& output_dir) :
        algorithms::BinSerializer(path_delim, aggregate_delim),
        path_id(0), output_dir_(output_dir)
    {
        std::filesystem::create_directories(output_dir);

        paths_ = create_writer(output_dir / "paths.parquet", Path::GetSchema());
        path_bins_ = create_writer(output_dir / "path_bins.parquet", PathBin::GetSchema());
        path_bin_ranges_ = create_writer(output_dir / "path_bin_ranges.parquet", PathBinRange::GetSchema());
        links_ = create_writer(output_dir / "links.parquet", Link::GetSchema());
    }

    void write_header(const uint64_t pangenome_length, const uint64_t bin_width) override {
        std::ofstream out(output_dir_ / "metadata.json");

        out << "{\"odgi_parquet_version\": " << ODGI_PARQUET_VERSION << ",";
        out << "\"bin_width\": " << bin_width << ",";
        out << "\"pangenome_length\": " << pangenome_length << "}" << std::endl;
    }

    void write_seq(const uint64_t& bin_id, const std::string& seq) override {}

    void write_path(const std::string& name, const link_vec_t& links, const bin_map_t &bins) override {
        paths_ << Path{path_id, name, get_path_prefix(name), get_path_suffix(name)};

        for (auto& bin: bins) {
            PathBin path_bin{path_id, bin.first, bin.second};
            path_bins_ << path_bin;
            for (const auto& range: path_bin.ranges())
                path_bin_ranges_ << PathBinRange{path_id, bin.first, range.first, range.second};
        }

        for (auto& link: links)
            links_ << Link{path_id, link.first, link.second};

        ++path_id;
    }
};

int main_bin(int argc, char** argv) {

    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi bin";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("binning of path information in the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> fa_out_file(parser, "FILE", "store the pangenome sequence in FASTA format in this file", {'f', "fasta"});
    args::ValueFlag<std::string> path_delim(parser, "path-delim", "annotate rows by prefix and suffix of this delimiter", {'D', "path-delim"});
    args::Flag output_json(parser, "write-json", "write JSON format output including additional path positional information", {'j', "json"});
    args::ValueFlag<std::string> parquet_out_dir(parser, "DIR", "write Parquet format output to this directory", {'p', "parquet"});
    args::Flag aggregate_delim(parser, "aggregate-delim", "aggregate on path prefix delimiter", {'a', "aggregate-delim"});
    args::ValueFlag<uint64_t> num_bins(parser, "N", "number of bins", {'n', "num-bins"});
    args::ValueFlag<uint64_t> bin_width(parser, "bp", "width of each bin in basepairs along the graph vector", {'w', "bin-width"});
    args::Flag write_seqs_not(parser, "write-seqs-not", "don't write out the sequences for each bin", {'s', "no-seqs"});
    args::Flag drop_gap_links(parser, "drop-gap-links", "don't include gap links in the output", {'g', "no-gap-links"});
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
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }

    if (args::get(num_bins) + args::get(bin_width) == 0) {
        std::cerr << "[odgi bin] error: a bin width or a bin count is required" << std::endl;
        return 1;
    }

    auto write_fasta = [&](const std::string& nuc_seq) {
        if (fa_out_file) {
            std::ofstream out(args::get(fa_out_file));
            std::string fa_out_name = args::get(fa_out_file).c_str();
            std::regex regex("/");
            auto token_it = std::sregex_token_iterator(fa_out_name.begin(), fa_out_name.end(), regex, -1);
            std::vector<std::string> splitted(token_it, std::sregex_token_iterator());
            fa_out_name = splitted[splitted.size() - 1];
            // Write header
            out << ">" << fa_out_name << std::endl;
            // Write the actual sequences, 80 nucleotides per line
            for (unsigned i = 0; i < nuc_seq.length(); i += 80) {
                std:: string sub_nuc_seq = nuc_seq.substr(i, 80);
                out << sub_nuc_seq << std::endl;
            }
        }
    };

    std::string delim = args::get(path_delim);
    bool agg_delim = args::get(aggregate_delim);

    std::shared_ptr<algorithms::BinSerializer> serializer;
    bool skip_seqs = args::get(write_seqs_not) || fa_out_file;
    if (args::get(output_json)) {
        serializer = std::make_shared<JsonSerializer>(delim, agg_delim, !skip_seqs);
    } else if (parquet_out_dir) {
        serializer = std::make_shared<ParquetSerializer>(delim, agg_delim, args::get(parquet_out_dir));
    } else {
        serializer = std::make_shared<TsvSerializer>(delim, agg_delim);
    }

    algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                              serializer, write_fasta,
                              args::get(num_bins), args::get(bin_width), args::get(drop_gap_links));
    return 0;
}

static Subcommand odgi_bin("bin", "bin path information across the graph",
                              PIPELINE, 3, main_bin);


}
