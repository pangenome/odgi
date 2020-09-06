#include "layout.hpp"

namespace odgi {
namespace algorithms {
namespace layout {

using namespace handlegraph;

void to_tsv(std::ostream &out, const std::vector<double> &X, const std::vector<double> &Y, const HandleGraph &graph) {
    uint64_t n = graph.get_node_count() * 2;
    out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    out << "idx" << "\t" << "X" << "\t" << "Y" << std::endl;
    for (uint64_t i = 0; i < n; ++i) {
        out << i << "\t" << X[i] << "\t" << Y[i] << std::endl;
    }
}

Layout::Layout(const std::vector<double> &X, const std::vector<double> &Y) {
    std::vector<uint64_t> vals;
    //assert(X.size)( == Y.size());
    for (auto& v : X) { min_value = std::min(v, min_value); }
    for (auto& v : Y) { min_value = std::min(v, min_value); }
    conv_t x, y;
    for (uint64_t i = 0; i < X.size(); ++i) {
        x.d = X[i] - min_value;
        y.d = Y[i] - min_value;
        vals.push_back(x.i);
        vals.push_back(y.i);
    }
    sdsl::util::assign(xy, sdsl::enc_vector<>(vals));
}

void Layout::serialize(std::ostream& out) {
    sdsl::write_member(min_value, out);
    xy.serialize(out);
}

void Layout::load(std::istream& in) {
    sdsl::read_member(min_value, in);
    xy.load(in);
}

xy_d_t Layout::coords(const handle_t& handle) {
    uint64_t idx = 2 * number_bool_packing::unpack_number(handle);
    conv_t x, y;
    x.i = xy[idx];
    y.i = xy[idx+1];
    return { x.d + min_value, y.d + min_value };
}

}
}
}
