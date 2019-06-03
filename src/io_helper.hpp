#ifndef gfa_io_helper
#define gfa_io_helper
#include "gfakluge.hpp"
#include "odgi.hpp"
#include <handlegraph/handle_graph.hpp>
#include "dgraph.hpp"
#include "bgraph.hpp"
#include "btypes.hpp"
#include <unordered_set>


namespace dg{

using namespace handlegraph;

struct id_emitter_factory{
    std::unordered_set<string> seen_identifiers;
    nid_t current_id = 0;

    id_emitter_factory(){
        current_id = 0;
        seen_identifiers.reserve(1000);
    };
    nid_t emit_id(string identifier){
        if (!seen_identifiers.count(identifier)){
            return ++current_id;
        }
        else{
            cerr << "Error: duplicated IDs in GFA file: " << identifier << "." << endl;
            exit(1);
        }
    };
};

inline void gfa_to_beta(const char* filename, betagraph::BGraph* bg){
    gfak::GFAKluge gg;
    gg.parse_gfa_file(string(filename));
    
};

inline void dank_to_gfa_stream(SuccinctDynamicSequenceGraph* sd, std::ostream& os) {
    gfak::GFAKluge gg;
    gg.set_version(2.0);

    std::function<bool(const dg::handle_t&)> node_handle_to_kluge = [&](const dg::handle_t& h){
        gfak::sequence_elem s;
        s.name = std::to_string(sd->get_id(h));
        s.length = sd->get_length(h);
        s.sequence = sd->get_sequence(h);
        gg.write_element(os, s);
        return true;
    };
    sd->for_each_handle(node_handle_to_kluge);

    std::function<bool(const edge_t&)> edge_to_kluge = [&](const dg::edge_t& e){
        gfak::edge_elem e_elem;
        e_elem.source_name = std::to_string(sd->get_id(e.first));
        e_elem.sink_name = std::to_string(sd->get_id(e.second));
        e_elem.source_orientation_forward = true;
        e_elem.sink_orientation_forward = true;
        e_elem.type = 1;
        e_elem.ends.set(0,1);
        e_elem.ends.set(1,1);
        e_elem.ends.set(2,0);
        e_elem.ends.set(3,0);
        e_elem.alignment = "0M";
        e_elem.id = "*";
        gg.write_element(os, e_elem);
        return true;
    };

    sd->for_each_edge(edge_to_kluge);
};

inline void dank_from_gfa_file(char* filename, SuccinctDynamicSequenceGraph* sd) {

    // We can't do ID renaming unless we put all S lines in before our edges,
    // and even then it's a little dicey. Ideally, this would get pushed all
    // the way down to GFAK or it'd be used with line-by-line input.
    //id_emitter_factory id_fac;
    gfak::GFAKluge gg;
    gg.parse_gfa_file(std::string(filename));

    // Load up our nodes
    map<string, gfak::sequence_elem, gfak::custom_key> seqs = gg.get_name_to_seq();
    map<string, std::vector<gfak::edge_elem>> edges = gg.get_seq_to_edges();
    for (auto s : seqs){
        //nid_t id = id_fac.emit_id(s.first);

        for (auto e : edges[s.first]){

        }

    };

};
};

#endif
