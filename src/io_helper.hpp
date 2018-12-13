#ifndef gfa_io_helper
#define gfa_io_helper
#include "gfakluge.hpp"
#include "handle_types.hpp"
#include <unordered_set>


struct id_emitter_factory{
    std::unordered_set<string> seen_identifiers;
    id_t current_id = 0;

    id_emitter_factory(){
        current_id = 0;
        seen_identifiers.reserce(1000);
    };
    id_t emit_id(string identifier){
        if (!seen_identifiers.count(identifier)){
            return ++current_id;
        }
        else{
            cerr << "Error: duplicated IDs in GFA file: " << identifier << "." << endl;
            exit(1);
        }
    };
};


inline void dank_to_gfa_stream(SuccinctDynamicSequenceGraph* sd, std::ostream& os) const {
    
};

inline void dank_from_gfa_file(char* filename, SuccinctDynamicSequenceGraph* sd) const{

    id_emitter_factory id_fac();
    gfak::GFAKluge gg;
    gg.parse_gfa_file(std::string(filename));


    
    // Load up our nodes
    map<string, gfak::sequence_elem,custom_key> seqs = gg.get_name_to_seq();
    map<string, gfak::edge_elem> edges = gg.get_seq_to_edges();
    for (auto s : seqs){
        id_t id = id_fac.emit_id(s.first);
    }
    int s_counter = 0;
    for (auto s : seqs){
        
    }

};



#endif