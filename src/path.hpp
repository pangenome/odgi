#ifndef dank_path_hpp
#define dank_path_hpp

#include "dynamic.hpp"
#include "dynamic_types.hpp"
#include "handle_types.hpp"
#include "handle_helper.hpp"
#include "dna.hpp"

namespace dankgraph {

struct step_t {
    id_t id = 0;
    bool strand = 0;
    std::string seq;
};

class path_t {
public:
    /// the path name
    std::string name;
    /// constructor
    path_t(const std::string& n) { name = n; }
    /// append a step across the given id with the orientation given by strand
    void append_occurrence(const id_t& id, bool strand);
    /// remove all elements
    void clear(void);
    /// the number of steps in the path
    uint64_t occurrence_count(void) const;
    /// construct a step object that describes the given step, which may include non-graph sequence
    step_t get_occurrence(uint64_t rank) const;
    /// unlink the occurrence from the graph handle, storing the sequence in the path itself
    void unlink_occurrence(uint64_t rank, const std::string& seq);
    /// relink the occurrence to the graph, using the given handle as target and sequence to verify correctness
    void link_occurrence(uint64_t rank, const handle_t& handle, const std::string& seq);
    /// replace the occurrence with a series of steps (used in divide_handle)
    void replace_occurrence(uint64_t rank, const std::vector<handle_t>& handles);
private:
    /// store the ids in the path
    /// zeros indicate privately stored sequences in unlinked occurrences
    //wt_str ids_wt;
    dyn::packed_vector ids_pv;
    /// the strand of each step
    suc_bv strands_wt;
    /// mark the unlinked sequences in ids_pv
    gap_bv unlinked_gbv;
    /// sequence that is in this path, but not represented in the graph
    /// for instance, after the removal of nodes from the graph
    wt_str seq_wt;
};

inline void path_t::clear(void) {
    //wt_str null_wt;
    dyn::packed_vector null_pv;
    wt_str null_wt;
    suc_bv null_bv;
    gap_bv null_gbv;
    ids_pv = null_pv;
    strands_wt = null_bv;
    unlinked_gbv = null_gbv;
    seq_wt = null_wt;
}

inline void path_t::append_occurrence(const id_t& id, bool strand) {
    ids_pv.push_back(id);
    strands_wt.push_back(strand);
    unlinked_gbv.push_back(0);
}

inline uint64_t path_t::occurrence_count(void) const {
    return ids_pv.size();
}

inline step_t path_t::get_occurrence(uint64_t rank) const {
    step_t result;
    result.id = ids_pv.at(rank);
    result.strand = strands_wt.at(rank);
    if (!result.id) {
        uint64_t x = unlinked_gbv.rank1(rank);
        for (uint64_t i = seq_wt.select(x, 0)+1; ; ++i) {
            char c = seq_wt.at(i);
            if (!c) break;
            result.seq += c;
        }
    }
    return result;
}

inline void path_t::unlink_occurrence(uint64_t rank, const std::string& seq) {
    // set the path step id to 0
    ids_pv.remove(rank);
    ids_pv.insert(rank, 0);
    unlinked_gbv.set(rank);
    uint64_t null_rank = unlinked_gbv.rank1(rank);
    // insert the sequence in the right place in seq_wt
    if (seq_wt.size() == 0) {
        // at the end
        assert(null_rank == 0);
        seq_wt.push_back(0); // start and end with 0s
        for (auto c : seq) seq_wt.push_back(c);
        seq_wt.push_back(0);
    } else {
        // at a particular place
        uint64_t i = seq_wt.select(null_rank, 0)+1;
        seq_wt.insert(i, 0); // insert the 0
        // write in reverse
        for (uint64_t j = seq.size()-1; j >= 0; --j) {
            seq_wt.insert(i, seq.at(j));
        }
    }
    // CAUTION we've appended the seq in its natural orientation in the graph
    // and the orientation is maintained in the strand_bv
    // we have to refer to this when e.g. serializing the path
}

inline void path_t::link_occurrence(uint64_t rank, const handle_t& handle, const std::string& seq) {
    // get the original step to verify validity of new step sequence and set orientation
    step_t step = get_occurrence(rank);
    // assert that it's unlinked
    assert(step.id == 0);
    // verify the sequence, assume that seq is that of the new handle
    assert(step.seq == (handle_helper::unpack_bit(handle) == step.strand ? seq : reverse_complement(seq)));
    // set the path step id to that of the handle
    ids_pv.remove(rank);
    ids_pv.insert(rank, handle_helper::unpack_number(handle));
    // set the strand flag to that of the handle
    strands_wt.remove(rank);
    strands_wt.insert(rank, handle_helper::unpack_bit(handle));
    // remove the sequence from seq_wt
    uint64_t i = seq_wt.select(unlinked_gbv.rank1(rank), 0)+1;
    while (seq_wt.at(i) != 0) {
        seq_wt.remove(i);
    }
    seq_wt.remove(i); // remove trailing delimiter 0
    // remove the unlinked marker bit
    unlinked_gbv.remove(rank);
}

inline void path_t::replace_occurrence(uint64_t rank, const std::vector<handle_t>& handles) {
    // delete the step
    ids_pv.remove(rank);
    strands_wt.remove(rank);
    unlinked_gbv.remove(rank);
    // insert the new steps in reverse order
    for (uint64_t i = handles.size()-1; i >= 0; --i) {
        auto& handle = handles.at(i);
        ids_pv.insert(rank, handle_helper::unpack_number(handle));
        strands_wt.insert(rank, handle_helper::unpack_bit(handle));
        unlinked_gbv.insert(rank, 0);
    }
}

}

#endif
