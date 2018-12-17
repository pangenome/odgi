#include "dgraph.hpp"
#include "dynamic_structs.hpp"

#include <vector>
#include <map>
#include <unordered_set>
#include <deque>

using namespace dankgraph;

void test_dynamic_vector() {
    
    enum vec_op_t {SET = 0, GET = 1, APPEND = 2, POP = 3};
    std::random_device rd;
    std::default_random_engine prng(rd());
    std::uniform_int_distribution<int> op_distr(0, 3);
    
    int num_runs = 10000;
    int num_ops = 300;
    int gets_per_op = 5;
    int sets_per_op = 5;
    int appends_per_op = 3;
    int pops_per_op = 1;
    
    for (size_t i = 0; i < num_runs; i++) {
        
        uint64_t next_val = 0;
        
        std::vector<uint64_t> std_vec;
        SuccinctDynamicVector dyn_vec;
        
        for (size_t j = 0; j < num_ops; j++) {
            
            vec_op_t op = (vec_op_t) op_distr(prng);
            switch (op) {
                case SET:
                    if (!std_vec.empty()) {
                        for (size_t k = 0; k < sets_per_op; k++) {
                            size_t idx = prng() % dyn_vec.size();
                            std_vec[idx] = next_val;
                            dyn_vec.set(idx, next_val);
                            next_val++;
                        }
                    }
                    
                    break;
                    
                case GET:
                    if (!std_vec.empty()) {
                        for (size_t k = 0; k < gets_per_op; k++) {
                            size_t idx = prng() % dyn_vec.size();
                            assert(std_vec[idx] == dyn_vec.get(idx));
                            next_val++;
                        }
                    }
                    
                    break;
                    
                case APPEND:
                    for (size_t k = 0; k < appends_per_op; k++) {
                        std_vec.push_back(next_val);
                        dyn_vec.append(next_val);
                        next_val++;
                    }
                    
                    break;
                    
                case POP:
                    if (!std_vec.empty()) {
                        for (size_t k = 0; k < pops_per_op; k++) {
                            std_vec.pop_back();
                            dyn_vec.pop();
                        }
                    }
                        
                    break;
                    
                default:
                    break;
            }
                        
            assert(std_vec.empty() == dyn_vec.empty());
            assert(std_vec.size() == dyn_vec.size());
        }
    }
}

void test_deque() {
    
    enum deque_op_t {SET = 0, GET = 1, APPEND_LEFT = 2, POP_LEFT = 3, APPEND_RIGHT = 4, POP_RIGHT = 5};
    std::random_device rd;
    std::default_random_engine prng(rd());
    std::uniform_int_distribution<int> op_distr(0, 5);
    
    int num_runs = 10000;
    int num_ops = 300;
    int gets_per_op = 5;
    int sets_per_op = 5;
    int appends_per_op = 3;
    int pops_per_op = 1;
    
    for (size_t i = 0; i < num_runs; i++) {
        
        uint64_t next_val = 0;
        
        std::deque<uint64_t> std_deq;
        SuccinctDeque suc_deq;
        
        for (size_t j = 0; j < num_ops; j++) {
            
            deque_op_t op = (deque_op_t) op_distr(prng);
            switch (op) {
                case SET:
                    if (!std_deq.empty()) {
                        for (size_t k = 0; k < sets_per_op; k++) {
                            size_t idx = prng() % std_deq.size();
                            std_deq[idx] = next_val;
                            suc_deq.set(idx, next_val);
                            next_val++;
                        }
                    }
                    
                    break;
                    
                case GET:
                    if (!std_deq.empty()) {
                        for (size_t k = 0; k < gets_per_op; k++) {
                            size_t idx = prng() % std_deq.size();
                            assert(std_deq[idx] == suc_deq.get(idx));
                            next_val++;
                        }
                    }
                    
                    break;
                    
                case APPEND_LEFT:
                    for (size_t k = 0; k < appends_per_op; k++) {
                        std_deq.push_front(next_val);
                        suc_deq.append_front(next_val);
                        next_val++;
                    }
                    
                    break;
                    
                case POP_LEFT:
                    for (size_t k = 0; k < pops_per_op && !std_deq.empty(); k++) {
                        std_deq.pop_front();
                        suc_deq.pop_front();
                    }
                    
                    break;
                    
                case APPEND_RIGHT:
                    for (size_t k = 0; k < appends_per_op; k++) {
                        std_deq.push_back(next_val);
                        suc_deq.append_back(next_val);
                        next_val++;
                    }
                    
                    break;
                    
                case POP_RIGHT:
                    for (size_t k = 0; k < pops_per_op && !std_deq.empty(); k++) {
                        std_deq.pop_back();
                        suc_deq.pop_back();
                    }
                    
                    break;
                    
                default:
                    break;
            }
            
            assert(std_deq.empty() == suc_deq.empty());
            assert(std_deq.size() == suc_deq.size());
        }
    }
}

void test_splay_tree() {
    
    enum tree_op_t {INSERT = 0, ERASE = 1, ACCESS = 2};
    std::random_device rd;
    std::default_random_engine prng(rd());
    std::uniform_int_distribution<int> op_distr(0, 2);
    
    int num_runs = 10000;
    int num_ops = 10;
    
    int inserts_per_op = 3;
    int erases_per_op = 1;
    int accesses_per_op = 5;
    
    size_t value_bank_size = num_ops * inserts_per_op + 1;
    size_t max_value = 10 * value_bank_size;
    
    for (size_t i = 0; i < num_runs; i++) {
        
        std::vector<size_t> bank(max_value);
        for (size_t j = 0; j < bank.size(); j++) {
            // only even so we can use lower bound on odds
            bank[j] = j * 2;
        }
        std::shuffle(bank.begin(), bank.end(), prng);
        bank.resize(value_bank_size);
        size_t bank_idx = 0;
        size_t erase_idx = 0;
        
        std::map<size_t, size_t> std_map;
        SuccinctSplayTree suc_map;
        
        std::unordered_set<size_t> erased;
        
        for (size_t j = 0; j < num_ops; j++) {
            
            tree_op_t op = (tree_op_t) op_distr(prng);
            switch (op) {
                case INSERT:
                    for (size_t k = 0; k < inserts_per_op; k++) {
                        std_map[bank[bank_idx]] = bank_idx;
                        suc_map.insert(bank[bank_idx], bank_idx);
                        bank_idx++;
                    }
                    
                    break;
                    
                case ERASE:
                    if (!std_map.empty()) {
                        for (size_t k = 0; k < erases_per_op; k++) {
                            std_map.erase(bank[erase_idx]);
                            suc_map.erase(bank[erase_idx]);
                            erased.insert(bank[erase_idx]);
                            erase_idx++;
                        }
                    }
                    
                    break;
                    
                case ACCESS:
                    if (!std_map.empty()) {
                        for (size_t k = 0; k < accesses_per_op; k++) {
                            size_t access_val = bank[prng() % bank_idx];
                            size_t x = suc_map.find(access_val);
                            auto xi = std_map.find(access_val);
                            
                            
                            
                            if (erased.count(access_val)) {
                                assert(x == 0);
                            }
                            else {
                                assert(suc_map.get_key(x) == xi->first);
                                assert(suc_map.get_value(x) == xi->second);
                            }
                            
                            
                            size_t y = suc_map.first_lower(access_val);
                            auto yi = std_map.lower_bound(access_val);
                            
                            if (!erased.count(access_val)) {
                                assert(suc_map.get_key(y) == yi->first);
                                assert(suc_map.get_value(y) == yi->second);
                            }
                            
                            size_t z = suc_map.first_lower(access_val + 1);
                            auto zi = std_map.lower_bound(access_val + 1);
                            
                            if (!erased.count(access_val)) {
                                zi--;
                                assert(suc_map.get_key(z) == zi->first);
                                assert(suc_map.get_value(z) == zi->second);
                            }
                            
                            if (x != 0) {
                                size_t w = suc_map.next(x);
                                auto wi = xi;
                                wi++;
                                if (wi == std_map.end()) {
                                    assert(w == 0);
                                }
                                else {
                                    assert(suc_map.get_key(w) == wi->first);
                                    assert(suc_map.get_value(w) == wi->second);
                                }
                            }
                        }
                    }
                    
                    break;
                    
                default:
                    break;
            }
            
            assert(std_map.empty() == suc_map.empty());
            assert(std_map.size() == suc_map.size());
        }
    }
    
}

void test_succinct_dynamic_graph() {
    
    {
        // creating and accessing a single node
        SuccinctDynamicSequenceGraph graph;
        handle_t h = graph.create_handle("ATG", 2);
        assert(graph.get_sequence(h) == "ATG");
        assert(graph.get_sequence(graph.flip(h)) == "CAT");
        assert(graph.get_length(h) == 3);
        
        assert(graph.get_handle(graph.get_id(h)) == h);
        assert(!graph.get_is_reverse(h));
        assert(graph.get_is_reverse(graph.flip(h)));
        
        assert(graph.node_size() == 1);
        assert(graph.min_node_id() == graph.get_id(h));
        assert(graph.max_node_id() == graph.get_id(h));
        
        graph.follow_edges(h, true, [](const handle_t& prev) {
            assert(false);
            return true;
        });
        graph.follow_edges(h, false, [](const handle_t& next) {
            assert(false);
            return true;
        });
        
        // creating and accessing a node at the beginning of ID space
        
        handle_t h2 = graph.create_handle("CT", 1);
        assert(graph.get_sequence(h2) == "CT");
        assert(graph.get_sequence(graph.flip(h2)) == "AG");
        assert(graph.get_length(h2) == 2);
        
        assert(graph.get_handle(graph.get_id(h2)) == h2);
        
        assert(graph.node_size() == 2);
        assert(graph.min_node_id() == graph.get_id(h2));
        assert(graph.max_node_id() == graph.get_id(h));
        
        graph.follow_edges(h2, true, [](const handle_t& prev) {
            assert(false);
            return true;
        });
        graph.follow_edges(h2, false, [](const handle_t& next) {
            assert(false);
            return true;
        });
        
        // creating and accessing a node at the end of ID space
        
        handle_t h3 = graph.create_handle("GAC", 4);
        assert(graph.get_sequence(h3) == "GAC");
        assert(graph.get_sequence(graph.flip(h3)) == "GTC");
        assert(graph.get_length(h3) == 3);
        
        assert(graph.get_handle(graph.get_id(h3)) == h3);
        
        assert(graph.node_size() == 3);
        assert(graph.min_node_id() == graph.get_id(h2));
        assert(graph.max_node_id() == graph.get_id(h3));
        
        graph.follow_edges(h3, true, [](const handle_t& prev) {
            assert(false);
            return true;
        });
        graph.follow_edges(h3, false, [](const handle_t& next) {
            assert(false);
            return true;
        });
        
        // creating and accessing in the middle of ID space
        
        handle_t h4 = graph.create_handle("T", 3);
        assert(graph.get_sequence(h4) == "T");
        assert(graph.get_sequence(graph.flip(h4)) == "A");
        assert(graph.get_length(h4) == 1);
        
        assert(graph.get_handle(graph.get_id(h4)) == h4);
        
        assert(graph.node_size() == 4);
        assert(graph.min_node_id() == graph.get_id(h2));
        assert(graph.max_node_id() == graph.get_id(h3));
        
        graph.follow_edges(h4, true, [](const handle_t& prev) {
            assert(false);
            return true;
        });
        graph.follow_edges(h4, false, [](const handle_t& next) {
            assert(false);
            return true;
        });
        
        // add an edge and traverse it
        
        graph.create_edge(h, h2);
        
        bool found1 = false, found2 = false, found3 = false, found4 = false;
        int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
        graph.follow_edges(h, false, [&](const handle_t& next) {
            if (next == h2) {
                found1 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(h2, true, [&](const handle_t& prev) {
            if (prev == h) {
                found2 = true;
            }
            count2++;
            return true;
        });
        graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
            if (prev == graph.flip(h2)) {
                found3 = true;
            }
            count3++;
            return true;
        });
        graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found4 = true;
            }
            count4++;
            return true;
        });
        assert(count1 == 1);
        assert(count2 == 1);
        assert(count3 == 1);
        assert(count4 == 1);
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        
        count1 = count2 = count3 = count4 = 0;
        found1 = found2 = found3 = found4 = false;
        
        // add another edge with a reversal and traverse it
        
        graph.create_edge(h, graph.flip(h3));
        
        bool found5 = false, found6 = false, found7 = false, found8 = false;
        int count5 = 0, count6 = 0;
        graph.follow_edges(h, false, [&](const handle_t& next) {
            if (next == h2) {
                found1 = true;
            }
            else if (next == graph.flip(h3)) {
                found2 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
            if (prev == graph.flip(h2)) {
                found3 = true;
            }
            else if (prev == h3) {
                found4 = true;
            }
            count2++;
            return true;
        });
        graph.follow_edges(h2, true, [&](const handle_t& prev) {
            if (prev == h) {
                found5 = true;
            }
            count3++;
            return true;
        });
        graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found6 = true;
            }
            count4++;
            return true;
        });
        graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
            if (prev == h) {
                found7 = true;
            }
            count5++;
            return true;
        });
        graph.follow_edges(h3, false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found8 = true;
            }
            count6++;
            return true;
        });
        assert(count1 == 2);
        assert(count2 == 2);
        assert(count3 == 1);
        assert(count4 == 1);
        assert(count5 == 1);
        assert(count6 == 1);
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        assert(found5);
        assert(found6);
        assert(found7);
        assert(found8);
        
        count1 = count2 = count3 = count4 = count5 = count6 = 0;
        found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
        
        // add a reversing self-loop
        
        graph.create_edge(h4, graph.flip(h4));
        graph.follow_edges(h4, false, [&](const handle_t& next) {
            if (next == graph.flip(h4)) {
                found1 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
            if (prev == h4) {
                found2 = true;
            }
            count2++;
            return true;
        });
        assert(count1 == 1);
        assert(count2 == 1);
        assert(found1);
        assert(found2);
        
        count1 = count2 = 0;
        found1 = found2 = false;
        
        // create new edges and then delete them, repeat the previous checks
        
        graph.create_edge(h, graph.flip(h4));
        graph.create_edge(graph.flip(h3), h4);
        
        graph.destroy_edge(h, graph.flip(h4));
        graph.destroy_edge(graph.flip(h3), h4);
        
        graph.follow_edges(h, false, [&](const handle_t& next) {
            if (next == h2) {
                found1 = true;
            }
            else if (next == graph.flip(h3)) {
                found2 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
            if (prev == graph.flip(h2)) {
                found3 = true;
            }
            else if (prev == h3) {
                found4 = true;
            }
            count2++;
            return true;
        });
        graph.follow_edges(h2, true, [&](const handle_t& prev) {
            if (prev == h) {
                found5 = true;
            }
            count3++;
            return true;
        });
        graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found6 = true;
            }
            count4++;
            return true;
        });
        graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
            if (prev == h) {
                found7 = true;
            }
            count5++;
            return true;
        });
        graph.follow_edges(h3, false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found8 = true;
            }
            count6++;
            return true;
        });
        assert(count1 == 2);
        assert(count2 == 2);
        assert(count3 == 1);
        assert(count4 == 1);
        assert(count5 == 1);
        assert(count6 == 1);
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        assert(found5);
        assert(found6);
        assert(found7);
        assert(found8);
        
        count1 = count2 = count3 = count4 = count5 = count6 = 0;
        found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
        
        graph.follow_edges(h4, false, [&](const handle_t& next) {
            if (next == graph.flip(h4)) {
                found1 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
            if (prev == h4) {
                found2 = true;
            }
            count2++;
            return true;
        });
        assert(count1 == 1);
        assert(count2 == 1);
        assert(found1);
        assert(found2);
        
        count1 = count2 = 0;
        found1 = found2 = false;
        
        // add and remove a node with edges, then repeat previous checks
        
        handle_t h5 = graph.create_handle("GGACC");
        
        // make some edges to ensure that this is difficult
        graph.create_edge(h, h5);
        graph.create_edge(h5, h);
        graph.create_edge(graph.flip(h5), h2);
        graph.create_edge(h3, graph.flip(h5));
        graph.create_edge(h3, h5);
        graph.create_edge(h5, h4);
        
        graph.destroy_handle(h5);
        
        graph.follow_edges(h, false, [&](const handle_t& next) {
            if (next == h2) {
                found1 = true;
            }
            else if (next == graph.flip(h3)) {
                found2 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
            if (prev == graph.flip(h2)) {
                found3 = true;
            }
            else if (prev == h3) {
                found4 = true;
            }
            count2++;
            return true;
        });
        graph.follow_edges(h2, true, [&](const handle_t& prev) {
            if (prev == h) {
                found5 = true;
            }
            count3++;
            return true;
        });
        graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found6 = true;
            }
            count4++;
            return true;
        });
        graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
            if (prev == h) {
                found7 = true;
            }
            count5++;
            return true;
        });
        graph.follow_edges(h3, false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found8 = true;
            }
            count6++;
            return true;
        });
        assert(count1 == 2);
        assert(count2 == 2);
        assert(count3 == 1);
        assert(count4 == 1);
        assert(count5 == 1);
        assert(count6 == 1);
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        assert(found5);
        assert(found6);
        assert(found7);
        assert(found8);
        
        count1 = count2 = count3 = count4 = count5 = count6 = 0;
        found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
        
        graph.follow_edges(h4, false, [&](const handle_t& next) {
            if (next == graph.flip(h4)) {
                found1 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
            if (prev == h4) {
                found2 = true;
            }
            count2++;
            return true;
        });
        assert(count1 == 1);
        assert(count2 == 1);
        assert(found1);
        assert(found2);
        
        count1 = count2 = 0;
        found1 = found2 = false;
        
        // do a few swaps and then repeat the topology checks
        
        graph.swap_handles(h, h2);
        graph.swap_handles(h2, h3);
        
        graph.follow_edges(h, false, [&](const handle_t& next) {
            if (next == h2) {
                found1 = true;
            }
            else if (next == graph.flip(h3)) {
                found2 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
            if (prev == graph.flip(h2)) {
                found3 = true;
            }
            else if (prev == h3) {
                found4 = true;
            }
            count2++;
            return true;
        });
        graph.follow_edges(h2, true, [&](const handle_t& prev) {
            if (prev == h) {
                found5 = true;
            }
            count3++;
            return true;
        });
        graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found6 = true;
            }
            count4++;
            return true;
        });
        graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
            if (prev == h) {
                found7 = true;
            }
            count5++;
            return true;
        });
        graph.follow_edges(h3, false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found8 = true;
            }
            count6++;
            return true;
        });
        assert(count1 == 2);
        assert(count2 == 2);
        assert(count3 == 1);
        assert(count4 == 1);
        assert(count5 == 1);
        assert(count6 == 1);
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        assert(found5);
        assert(found6);
        assert(found7);
        assert(found8);
        
        count1 = count2 = count3 = count4 = count5 = count6 = 0;
        found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
        
        graph.follow_edges(h4, false, [&](const handle_t& next) {
            if (next == graph.flip(h4)) {
                found1 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
            if (prev == h4) {
                found2 = true;
            }
            count2++;
            return true;
        });
        assert(count1 == 1);
        assert(count2 == 1);
        assert(found1);
        assert(found2);
        
        count1 = count2 = 0;
        found1 = found2 = false;
        
        // make sure for each handle visits every handle
        
        graph.for_each_handle([&](const handle_t& handle) {
            if (handle == h) {
                found1 = true;
            }
            else if (handle == h2) {
                found2 = true;
            }
            else if (handle == h3) {
                found3 = true;
            }
            else if (handle == h4) {
                found4 = true;
            }
            else {
                assert(false);
            }
            return true;
        });
        
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        
        found1 = found2 = found3 = found4 = false;
        
        // reverse a node in place
        
        // make sure the sequence reverse complemented correctly
        
        auto check_rev_comp = [](const std::string& seq1, const std::string& seq2) {
            assert(seq1.size() == seq2.size());
            auto it = seq1.begin();
            auto rit = seq2.rbegin();
            for (; it != seq1.end(); it++) {
                if (*it == 'A') {
                    assert(*rit == 'T');
                }
                else if (*it == 'C') {
                    assert(*rit == 'G');
                }
                else if (*it == 'G') {
                    assert(*rit == 'C');
                }
                else if (*it == 'T') {
                    assert(*rit == 'A');
                }
                else if (*it == 'N') {
                    assert(*rit == 'N');
                }
                else {
                    assert(false);
                }
                
                rit++;
            }
        };
        
        std::string seq1 = graph.get_sequence(h);
        h = graph.apply_orientation(graph.flip(h));
        std::string rev_seq1 = graph.get_sequence(h);
        check_rev_comp(seq1, rev_seq1);
        
        // check that the edges are still what we expect
        
        int count7 = 0, count8 = 0;
        graph.follow_edges(h, false, [&](const handle_t& next) {
            count1++;
            return true;
        });
        graph.follow_edges(h, true, [&](const handle_t& prev) {
            if (prev == graph.flip(h2)) {
                found1 = true;
            }
            else if (prev == h3) {
                found2 = true;
            }
            count2++;
            return true;
        });
        graph.follow_edges(graph.flip(h), true, [&](const handle_t& next) {
            count3++;
            return true;
        });
        graph.follow_edges(graph.flip(h), false, [&](const handle_t& prev) {
            if (prev == h2) {
                found3 = true;
            }
            else if (prev == graph.flip(h3)) {
                found4 = true;
            }
            count4++;
            return true;
        });
        graph.follow_edges(h2, true, [&](const handle_t& prev) {
            if (prev == graph.flip(h)) {
                found5 = true;
            }
            count5++;
            return true;
        });
        graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
            if (next == h) {
                found6 = true;
            }
            count6++;
            return true;
        });
        graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
            if (prev == graph.flip(h)) {
                found7 = true;
            }
            count7++;
            return true;
        });
        graph.follow_edges(h3, false, [&](const handle_t& next) {
            if (next == h) {
                found8 = true;
            }
            count8++;
            return true;
        });
        assert(count1 == 0);
        assert(count2 == 2);
        assert(count3 == 0);
        assert(count4 == 2);
        assert(count5 == 1);
        assert(count6 == 1);
        assert(count7 == 1);
        assert(count8 == 1);
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        assert(found5);
        assert(found6);
        assert(found7);
        assert(found8);
        
        count1 = count2 = count3 = count4 = count5 = count6 = count7 = count8 = 0;
        found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
        
        // and now switch it back to the same orientation and repeat the topology checks
        
        h = graph.apply_orientation(graph.flip(h));
        
        graph.swap_handles(h, h2);
        graph.swap_handles(h2, h3);
        
        graph.follow_edges(h, false, [&](const handle_t& next) {
            if (next == h2) {
                found1 = true;
            }
            else if (next == graph.flip(h3)) {
                found2 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h), true, [&](const handle_t& prev) {
            if (prev == graph.flip(h2)) {
                found3 = true;
            }
            else if (prev == h3) {
                found4 = true;
            }
            count2++;
            return true;
        });
        graph.follow_edges(h2, true, [&](const handle_t& prev) {
            if (prev == h) {
                found5 = true;
            }
            count3++;
            return true;
        });
        graph.follow_edges(graph.flip(h2), false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found6 = true;
            }
            count4++;
            return true;
        });
        graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
            if (prev == h) {
                found7 = true;
            }
            count5++;
            return true;
        });
        graph.follow_edges(h3, false, [&](const handle_t& next) {
            if (next == graph.flip(h)) {
                found8 = true;
            }
            count6++;
            return true;
        });
        assert(count1 == 2);
        assert(count2 == 2);
        assert(count3 == 1);
        assert(count4 == 1);
        assert(count5 == 1);
        assert(count6 == 1);
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        assert(found5);
        assert(found6);
        assert(found7);
        assert(found8);
        
        count1 = count2 = count3 = count4 = count5 = count6 = 0;
        found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = false;
        
        graph.follow_edges(h4, false, [&](const handle_t& next) {
            if (next == graph.flip(h4)) {
                found1 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(graph.flip(h4), true, [&](const handle_t& prev) {
            if (prev == h4) {
                found2 = true;
            }
            count2++;
            return true;
        });
        assert(count1 == 1);
        assert(count2 == 1);
        assert(found1);
        assert(found2);
        
        count1 = count2 = 0;
        found1 = found2 = false;
        
        // divide a node into multiple parts
        
        std::vector<size_t> offsets{1, 2};
        std::vector<handle_t> parts = graph.divide_handle(h, offsets);
        
        assert(parts.size() == 3);
        
        assert(graph.get_sequence(parts[0]) == "A");
        assert(graph.get_length(parts[0]) == 1);
        assert(graph.get_sequence(parts[1]) == "T");
        assert(graph.get_length(parts[1]) == 1);
        assert(graph.get_sequence(parts[2]) == "G");
        assert(graph.get_length(parts[2]) == 1);
        
        int count9 = 0, count10 = 0, count11 = 0, count12 = 0;
        bool found9 = false, found10 = false, found11 = false, found12 = false, found13 = false, found14 = false;
        graph.follow_edges(parts[0], false, [&](const handle_t& next) {
            if (next == parts[1]) {
                found1 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(parts[0], true, [&](const handle_t& prev) {
            count2++;
            return true;
        });
        graph.follow_edges(graph.flip(parts[0]), true, [&](const handle_t& prev) {
            if (prev == graph.flip(parts[1])) {
                found2 = true;
            }
            count3++;
            return true;
        });
        graph.follow_edges(graph.flip(parts[0]), false, [&](const handle_t& next) {
            count4++;
            return true;
        });
        
        graph.follow_edges(parts[1], false, [&](const handle_t& next) {
            if (next == parts[2]) {
                found3 = true;
            }
            count5++;
            return true;
        });
        graph.follow_edges(parts[1], true, [&](const handle_t& prev) {
            if (prev == parts[0]) {
                found4 = true;
            }
            count6++;
            return true;
        });
        graph.follow_edges(graph.flip(parts[1]), true, [&](const handle_t& prev) {
            if (prev == graph.flip(parts[2])) {
                found5 = true;
            }
            count7++;
            return true;
        });
        graph.follow_edges(graph.flip(parts[1]), false, [&](const handle_t& next) {
            if (next == graph.flip(parts[0])) {
                found6 = true;
            }
            count8++;
            return true;
        });
        
        graph.follow_edges(parts[2], false, [&](const handle_t& next) {
            if (next == h2) {
                found7 = true;
            }
            else if (next == graph.flip(h3)) {
                found8 = true;
            }
            count9++;
            return true;
        });
        graph.follow_edges(parts[2], true, [&](const handle_t& prev) {
            if (prev == parts[1]) {
                found9 = true;
            }
            count10++;
            return true;
        });
        graph.follow_edges(graph.flip(parts[2]), true, [&](const handle_t& prev) {
            if (prev == graph.flip(h2)) {
                found10 = true;
            }
            else if (prev == h3) {
                found11 = true;
            }
            count11++;
            return true;
        });
        graph.follow_edges(graph.flip(parts[2]), false, [&](const handle_t& next) {
            if (next == graph.flip(parts[1])) {
                found12 = true;
            }
            count12++;
            return true;
        });
        graph.follow_edges(graph.flip(h3), true, [&](const handle_t& prev) {
            if (prev == parts[2]) {
                found13 = true;
            }
            return true;
        });
        graph.follow_edges(h2, true, [&](const handle_t& prev) {
            if (prev == parts[2]) {
                found14 = true;
            }
            return true;
        });
        
        assert(count1 == 1);
        assert(count2 == 0);
        assert(count3 == 1);
        assert(count4 == 0);
        assert(count5 == 1);
        assert(count6 == 1);
        assert(count7 == 1);
        assert(count8 == 1);
        assert(count9 == 2);
        assert(count10 == 1);
        assert(count11 == 2);
        assert(count12 == 1);
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        assert(found5);
        assert(found6);
        assert(found7);
        assert(found8);
        assert(found9);
        assert(found10);
        assert(found11);
        assert(found12);
        assert(found13);
        assert(found14);
        
        count1 = count2 = count3 = count4 = count5 = count6 = count7 = count8 = count9 = count10 = count11 = count12 = 0;
        found1 = found2 = found3 = found4 = found5 = found6 = found7 = found8 = found9 = found10 = found11 = found12 = false;
        
        // make a division on the reverse strand

        std::vector<handle_t> rev_parts = graph.divide_handle(graph.flip(h3), {1});
        
        assert(graph.get_sequence(rev_parts[0]) == "G");
        assert(graph.get_length(rev_parts[0]) == 1);
        assert(graph.get_is_reverse(rev_parts[0]));
        assert(graph.get_sequence(rev_parts[1]) == "TC");
        assert(graph.get_length(rev_parts[1]) == 2);
        assert(graph.get_is_reverse(rev_parts[1]));
        
        graph.follow_edges(rev_parts[0], false, [&](const handle_t& next) {
            if (next == rev_parts[1]) {
                found1 = true;
            }
            count1++;
            return true;
        });
        graph.follow_edges(rev_parts[1], true, [&](const handle_t& prev) {
            if (prev == rev_parts[0]) {
                found2 = true;
            }
            count2++;
            return true;
        });
        graph.follow_edges(graph.flip(rev_parts[1]), false, [&](const handle_t& next) {
            if (next == graph.flip(rev_parts[0])) {
                found3 = true;
            }
            count3++;
            return true;
        });
        graph.follow_edges(graph.flip(rev_parts[0]), true, [&](const handle_t& prev) {
            if (prev == graph.flip(rev_parts[1])) {
                found4 = true;
            }
            count4++;
            return true;
        });
        graph.follow_edges(rev_parts[0], true, [&](const handle_t& prev) {
            if (prev == parts[2]) {
                found5 = true;
            }
            count5++;
            return true;
        });
        graph.follow_edges(rev_parts[1], false, [&](const handle_t& next) {
            count6++;
            return true;
        });
        
        assert(count1 == 1);
        assert(count2 == 1);
        assert(count3 == 1);
        assert(count4 == 1);
        assert(count5 == 1);
        assert(count6 == 0);
        assert(found1);
        assert(found2);
        assert(found3);
        assert(found4);
        assert(found5);
    }
}

int main(void){
    
    test_succinct_dynamic_graph();
    test_dynamic_vector();
    test_deque();
    test_splay_tree();
    
    return 0;
}
