#include "dgraph.hpp"
#include "dynamic_structs.hpp"

#include <vector>
#include <map>
#include <unordered_set>

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
                        std::cerr << "inserting " << bank[bank_idx] << std::endl;
                        std_map[bank[bank_idx]] = bank_idx;
                        suc_map.insert(bank[bank_idx], bank_idx);
                        std::cerr << "inserting finished" << std::endl;
                        suc_map.print_topology(std::cerr);
                        bank_idx++;
                    }
                    
                    break;
                    
                case ERASE:
                    if (!std_map.empty()) {
                        for (size_t k = 0; k < erases_per_op; k++) {
                            std::cerr << "erasing " << bank[erase_idx] << std::endl;
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
                            std::cerr << "accessing " << access_val << std::endl;
                            size_t x = suc_map.find(access_val);
                            auto xi = std_map.find(access_val);
                            
                            if (erased.count(access_val)) {
                                assert(x == 0);
                            }
                            else {
                                assert(suc_map.get_key(x) == xi->first);
                                assert(suc_map.get_value(x) == xi->second);
                            }
                            
                            size_t y = suc_map.lower_bound(access_val);
                            auto yi = std_map.lower_bound(access_val);
                            
                            if (!erased.count(access_val)) {
                                assert(suc_map.get_key(y) == yi->first);
                                assert(suc_map.get_value(y) == yi->second);
                            }
                            
                            size_t z = suc_map.lower_bound(access_val + 1);
                            auto zi = std_map.lower_bound(access_val + 1);
                            
                            if (!erased.count(access_val)) {
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

int main(void){
    
    test_dynamic_vector();
    test_splay_tree();
    
    return 0;
}
