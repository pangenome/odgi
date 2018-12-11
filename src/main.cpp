#include "dgraph.hpp"
#include "dynamic_structs.hpp"

void test_dynamic_vector() {
    
    enum vec_op_t {SET = 0, GET = 1, APPEND = 2, POP = 3};
    std::default_random_engine prng(std::random_device()());
    std::uniform_int_distribution<int> op_distr(0, 3);
    
    int num_runs = 10000;
    int num_ops = 500;
    int mods_per_op = 5;
    
    
    for (size_t i = 0; i < num_runs; i++) {
        
        uint64_t next_val = 0;
        
        vector<uint64_t> std_vec;
        SuccinctDynamicVector dyn_vec;
        
        for (size_t j = 0; j < num_ops; j++) {
            
            vec_op_t op = (vec_op_t) op_distr(prng);
            switch (op) {
                case SET:
                    
                    for (size_t k = 0; k < mods_per_op; k++) {
                        size_t idx = prng() % dyn_vec.size();
                        std_vec[idx] = next_val;
                        dyn_vec[idx] = next_val;
                        next_val++;
                    }
                    
                    break;
                    
                case GET:
                    
                    for (size_t k = 0; k < mods_per_op; k++) {
                        size_t idx = prng() % dyn_vec.size();
                        assert(std_vec[idx] == dyn_vec[idx]);
                        next_val++;
                    }
                    
                    break;
                    
                case APPEND:
                    
                    std_vec.push_back(next_val);
                    dyn_vec.append(next_val);
                    next_val++;
                    
                    break;
                    
                case POP:
                    
                    std_vec.pop_back();
                    dyn_vec.pop();
                    
                    break;
                    
                default:
                    break;
            }
            
            assert(std_vec.empty() == dyn_vec.empty());
            assert(std_vec.size() == dyn_vec.size());
        }
    }
}

int main(void){
    
    test_dynamic_vector();
    
    return 0;
}
