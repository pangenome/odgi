#include "odgi.hpp"
#include "dynamic_structs.hpp"

using namespace dg;

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

int main(void){
    
    test_dynamic_vector();
    
    return 0;
}
