#include "threads.hpp"

namespace odgi {

int get_thread_count(void) {
    int thread_count = 1;
#pragma omp parallel
    {
#pragma omp master
        thread_count = omp_get_num_threads();
    }
    return thread_count;
}

}
