#include "diffpriv.hpp"

namespace odgi {
namespace algorithms {

void diff_priv(
    const PathHandleGraph& input,
    PathHandleGraph& priv,
    double target_coverage,
    uint64_t bp_limit) {

    // algorithm
    // we randomly sample a starting node (todo: step)
    // we collect all steps on the node
    // we then sample a direction
    // and look for potential extensions of open path intervals

}

}
}
