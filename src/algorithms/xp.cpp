#include "xp.hpp"

namespace xp {

    XP::~XP(void) {
        // Clean up any created XGPaths
        while (!paths.empty()) {
            delete paths.back();
            paths.pop_back();
        }
    }

}