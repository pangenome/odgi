#include "perfect_neighbors.hpp"

namespace odgi {
namespace algorithms {

//#define debug_perfect_neighbors

/// Return true if nodes share all paths and the mappings they share in these paths
/// are adjacent, in the specified relative order and orientation.
bool nodes_are_perfect_path_neighbors(const PathHandleGraph& graph, handle_t left_handle, handle_t right_handle) {
    
#ifdef debug_perfect_neighbors
    std::cerr << "Check if " << graph.get_id(left_handle) << (graph.get_is_reverse(left_handle) ? "-" : "+") << " and "
              << graph.get_id(right_handle) << (graph.get_is_reverse(right_handle) ? "-" : "+") << " are perfect path neighbors" << std::endl;
#endif

    // Set this false if we find an impermissible step
    bool ok = true;

    // Count the number of permissible steps on the next node we find
    size_t expected_next = 0;

    graph.for_each_step_on_handle(left_handle, [&](const step_handle_t& here) {
        // For each path step on the left

        // We need to work out if the path traverses this handle backward.
        bool step_is_to_reverse_of_handle = (graph.get_handle_of_step(here) != left_handle);

#ifdef debug_perfect_neighbors
        std::cerr << "Consider visit of path " << graph.get_path_name(graph.get_path_handle_of_step(here))
                  << " to " << (step_is_to_reverse_of_handle ? "reverse" : "forward") << " orientation of handle" << std::endl;
#endif

        if (!(step_is_to_reverse_of_handle ? graph.has_previous_step(here) : graph.has_next_step(here))) {
            // If there's no visit to the right of this handle, it can't be to the right next place
#ifdef debug_perfect_neighbors
            std::cerr << "Path stops here so no subsequent handle is a perfect path neighbor" << std::endl;
#endif
            ok = false;
            return false;
        }

        // Walk along the path whatever direction is forward relative to our left handle.
        step_handle_t step_to_right = step_is_to_reverse_of_handle ? graph.get_previous_step(here) : graph.get_next_step(here);
        handle_t handle_to_right = graph.get_handle_of_step(step_to_right);
        if (step_is_to_reverse_of_handle) {
            // If we had to walk back along the reverse strand of the path, we have to flip where we ended up.
            handle_to_right = graph.flip(handle_to_right);
        }

        if (handle_to_right != right_handle) {
            // It goes to the wrong next place

#ifdef debug_perfect_neighbors
            std::cerr << "Path goes to the wrong place ("
                << graph.get_id(handle_to_right) << (graph.get_is_reverse(handle_to_right) ? "-" : "+")
                << ") and so these nodes are not perfect path neighbors"  << std::endl;
#endif
            ok = false;
            return false;
        }

        // Otherwise, record a step that is allowed to exist on the next handle
        expected_next++;
        return true;
    });

    if (!ok) {
        // We found a bad step, or the path stopped.
        return false;
    }

    // Now count up the path steps on the right handle
    size_t observed_next = 0;
    graph.for_each_step_on_handle(right_handle, [&](const step_handle_t& ignored) {
#ifdef debug_perfect_neighbors
        std::cerr << "Next node has path " << graph.get_path_name(graph.get_path_handle_of_step(ignored)) << std::endl;
#endif
        observed_next++;
    });

#ifdef debug_perfect_neighbors
    if (observed_next != expected_next) {
        std::cerr << "Next node has " << observed_next << " path visits but should have " << expected_next << std::endl;
    }
#endif

    // If there are any steps on the right node that weren't accounted for on
    // the left node, fail. Otherwise, succeed.
    return observed_next == expected_next;

}

}
}
