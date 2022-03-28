#include "bubble.hpp"

namespace odgi {

namespace algorithms {

uint64_t hash_vec(std::vector<uint64_t> const& vec) {
    uint64_t seed = vec.size();
    for(auto& i : vec) {
        seed ^= i + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

size_t hash_steps_on_handle(const PathHandleGraph& graph,
                            const handle_t& handle) {
    std::vector<uint64_t> v;
    graph.for_each_step_on_handle(
        handle,
        [&](const step_handle_t& step) {
            v.push_back(as_integer(graph.get_path_handle_of_step(step)));
        });
    std::sort(v.begin(), v.end());
    return hash_vec(v);
}

void for_each_bubble(const PathHandleGraph& graph,
                     const step_index_t& step_index,
                     const path_handle_t& path,
                     const std::function<void(const bubble_t& bubble)>& func) {

    std::vector<step_hash_t> step_hashes;
    uint64_t pos = 0;
    graph.for_each_step_in_path(
        path, [&](const step_handle_t& step) {
            auto handle = graph.get_handle_of_step(step);
            auto depth = graph.get_step_count(handle);
            auto length = graph.get_length(handle);
            step_hashes.push_back(
                { step,
                  pos,
                  length,
                  depth,
                  hash_steps_on_handle(
                      graph,
                      handle)
                });
            pos += length;
        });
    auto j = step_hashes.begin();
    std::vector<bubble_t> bubbles;
    uint64_t max_bubble = 1e6;
    for (auto i = step_hashes.begin(); i != step_hashes.end(); ++i) {
        auto n = i;
        uint64_t distance = 0;
        while (++n != step_hashes.end() && n != j && distance < max_bubble) {
            distance += n->length;
            if (n->hash == i->hash) {
                //bubbles.push_back({&*i, &*n});
                i->other_tail = &*n;
                n->other_head = &*i;
                j = n;
                break;
            }
        }
    }

    /*
    for (auto& bubble : bubbles) {
        std::cout << bubble.head->pos << "@" << bubble.head->depth << " "
                  << bubble.tail->pos << "@" << bubble.tail->depth << std::endl;
    }
    */

    for (auto i = step_hashes.begin(); i != step_hashes.end(); ++i) {
        if (i->other_tail == nullptr || !i->other_tail->is_head_tail()) { // tail is filtered
            uint64_t distance = 0;
            //auto n = (std::vector<step_hash_t>::iterator)i->other_tail;
            auto n = i;
            while (++n != step_hashes.end() && distance < max_bubble) {
                distance += n->length;
                if (n->hash == i->hash
                    && n->other_head == nullptr
                    && n->other_tail != nullptr) {
                    bubbles.push_back({&*i, &*n});
                    i->other_tail = &*n;
                    n->other_head = &*i;
                    break;
                }
            }
        } else {
            auto n = i->other_tail;
            bubbles.push_back({&*i, &*n});
        }
    }

    for (auto i = step_hashes.rbegin(); i != step_hashes.rend(); ++i) {
        if (i->other_head == nullptr || !i->other_head->is_head_tail()) { // tail is filtered
            uint64_t distance = 0;
            //auto n = (std::vector<step_hash_t>::iterator)i->other_tail;
            auto n = i;
            while (++n != step_hashes.rend() && distance < max_bubble) {
                distance += n->length;
                if (n->hash == i->hash
                    && n->other_tail == nullptr
                    && n->other_head != nullptr) {
                    bubbles.push_back({&*n, &*i});
                    i->other_head = &*n;
                    n->other_tail = &*i;
                    break;
                }
            }
        } else {
            auto n = i->other_head;
            bubbles.push_back({&*n, &*i});
        }
    }

    /*
    std::vector<bubble_chain_t> chains;
    std::vector<bool> seen(step_hashes.size());
    for (auto i = step_hashes.begin(); i != step_hashes.end(); ++i) {
        chains.emplace_back();
        auto& c = chains.back().bubbles;
        // follow the chain to the end
        // 
    }
    */
    
    // to filter the bubbles to keep those whose head is the tail of another and vice-versa
    /*
    bubbles.erase(
        std::remove_if(
            bubbles.begin(),
            bubbles.end(),
            [](const bubble_t& bubble) {
                return !(bubble.head->is_head_tail()
                         && bubble.tail->is_head_tail());
            }),
        bubbles.end());
    */

    //if (false) {
    // now draw chains
    /*
    std::cout << "digraph graphname {\n"
              << "node [shape=plaintext];\n";
        //<< "rankdir=LR;\n";
    for (auto& bubble : bubbles) {
        std::cout << "\"" << bubble.head->pos << "\" [label=\"" << bubble.head->pos
                  << "@" << bubble.head->depth << "x\"];\n";
    }
    for (auto& bubble : bubbles) {
        std::cout << "\"" << bubble.head->pos << "\" -> \"" << bubble.tail->pos << "\";\n";
        //std::cout << "bubble " << bubble.head->pos << " " << bubble.head->depth << " -> " << bubble.tail->pos << " " << bubble.tail->depth << std::endl;
        func(bubble);
    }
    std::cout << "}\n";
    */
    //}
    std::cout << "digraph graphname {\n"
              << "node [shape=plaintext];\n";
    for (auto i = step_hashes.begin(); i != step_hashes.end(); ++i) {
        std::cout << "\"" << i->pos << "\" [label=\"" << i->pos
                  << "@" << i->depth << "x\"];\n";
    }
    for (auto i = step_hashes.begin(); i != step_hashes.end(); ++i) {
        std::cout << "\"" << (i->other_head != nullptr ? std::to_string(i->other_head->pos) : "null" )
                  << "\" -> \"" << i->pos << "\";\n";
        std::cout << "\"" << i->pos << "\" -> \""
                  << (i->other_tail != nullptr ? std::to_string(i->other_tail->pos) : "null"  ) << "\";\n";
    }
    std::cout << "}\n";
    
}

}

}
