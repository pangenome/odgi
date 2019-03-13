#ifndef __odgi_split_hpp
#define __odgi_split_hpp

#include <string>
#include <sstream>
#include <vector>
#include <iterator>

namespace odgi {

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim);

}

#endif
