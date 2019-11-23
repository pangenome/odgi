#pragma once

/**
 * \file temp_file.hpp
 *
 * Defines a static structure to handle temp files and clean them up at exit
 */

#include <iostream>
#include <mutex>
#include <string>
#include <set>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <dirent.h>

namespace odgi {
namespace algorithms {

/**
 * Temporary files. Create with create() and remove with remove(). All
 * temporary files will be deleted when the program exits normally or with
 * std::exit(). The files will be created in a directory determined from
 * environment variables, though this can be overridden with set_dir().
 * The interface is thread-safe.
 */
namespace temp_file {

/// Create a temporary file starting with the given base name
std::string create(const std::string& base);

/// Create a temporary file
std::string create();

/// Remove a temporary file
void remove(const std::string& filename);

/// Set a temp dir, overriding system defaults and environment variables.
void set_dir(const std::string& new_temp_dir);

/// Get the current temp dir
std::string get_dir();

} // namespace temp_file

}
}
