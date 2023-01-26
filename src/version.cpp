// modified from https://github.com/vgteam/vg/blob/master/src/version.cpp

#include "version.hpp"

// Get the git version macro from the build system
#include "../include/odgi_git_version.hpp"

#include <iostream>
#include <sstream>

// If the ODGI_GIT_VERSION deosn't exist at all, define a placeholder
// This lets us be somewhat robust to undeterminable versions
#ifndef ODGI_GIT_VERSION
    #define ODGI_GIT_VERSION "not-from-git"
#endif

// Define a way to quote macro values.
// See https://stackoverflow.com/a/196093
//#define QUOTE(arg) #arg
// We need another level to get the macro's value and not its name.
//#define STR(macro) QUOTE(macro)

namespace odgi {

using namespace std;

// Define all the strings as the macros' values
const string Version::VERSION = ODGI_GIT_VERSION;

// Keep the list of codenames.
// Add new codenames here
const unordered_map<string, string> Version::codenames = {
    {"v0.1.0", "initial release"},
    {"v0.2.0", "point release often"},
    {"v0.3.0", "quarantine lockdown magic"},
    {"v0.4.0", "edgy"},
    {"v0.4.1", "back to old ABI"},
    {"v0.5.0", "fastify everything"},
    {"v0.5.1", "Phoenix"},
    {"0.6", "Domani"},
    {"v0.6.1", "Froggi"},
    {"v0.6.2", "Auff"},
	{"v0.6.3", "Pulizia"},
    {"v0.7.0", "Presente"},
    {"v0.7.1", "Pasticcione"},
    {"v0.7.2", "Radler"},
    {"v0.7.3", "Fissaggio"},
    {"v0.8.0", "Nascondino"},
    {"v0.8.1", "Piccino"}
    // Add more codenames here
};

string Version::get_version() {
    return VERSION;
}

string Version::get_release() {
    auto dash = VERSION.find('-');
    if (dash == -1) {
        // Pure tag versions have no dash
        return VERSION;
    } else {
        // Otherwise it is tag-count-hash and the tag describes the release
        return VERSION.substr(0, dash);
    }
}

string Version::get_codename() {
    auto release = get_release();

    auto found = codenames.find(release);

    if (found == codenames.end()) {
        // No known codename for this release.
        // Return an empty string so we can just not show it.
        return "";
    } else {
        // We have a known codename!
        return found->second;
    }
}

string Version::get_short() {
    stringstream s;
    s << VERSION;

    auto codename = get_codename();
    if (!codename.empty()) {
        // Add the codename if we have one
        s << " \"" << codename << "\"";
    }

    return s.str();
}

}
