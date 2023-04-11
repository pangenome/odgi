INC_DIR=$1

# Go to the directory where the script is
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"

GIT_VERSION=$(git describe --always --tags --long)

echo "#define ODGI_GIT_VERSION" \"$GIT_VERSION\" > $1/odgi_git_version.hpp.tmp
diff $1/odgi_git_version.hpp.tmp $1/odgi_git_version.hpp >/dev/null || cp $1/odgi_git_version.hpp.tmp $1/odgi_git_version.hpp
rm -f $1/odgi_git_version.hpp.tmp
