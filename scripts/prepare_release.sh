NEW_VER=$1 # In format 0.0.0
NEW_CODENAME=$2 # Simple string, could be MAD DOG

if [ -z "$1" ]
  then
    echo "No new version supplied."
    exit 1
fi

if [ -z "$2" ]
  then
    echo "No new codename supplied."
    exit 1
fi

# get the line of the latest release version
# assumes we already have build odgi
# ../src/version.cpp
CUR_VER=$(grep "{\"v" ../src/version.cpp | tail -n 1)
LINE=$(grep -n "$CUR_VER" ../src/version.cpp | awk -F: '{print $1}')
sed -i "s|$CUR_VER|$CUR_VER,|g" ../src/version.cpp

LINE=$(($LINE + 1))
sed -i "$LINE i \ \ \ \ {\"v$NEW_VER\", \"$NEW_CODENAME\"}" ../src/version.cpp

# TODO generate updated documentation