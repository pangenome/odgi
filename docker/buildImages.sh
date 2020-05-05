#!/bin/sh
user=$1
repository=$2
if [ -z $repository ]
then
    repository=localhost
fi
docker build --target binary -t ${repository}/${user}/odgi:latest .
docker build --target debug -t ${repository}/${user}/odgi-debug:latest .
docker build --target python3 -t ${repository}/${user}/odgi-python3:latest .
