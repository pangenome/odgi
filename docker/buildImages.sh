#!/bin/sh
user=$1
docker build --target binary -t ${user}/odgi:latest .
docker build --target python3 -t ${user}/odgi-python3:latest .
