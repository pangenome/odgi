#!/bin/bash

# path to the ODGI executable
OG=$1
# path to the ODGI test folder
TEST=$2

echo " [binary_tester::untangle] INFO: Testing with default settings."
diff -u "$TEST"/binary/untangle/default <("$OG" untangle -i "$TEST"/overlap.gfa)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::untangle] SUCCESS: Testing with default settings."
else
    echo " [binary_tester::untangle] FAILED: Testing with default settings."
    exit 1
fi