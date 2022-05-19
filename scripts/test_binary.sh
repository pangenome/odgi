#!/bin/bash

# path to the ODGI executable
OG=$1
# path to the ODGI test folder
TEST=$2
# path to scripts folder
SC=$3

if [[ ! $# -eq 3 ]] ; then
    echo "[binary_tester] ERROR: The binary tests need exactly three arguments provided in the following order: 1. The path to the ODGI executable. 2. The path to the ODGI test folder. 3. The path to the scripts folder."
    exit 1
fi


echo "[binary_tester] INFO: Path to ODGI executable: ""$OG"
echo "[binary_tester] INFO: Path to ODGI test folder: ""$TEST"

echo "[binary_tester] INFO: Running binary tests of odgi position."
bash "$SC"/position.sh "$OG" "$TEST"
ret=$?
if [[ $ret -eq 0 ]]; then
    echo "[binary_tester] SUCCESS: All binary tests for odgi position passed."
else
    echo "[binary_tester] FAILED: At least one binary test for odgi position failed."
    exit 1
fi

echo "[binary_tester] INFO: Running binary tests of odgi degree."
bash "$SC"/degree.sh "$OG" "$TEST"
ret=$?
if [[ $ret -eq 0 ]]; then
    echo "[binary_tester] SUCCESS: All binary tests for odgi degree passed."
else
    echo "[binary_tester] FAILED: At least one binary test for odgi degree failed."
    exit 1
fi

echo "[binary_tester] INFO: Running binary tests of odgi untangle."
bash "$SC"/untangle.sh "$OG" "$TEST"
ret=$?
if [[ $ret -eq 0 ]]; then
    echo "[binary_tester] SUCCESS: All binary tests for odgi untangle passed."
else
    echo "[binary_tester] FAILED: At least one binary test for odgi untangle failed."
    exit 1
fi
