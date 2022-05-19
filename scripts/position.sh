#!/bin/bash

# path to the ODGI executable
OG=$1
# path to the ODGI test folder
TEST=$2

# echo " [binary_tester::position] INFO: Path to ODGI executable: ""$OG"
# echo " [binary_tester::position] INFO: Path to ODGI test folder: ""$TEST"

echo " [binary_tester::position] INFO: Testing path to node mapping."
diff -u "$TEST"/binary/position/path_node_mapping <("$OG" position -i "$TEST"/k.gfa -p y,10 -v)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing path to node mapping."
else
    echo " [binary_tester::position] FAILED: Testing path to node mapping."
    exit 1
fi

echo " [binary_tester::position] INFO: Testing path to node mapping with reference."
diff -u "$TEST"/binary/position/path_node_mapping_ref <("$OG" position -i "$TEST"/k.gfa -p y,10 -r x)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing path to node mapping with reference."
else
    echo " [binary_tester::position] FAILED: Testing path to node mapping with reference."
    exit 1
fi

echo " [binary_tester::position] INFO: Testing node to node mapping."
diff -u "$TEST"/binary/position/node_node_mapping <("$OG" position -i "$TEST"/k.gfa -g 6)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing node to node mapping."
else
    echo " [binary_tester::position] FAILED: Testing node to node mapping."
    exit 1
fi

echo " [binary_tester::position] INFO: Testing node to node mapping with offset."
diff -u "$TEST"/binary/position/node_node_mapping_offset <("$OG" position -i "$TEST"/k.gfa -g 6,2)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing node to node mapping with offset."
else
    echo " [binary_tester::position] FAILED: Testing node to node mapping with offset."
    exit 1
fi

echo " [binary_tester::position] INFO: Testing node to node mapping with ref."
diff -u "$TEST"/binary/position/node_node_mapping_ref <("$OG" position -i "$TEST"/k.gfa -g 4 -r x)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing node to node mapping with ref."
else
    echo " [binary_tester::position] FAILED: Testing node to node mapping with ref."
    exit 1
fi

echo " [binary_tester::position] INFO: Testing path to path mapping."
diff -u "$TEST"/binary/position/path_path_mapping_1 <("$OG" position -i "$TEST"/overlap.gfa -r target -p query3,0)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing path to path mapping 1."
else
    echo " [binary_tester::position] FAILED: Testing path to path mapping 1."
    exit 1
fi

diff -u "$TEST"/binary/position/path_path_mapping_2 <("$OG" position -i "$TEST"/overlap.gfa -r target -p query3,1)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing path to path mapping 2."
else
    echo " [binary_tester::position] FAILED: Testing path to path mapping 2."
    exit 1
fi

diff -u "$TEST"/binary/position/path_path_mapping_3 <("$OG" position -i "$TEST"/overlap.gfa -r target -p query3,2)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing path to path mapping 3."
else
    echo " [binary_tester::position] FAILED: Testing path to path mapping 3."
    exit 1
fi

diff -u "$TEST"/binary/position/path_path_mapping_4 <("$OG" position -i "$TEST"/overlap.gfa -r target -p query3,5)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing path to path mapping 4."
else
    echo " [binary_tester::position] FAILED: Testing path to path mapping 4."
    exit 1
fi

echo " [binary_tester::position] INFO: Testing path to path mapping with jaccard."
diff -u "$TEST"/binary/position/path_path_mapping_jaccard <("$OG" position -i "$TEST"/overlap.gfa -r target -p query1,5 -w 2)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing path to path mapping with jaccard."
else
    echo " [binary_tester::position] FAILED: Testing path to path mapping with jaccard."
    exit 1
fi

echo " [binary_tester::position] INFO: Testing GFF lifting for Bandage."
diff -u "$TEST"/binary/position/gff <("$OG" position -i "$TEST"/overlap.gfa -E "$TEST"/overlap.gtf)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::position] SUCCESS: Testing GFF lifting for Bandage."
else
    echo " [binary_tester::position] FAILED: Testing GFF lifting for Bandage."
    exit 1
fi