#!/bin/bash

# path to the ODGI executable
OG=$1
# path to the ODGI test folder
TEST=$2

# echo " [binary_tester::degree] INFO: Path to ODGI executable: ""$OG"
# echo " [binary_tester::degree] INFO: Path to ODGI test folder: ""$TEST"

echo " [binary_tester::degree] INFO: Testing with default settings."
diff -u "$TEST"/binary/degree/default <("$OG" degree -i "$TEST"/overlap.gfa)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing with default settings."
else
    echo " [binary_tester::degree] FAILED: Testing with default settings."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing subset paths."
diff -u "$TEST"/binary/degree/subset_paths <("$OG" degree -i "$TEST"/overlap.gfa -s "$TEST"/binary/degree/paths)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing subset paths."
else
    echo " [binary_tester::degree] FAILED: Testing subset paths."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing path."
diff -u "$TEST"/binary/degree/path <("$OG" degree -i "$TEST"/overlap.gfa -r target)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing path."
else
    echo " [binary_tester::degree] FAILED: Testing path."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing paths."
diff -u "$TEST"/binary/degree/paths_ <("$OG" degree -i "$TEST"/overlap.gfa -R "$TEST"/binary/degree/paths)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing paths."
else
    echo " [binary_tester::degree] FAILED: Testing paths."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing graph position."
diff -u "$TEST"/binary/degree/graph_pos <("$OG" degree -i "$TEST"/overlap.gfa -g 8,2)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing graph position."
else
    echo " [binary_tester::degree] FAILED: Testing graph position."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing graph positions."
diff -u "$TEST"/binary/degree/graph_pos_file_ <("$OG" degree -i "$TEST"/overlap.gfa -G "$TEST"/binary/degree/graph_pos_file)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing graph positions."
else
    echo " [binary_tester::degree] FAILED: Testing graph positions."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing path position."
diff -u "$TEST"/binary/degree/path_pos <("$OG" degree -i "$TEST"/overlap.gfa -p target,3,+)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing path position."
else
    echo " [binary_tester::degree] FAILED: Testing path position."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing path positions."
diff -u "$TEST"/binary/degree/path_pos_file_ <("$OG" degree -i "$TEST"/overlap.gfa -F "$TEST"/binary/degree/path_pos_file)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: path positions."
else
    echo " [binary_tester::degree] FAILED: path positions."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing bed input."
diff -u "$TEST"/binary/degree/bed_input <("$OG" degree -i "$TEST"/overlap.gfa -b "$TEST"/binary/degree/bed)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: bed input."
else
    echo " [binary_tester::degree] FAILED: bed input."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing graph degree table."
diff -u "$TEST"/binary/degree/graph_degree_table <("$OG" degree -i "$TEST"/overlap.gfa -d)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing graph degree table."
else
    echo " [binary_tester::degree] FAILED: Testing graph degree table."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing graph degree vec."
diff -u "$TEST"/binary/degree/graph_degree_vec <("$OG" degree -i "$TEST"/overlap.gfa -v)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing graph degree vec."
else
    echo " [binary_tester::degree] FAILED: Testing graph degree vec."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing path degree."
diff -u "$TEST"/binary/degree/path_degree <("$OG" degree -i "$TEST"/overlap.gfa -D)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing path degree."
else
    echo " [binary_tester::degree] FAILED: Testing path degree."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing self degree."
diff -u "$TEST"/binary/degree/self_degree <("$OG" degree -i "$TEST"/overlap.gfa -a)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing self degree."
else
    echo " [binary_tester::degree] FAILED: Testing self degree."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing summarize."
diff -u "$TEST"/binary/degree/summarize <("$OG" degree -i "$TEST"/overlap.gfa -S)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing summarize."
else
    echo " [binary_tester::degree] FAILED: Testing summarize."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing windows in."
diff -u "$TEST"/binary/degree/windows_in <("$OG" degree -i "$TEST"/overlap.gfa -w 10:0:5)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing windows in."
else
    echo " [binary_tester::degree] FAILED: Testing windows in."
    exit 1
fi

echo " [binary_tester::degree] INFO: Testing windows out."
diff -u "$TEST"/binary/degree/windows_out <("$OG" degree -i "$TEST"/overlap.gfa -W 10:0:5)
ret=$?
if [[ $ret -eq 0 ]]; then
    echo " [binary_tester::degree] SUCCESS: Testing windows out."
else
    echo " [binary_tester::degree] FAILED: Testing windows out."
    exit 1
fi