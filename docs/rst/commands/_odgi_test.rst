.. _odgi test:

#########
odgi test
#########

Run ODGI unit tests.

SYNOPSIS
========

**odgi test** [<TEST NAME|PATTERN|TAGS> …] [*OPTION*]…

DESCRIPTION
===========

The odgi test command starts all unit tests that are implemented in
odgi. For targeted testing, a subset of tests can be selected. odgi
test depends on `Catch2 <https://github.com/catchorg/Catch2>`__. In
the default setting, all results are printed to stdout.

OPTIONS
=======

Testing Options
---------------

| **-l, --list-tests**
| List all test cases. If a pattern was specified, all matching test
  cases are listed.

| **-t, --list-tags**
| List all tags. If a pattern was specified, all matching tags are
  listed.

| **-s, --success**
| Include successful tests in output.

| **-b, --break**
| Break into debugger mode upon failed test.

| **-e, --nothrow**
| Skip exception tests.

| **-i, --invisibles**
| Show invisibles like tabs or newlines.

| **-o, --out**\ =\ *FILE*
| Write all output to *FILE*.

| **-r, --reporter**\ =\ *STRING*
| Reporter to use. Default is console.

| **-n, --name**\ =\ *STRING*
| Suite name.

| **-a, --abort**
| Abort at first failure.

| **-x, --abortx**\ =\ *N*
| Abort after *N* failures.

| **-w, --warn**\ =\ *STRING*
| Enable warnings.

| **-d, --durations**\ =\ *yes|no*
| Show test durations. Default is *no*.

| **-f, --input-file**\ =\ *FILE*
| Load test names from a file.

| **-#, --filenames-as-tags**
| Adds a tag for the file name.

| **-c, --section**\ =\ *STRING*
| Specify the section to run the tests on.

| **-v, --verosity**\ =\ *quiet|normal|high*
| Set output verbosity. Default is *normal*.

| **–list-test-names-only**
| List all test cases names only. If a pattern was specified, all
  matching test cases are listed.

| **–list-reporters**
| List all reporters.

| **–order**\ =\ *decl|lex|rand*
| Test case order. Default ist *decl*.

| **–rng-seed**\ =\ *time|number*
| Set a specific seed for random numbers.

| **–use-color**\ =\ *yes|no*
| Should the output be colorized? Default is *yes*.

| **–libidentify**
| Report name and version according to libidentify.

| **–wait-for-keypress**\ =\ *start|exit|both*
| Waits for a keypress before *start|exit|both*.

| **–benchmark-resolution-multiple**
| Multiple of clock resolution to run benchmarks.

Program Information
-------------------

| **-?, -h, --help**
| Print a help message for **odgi test**.

..
	EXIT STATUS
	===========
	
	| **0**
	| Success.
	
	| **1**
	| Failure (syntax or usage error; parameter error; file processing
	  failure; unexpected error).
	
	BUGS
	====
	
	Refer to the **odgi** issue tracker at
	https://github.com/pangenome/odgi/issues.
