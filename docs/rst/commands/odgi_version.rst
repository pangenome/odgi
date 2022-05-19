.. _odgi version:

#########
odgi version
#########

Print the version of ODGI to stdout.

SYNOPSIS
========

**odgi version** [*OPTION*]…

DESCRIPTION
===========

The odgi version command prints the current git version with tags and
codename to stdout (like *v-44-g89d022b “back to old ABI”*). Optionally,
only the release, version or codename can be printed.

OPTIONS
=======

Version Options
---------------

| **-v, --version**\ =
| Print only the version (like *v0.4.0-44-g89d022b*).

| **-c, --codename**
| Print only the codename (like *back to old ABI*).

| **-r, --release**
| Print only the release (like *v0.4.0*).

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi version**.

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
