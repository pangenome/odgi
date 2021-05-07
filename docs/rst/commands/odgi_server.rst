.. _odgi server:

#########
odgi server
#########

start a HTTP server with a given index file to query a
pangenome position

SYNOPSIS
========

**odgi server** [**-i, –idx**\ =\ *FILE*] [**-p, –port**\ =\ *N*]
[*OPTION*]…

DESCRIPTION
===========

| The odgi server(1) command starts an HTTP server with a given path
  index as input. The idea is that we can go from **path:position** →
  **pangenome:position** via GET requests to the HTTP server. The server
  headers do not block cross origin requests. Example GET request:
  **http://localost:3000/path_name/nucleotide_position**.
| The required path index can be created with :ref:`odgi pathindex`. Going from
  **path:position** → **pangenome:position** is important when
  navigating large graphs in an interactive manner like in the
  `Pantograph <https://graph-genome.github.io/>`__ project. All input
  and output positions are 1-based. If no IP address is specified, the
  server will run on localhost.

OPTIONS
=======

Graph Files IO
--------------

| **-i, –idx**\ =\ *FILE*
| File containing the succinct variation graph index to host in a HTTP
  server. The file name usually ends with *.xp*.

HTTP Options
------------

| **-p, –port**\ =\ *N*
| Run the server under this port.

| **-a, –ip**\ =\ *IP*
| Run the server under this IP address. If not specified, *IP* will be
  *localhost*.

Program Information
-------------------

| **-h, –help**
| Print a help message for **odgi server**.

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
