#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file shows some code examples of how to use the python bindings to odgi.
# See the odgi documentation for more information on the python bindings.

import sys
if len(sys.argv) == 1:
    print("usage: PYTHONPATH=../lib", sys.argv[0], "graph.og")
    print("displays information about the graph and iterates over it")
    print("to build an example input graph, use `odgi build -g graph.gfa -o graph.og`")
    exit(1)

import odgi

g = odgi.graph()
g.load(sys.argv[1])

# the number of nodes and edges is given
print("node count:", g.get_node_count())
print("path count:", g.get_path_count())

# iterate over all nodes and sum the sequence length
# we have to use a generic counter object to count
# because in python you can't assign inside lambdas 🤦
# and there are no true closures over local variables 😭
class counter:
    total = 0
    def add(self, l):
        self.total += l

s = counter()
g.for_each_handle(lambda h: s.add(g.get_length(h)))
print("graph length:", s.total)

# iterate over the handles and use follow_edges to list the edges
def show_edge(a, b):
    print(g.get_id(a), "-->", g.get_id(b))

def display_node_edges(h):
    print("node", g.get_id(h))
    g.follow_edges(
        h, False,
        lambda n:
        show_edge(h, n))
    g.follow_edges(
        h, True,
        lambda n:
        show_edge(n, h))

# displays all the edges twice, once for each of their ends
g.for_each_handle(display_node_edges)

# count all the edges
e = counter()
def handle_edges(h):
    g.follow_edges(h, True, lambda x: e.add(1))

g.for_each_handle(handle_edges)
print("edge count", e.total)
print("edges/node", float(e.total) / float(g.get_node_count()))
print("bp/edge", float(s.total) / float(e.total))

# getting path names from path handles
path_names = []
g.for_each_path_handle(lambda p: path_names.append(g.get_path_name(p)))
print(path_names)

# getting the path handle from a path name
print([g.get_path_handle(name) for name in path_names])

# print a string representing each step
# this shows the path and orientation relative to the handle
def step_str(step):
    #id_str = str(g.get_id(g.get_handle_of_step(step)))
    path_str = g.get_path_name(g.get_path_handle_of_step(step))
    dir_str = "+" if not g.get_is_reverse(g.get_handle_of_step(step)) else "-"
    return path_str + dir_str

# getting the steps on each handle
def show_steps(handle):
    steps = []
    g.for_each_step_on_handle(handle, lambda step: steps.append(step))
    print("node:", g.get_id(handle), "steps:", " ".join([step_str(s) for s in steps]))

g.for_each_handle(show_steps)
