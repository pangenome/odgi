
# This script convert an adjacency lists from https://github.com/jxz12/s_gd2 in a GFA file. The file is printed in
# the standard output.

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--adjacency-list',type=str, help='text file with the adjacency list', required=True)

args = parser.parse_args()


node_id_to_neighborsIds_dict = dict()

with open(args.adjacency_list) as f:
	for line in f:
		node_i, node_j = [int(x) for x in line.strip('\n').split(' ')]

		if node_i not in node_id_to_neighborsIds_dict:
			node_id_to_neighborsIds_dict[node_i] = set()
		node_id_to_neighborsIds_dict[node_i].add(node_j)


nucleotides = 'ATCG'

node_id_to_write_set = set(node_id_to_neighborsIds_dict.keys())
for x_set in node_id_to_neighborsIds_dict.values():
	node_id_to_write_set.update(x_set)


print('\t'.join(['H', 'VN:Z:1.0']))
for node_id in sorted(node_id_to_write_set):
	print('\t'.join(['S', str(node_id), nucleotides[(node_id - 1) % 4], 'DP:i:0', 'RC:i:0']))

	if node_id in node_id_to_neighborsIds_dict:
		for neighbor_id in sorted(node_id_to_neighborsIds_dict[node_id]):
			print('\t'.join(['L', str(node_id), '+', str(neighbor_id), '+', '0M']))