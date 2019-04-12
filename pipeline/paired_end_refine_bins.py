#!/usr/bin/env python2.7

# Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan
#
# This file is part of Autometa.
#
# Autometa is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Autometa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Autometa. If not, see <http://www.gnu.org/licenses/>.

# This program uses paired-end read information to refine bins produced by Autometa

from __future__ import division
import numpy as np
import math
import argparse
import os
import pandas as pd

# See https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
def weighted_av_and_stdev(values, weights):
	value_array = np.asarray(values)
	weight_array = np.asarray(weights)
	average = np.average(value_array, weights=weight_array)
	variance = np.average((value_array - average)**2, weights=weight_array)
	return(average, math.sqrt(variance))

def bfs(graph,start,bin_designation):
	# Expects graph to be a dictionary of dictionaries, and start to be a list of starting nodes
	# Note starting nodes are bare contig names, and they are appended with both 's' and 'e' in the graph

	# keep track of all visited nodes
	explored = []
	# keep track of nodes to be checked
	queue = []
	for contig in start:
		queue.append(contig + 's')
		queue.append(contig + 'e')

	# keep looping until there are no nodes still to be checked
	while queue:
		# pop shallowest node (first node) from queue
		node = queue.pop(0)
		if node not in explored:
			# add node to list of checked nodes
			explored.append(node)
			if node in graph:
				neighbors = []
				for neighbor in graph[node]:
					if neighbor in bin_sets[bin_designation]:
						neighbors.append(neighbor)
					elif node[:-1] == neighbor[:-1]:
						# We keep any connection between s and e of the same contig
						neighbors.append(neighbor)
					elif neighbor in repeat_contigs:
						continue
					elif bin_designation in bin_stats:
						# Work out z score
						if bin_stats[bin_designation]['connections_stdev'] == 0:
							continue
						else:
							z_score = (graph[node][neighbor] - bin_stats[bin_designation]['connections_mean']) / bin_stats[bin_designation]['connections_stdev']
							if z_score <= 3:
								neighbors.append(neighbor)
			else:
				neighbors = []

			# add neighbors of node to queue
			for neighbor in neighbors:
				queue.append(neighbor)

	return set(explored)

parser = argparse.ArgumentParser(description='Script to refine bins made by run_autometa.py or ML_recruitment.py using information from paired-end read alignment')
parser.add_argument('-b', '--bin_table', metavar='<bin.tab>', help='path to the output from either run_autometa.py or ML_recruitment.py', required=True)
parser.add_argument('-c', '--column', metavar='<bin column name>', help='the name of the column to use for binning purposes', default='cluster')
parser.add_argument('-o', '--output_dir', metavar='<dir>', help='path to the directory where output files will go', default='.')
parser.add_argument('-g', '--graph_file', metavar='cytoscape.connections.tab', help='Output from cytoscapeviz.pl (http://madsalbertsen.github.io/multi-metagenome/), it is recommended that you use the error corrected reads made by metaSPAdes to align to your contig, and use the following cytoscapeviz.pl parameters -f 2 -e 500 -m 3000 -a <read size>', required=True)

args = vars(parser.parse_args())

bin_table_path = args['bin_table']
cluster_column_heading = args['column']
output_dir = args['output_dir']
graph_file_path = args['graph_file']

# Check paths exist
if not os.path.isfile(bin_table_path):
	print('Error! Could not find a bin table at the following path: ' + bin_table_path)
	exit(1)

if not os.path.isfile(graph_file_path):
	print('Error! Cannot find a graph file at the following path: ' + graph_file_path)
	exit(1)

# Make output directory if it doesn't exist
if not os.path.isdir(output_dir):
	os.makedirs(output_dir)

master_table = pd.read_table(bin_table_path)

# Format check for the table
columns_to_check = [cluster_column_heading, 'contig', 'length', 'cov', 'gc', 'single_copy_PFAMs']

for column in columns_to_check:
	if column not in master_table.columns:
		print('Error! Could not find a column called ' + column + ' in table ' + bin_table_path)
		exit(1)

# First we make sets that contain the contigs in each bin
bin_sets = dict() # Dictionary of sets
bin_lookup = dict() # Dictionary of contigs, holds bin
bin_lengths = dict() # Dictionary keyed by bin, holds lists of contig lengths
bin_gc = dict() # Dictionary keyed by bin, holds lists of contig gc percents
bin_cov = dict() # Dictionary keyed by bin, holds lists of contig coverages
for i,row in master_table.iterrows():
	contig = row['contig']
	bin_name = row[cluster_column_heading]
	bin_lookup[contig] = bin_name
	length = row['length']
	gc = row['gc']
	cov = row['cov']
	if bin_name in bin_sets:
		bin_sets[bin_name].add(contig)
		bin_lengths[bin_name].append(length)
		bin_gc[bin_name].append(gc)
		bin_cov[bin_name].append(cov)
	else:
		bin_sets[bin_name] = set([contig])
		bin_lengths[bin_name] = [ length ]
		bin_gc[bin_name] = [ gc ]
		bin_cov[bin_name] = [ cov ]

bin_stats = dict() # Keyed by bin_name, holds dictionarys

# Calculate weighted averages and stdevs for cov and gc for each bin
for bin_name in bin_lengths:
	weights = list()
	for length in bin_lengths[bin_name]:
		weights.append(length / sum(bin_lengths[bin_name]))
	( weighted_av_gc, gc_stdev ) = weighted_av_and_stdev(bin_gc[bin_name], weights)
	( weighted_av_cov, cov_stdev ) = weighted_av_and_stdev(bin_cov[bin_name], weights)
	bin_stats[bin_name] = { 'weighted_av_gc': weighted_av_gc, 'gc_sdpc': (gc_stdev / weighted_av_gc) * 100, 'weighted_av_cov': weighted_av_cov, 'cov_sdpc': (cov_stdev / weighted_av_cov) * 100 }

# Now make graph data structure
connection_graph = dict() # Dictionary of dictionaries (stores number of read connections)
with open(graph_file_path) as graph_input:
	for i,line in enumerate(graph_input):
		if i > 0:
			line_list = line.rstrip().split('\t')
			contig1 = line_list[0]
			contig2 = line_list[2]
			read_connections = int(line_list[3])

			# Adjust read_connections if it is an intracontig connection
			if contig1[:-1] == contig2[:-1] and len(line_list) > 5:
				read_connections = int(line_list[6])

			if contig1 in connection_graph:
				connection_graph[contig1][contig2] = read_connections
			else:
				connection_graph[contig1] = { contig2: read_connections }

			if contig2 in connection_graph:
				connection_graph[contig2][contig1] = read_connections
			else:
				connection_graph[contig2] = { contig1: read_connections }

# Now we work out which contigs are repeats (defined as having more than one two connections)
repeat_contigs = set()
for source_contig in connection_graph:
	if len(connection_graph[source_contig]) > 2:
		repeat_contigs.add(source_contig)

# Now for each bin we want to calculate the mean number of intrabin connections, and their standard deviation
for bin_name in bin_sets:
	read_connections_list = []
	for source_contig in bin_sets[bin_name]:
		for source_contig_name in [ source_contig + 's', source_contig + 'e' ]:
			if source_contig_name in connection_graph:
				for target_contig in connection_graph[source_contig_name]:
					if target_contig[:-1] in bin_sets[bin_name]:
						read_connections_list.append(connection_graph[source_contig_name][target_contig])
	if len(read_connections_list) > 0:
		bin_stats[bin_name]['connections_mean'] = np.mean(read_connections_list)
		bin_stats[bin_name]['connections_stdev'] = np.std(read_connections_list)

# We do a BFS, but only count connections which have a z score <= 3 with respect to the bin of interest
# In other words, we only follow a connection if it is not an outlier in terms of number of reads
# Z score: (https://medium.com/datadriveninvestor/finding-outliers-in-dataset-using-python-efc3fce6ce32)

bfs_sets = dict() # Dictionary, holds sets of contigs by bin

for bin_name in bin_sets:
	contig_list = list(bin_sets[bin_name])
	bfs_sets[bin_name] = bfs(connection_graph, contig_list, bin_name)

# Measure how heterogenous each bfs_set is
for bin_name in bfs_sets:
	bin_dict = dict()
	for contig in bfs_sets[bin_name]:
		contig_bin = bin_lookup[contig[:-1]]
		if contig_bin in bin_dict:
			bin_dict[contig_bin] += 1
		else:
			bin_dict[contig_bin] = 1
	output_string_list = list()
	for contig_bin in bin_dict:
		percent = (bin_dict[contig_bin] / len(bfs_sets[bin_name])) * 100
		output_string = str(percent) + '% ' + contig_bin
		output_string_list.append(output_string)

	output_string = bin_name + ' BFS set: ' + ','.join(output_string_list)
	print(output_string)

# Now for each bin, write a subset of the connection graph
filehandles = dict() # binned by bin name
for bin_name in bfs_sets:
	filepath = os.path.join(output_dir, bin_name + '.connections.tab')
	filehandles[bin_name] = open(filepath, 'w')
#pdb.set_trace()
with open(graph_file_path) as graph_input:
	for i,line in enumerate(graph_input):
		if i == 0:
			for bin_name in filehandles:
				filehandles[bin_name].write(line)
		else:
			line_list = line.rstrip().split('\t')
			contig1 = line_list[0]
			contig2 = line_list[2]

			for bin_name in bfs_sets:
				if (contig1 in bfs_sets[bin_name]) or (contig2 in bfs_sets[bin_name]):
					filehandles[bin_name].write(line)

for bin_name in filehandles:
	filehandles[bin_name].close()
