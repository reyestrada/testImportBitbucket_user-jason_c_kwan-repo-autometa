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
from itertools import combinations
import operator
from scipy import stats
import csv
import numbers
import pdb
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
					elif node in repeat_contigs:
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

def shouldMerge(donor_bin, host_bin):
	# First we make arrays of coverages and gc for each bin, for t-test
	donor_bin_coverages = list()
	donor_bin_bhtsne_x = list()
	donor_bin_bhtsne_y = list()
	host_bin_coverages = list()
	host_bin_bhtsne_x = list()
	host_bin_bhtsne_y = list()

	for i,row in master_table.iterrows():
		contig = row['contig']
		cov = row['cov']
		bh_tsne_x = row['bh_tsne_x']
		bh_tsne_y = row['bh_tsne_y']

		if contig in bfs_sets_simple[donor_bin]:
			donor_bin_coverages.append(cov)
			donor_bin_bhtsne_x.append(bh_tsne_x)
			donor_bin_bhtsne_y.append(bh_tsne_y)

		if contig in bfs_sets_simple[host_bin]:
			host_bin_coverages.append(cov)
			host_bin_bhtsne_x.append(bh_tsne_x)
			host_bin_bhtsne_y.append(bh_tsne_y)

	( cov_t_statistic, cov_p_value ) = stats.ttest_ind(donor_bin_coverages, host_bin_coverages, equal_var=False)
	( bhtsne_x_statistic, bhtsne_x_p_value ) = stats.ttest_ind(donor_bin_bhtsne_x, host_bin_bhtsne_x, equal_var=False)
	( bhtsne_y_statistic, bhtsne_y_p_value ) = stats.ttest_ind(donor_bin_bhtsne_y, host_bin_bhtsne_y, equal_var=False)

	if cov_p_value < 0.05 or bhtsne_x_p_value < 0.05 or bhtsne_y_p_value < 0.05:
		return False

	# Now we check whether merging would give a bin >95% pure
	combined_contig_set = bfs_sets_simple[donor_bin].union(bfs_sets_simple[host_bin])

	pfam_counts = dict()
	for i,row in master_table.iterrows():
		contig = row['contig']
		if contig + 's' in combined_contig_set or contig + 'e' in combined_contig_set:
			if row['single_copy_PFAMs'] == 'NA' or (isinstance(row['single_copy_PFAMs'], numbers.Number) and math.isnan(row['single_copy_PFAMs'])):
				continue
			pfam_list = row['single_copy_PFAMs'].split(',')
			for pfam in pfam_list:
				if pfam in pfam_counts:
					pfam_counts[pfam] += 1
				else:
					pfam_counts[pfam] = 1
	number_unique_markers = 0
	for pfam in pfam_counts:
		if pfam_counts[pfam] == 1:
			number_unique_markers += 1
	number_markers_found = len(pfam_counts.keys())
	if number_markers_found == 0:
		purity = 100
	else:
		purity = (number_unique_markers / number_markers_found) * 100
	if purity > 95.0:
		return True
	else:
		return False


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
columns_to_check = [cluster_column_heading, 'contig', 'length', 'cov', 'gc', 'bh_tsne_x', 'bh_tsne_y', 'single_copy_PFAMs']

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
contig_lengths = dict()
for i,row in master_table.iterrows():
	contig = row['contig']
	bin_name = row[cluster_column_heading]
	bin_lookup[contig] = bin_name
	length = row['length']
	gc = row['gc']
	cov = row['cov']
	contig_lengths[contig] = length
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

# Make a new dataset with simple contig names (i.e. without 's' and 'e' on the end)
bfs_sets_simple = dict()

for bin_name in bfs_sets:
	if bin_name == 'unclustered': # We assume all other bins have better claim to contigs in unclustered where there is overlap
		continue
	bfs_sets_simple[bin_name] = set()
	for contig_name in bfs_sets[bin_name]:
		bfs_sets_simple[bin_name].add(contig_name[:-1])

# Now we merge bins, if there is >10% of contigs overlapping, and the coverage, bh_tsne_x and bh_tsne_y is not significantly
# different (t-test), and only if the combination is >95% pure.
notFinished = True
while(notFinished):
	ordered_bin_list = bfs_sets_simple.keys()
	overlap = dict() # Dictionary of dictionaries. Keyed by bin, holds dictionary of the number of contigs shared with other bins
	for combination_tuple in combinations(ordered_bin_list, 2):
		number_shared_contigs = len(bfs_sets_simple[combination_tuple[0]].intersection(bfs_sets_simple[combination_tuple[1]]))
		if number_shared_contigs == 0:
			continue
		if combination_tuple[0] in overlap:
			overlap[combination_tuple[0]][combination_tuple[1]] = number_shared_contigs
		else:
			overlap[combination_tuple[0]] = { combination_tuple[1]: number_shared_contigs }

		if combination_tuple[1] in overlap:
			overlap[combination_tuple[1]][combination_tuple[0]] = number_shared_contigs
		else:
			overlap[combination_tuple[1]] = { combination_tuple[0]: number_shared_contigs }

	overlap_percents = dict()
	for source_bin in overlap:
		source_bin_length = 0
		for contig in bfs_sets_simple[source_bin]:
			source_bin_length += contig_lengths[contig]
		overlap_percents[source_bin] = {}
		for target_bin in overlap[source_bin]:
			target_overlap_length = 0
			for contig in bfs_sets_simple[source_bin].intersection(bfs_sets_simple[target_bin]):
				target_overlap_length += contig_lengths[contig]
			percent = (target_overlap_length / source_bin_length) * 100
			overlap_percents[source_bin][target_bin] = percent

	# Order the bins in descending order of highest overlap, to go through and see which ones should be 
	# merged
	biggest_overlaps = dict()
	biggest_overlap_target_bins = dict()
	for bin_name in overlap:
		largest_local_overlap = 0
		largest_target_bin = None
		for target_bin in overlap_percents[bin_name]:
			if overlap_percents[bin_name][target_bin] > largest_local_overlap:
				largest_local_overlap = overlap_percents[bin_name][target_bin]
				largest_target_bin = target_bin
		biggest_overlaps[bin_name] = largest_local_overlap
		biggest_overlap_target_bins[bin_name] = largest_target_bin

	sorted_bin_tuples = sorted(biggest_overlaps.items(), key=operator.itemgetter(1), reverse=True)

	# Go through sorted_bin_tuples, testing whether the overlapping pairs should be merged
	# If so, we stop and continue the while loop after merging
	merged = 0
	for bin_tuple in sorted_bin_tuples:
		source_bin = bin_tuple[0]
		target_bin = biggest_overlap_target_bins[source_bin]
		if bin_tuple[1] > 10 and shouldMerge(source_bin, target_bin):
			print('Merging ' + source_bin + ' into ' + target_bin)
			bfs_sets_simple[target_bin] = bfs_sets_simple[target_bin].union(bfs_sets_simple[source_bin])
			del bfs_sets_simple[source_bin]
			merged += 1
			break

	if merged == 0:
		notFinished = False

# Now we identify contigs that are shared between more than one bfs set, and delete them from the sets
contigs_to_delete = set()
ordered_bin_list = bfs_sets_simple.keys()
for combination_tuple in combinations(ordered_bin_list, 2):
	shared_contigs = bfs_sets_simple[combination_tuple[0]].intersection(bfs_sets_simple[combination_tuple[1]])
	for contig in shared_contigs:
		contigs_to_delete.add(contig)

for contig in contigs_to_delete:
	for bin_name in bfs_sets_simple:
		if contig in bfs_sets_simple[bin_name]:
			bfs_sets_simple[bin_name].remove(contig)

# Make lookup data structure
paired_end_refined_bins = dict() # Keyed by contig
for bin_name in bfs_sets_simple:
	for contig in bfs_sets_simple[bin_name]:
		paired_end_refined_bins[contig] = bin_name

# Now we write a new column to the master table
new_column = list()
for i,row in master_table.iterrows():
	contig = row['contig']
	if contig in paired_end_refined_bins:
		new_column.append(paired_end_refined_bins[contig])
	else:
		new_column.append('unclustered')

master_table['paired_end_refined_bin'] = pd.Series(new_column, index = master_table.index)

output_table_path = os.path.join(output_dir, 'paired_end_refine_bins_output.tab')

master_table.to_csv(path_or_buf=output_table_path, sep='\t', index=False, quoting=csv.QUOTE_NONE)

# Now output subset graphs using the new bins
filehandles = dict()
for bin_name in bfs_sets_simple:
	file_path = os.path.join(output_dir, bin_name + '.connections.tab')
	filehandles[bin_name] = open(file_path, 'w')

with open(graph_file_path) as graph_input:
	for line in graph_input:
		if i == 0:
			for bin_name in filehandles:
				filehandles[bin_name].write(line)
		else:
			line_list = line.rstrip().split('\t')
			contig1 = line_list[0][:-1]
			contig2 = line_list[2][:-1]
			if contig1 in paired_end_refined_bins and contig2 in paired_end_refined_bins:
				if paired_end_refined_bins[contig1] == paired_end_refined_bins[contig2]:
					bin_name = paired_end_refined_bins[contig1]
					filehandles[bin_name].write(line)

for bin_name in filehandles:
	filehandles[bin_name].close()