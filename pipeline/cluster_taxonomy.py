#!/usr/bin/env python2.7

# Copyright 2018 Ian J. Miller, Evan R. Rees, Izaak Miller, Jason C. Kwan
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

# Program that determines the taxonomy of clusters called by dbscan
# Uses contig taxonomy information in the same way that add_contig_taxonomy.py uses protein taxonomy
# (except this program uses contig length as weighting rather than 1 protein 1 vote)
# Algorithm:
# In descending order of species, genus, family, order, class, phylum, superkingdom:
#    In descending order of votes:
#        If classification shares common ancestry with majority of other proteins, accept result
#    If no result, move up to next taxonomic level

import argparse
import os
import pprint
import sys
import subprocess

from time import strftime
from tqdm import tqdm


pp = pprint.PrettyPrinter(indent=4)

rank_priority = [
	'species',
	'genus',
	'family',
	'order',
	'class',
	'phylum',
	'superkingdom',
	'root',
]

canonical_ranks = {
	'superkingdom':1,
	'phylum':1,
	'class':1,
	'order':1,
	'family':1,
	'genus':1,
	'species':1,
}

def parse_names(names_dmp_path):
    names = {}
    print(strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid names')
    wc_output = subprocess.check_output(['wc', '-l', names_dmp_path])
    wc_list = wc_output.split()
    number_of_lines = int(wc_list[0])
    with open(names_dmp_path) as names_dmp:
        for line in tqdm(names_dmp, total=number_of_lines):
            taxid, name, _, classification = line.strip('\t|\n').split('\t|\t')[:4]
            taxid = int(taxid)
            # Only add scientific name entries
            scientific = classification == 'scientific name'
            if scientific:
                # line_list[1] = line_list[1].replace(' ', '_')
                names.update({taxid:name})
    return names

def parse_nodes(nodes_dmp_path):
    print(strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid nodes')
    wc_output = subprocess.check_output(['wc', '-l', nodes_dmp_path])
    wc_list = wc_output.split()
    number_of_lines = int(wc_list[0])
    nodes_dmp = open(nodes_dmp_path)
    root_line = nodes_dmp.readline()
    nodes = {}
    nodes.update({1:{'parent':1, 'rank':'root'}})
    for line in tqdm(nodes_dmp, total=number_of_lines):
        child, parent, rank = line.split('\t|\t')[:3]
        parent, child = map(int,[parent, child])
        nodes.update({child:{'parent':parent, 'rank':rank}})
    nodes_dmp.close()
    return nodes

def isConsistentWithOtherContigs(taxid, rank, clusters_taxids, nodes_dict):
    """
    Function that determines for a given taxid, whether the majority of proteins
    in a contig, with rank equal to or above the given rank, are common
    ancestors of the taxid.  If the majority are, this function returns True,
    otherwise it returns False
    clusters_taxids = {ctg:{canonical_rank:{taxid:num_hits,...},...},ctg2:{...},...}
    """
    # First we make a modified rank_priority list that only includes the current rank and above
    rank_index = rank_priority.index(rank)
    ranks_to_consider = rank_priority[rank_index:]

    # Now we total up the consistent and inconsistent ORFs
    consistent = 0
    inconsistent = 0

    for rankName in ranks_to_consider:
        if rankName not in clusters_taxids:
            continue
        for ctg_lca in clusters_taxids[rankName]:
            if isCommonAncestor(ctg_lca, taxid, nodes_dict):
                consistent += clusters_taxids[rankName][ctg_lca]
            else:
                inconsistent += clusters_taxids[rankName][ctg_lca]

    if consistent > inconsistent:
        return True
    else:
        return False

def isCommonAncestor(parent_taxid, child_taxid, nodes_dict):
    ancestor_taxid = child_taxid
    while ancestor_taxid != 1:
        if parent_taxid == ancestor_taxid:
            return True
        ancestor_taxid = nodes_dict[ancestor_taxid]['parent']
    return False

def lowest_majority(clusters_taxids, nodes_dict):
    # clusters_taxids = {canonical_rank:{taxid:num_hits, taxid2:#,...},rank2:{...},...}
    taxid_totals = {}

    for rank in rank_priority:
        if rank not in clusters_taxids:
            continue

        rank_index = rank_priority.index(rank)
        ranks_to_consider = rank_priority[rank_index:]

        for taxid in clusters_taxids[rank]:
            # Make a dictionary to total the number of canonical ranks hit
            # while traversing the path so that we can add 'unclassified' to
            # any that don't exist. Later we need to make sure that
            # 'unclassified' doesn't ever win
            ranks_in_path = {rank_to_consider:0 for rank_to_consider in ranks_to_consider}

            # We need to add to taxid_totals for each taxid in the tax_path
            current_taxid = taxid
            current_rank = rank
            while current_taxid != 1:
                if current_rank not in set(rank_priority):
                    current_taxid = nodes_dict[current_taxid]['parent']
                    current_rank = nodes_dict[current_taxid]['rank']
                    continue

                ranks_in_path[current_rank] += 1

                if current_rank not in taxid_totals:
                    taxid_totals.update({current_rank:{current_taxid:1}})
                    current_taxid = nodes_dict[current_taxid]['parent']
                    current_rank = nodes_dict[current_taxid]['rank']
                    continue

                if current_taxid in taxid_totals[current_rank]:
                    taxid_totals[current_rank][current_taxid] += 1
                else:
                    taxid_totals[current_rank][current_taxid] = 1

                current_taxid = nodes_dict[current_taxid]['parent']
                current_rank = nodes_dict[current_taxid]['rank']

            # Now go through ranks_in_path. Where total = 0, add 'unclassified'
            for rank_to_consider in ranks_to_consider:
                if ranks_in_path[rank_to_consider] == 0:
                    if rank_to_consider not in taxid_totals:
                        taxid_totals[rank_to_consider] = {'unclassified':1}
                        continue
                    if 'unclassified' in taxid_totals[rank_to_consider]:
                        taxid_totals[rank_to_consider]['unclassified'] += 1
                    else:
                        taxid_totals[rank_to_consider]['unclassified'] = 1
    # If there are any gaps in the taxonomy paths for any of the contigs in the cluster,
    # we need to add 'unclassified' to the relevant canonical taxonomic rank.
    # However, we must never allow 'unclassified' to win! (That just won't really tell us anything)
    # Now we need to determine which is the first level to have a majority
    for rank in rank_priority:
        total_votes = 0
        taxid_leader = None
        taxid_leader_votes = 0
        if not rank in taxid_totals:
            continue
        for taxid in taxid_totals[rank]:
            taxid_votes = taxid_totals[rank][taxid]
            total_votes += taxid_votes
            if taxid_votes > taxid_leader_votes:
                taxid_leader = taxid
                taxid_leader_votes = taxid_votes
        majority_threshold = float(total_votes)/2
        if taxid_leader_votes > majority_threshold and taxid_leader != 'unclassified':
            return taxid_leader
    # Just in case
    return 1

def parse_taxonomy(infpath, nodes_dict):

	print strftime("%Y-%m-%d %H:%M:%S") + ' Parsing taxonomy table'

	# Work out number of lines in file
	wc_output = subprocess.check_output(['wc', '-l', infpath])
	wc_list = wc_output.split()
	number_of_lines = int(wc_list[0])

	# Determine contig, length and cluster indexes
	fh = open(infpath, 'r')
	header_cols = fh.readline().split('\t')
	assert len(header_cols) == len(set(header_cols)), 'DuplicateColumnsInTable'
	cols = ['contig', 'length', 'taxid', 'cluster']
	if single_genome_mode:
		cols.remove('cluster')
	indices = {}
	for col in cols:
		try:
			indices.update({col:header_cols.index(col)})
		except ValueError as err:
			print('ColumnNotFound: {} not in {}'.format(col, infpath))
			sys.exit(2)

	clusters = {}
	taxid_index = indices['taxid']
	for line in tqdm(fh, total=number_of_lines):
		line_list = line.rstrip('\n').split('\t')

		if single_genome_mode:
			cluster = 'unclustered'
		else:
			cluster_index = indices['cluster']
			cluster = line_list[cluster_index]

		taxid = int(line_list[taxid_index])

		taxRank = None
		if taxid == 1:
			taxRank = 'root'
		elif taxid < 0:
			# Unclassified taxid == -1
			taxRank = 'no rank'
		elif taxid not in nodes_dict:
			# This happens sometimes when the taxid database and the NR database are not in sync
			print('Warning: {} not found in nodes.dmp.'.format(taxid))
			print('Consider updating your databases')
			taxid = 1
			taxRank = 'root'
		else:
			taxRank = nodes_dict[taxid]['rank']

		# Now get the taxid of the next canonical rank (if applicable)
		if taxRank == 'no rank':
			taxid = 1
			taxRank = 'root'

		if taxid != 1:
			while taxRank not in rank_priority:
				taxid = nodes_dict[taxid]['parent']
				taxRank = nodes_dict[taxid]['rank']

		# Keep running total of taxids for each cluster
		if cluster not in clusters:
			clusters.update({cluster:{taxRank:{taxid:1}}})
			continue

		if taxRank not in clusters[cluster]:
			clusters[cluster].update({taxRank:{taxid:1}})
			continue

		if taxid not in clusters[cluster][taxRank]:
			clusters[cluster][taxRank].update({taxid:1})
		else:
			clusters[cluster][taxRank][taxid] += 1

	fh.close()
	return clusters

def rank_taxids(clusters_taxids, nodes_dict):
	print strftime("%Y-%m-%d %H:%M:%S") + ' Ranking taxids'
	top_taxids = {}
	total_clusters = len(clusters_taxids)

	for cluster in tqdm(clusters_taxids, total=total_clusters):
		acceptedTaxid = None
		for rank in rank_priority:
			if acceptedTaxid:
				break
			# Order in descending order of votes
			if rank in clusters_taxids[cluster]:
				ordered_taxids = sorted(clusters_taxids[cluster][rank], key=clusters_taxids[cluster][rank].__getitem__, reverse=True)
				#sys.exit()
				for taxid in ordered_taxids:
					if isConsistentWithOtherContigs(taxid, rank, clusters_taxids[cluster], nodes_dict):
						acceptedTaxid = taxid
						break

		# If acceptedTaxid is still None at this point, there was some kind of draw, so we need to find the lowest taxonomic level where there is a
		# majority
		if not acceptedTaxid:
			acceptedTaxid = lowest_majority(clusters_taxids[cluster], nodes_dict)

		top_taxids[cluster] = acceptedTaxid

	return top_taxids

def resolve_taxon_paths(ctg2taxid, names_dict, nodes_dict):
    print(strftime("%Y-%m-%d %H:%M:%S") + ' Resolving taxon paths')
    contig_paths = {}
    # {contig:{rank1:name1,rank2,name2},contig2:{rank1:name1,rank2:name2,...},...}
    n_contigs = len(ctg2taxid)
    for contig in tqdm(ctg2taxid, total=n_contigs):
        taxid = ctg2taxid[contig]
        if taxid == 1:
            contig_paths.update({contig:{'root':names_dict[taxid]}})
        while taxid != 1:
            current_rank = nodes_dict[taxid]['rank']
            if current_rank not in set(rank_priority):
                taxid = nodes_dict[taxid]['parent']
                continue

            name = names_dict[taxid]

            if contig not in contig_paths:
                contig_paths.update({contig:{current_rank:name}})
            else:
                contig_paths[contig][current_rank] = name
            taxid = nodes_dict[taxid]['parent']

        for rank in rank_priority:
            if rank not in contig_paths[contig]:
                contig_paths[contig][rank] = 'unclassified'

    for contig in contig_paths:
        contig_paths[contig].pop('root')
        contig_paths[contig]['taxid'] = ctg2taxid[contig]

    return contig_paths

def write_cluster_taxonomy(taxids, clusters_taxids, outfpath):
	print strftime("%Y-%m-%d %H:%M:%S") + ' Writing table'

	output_table = open(outfpath, 'w')
	header_cols = [
		'cluster',
		'kingdom',
		'phylum',
		'class',
		'order',
		'family',
		'genus',
		'species',
		'taxid',
	]
	if single_genome_mode:
		header_cols.remove('cluster')

	header = '\t'.join(header_cols) + '\n'
	output_table.write(header)

	for cluster in taxids:
		cols = map(str,[
			cluster,
			taxids[cluster]['superkingdom'],
			taxids[cluster]['phylum'],
			taxids[cluster]['class'],
			taxids[cluster]['order'],
			taxids[cluster]['family'],
			taxids[cluster]['genus'],
			taxids[cluster]['species'],
			clusters_taxids[cluster],
		])
		if single_genome_mode:
			cols.remove(cluster)
		line = '\t'.join(cols) + '\n'
		output_table.write(line)

	output_table.close()

parser = argparse.ArgumentParser(
	description='Summarize the taxonomy of clusters in a table.'
	' Uses same modified majority voting algorithm as contig '
	'taxonomy addition.'
)
parser.add_argument(
	'-i',
	'--contig_tab',
	help='</path/to/autometa_output.tsv>',
	required=True,
)
parser.add_argument(
	'-db',
	'--taxdump',
	help='</path/to/taxdump/dir>'
	' (Needs names.dmp & nodes.dmp downloaded from NCBI taxonomy)',
	required=True,
)
parser.add_argument('-o', '--out', help='</path/to/output.tsv>', required=True)
parser.add_argument(
	'-c',
	'--cluster_column',
	help='Cluster column name',
	default='cluster',
)
parser.add_argument(
	'-s',
	'--single_genome',
	help='Specifies single genome mode',
	action='store_true',
)

args = vars(parser.parse_args())

contig_table_path = args['contig_tab']
taxdump_dir_path = args['taxdump']
output_file_path = args['out']
cluster_column_name = args['cluster_column']
single_genome_mode = args['single_genome']

# Process NCBI taxdump files
names_dmp_path = os.path.join(taxdump_dir_path, 'names.dmp')
nodes_dmp_path = os.path.join(taxdump_dir_path, 'nodes.dmp')

names = parse_names(names_dmp_path)
nodes = parse_nodes(nodes_dmp_path)

clusters_classifications = parse_taxonomy(contig_table_path, nodes)

cluster_taxids = rank_taxids(clusters_classifications, nodes)

resolved_taxids = resolve_taxon_paths(cluster_taxids, names, nodes)

write_cluster_taxonomy(resolved_taxids, cluster_taxids, output_file_path)

print('Written: {}'.format(output_file_path))
