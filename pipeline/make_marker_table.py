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

import pandas as pd
import argparse
import subprocess
import os

#argument parser
parser = argparse.ArgumentParser(description='Script to tabulate single copy markers \
	from a metagenome assembly. Dependencies: prodigal v2.6.2 (from "GoogleImport" branch), hmmscan (hmmer 3.1b2)')
parser.add_argument('-a','--assembly', help='Input assembly file', required=True)
parser.add_argument('-p','--processors', help='Number of processors to use for hmmscan', default=1)
parser.add_argument('-c','--cutoffs', help='Bacterial single copy hmm cutoffs as defined by Rinke et al. Default path is home directory.', default="~/Bacteria_single_copy_cutoffs.txt")
parser.add_argument('-m','--hmm', help='Bacteria_single_copy_cutoffs.hmm. Default path is home directory.', default="~/Bacteria_single_copy.hmm")
parser.add_argument('-o','--out', help='outfile.tab, three column table with contig, single copy PFAMS, and # of markers', required=False)
args = vars(parser.parse_args())

assembly = os.path.abspath(args['assembly'])

def get_contig_list(path_to_assembly):
    #Get list of spades contigs
    assembly_handle = open(path_to_assembly,"rU")
    contig_name_list = []
    for line in assembly_handle:
        if '>' in line:
            contig_name = line.rstrip("\n").split()[0][1:]
            contig_name_list.append(contig_name)
    assembly_handle.close()
    return contig_name_list

def run_prodigal(path_to_assembly):
	assembly_filename = os.path.splitext(os.path.basename(path_to_assembly))[0]
	orfs_fpath = os.path.join(output_dir, assembly_filename+'.orfs.faa')
	prodigal_txt = os.path.join(output_dir, assembly_filename+'.txt')
	if os.path.isfile(orfs_fpath):
		print('{} already exists...continuing'.format(orfs_fpath))
		return orfs_fpath
	#When "shell = True", need to give one string, not a list
	cmd = ' '.join([
		'prodigal',
		'-i', path_to_assembly,
		'-a', orfs_fpath,
		'-p meta',
		'-m',
		'-o', prodigal_txt
	])
	subprocess.call(cmd, shell = True)
	return orfs_fpath

def run_hmmscan(orfs_fp,hmmdb):
	hmm_outfpath = orfs_fp+'.hmm.tbl'
	cmd = ' '.join(map(str,[
		'hmmscan',
		'--cpu',args['processors'],
		'--tblout',hmm_outfpath,
		hmmdb,
		orfs_fp,
	]))
	subprocess.call(cmd, shell = True)
	return hmm_outfpath

output_dir = os.path.dirname(os.path.realpath(args['out']))

orfs_fpath = run_prodigal(assembly)
hmm_table_path = run_hmmscan(orfs_fpath, args['hmm'])

hmm_table = pd.read_csv(
	hmm_table_path,
	sep='\s+',
	usecols=[1, 2, 5],
	skiprows=3,
	header=None,
	index_col=False,
	engine='python',
)
cutoffs_table = pd.read_csv(
	args['cutoffs'],
	sep='\s',
	engine='python',
	header=None,
)

#Search for contigs/ORFs that contain single copy PFAM domains that pass cutoff
#for loop to search for PFAM domains in hmm table column 1:
contig_ORFs_that_pass_cutoffs = {}
for index,PFAM_cutoffs_id in enumerate(cutoffs_table[0]):
	for count,PFAM_hmm_scan_id in enumerate(hmm_table[1]):
		#If the PFAMs domain match (cutoffs ID names format are PFAMXXXXX, not PFAMXXXXX.1)
		if str(PFAM_cutoffs_id) in str(PFAM_hmm_scan_id):
			#and hmm score value > cutoff value,
			if float(hmm_table[5][count]) > float(cutoffs_table[1][index]):
				#make dictionary of lists for PFAMs, where PFAM is key and contigs populate the list
				if PFAM_hmm_scan_id not in contig_ORFs_that_pass_cutoffs:
					contig_ORFs_that_pass_cutoffs[PFAM_hmm_scan_id] = []
					contig_ORFs_that_pass_cutoffs[PFAM_hmm_scan_id].append(hmm_table[2][count])
				else:
					contig_ORFs_that_pass_cutoffs[PFAM_hmm_scan_id].append(hmm_table[2][count])

#Make a dictionary of dictionaries with contigs that have single copy genes (list PFAMS, for final table count length of list),
#contig length, contig GC, contig len, passecd PFAM domains
#write out tab-delimited table

header = 'contig\tsingle_copy_PFAMs\tnum_single_copies\n'
lines = header
contig_dictionary = {}
for count,contig in enumerate(get_contig_list(assembly)):
	contig_dictionary[contig] = {}
	contig_dictionary[contig]['single_copy_PFAMs'] = []
	contig_dictionary[contig]['num_single_copies'] = 0
	for PFAM_key,contigs in contig_ORFs_that_pass_cutoffs.items():
		for item in contigs:
			if str(contig) in item:
				contig_dictionary[contig]['single_copy_PFAMs'].append(PFAM_key)
	num_copies = len(contig_dictionary[contig]['single_copy_PFAMs'])
	if num_copies > 0:
		pfam_list = ','.join(contig_dictionary[contig]['single_copy_PFAMs'])
		lines += '\t'.join(map(str,[contig,pfam_list,num_copies]))+'\n'
	else:
		lines += '\t'.join(map(str,[contig,'NA',num_copies]))+'\n'

if args['out'] != None:
	outfile_handle = args['out']
else:
	outfile_handle = assembly+".marker.tab"

with open(outfile_handle, 'w') as outfile:
	outfile.write(lines)

print("\nDone!")
