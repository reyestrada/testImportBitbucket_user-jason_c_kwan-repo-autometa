#!/usr/bin/env python2

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

import sys
import subprocess
import time
import logging
import os

import pandas as pd

from multiprocessing import cpu_count
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

from Bio import SeqIO

def init_logger(autom_path, db_path, out_path):
	#logger
	logger = logging.getLogger(out_path + '/run_autometa.py')
	hdlr = logging.FileHandler(out_path + '/run_autometa.log')
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	logger.addHandler(hdlr)
	logger.setLevel(logging.DEBUG)

	#Check user CPUs
	# NOTE: As far as I can tell this is not being used in this script
	#user_CPU_number = cpu_count()

	# Output current branch and commit
	branch_command = "git -C " + autom_path + " branch | grep \* | sed 's/^..//'"
	branch = subprocess.Popen(branch_command, shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()

	commit_command = 'git -C ' + autom_path + ' rev-parse --short HEAD'
	commit = subprocess.Popen(commit_command, shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
	logger.info('Autometa branch: {}'.format(branch))
	logger.info('Autometa commit: {}'.format(commit))
	#Check programs
	progs = ['md5sum','gzip','tar','gunzip','wget']
	for prog in progs:
		prog_version = subprocess.check_output([prog,'--version']).split('\n')[0]
		logger.info('{}'.format(prog_version))
	#Check 3rd party dependencies
	dmnd_v = subprocess.check_output(['diamond','version']).strip()
	logger.info('{}'.format(dmnd_v))
	# hmmscan_v = subprocess.check_output(['hmmpress','-h']).split('\n')[1].replace('# ','')
	# logger.info('{}'.format(hmmscan_v))
	prodigal_v = subprocess.Popen("prodigal -v", stderr=subprocess.PIPE, shell=True).communicate()[1].replace('\n','')
	logger.info('{}'.format(prodigal_v))
	logger.info('DB Dir: {}'.format(db_path))
	db_fpaths = os.listdir(db_path)
	for fpath in db_fpaths:
		logger.info('DB (fname, size): {} {}'.format(fpath, os.stat(db_path+'/'+fpath).st_size))
	return logger

def run_command(command_string, stdout_path = None):
	# Function that checks if a command ran properly. If it didn't, then print an error message then quit
	logger.info('run_autometa.py, run_command: ' + command_string)
	if stdout_path:
		f = open(stdout_path, 'w')
		exit_code = subprocess.call(command_string, stdout=f, shell=True)
		f.close()
	else:
		exit_code = subprocess.call(command_string, shell=True)

	if exit_code != 0:
		print('run_autometa.py: Error, the command:')
		print(command_string)
		print('failed, with exit code ' + str(exit_code))
		exit(1)

def run_command_quiet(command_string):
	exit_code = subprocess.call(command_string, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

	if exit_code !=0:
		print('run_autometa.py: Error, the command:')
		print(command_string)
		print('failed, with exit code ' + str(exit_code))
		exit(1)

# def cythonize_lca_functions():
# 	logger.info("{}/lca_functions.so not found, cythonizing lca_function.pyx for make_taxonomy_table.py".format(pipeline_path))
# 	current_dir = os.getcwd()
# 	os.chdir(pipeline_path)
# 	run_command("python setup_lca_functions.py build_ext --inplace")
# 	os.chdir(current_dir)

def run_make_taxonomy_tab(fasta, length_cutoff):
	"""Runs make_taxonomy_table.py and directs output to taxonomy.tab for run_autometa.py"""
	# Note we don't have to supply the cov_table here because earlier in this script we already run make_contig_table.py
	outfpath = os.path.join(output_dir,'taxonomy.tab')
	cmd = ' '.join(map(str,[
					os.path.join(pipeline_path,'make_taxonomy_table.py'),
					'-a',fasta,
					'-db',db_dir_path,
					'-p',processors,
					'-l',length_cutoff,
					'-o',output_dir
				]))
	if cov_table:
		cmd += ' '.join(['', '-v', cov_table])
	run_command(cmd)
	return outfpath

def length_trim(fasta, length_cutoff, outfpath=None):
	assert type(length_cutoff) == int, 'length_cutoff must be an integer value'

	if fasta.endswith('.gz'):
		fasta_stripped = fasta.rstrip('.gz')
		outfname, ext = os.path.splitext(os.path.basename(fasta_stripped))
	#will need to update path of this perl script
	else:
		outfname, ext = os.path.splitext(os.path.basename(fasta))

	outfname += ".filtered" + ext
	output_path = os.path.join(output_dir, outfname)
	cmd = ' '.join(map(str,[
					os.path.join(pipeline_path,'fasta_length_trim.py'),
					fasta,
					length_cutoff,
					output_path]))
	run_command(cmd)
	return output_path

def make_cov_table(asm_fpath, reads_fpath, dirpath, proc=1):
	asm_base, _ = os.path.splitext(os.path.basename(asm_fpath))
	cov_tab_fpath = os.path.join(dirpath, asm_base+'.coverage.tab')
	if os.path.isfile(cov_tab_fpath):
		return cov_tab_fpath
	cmd = ' '.join(map(str,[
		os.path.join(pipeline_path,'calculate_read_coverage.py'),
		'-a',asm_fpath,
		'-i',reads_fpath,
		'-p',proc,
		'-o',dirpath,
	]))
	run_command(cmd)
	return cov_tab_fpath

def make_contig_table(fasta, cov_tab_fpath=None):
	# Fasta is an absolute path
	outfname, _ = os.path.splitext(os.path.basename(fasta))
	outfname += '.tab'
	outfpath = os.path.join(output_dir, outfname)
	cmd = [
		os.path.join(pipeline_path,'make_contig_table.py'),
		'-a',fasta,
		'-o',outfpath,
	]
	if cov_tab_fpath:
		cmd.extend(['-c', cov_tab_fpath])
	cmd = ' '.join(cmd)
	run_command(cmd)
	return outfpath

def make_marker_table(fasta, all_orfs=None, domain='bacteria'):
	if domain == 'bacteria':
		hmms = 'Bacteria_single_copy.hmm'
		cutoffs = 'Bacteria_single_copy_cutoffs.txt'
	elif domain == 'archaea':
		hmms = 'Archaea_single_copy.hmm'
		cutoffs = 'Archaea_single_copy_cutoffs.txt'

	marker_dir = os.path.join(autometa_path, 'single-copy_markers')
	hmm_marker_path = os.path.join(marker_dir, hmms)
	hmm_cutoffs_path = os.path.join(marker_dir, cutoffs)
	#need to add processors to this script
	outfname, _ = os.path.splitext(os.path.basename(fasta))
	outfname += '.marker.tab'
	outfpath = os.path.join(output_dir, outfname)
	if os.path.isfile(outfpath):
		print "{} file already exists!".format(outfpath)
		print "Continuing to next step..."
		logger.info('{} file already exists!'.format(outfpath))
		logger.info('Continuing to next step...')
		return outfpath
	if all_orfs and os.path.isfile(all_orfs):
		# Retrieve orfs corresponding to fasta to skip prodigal step
		# 1. Retrieve contigs that need to be in file
		fasta_seqs = {seq.id for seq in SeqIO.parse(fasta, 'fasta')}
		# 2. Retrieve all orfs and subset
		orfs = (orf for orf in SeqIO.parse(all_orfs, 'fasta')
			if orf.id.rsplit('_',1)[0] in fasta_seqs)
		orfs_fname = os.path.splitext(os.path.basename(fasta))[0]+'.orfs.faa'
		orfs_fpath = os.path.join(output_dir, orfs_fname)
		SeqIO.write(orfs, orfs_fpath, 'fasta')

	print("Making marker tab w/prodigal & hmmscan")
	logger.info('Making {}: Running prodigal and hmmscan'.format(outfname))
	cmd = ' '.join(map(str,[
		os.path.join(pipeline_path,'make_marker_table.py'),
		'-a',fasta,
		'-m',hmm_marker_path,
		'-c',hmm_cutoffs_path,
		'-o',outfpath,
		'-p',processors,
	]))
	run_command_quiet(cmd)
	return outfpath

def recursive_dbscan(input_table, fasta_fp, domain):
	fname = '{}_recursive_dbscan_output.tab'.format(domain)
	dbscan_outfpath = os.path.join(output_dir,fname)
	kmer_fpath = os.path.join(output_dir,'k-mer_matrix')
	cmd = ' '.join([
		os.path.join(pipeline_path,'recursive_dbscan.py'),
		'-t',input_table,
		'-a',fasta_fp,
		'-d',output_dir,
		'-k',domain,
		'-o',dbscan_outfpath
	])
	run_command(cmd)
	return dbscan_outfpath, kmer_fpath

def combine_tables(table1_path, table2_path, outfname):
	comb_table_path = os.path.join(output_dir,outfname)
	# Note: in this sub we assume that the tables have column 1 in common
	# Store lines of table 2, keyed by the value of the first column

	# First make data structure from table 2
	table2_lines = dict()
	with open(table2_path) as table2:
		for i, line in enumerate(table2):
			line_list = line.rstrip().split('\t')
			contig = line_list.pop(0)
			if i == 0:
				table_2_header = '\t'.join(line_list)
			else:
				table2_lines[contig] = '\t'.join(line_list)

	comb_table = open(comb_table_path, 'w')
	with open(table1_path) as table1:
		for i, line in enumerate(table1):
			line_list = line.rstrip().split('\t')
			contig = line_list.pop(0)
			if i == 0:
				new_header = line.rstrip()+'\t'+table_2_header+'\n'
				comb_table.write(new_header)
			else:
				# We have to check whether the line exists in table 2
				if contig in table2_lines:
					new_line = line.rstrip()+'\t'+table2_lines[contig]+'\n'
					comb_table.write(new_line)

	comb_table.close()
	return comb_table_path

def ML_recruitment(input_table, matrix):
	outfpath = os.path.join(output_dir,'ML_recruitment_output.tab')
	cmd = ' '.join(map(str,[
					os.path.join(pipeline_path,'ML_recruitment.py'),
					'-t', input_table,
					'-p', processors,
					'-m', matrix,
					'-o', outfpath,
					'-r'
				]))
	run_command(cmd)
	return outfpath

def cami_format(infpath, out_dpath):
	fname = os.path.splitext(os.path.basename(infpath))[0]
	master_output = os.path.join(out_dpath, fname+'.binning')
	if not os.path.exists(master_output):
		version = '@Version:0.9.0'
		#@Version:CAMI2BinningFileFormat
		sample_id = '@SampleID:{}.autometaRun'.format(os.path.basename(infpath))
		#@SampleID:SAMPLEID
		#@@SEQUENCEID\tBINID
		with open(master_output, 'w') as outfile:
			lines = '\n'.join([version, sample_id, '@@SEQUENCEID\tBINID\n'])
			outfile.write(lines)
	cols = ['contig', 'cluster']
	df = pd.read_csv(infpath, sep='\t', usecols=cols, index_col='contig')
	df.to_csv(master_output, mode='a', sep='\t', header=False)
	return master_output

def rename_bins(results, outfpath, assembly_fp):
	if not os.path.exists(outfpath):
		bname = os.path.basename(assembly_fp)
		version = '@Version:0.9.0'
		sample_id = '@SampleID:AutometaRun-{}'.format(bname)
		cols = '@@SEQUENCEID\tBINID'
		header = '\n'.join([version,sample_id,cols])+'\n'
		with open(outfpath, 'w') as outfile:
			outfile.write(header)
		bin_num = 0
	else:
		bin_num = float('-inf')
		with open(outfpath) as fh:
			# Skip header lines
			for _ in range(3):
				fh.readline()
			# Read through to get current bin number
			for line in fh:
				ctg, bin = line.strip().split('\t')
				current_bin_num = int(bin.split('_')[1])
				if current_bin_num > bin_num:
					bin_num = current_bin_num
	# Construct bins dict
	for result_fp in results:
		df = pd.read_csv(
			result_fp,
			sep='\t',
			names=['contig','binid'],
			skiprows=3,
			header=None,
			index_col='contig',
		)
		bins = dict(list(df.groupby('binid')))
		for binid,dff in bins.items():
			if binid == 'unclustered':
				dff['binid'] = 'bin_0'
			else:
				bin_num += 1
				dff['binid'] = 'bin_{}'.format(bin_num)
			dff.to_csv(outfpath, mode='a', sep='\t', header=False)
	return outfpath

pipeline_path = sys.path[0]
autometa_path = os.path.dirname(pipeline_path)

#argument parser
parser = ArgumentParser(description="Script to run the Autometa pipeline.",\
 epilog="Please do not forget to cite us. Thank you for using Autometa!",\
  formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-a', '--assembly', metavar='<assembly.fasta>', help='Path to metagenomic assembly fasta', required=True)
parser.add_argument('-i', '--i_reads', metavar='<reads.fastq.gz>', help='</path/to/interleaved/reads.fq>', required=False)
parser.add_argument('-p', '--processors', metavar='<int>', help='Number of processors to use', type=int, default=1)
parser.add_argument('-l', '--length_cutoff', metavar='<int>', help='Contig length cutoff to consider for binning in bp', default=10000, type=int)
parser.add_argument('-c', '--completeness_cutoff', metavar='<float>', help='Completeness cutoff (in percent) to use for accepting clusters', type=float, default=20.0)
parser.add_argument('-k', '--kingdom', metavar='<archaea|bacteria>', help='Kingdom to consider',\
choices=['bacteria','archaea'], default = 'bacteria')
parser.add_argument('-t', '--taxonomy_table', metavar='<taxonomy.tab>', help='Path to output of make_taxonomy_table.py')
parser.add_argument('-o', '--output_dir', metavar='<dir>', help='Path to directory to store all output files', default = '.')
parser.add_argument('-r', '--ML_recruitment', help='Use ML to further recruit unclassified contigs', action='store_true')
parser.add_argument('-m', '--maketaxtable', action='store_true',\
help='runs make_taxonomy_table.py before performing autometa binning. Must specify databases directory (-db)')
parser.add_argument('-db', '--db_dir', metavar='<dir>', help="Path to directory with taxdump files. If this doesn't exist, the files will be automatically downloaded", required=False, default=os.path.join(autometa_path,'databases'))
parser.add_argument('-v', '--cov_table', metavar='<coverage.tab>', help="Path to coverage table made by calculate_read_coverage.py. If this is not specified then coverage information will be extracted from contig names (SPAdes format)", required=False)

args = vars(parser.parse_args())

length_cutoff = args['length_cutoff']
fasta_assembly = os.path.abspath(args['assembly'])
i_reads = args['i_reads']
processors = args['processors']
cluster_completeness = args['completeness_cutoff']
kingdom = args['kingdom'].lower()
taxonomy_table_path = args['taxonomy_table']
output_dir = os.path.abspath(args['output_dir'])
do_ML_recruitment = args['ML_recruitment']
make_tax_table = args['maketaxtable']
db_dir_path = os.path.abspath(args['db_dir'])
cov_table = args['cov_table']

# Make output directory if it doesn't exist
if not os.path.isdir(output_dir):
	os.mkdir(output_dir)

logger = init_logger(autometa_path, db_dir_path, output_dir)

#check if fasta in path
if not os.path.isfile(fasta_assembly):
	print "Could not find {}...".format(fasta_assembly)
	logger.debug('Could not find {}...'.format(fasta_assembly))
	exit(1)

#if make_tax_table specified but taxonomy_table_path not defined
if make_tax_table and not taxonomy_table_path:
	taxonomy_table_path = os.path.join(output_dir,'taxonomy.tab')

#If coverage table is given, it must exist
if cov_table and not os.path.isfile(cov_table):
	print('File path not found {}'.format(cov_table))
	exit(1)

#what input variables were and when you ran it (report fill path based on argparse)
logger.info('Input command: {}'.format(' '.join(sys.argv)))

start_time = time.time()
FNULL = open(os.devnull, 'w')

#run length trim and store output name
if fasta_assembly.endswith('.gz'):
	fstripped = fasta_assembly.rstrip('.gz')
	fname, ext = os.path.splitext(os.path.basename(fstripped))
else:
	fname, ext = os.path.splitext(os.path.basename(fasta_assembly))

filtered_asm_fpath = os.path.join(output_dir,fname)
filtered_asm_fpath += ".filtered"+ext
if not os.path.isfile(filtered_asm_fpath):
	filtered_assembly = length_trim(fasta_assembly, length_cutoff)
else:
	filtered_assembly = filtered_asm_fpath

# make_cov_table(asm_fpath, reads_fpath, proc=1, dirpath=output_dir)
if i_reads and not cov_table:
	i_reads = os.path.abspath(i_reads)
	cov_table = make_cov_table(filtered_assembly, i_reads, output_dir, processors)

if cov_table:
	contig_table = make_contig_table(filtered_assembly, cov_table)
else:
	contig_table = make_contig_table(filtered_assembly)

# Combine table if taxonomy should not be taken into account
if taxonomy_table_path and not make_tax_table:
	marker_tab_path = make_marker_table(filtered_assembly, domain=kingdom)
	combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path, 'combined_contig_info.tab')
elif taxonomy_table_path and make_tax_table:
	if not os.path.isfile(taxonomy_table_path):
		print "Could not find {}, running make_taxonomy_table.py".format(taxonomy_table_path)
		logger.debug('Could not find {}, running make_taxonomy_table.py'.format(taxonomy_table_path))
		# if not os.path.isfile(pipeline_path+"/lca_functions.so"):
		# 	cythonize_lca_functions()
		taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
		# combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
	elif os.path.isfile(taxonomy_table_path) and os.stat(taxonomy_table_path).st_size == 0:
		print "{} file is empty, running make_taxonomy_table.py".format(taxonomy_table_path)
		logger.debug('{} file is empty, running make_taxonomy_table.py'.format(taxonomy_table_path))
		# if not os.path.isfile(pipeline_path+"/lca_functions.so"):
		# 	cythonize_lca_functions()
		taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
		# combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
	else:
		print "{} already exists, not performing make_taxonomy_table.py".format(taxonomy_table_path)
		# combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
elif not taxonomy_table_path and make_tax_table:
	# if not os.path.isfile(pipeline_path+"/lca_functions.so"):
	# 	cythonize_lca_functions()
	taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
	# combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
else:
	marker_tab_path = make_marker_table(filtered_assembly, domain=kingdom)
	combined_table_path = combine_tables(contig_table, marker_tab_path, 'combined_contig_info.tab')

# If we carried out make_taxonomy_table.py, we now should use the appropriate kingdom bin instead of the raw
# input fasta
if make_tax_table:
	all_results = []
	# First, check that the expected kingdom bin is there
	all_orfs_fname = os.path.splitext(os.path.basename(filtered_assembly))[0]
	all_orfs_fname += '.orfs.faa'
	all_orfs_fpath = os.path.join(output_dir, all_orfs_fname)
	for kingdom in ['archaea','bacteria']:
		# taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
		kingdom_fpath = os.path.join(output_dir, kingdom.title()+'.fasta')
		if (not os.path.isfile(kingdom_fpath)) or os.stat(kingdom_fpath).st_size == 0:
			print('Kingdom {} is either not there or empty'.format(kingdom.title()))
			continue
		print('Binning {}'.format(kingdom.title()))
		marker_tab_path = os.path.join(output_dir, kingdom+'.marker.tab')
		if not os.path.exists(marker_tab_path):
			marker_tab_path = make_marker_table(
				fasta=kingdom_fpath,
				all_orfs=all_orfs_fpath,
				domain=kingdom,
			)
		combined_table_fname = kingdom+'_contig_info.tab'
		combined_table_path = os.path.join(output_dir, combined_table_fname)
		if not os.path.exists(combined_table_path):
			combined_table_path = combine_tables(
				table1_path=taxonomy_table_path,
				table2_path=marker_tab_path,
				outfname=combined_table_fname,
			)
		binning_fname = '{}_recursive_dbscan_output.tab'.format(kingdom)
		binning_outfpath = os.path.join(output_dir, binning_fname)
		if not os.path.exists(binning_outfpath):
			binning_outfpath, matrix_file = recursive_dbscan(
				input_table=combined_table_path,
				fasta_fp=kingdom_fpath,
				domain=kingdom,
			)

		if do_ML_recruitment:
			binning_outfpath = ML_recruitment(binning_outfpath, matrix_file)

		cami_formatted_results = cami_format(binning_outfpath, output_dir)
		elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))
		all_results.append(cami_formatted_results)
		print('Wrote Binning results to {}'.format(cami_formatted_results))
	master_outfpath = rename_bins(
		results=all_results,
		outfpath=os.path.join(output_dir,'autometa_cami2.binning'),
		assembly_fp=fasta_assembly,
	)

print "Done!"
print "Elapsed time is {} (HH:MM:SS)".format(elapsed_time)
logger.info('Done!')
logger.info('Elapsed time is {} (HH:MM:SS)'.format(elapsed_time))
logger.info('Binning: {}'.format((os.path.abspath(master_outfpath))))
FNULL.close()
