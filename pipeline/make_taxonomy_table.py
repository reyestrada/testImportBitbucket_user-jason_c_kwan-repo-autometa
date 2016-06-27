#!/usr/bin/env python

import subprocess
import argparse
import os

#argument parser
parser = argparse.ArgumentParser(description="Script to generate the contig taxonomy table.")
parser.add_argument('-a','--assembly', help='assembly.fasta', required=True)
parser.add_argument('-p','--processors', help='assembly.fasta', default=1)
parser.add_argument('-l','--length_cutoff', help='Contig length cutoff to consider for binning.\
 Default is 10,000 bp.', default=10000, type = int)
args = vars(parser.parse_args())

num_processors = args['processors']
length_cutoff = args['length_cutoff']
fasta_assembly = args['assembly']


def length_trim(fasta,length_cutoff):
	#Trim the length of fasta file
	outfile_name = str(args['assembly'].split("/")[-1].split(".")[0]) + "_filtered.fasta"
	subprocess.call("{}fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta, length_cutoff, outfile_name), shell = True)
	return outfile_name

def run_prodigal(path_to_assembly):
	#When "shell = True", need to give one string, not a list
	subprocess.call(" ".join(['prodigal ','-i ' + path_to_assembly, '-a ' + path_to_assembly.split(".")[0] +\
	 '.orfs.faa','-p meta' '-m', '-o ' + path_to_assembly.split(".")[0] + '.txt']), shell = True)

def run_diamond(diamond_path, prodigal_output, diamond_database_path, num_processors, prodigal_daa):
	view_output = prodigal_output + ".tab"
	subprocess.call("{} blastp --query {}.faa --db {} --evalue 1e-5 --max-target-seqs 200 -p {} --daa {}"\
		.format(diamond_path,prodigal_output, diamond_database_path, num_processors, prodigal_daa), shell = True)
	subprocess.call("{} view -a {} -f tab -o {}".format(diamond_path,prodigal_daa, view_output), shell = True)
	return view_output
	#return  view_output
	#might want to change name of outputfile

def run_blast2lca(input_file):
	output = input_file.rstrip(".tab") + ".lca"
	subprocess.call("/home/ijmiller/go/bin/blast2lca -savemem -dict /home/jkwan/blast2lca_taxdb/gi_taxid.bin -nodes /home/jkwan/blast2lca_taxdb/nodes.dmp -names /home/jkwan/blast2lca_taxdb/names.dmp {} > {}"\
		.format(input_file, output), shell = True) 
	return output

def run_taxonomy(add_contig_path, contig_table_path, tax_table_path, taxdump_dir_path): #Have to update this
	subprocess.call("{} {} {} {} taxonomy.tab".format(add_contig_path, contig_table_path,tax_table_path, taxdump_dir_path), shell = True)
	#contig_table_path, tax_table_path, taxdump_dir_path, output_file_path

diamond_path = subprocess.check_output('find ~ -name "diamond"', shell=True).rstrip("\n")
taxdump_dir_path = '/home/jkwan/blast2lca_taxdb'
prodigal_output = args['assembly'].rstrip(".fasta") + "_filtered.orfs"
prodigal_daa = prodigal_output + ".daa"
autometa_path = subprocess.check_output('find ~ -name "autometa"', shell=True).rstrip("\n")
pipeline_path = autometa_path + "/pipeline/"
diamond_database_path = subprocess.check_output('find /mnt/not_backed_up/nr_diamond/ -name "nr.dmnd"', shell=True).strip("\n")
add_contig_path = subprocess.check_output('find ~ -name "add_contig_taxonomy.py"', shell=True).rstrip("\n")

if not os.path.isfile(prodigal_output + ".faa"):
	print "Prodigal output not found. Running prodigal..."
	#Check for file and if it doesn't exist run make_marker_table
	filtered_assembly = length_trim(fasta_assembly, length_cutoff)
	run_prodigal(filtered_assembly)


if not os.path.isfile(prodigal_output + ".tab"):
	print "Running diamond blast... "
	#diamond_output = 
	diamond_output = run_diamond(diamond_path, prodigal_output, diamond_database_path, num_processors, prodigal_daa)
else:
	diamond_output = prodigal_output + ".tab"

blast2lca_output = run_blast2lca(diamond_output)
print "Running add_contig_taxonomy.py... "
run_taxonomy(add_contig_path, diamond_output, blast2lca_output, taxdump_dir_path)
print "Done!"