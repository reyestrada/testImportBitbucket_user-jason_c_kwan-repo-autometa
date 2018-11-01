#!/usr/bin/env python

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

import argparse
import subprocess
import os
#import multiprocessing

#Input: fastq (or fastq.gz?) reads, metagenomic assembly contigs as fasta file
#Output: "contig\tread_coverage"

#: Dependencies: bowtie2 (2.2.5), samtools (0.1.18),
#: Bedtools (genomeCoverageBed: v2.17.0)
#: contig_coverage_from_bedtools.pl
#: fasta_length_table.pl

#1. Align the reads to your assembly with bowtie2.
#2. Convert the SAM file to a sorted BAM file, and create a BAM index.
#3. Tabulate the average coverage of each contig.

#num_cpus = multiprocessing.cpu_count()
#use_cpus = int(round((num_cpus / 2)))

def run_cmd(command_string, stdout_path = None):
    # Function checks if command ran properly. If not print error message & quit
    print('calculate_read_coverage.py\ncommand: ' + command_string)
    if stdout_path:
        f = open(stdout_path, 'w')
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        print('calculate_read_coverage.py\nError!\ncommand:')
        print(command_string)
        print('Failed!\nExit code: ' + str(exit_code))
        exit(1)

def run_bowtie2(asm_fpath, f_reads, r_reads, s_reads, num_cpus=1):
    # Check the number/type of reads we have as input
    if not (len(f_reads) or len(r_reads) or len(s_reads)):
        print ('Error! Specify at least one read file')
        quit()

    if not (len(f_reads) == len(r_reads)):
        print('Error! Specify the same number of forward and reverse read files')
        quit()

    if f_reads:
        f_read_path_str = ','.join(f_reads)
        r_read_path_str = ','.join(r_reads)

    if s_reads:
        s_read_path_str = ','.join(s_reads)

	#When "shell = True", need to give one string, not a list
    #Build bowtie2 database
    assembly_fname, _ = os.path.splitext(os.path.basename(asm_fpath))
    bowtie2_db = '{}/{}'.format(output_dir, assembly_fname)
    cmd = ' '.join(['bowtie2-build', asm_fpath, bowtie2_db])
    run_cmd(cmd)
    #Run bowtie2 alignment
    sam_fname = bowtie2_db + '.sam'
    bt2_cmd_list = [
        'bowtie2 -x', bowtie2_db,
        '-q --phred33 --very-sensitive --no-unal',
        '-p', str(num_cpus),
        '-S', sam_fname
        ]

    if f_reads:
        bt2_cmd_list.extend(['-1', f_read_path_str, '-2', r_read_path_str])

    if s_reads:
        bt2_cmd_list.extend(['-U', s_read_path_str])

    bt2_cmd = ' '.join(bt2_cmd_list)
    run_cmd(bt2_cmd)
    return sam_fname


#argument parser
parser = argparse.ArgumentParser(description='Script to tabulate paired-end or \
single read coverage using Samtools and Bedtools. Note: it is recommended that \
you quality filter your *.fastq or *fastq.gz files prior to calculating read \
coverage.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a','--assembly', metavar='<assembly.fasta>',
help='Path to input assembly file', required=True)
parser.add_argument('-p','--processors', metavar='<int>', type=int, default=1,
help='Number of processors')
parser.add_argument('-fwd','--fwd_reads', metavar='<PE_F_reads>',  nargs='*',
help='Paired (*.fastq|*.fastq.gz) forward reads (must be in same order as reverse list)')
parser.add_argument('-rev','--rev_reads', metavar='<PE_R_reads>', nargs='*',
help='Paired (*.fastq|*.fastq.gz) reverse reads (must be in same order as forward list)')
parser.add_argument('-S','--single_reads', metavar='<SE_reads>',
help='Single (*.fastq|*.fastq.gz) reads', nargs='*')
parser.add_argument('-o','--output_dir', metavar='<dir>', default=os.getcwd(),
help='Path to output directory')

args = vars(parser.parse_args())

#Check for dependencies in $PATH

#1. Align the comparison dataset reads to your assembly with bowtie2.
assembly_fpath = os.path.abspath(args['assembly'])
assembly_fname, _ = os.path.splitext(os.path.basename(assembly_fpath))
read_paths = {
    'fwd':args['fwd_reads'],
    'rev':args['rev_reads'],
    'single': args['single_reads']
}

output_dir = os.path.abspath(args['output_dir'])

if not os.path.isfile(assembly_fpath):
    print('Error! Assembly path not found.\nPath: {}'.format(assembly_fpath))
    exit(1)

# Check that at least one of the above lists is not empty
if not (read_paths['fwd'] or read_paths['rev'] or read_paths['single']):
    print('Error! Must provide -fwd/-rev (paired) and/or -S (single)')
    exit(1)

# Check that all read paths exist
for read_type in read_paths:
    for fpath in read_paths[read_type]:
        print('Read type: {} Path: {}'.format(read_type, fpath))
        if not os.path.isfile(fpath):
            print('Error! Cannot find the read file: ' + path)
            exit(1)

# Check that the output dir exists
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

sam_file = run_bowtie2(
    asm_fpath = assembly_fpath,
    f_reads = read_paths['fwd'],
    r_reads = read_paths['rev'],
    s_reads = read_paths['single'],
    num_cpus = args['processors']
)

#2. Convert the SAM file to a sorted BAM file, and create a BAM index.
sorted_bam_fh = '{}/{}.sort.bam'.format(output_dir, assembly_fname)
cmd = ' '.join(['samtools view -bS', sam_file, '| samtools sort -o', sorted_bam_fh])
run_cmd(cmd)

#Clean up the SAM file, which will be a lot larger than the sorted BAM file
os.remove(sam_file)

#3. Tabulate the average coverage of each contig.
#Run fasta_length_table.pl
contig_len_tab = '{}/{}.tab'.format(output_dir, assembly_fname)
cmd = ' '.join(['fasta_length_table.pl', assembly_fpath, contig_len_tab])
run_cmd(cmd)

#Run genomeCoverageBed
genomebed_fh = '{}/{}.txt'.format(output_dir, assembly_fname)
cmd = ' '.join(['genomeCoverageBed',
                '-ibam', sorted_bam_fh,
                '-g', contig_len_tab])
run_cmd(cmd, genomebed_fh)

#Build final table
outfile = '{}/coverage.tab'.format(output_dir)
cmd = ' '.join(['contig_coverage_from_bedtools.pl', genomebed_fh])
run_cmd(cmd, outfile)

#Future features: check for .fastq.gz, dependencies in path, default processors
