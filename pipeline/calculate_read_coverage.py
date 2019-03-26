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

#import multiprocessing

#Input: fastq (or fastq.gz?) reads, metagenomic assembly contigs as fasta file
#Output: "contig\tread_coverage"

#See: /home/ijmiller/slurm_scripts/2016/October/7-Oct-16_Bpac5Tips_coverage_calculation.sh
#Dependencies: bowtie2 (2.2.5), samtools (0.1.18), \
    #Bedtools (genomeCoverageBed: v2.17.0), contig_coverage_from_bedtools.pl,\
    #fasta_length_table.pl, python modules: suprocess, argparse

#1. Align the reads to your assembly with bowtie2.
#2. Convert the SAM file to a sorted BAM file, and create a BAM index.
#3. Tabulate the average coverage of each contig.

#num_cpus = multiprocessing.cpu_count()
#use_cpus = int(round((num_cpus / 2)))

import argparse
import subprocess
import os

def run_command(command_string, stdout_path = None):
    # Function that checks if a command ran properly. If it didn't, then print an error message then quit
    print('calculate_read_coverage.py, run_command: ' + command_string)
    if stdout_path:
        f = open(stdout_path, 'w')
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        print('calculate_read_coverage.py: Error, the command:')
        print(command_string)
        print('failed, with exit code ' + str(exit_code))
        exit(1)

def run_bowtie2(asm_fpath, reads_fpath, output_dir, num_processors=1):
	#When "shell = True", need to give one string, not a list
    #Build bowtie2 database
    asm_bname, ext = os.path.splitext(os.path.basename(asm_fpath))
    bowtie2_db = os.path.join(output_dir, asm_bname)
    cmd = ' '.join(['bowtie2-build', asm_fpath, bowtie2_db])
    run_command(cmd)
    #Run bowtie2 alignment
    sam_outfile = bowtie2_db + '.sam'
    bt2_cmd = ' '.join(['bowtie2 -x',bowtie2_db,
                        '--interleaved',reads_fpath,
                        '-q',
                        '--phred33',
                        '--very-sensitive',
                        '--no-unal',
                        '-p',str(num_processors),
                        '-S',sam_outfile])

    run_command(bt2_cmd)
    return sam_outfile

def main():
    #argument parser
    parser = argparse.ArgumentParser(description='Script to tabulate paired-end or single read \
    coverage using Samtools and Bedtools. Note: it is recommended that you quality \
    filter your *.fastq or *fastq.gz files prior to calculating read coverage.')
    parser.add_argument('-a','--assembly', metavar='<assembly.fasta>', help='Path to input assembly file', required=True)
    parser.add_argument('-i','--ireads', metavar='<ireads.fastq.gz>', help='/path/to/interleaved/reads.fastq.gz', required=True)
    parser.add_argument('-p','--processors', metavar='<int>', help='Number of processors', type=int, default=1)
    parser.add_argument('-o','--outdir', metavar='<dir>', help='Path to output directory', default='.')
    args = vars(parser.parse_args())

    #Check for dependencies in $PATH
    assembly_fpath = os.path.abspath(args['assembly'])
    assembly_bname = os.path.basename(assembly_fpath)
    i_reads = os.path.abspath(args['i_reads'])
    outdir = os.path.abspath(args['outdir'])
    num_proc = str(args['processors'])
    # forward_read_path_list = args['forward_reads']
    # reverse_read_path_list = args['reverse_reads']
    # Check that read fpath exists
    for fpath in [i_reads, assembly_fpath]:
        if not os.path.isfile(fpath):
            print('File path not found {}'.format(fpath))
            exit(1)

    # Check that the output dir exists
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    #1. Align the comparison dataset reads to your assembly with bowtie2.
    samfpath = run_bowtie2(assembly_fpath, i_reads, outdir, num_proc)

    #2. Convert the SAM file to a sorted BAM file, and create a BAM index.
    sorted_bam_fpath = '{0}/{1}.sort.bam'.format(outdir, assembly_bname)
    cmd = ' '.join(['samtools','view','-bS', samfpath,'|',
                    'samtools','sort','-o',sorted_bam_fpath])
    run_command(cmd)

    #Clean up the SAM file, which will be a lot larger than the sorted BAM file
    os.remove(samfpath)

    #3. Tabulate the average coverage of each contig.

    #Run fasta_length_table.py
    ctg_len_tab_fpath = '{0}/{1}.tab'.format(outdir, assembly_bname)
    cmd = ' '.join(['fasta_length_table.py', assembly_fpath, ctg_len_tab_fpath])
    run_command(cmd)

    #Run genomeCoverageBed
    genome_bed_fpath = '{0}/{1}.txt'.format(outdir,assembly_bname)
    cmd = ' '.join(['genomeCoverageBed',
                    '-ibam',sorted_bam_fpath,
                    '-g',ctg_len_tab_fpath])
    run_command(cmd, genome_bed_fpath)

    #Build final table
    asm_base, ext = os.path.splitext(assembly_bname)
    outfile = '{0}/{1}.coverage.tab'.format(outdir, asm_base)
    cmd = ' '.join(['contig_coverage_from_bedtools.pl',genome_bed_fpath])
    run_command(cmd, outfile)

if __name__ == '__main__':
    main()
