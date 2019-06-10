#!/usr/bin/env python

import sys
import os

from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from math import log


def sizeof_fmt(num):
    # https://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
    unit_list = zip(['bytes', 'kB', 'MB', 'GB', 'TB', 'PB'], [0, 0, 1, 2, 2, 2])
    """Human friendly file size"""
    if num > 1:
        exponent = min(int(log(num, 1024)), len(unit_list) - 1)
        quotient = float(num) / 1024**exponent
        unit, num_decimals = unit_list[exponent]
        format_string = '{:.%sf} {}' % (num_decimals)
        return format_string.format(quotient, unit)
    if num == 0:
        return '0 bytes'
    if num == 1:
        return '1 byte'

parser = ArgumentParser(description='Splits fasta file into equal sized number of files',
formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-f', '--fasta', metavar='<file>', required=True,
help='Path to Fasta')
parser.add_argument('-c', '--num_files', metavar='<int>', default=4,
required=False, help='Number of approximately equal sized files to generate')
parser.add_argument('-fpd', '--files_per_dir', metavar='<int>', default=200,
required=False, help='Overflow for multiple directories.')
parser.add_argument('-o', '--output_dir', metavar='<dir>', default = os.getcwd(),
required=False, help='Path to store output files')
parser.add_argument('-r', '--report', action='store_true', default=False,
help='report file size if fasta split into specified number of files and exit')

args = vars(parser.parse_args())
fasta = args['fasta']
outdir = os.path.abspath(args['output_dir'])
outdirname = '{0}/splits_1'.format(outdir)
num_files = int(args['num_files'])
files_per_dir = int(args['files_per_dir']) - 1 #: python starts count at 0.
max_file_size = os.stat(fasta).st_size/num_files

split_info = ' '.join(map(str,['Fasta:',os.path.basename(fasta),
    '\nsplit into', num_files,'files',
    '\nMax files per directory:',files_per_dir+1,
    '\nMax file size:',sizeof_fmt(max_file_size)]))

if args['report']:
    print(split_info)
    exit()

print(split_info)

if not os.path.isdir(outdir):
    os.mkdir(outdir)
if not os.path.isdir(outdirname):
    os.mkdir(outdirname)

fname, ext = os.path.splitext(os.path.basename(fasta)) #: ext == '.ext'
fnum = '0'
fname = fname.replace('.','_') if '.' in fname else fname
outfpath = '{dir}/{base}_split.{num}{ext}'.format(dir=outdirname, base=fname, num=fnum, ext=ext)
outfile = open(outfpath, 'w')
print('generating file {}'.format(os.path.basename(outfpath)))
with open(fasta) as fh:
    for line in fh:
        if len(line) == 0 or line == '\n':
            #: Removing blank and newline character lines
            continue
        if line.startswith('>'):
            outfile.close()
            if os.stat(outfile.name).st_size >= max_file_size:
                fpath, fnum, ext = outfile.name.split('.')
                outdirpath = os.path.dirname(os.path.abspath(outfile.name))
                fn = os.path.basename(outfile.name).replace('_','.')
                fname = os.path.join(outdirpath,fn)
                os.rename(os.path.abspath(outfile.name), fname)
                if int(fnum) == files_per_dir:
                    print('File limit reached. Moving to next directory')
                    #If at max files per directory add directory & restart count
                    dirname, dirnum = os.path.basename(outdirpath).split('_')
                    outdirpath = '_'.join([dirname, str(int(dirnum)+1)])
                    outdirpath = os.path.join(outdir, outdirpath)
                    if not os.path.isdir(outdirpath):
                        os.mkdir(outdirpath)
                    fnum = '0'
                    fpath = os.path.join(outdirpath, os.path.basename(fpath))
                else:
                    fnum = str(int(fnum) + 1)
                outfpath = '.'.join([fpath, fnum, ext])
                outfile = open(outfpath, 'w')
                print('generating file {}'.format(os.path.basename(outfpath)))
        if outfile.closed:
            outfile = open(outfile.name, 'a')
        outfile.write(line)
outfile.close()
outdirpath = os.path.dirname(os.path.abspath(outfile.name))
fname = outdirpath + '/' + os.path.basename(outfile.name).replace('_','.')
os.rename(os.path.abspath(outfile.name), fname)
