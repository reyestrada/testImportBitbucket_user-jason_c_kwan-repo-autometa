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


# Given a fasta file construct a tab-delimited file of 'contig\tlength\n'


import sys
import os
import gzip
from Bio import SeqIO

def main():
    usage = 'fasta_length_table.py </path/to/fasta> </path/to/outfile>'
    if not len(sys.argv) >= 2:
        exit(usage)
    fasta = os.path.abspath(sys.argv[1])
    outfpath = os.path.abspath(sys.argv[2])
    outfile = open(outfpath, 'w')
    header = 'sequence\tlength\n'
    outfile.write(header)
    if fasta.endswith('.gz'):
        fh = gzip.open(fasta)
    else:
        fh = open(fasta)
    for s in SeqIO.parse(fh, 'fasta'):
        line = '{seq}\t{length}\n'.format(seq=s.id, length=len(s.seq))
        outfile.write(line)
    fh.close()
    outfile.close()
    print('written: {}'.format(outfpath))

if __name__ == '__main__':
    main()
