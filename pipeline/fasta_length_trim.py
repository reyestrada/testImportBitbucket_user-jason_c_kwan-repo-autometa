#!/usr/bin/env python
# Trims a fasta file based on length provided in args

import os
import sys
import gzip

from Bio import SeqIO


usage = 'Usage: fasta_length_trim.py </path/to/fasta[.gz]> <LengthCutoff> </path/to/outfile>'

if not len(sys.argv) >= 4:
    exit(usage)

infpath = os.path.abspath(sys.argv[1])
outfpath = os.path.abspath(sys.argv[3])

try:
    cutoff = int(sys.argv[2])
except ValueError as e:
    exit('{} must be an integer\n{}'.format(cutoff,usage))

if infpath.endswith('.gz'):
    fh = gzip.open(infpath)
else:
    fh = open(infpath)
# Parse seqs into dict
seqs = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
fh.close()
# Filter based on cutoff length
recs = [seqs[s] for s in seqs if len(seqs[s].seq) >= cutoff]
# Write out seqs
n_seqs = SeqIO.write(recs, outfpath, 'fasta')
# Print File information
print('wrote {} seqs to {}'.format(n_seqs, outfpath))
