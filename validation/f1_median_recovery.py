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

import sys
import os
import pandas as pd
import numpy as np

def summarize_f1_stats(f1_fh, max_bins):
    df = pd.read_table(f1_fh)
    df.dropna(inplace=True)
    df = df[df['ref_genome'] != 'misassembled']
    recovery = round(sum(df.F1/100.0)/max_bins, 4)
    median = round(np.median(df.F1/100.0), 4)
    return(os.path.basename(f1_fh), recovery, median)

if not len(sys.argv) >=3:
    exit('Usage: f1_median_recovery.py <F1.tab> <theoretical max no. bins> ...')

f1_tab = sys.argv[1]
try:
    theo_max_bins = int(sys.argv[2])
except:
    exit("Must provide an int value for the max no. bins")

fname, recovery, median = summarize_f1_stats(f1_tab, theo_max_bins)
print('Table: {}\t F1 recovery: {}\tmedian F1: {}'.format(fname, recovery, median))
