FROM continuumio/anaconda
MAINTAINER Jason C. Kwan "jason.kwan@wisc.edu"

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

# Bioinformatics progs. and system utils.
RUN apt-get update
RUN apt-get install -y \
  apt-utils \
  bedtools \
  bowtie2 \
  build-essential \
  hmmer \
  libatlas-base-dev \
  libncurses5-dev \
  libncursesw5-dev \
  libbz2-dev \
  liblzma-dev \
  prodigal \
  zlib1g-dev

# Compile diamond
RUN mkdir diamond \
  && cd diamond \
  && wget http://github.com/bbuchfink/diamond/releases/download/v0.9.14/diamond-linux64.tar.gz \
  && tar xvf diamond-linux64.tar.gz \
  && mv diamond /usr/bin \
  && cd ../ \
  && rm -rf /diamond

# Compile samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 \
  && tar -vxjf samtools-1.6.tar.bz2 \
  && rm samtools-1.6.tar.bz2 \
  && cd samtools-1.6 && ./configure --prefix=/usr/bin && make && make install

# Python Dependencies
RUN pip install tqdm joblib biopython cython

# Compile Barnes-hut Stochastic Neighbor Embedding
RUN git clone https://github.com/danielfrg/tsne.git \
  && cd tsne \
  && python setup.py install

# Install Autometa
RUN git clone https://bitbucket.org/jason_c_kwan/autometa \
  && cd autometa/pipeline \
  && python setup_lca_functions.py build_ext --inplace

# Prepare Databases.
RUN hmmpress -f autometa/single-copy_markers/Archaea_single_copy.hmm \
  && hmmpress -f autometa/single-copy_markers/Bacteria_single_copy.hmm

ENV PATH="/autometa/pipeline:${PATH}"
