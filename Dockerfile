FROM continuumio/miniconda
MAINTAINER Evan R. Rees "erees@wisc.edu"

# Copyright 2018 Ian J. Miller, Evan R. Rees, Izaak Miller, Jason C. Kwan
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

# Download and install dependencies
RUN conda install -c bioconda -y \
    prodigal \
    hmmer \
    bowtie2 \
    bedtools \
    diamond \
    tqdm \
    biopython \
    samtools \
    numpy \
    pandas \
    cython

RUN conda update -y conda
RUN apt-get update && apt-get install -y \
    zlib1g-dev \
    libc-dev \
    libatlas-base-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    build-essential

# Download and install Autometa
RUN git clone https://bitbucket.org/jason_c_kwan/autometa && \
    cd autometa && \
    git pull && \
    git checkout cami2 && \
    chmod 777 single-copy_markers/Bacteria_single_copy.hmm && \
    chmod 777 single-copy_markers/Archaea_single_copy.hmm && \
    hmmpress -f single-copy_markers/Bacteria_single_copy.hmm && \
    hmmpress -f single-copy_markers/Archaea_single_copy.hmm

RUN pip install tsne

ENV PATH="/autometa/pipeline:${PATH}"
