#!/usr/bin/env python

import editdistance
import pysam
import numpy as np
from Bio import SeqIO
from typing import Tuple, Dict
from collections import Counter, defaultdict
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.csgraph import connected_components


def bc_and_umi_given_read_name(read_name: str) -> Tuple[str, str]:
    """Split a read name into its cell barcode and UMI.

    Read names must have the form `{barcode}_{umi}#{original_read_name}`.  The original read name
    may itself contain `_` and `#`. A barcode or UMI may not contain `_`.

    Args:
        read_name: The full read name from the FASTQ or BAM.

    Returns:
        `(barcode, umi)`.
    """
    return read_name.split('#')[0].split('_')

def umi_given_read_name(read_name: str) -> str:
    return bc_and_umi_given_read_name(read_name)[1]

def get_connected_components(umis: list, max_dist: int = 1) -> Tuple[int, np.array]:
    """
    Builds adjacency matrix of umis given max_dist and calls connected_componenets

    returns
        :int:       n_vals
        :int_array: component_array
    """
    adj_mat = lil_matrix((len(umis), len(umis)), dtype=np.uint8)
    for i, umi_i in enumerate(umis):
        for j in range(i+1, len(umis)):
            umi_j = umis[j]
            if editdistance.distance(umi_i, umi_j) <= max_dist:
                adj_mat[i, j] = 1
                adj_mat[j, i] = 1
    return connected_components(adj_mat)

def get_umi_map_from_cntr(umi_cntr: Counter, max_dist: int = 1) -> dict:
    """
    Builds a dict from observed umis to connected component umi with max count.
    """
    umi_list = list(umi_cntr.keys())
    n_vals, component_array = get_connected_components(umi_list, max_dist=max_dist)
    component_max_umi = [None for _ in range(n_vals)]
    for umi, component in zip(umi_list, component_array):
        if component_max_umi[component] is None or umi_cntr[umi] > umi_cntr[component_max_umi[component]]:
            component_max_umi[component] = umi
    umi_map = {umi: component_max_umi[component] for umi, component in zip(umi_list, component_array)}
    return umi_map

def get_umi_maps_from_fastq_or_bam_file(
        fastq_fpath: str = None,
        bam_fpath: str = None,
        chrm: str = None,
        start: int = None,
        end: int = None,
        max_dist: int = 1
        ) -> Dict[str, Counter]:
    """
    Builds umi_map_given_bc for reads from (specified region of) a given file.

    returns
        :dict: umi_map_given_bc
    """

    # Select correct iterator and read name function for fastq/bam file
    if fastq_fpath is not None:
        assert bam_fpath is None, "Only one file allowed."
        read_iter = SeqIO.parse(open(fastq_fpath), 'fastq')
        get_read_name = lambda read: read.name
    else:
        assert bam_fpath is not None, "File must be specified."
        read_iter = pysam.AlignmentFile(bam_fpath).fetch(chrm, start, end)
        get_read_name = lambda read: read.qname

    # build umi_maps
    umi_cntr_given_bc = defaultdict(Counter)
    for read in read_iter:
        bc, umi = bc_and_umi_given_read_name(get_read_name(read))
        umi_cntr_given_bc[bc][umi] += 1
    umi_cntr_given_bc = dict(umi_cntr_given_bc)

    umi_map_given_bc = {}
    for bc, umi_cntr in umi_cntr_given_bc.items():
        umi_map_given_bc[bc] = get_umi_map_from_cntr(umi_cntr, max_dist=max_dist)

    return umi_map_given_bc

def make_fastq_with_corrected_umis(
        input_fq_fname: str, 
        output_fq_fname: str, 
        max_dist: int = 1
        ) -> None:
    """
    Outputs a fastq with corrected UMIs
    """
    umi_map_given_bc = get_umi_maps_from_fastq_or_bam_file(input_fq_fname, max_dist=max_dist)
    with open(output_fq_fname, 'w') as out:
        for rec in SeqIO.parse(open(input_fq_fname), 'fastq'):
            bc, umi = bc_and_umi_given_read_name(rec.id)
            orig_readname = str(rec.id).split('#')[1]
            new_umi = umi_map_given_bc[bc][umi]
            new_readname = f'{bc}_{new_umi}#{orig_readname}'
            rec.id = new_readname
            rec.name = new_readname
            SeqIO.write(rec, out, 'fastq')

