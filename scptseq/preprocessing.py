import logging
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import pysam
import yaml
import glob
import os
from collections import Counter, defaultdict
from statsmodels.distributions.empirical_distribution import ECDF
from . import umi_tools
from .constants import haplotypes
from .misc import parse_mutation, bc_from_fpath, TargetInfo

log = logging.getLogger(__name__)
plt.set_loglevel('critical')


def get_start_clip_and_last_rr(read):
    """Return number of bases to clip from beginning and last aligned read pos"""
    start_clip = next(i for i, (qq, rr) in enumerate(read.aligned_pairs) if rr)
    last_rr = next(rr for i, (qq, rr) in enumerate(read.aligned_pairs[::-1]) if rr) # last_rr is an integer position
    return start_clip, last_rr  


def get_splicing_junctions_from_ctrl(read, splicing_threshold=15):
    """Return candidate splice junctions observed in one control read.

    The catalogue of candidate splice junctions is built from a control sample, where cut-site
    deletions are not expected. The CIGAR string is not consulted. 

    Args:
        read: A pysam aligned read.
        splicing_threshold: A gap is recorded when its last reference position minus
            its first exceeds this value 

    Returns:
        Set of `(first_position, last_position)` tuples, one per gap recorded.
        Positions are 0-based, and soft-clipped ends plus the final aligned position
        are excluded.
    """
    splicing_junctions = set()

    chrm = read.reference_name
    query = read.query
    start_clip, last_rr = get_start_clip_and_last_rr(read)
    ins_seq = ''
    ins_after = None
    ins_active = False
    del_start = None
    del_active = False
    prev_rr = None
    for qq, rr in read.aligned_pairs:
        # skip clips at beginning and cend
        if rr is None and prev_rr is None:
            continue
        if rr is not None and rr >= last_rr:
            break

        # Process deletions and insertions after completed
        if del_active and qq is not None:
            if prev_rr - del_start > splicing_threshold:
                splicing_junctions.add((del_start, prev_rr))
            del_active = False
            del_start = None
        if ins_active and rr is not None:
            ins_seq = ''
            ins_after = None
            ins_active = False

        # Start/continue indels and find mismatches and seg sites
        if qq is None:
            # deletion or splicing
            if not del_active:
                del_active = True
                del_start = rr
        elif rr is None:
            # insertion
            qbase = query[qq-start_clip].upper()
            ins_seq += qbase
            if not ins_active:
                ins_active = True
                ins_after = prev_rr

        prev_rr = rr

    return splicing_junctions


def load_all_muts_with_info(fpath):
    with open(fpath) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
    all_muts = defaultdict(list)
    for read_name, alignment_ends, haplotype, seg_site_bases, muts in data:
        all_muts[haplotype].append((read_name, alignment_ends, muts))
    return dict(all_muts)


def get_haplotyped_mutations(read, genome, splicing_junction_str, seg_sites=[], skip_sites=[]):
    """Extract one read's mutations and the bases it shows at the segregating sites.

    Walks the aligned pairs, outputting mutation strings in the notation described by
    `misc.parse_mutation`.  Segregating sites are recorded as haplotype evidence rather than called
    as substitutions, and positions in `skip_sites` are ignored entirely.

    Args:
        read: A pysam aligned read.
        genome: Mapping of contig name to sequence, as returned by `SeqIO.to_dict`.
        splicing_junction_str: Maps a `(first_position, last_position)` gap to its
            mutation string, as built by `find_all_ref_splice_junctions`.
        seg_sites: Segregating-site positions for this gene.
        skip_sites: Reference positions to ignore.

    Returns:
        `(seg_site_bases, muts)` where `seg_site_bases` maps each covered
        segregating site to the base observed there, and `muts` is the list of
        mutation strings in reference order. Soft-clipped ends and the final aligned
        reference position are excluded from mutation calling.
    """
    muts = []
    seg_site_bases = {}

    chrm = read.reference_name
    query = read.query
    start_clip, last_rr = get_start_clip_and_last_rr(read)
    ins_seq = ''
    ins_after = None
    ins_active = False
    del_start = None
    del_active = False
    prev_rr = None
    for qq, rr in read.aligned_pairs:
        # skip clips at beginning and cend
        if rr is None and prev_rr is None and not ins_active:
            continue
        if rr is not None and rr >= last_rr:
            break

        # Process deletions and insertions after completed
        if del_active and qq is not None:
            if (del_start, prev_rr) not in splicing_junction_str:
                muts.append(f'D{del_start}-{prev_rr}')
            else:
                muts.append(splicing_junction_str[(del_start, prev_rr)])
            del_active = False
            del_start = None
        if ins_active and rr is not None:
            muts.append(f'I{ins_after}{ins_seq}')
            ins_seq = ''
            ins_after = None
            ins_active = False

        # Start/continue indels and find mismatches and seg sites
        if qq is None:
            # deletion or splicing
            if not del_active:
                del_active = True
                del_start = rr
        elif rr is None:
            # insertion
            qbase = query[qq-start_clip].upper()
            ins_seq += qbase
            if not ins_active:
                ins_active = True
                ins_after = prev_rr
        else:
            # aligned bases
            qbase = query[qq-start_clip].upper()
            rbase = genome[chrm][rr].upper()

            if rr in seg_sites:
                seg_site_bases[rr] = qbase
            elif rr in skip_sites:
                pass
            elif qbase != rbase:
                muts.append(f'{rbase}{rr}{qbase}')

        prev_rr = rr

    return seg_site_bases, muts


def find_all_ref_splice_junctions(arguments):
    """Haplotyped mutation statistics counting pipeline"""

    fig_dir = os.path.join(arguments.results_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    ### Load annotation 
    log.info('Loading annotation...')
    target_info = TargetInfo(arguments.target_info_file, arguments.gene_name)

    ### umi aware splice-junction determination from control sample
    # 
    # We find splicing junctions by finding all observed splice junctions in control cells that have at
    # least 2 umis in the same cell

    log.info('Finding splice junctions in control sample...')
    log.info('Finding control bam files...')
    ctrl_bam_files = glob.glob(os.path.join(arguments.control_bam_dir, '*bam'))
    ctrl_bam_files.sort()
    log.info(f'Found {len(ctrl_bam_files)} files')

    log.info('Building umi maps...')
    ctrl_bam_umi_map_given_bc = {ctrl_bam_fpath:
            umi_tools.get_umi_maps_from_fastq_or_bam_file(
                bam_fpath=ctrl_bam_fpath,
                chrm=target_info.gene_chrm,
                start=target_info.gene_start,
                end=target_info.gene_end)
            for ctrl_bam_fpath in ctrl_bam_files}

    log.info('Counting splice junction umis per cell...')
    # splicing_junctions_list_of_cntrs contains a list of one counter for each cell, where the counter
    # indicates for each splicing junction how many umis support that splicing junction
    splicing_junctions_list_of_cntrs = []
    splicing_junctions_cntr_given_bc_fpath = {}
    n_umis_given_bc_fpath = {}
    for i, bam_fpath in enumerate(ctrl_bam_files):
        if i % 100 == 0:
            log.info(f'  {i}/{len(ctrl_bam_files)}')
        umi_map_given_bc = ctrl_bam_umi_map_given_bc[bam_fpath]
        cell_junctions_given_umi = defaultdict(set)
        for read in pysam.AlignmentFile(bam_fpath).fetch(
                target_info.gene_chrm,
                target_info.gene_start,
                target_info.gene_end):
            if not read.query:
                continue
            bc, umi = umi_tools.bc_and_umi_given_read_name(read.qname)
            corrected_umi = umi_map_given_bc[bc][umi]
            cell_junctions_given_umi[corrected_umi].update(get_splicing_junctions_from_ctrl(read))
        n_umis_given_bc_fpath[bam_fpath] = len(cell_junctions_given_umi)

        cell_junctions_cntr = Counter()
        for umi_splicing_junctions in cell_junctions_given_umi.values():
            cell_junctions_cntr.update(umi_splicing_junctions)

        splicing_junctions_list_of_cntrs.append(cell_junctions_cntr)
        splicing_junctions_cntr_given_bc_fpath[bam_fpath] = cell_junctions_cntr


    all_possible_splicing_junctions = set()
    for cntr in splicing_junctions_list_of_cntrs:
        all_possible_splicing_junctions.update(cntr.keys())

    #### Accept splicing junction if it exists with at least 2 umis in one cell

    junction_cell_counts = Counter()
    junction_total_counts = Counter()
    for junction in all_possible_splicing_junctions:
        for cntr in splicing_junctions_list_of_cntrs:
            if junction in cntr:
                if cntr[junction] >= 2:  # at least 2 umis in one cell
                    junction_cell_counts[junction] += 1
                    junction_total_counts[junction] += cntr[junction]


    log.info(f'Top 50 of {len(junction_cell_counts)} observed splicing junctions:')
    log.info('Splicing junction\tIntron length\tCell count\tumi count')
    for tup, cell_count in junction_cell_counts.most_common()[:50]:
        intron_len = tup[1] - tup[0]
        log.info(f'{tup}\t{intron_len:,d}\t\t{cell_count:,d}\t\t{junction_total_counts[tup]:,d}')


    splicing_junctions = set(junction_cell_counts.keys())
    out_fpath = os.path.join(arguments.results_dir, f'{arguments.run_name}_splicing_junctions.yml')
    with open(out_fpath, 'w') as out:
        yaml.dump(splicing_junctions, out)
    log.info(f'Splice junctions saved to {out_fpath}')


    fig, ax = plt.subplots()
    ax.hist(n_umis_given_bc_fpath.values(), 50)
    ax.set_xlabel('umis per cell in control')
    out_fpath = os.path.join(fig_dir, f'{arguments.run_name}_umis_per_control_cell.pdf')
    fig.savefig(out_fpath)
    log.info(f'umis per control cell saved to {out_fpath}')


    #### Definition
    # We are here defining the `standard_splicing_junctions` as those which are present in >5% of umis
    # in at least 50% of cells.

    log.info('Finding standard splicing junctions (>5% umis in >50% cells)')
    junction_fracs_per_cell = defaultdict(list)
    for junction in splicing_junctions:
        for bam_fpath in ctrl_bam_files:
            n_umis = n_umis_given_bc_fpath[bam_fpath]
            n_junction_umis = splicing_junctions_cntr_given_bc_fpath[bam_fpath][junction]
            if n_umis == 0:
                assert n_junction_umis == 0, bam_fpath
            else:
                junction_fracs_per_cell[junction].append(float(n_junction_umis)/n_umis)

    log.info('Splicing junctions 99% confidence intervals of Fraction umis per cell:')
    standard_splicing_junctions = set()
    conf_intvl_99 = {}
    fig, ax = plt.subplots()
    for junction in splicing_junctions:
        ecdf = ECDF(junction_fracs_per_cell[junction])
        if ecdf.x[int(len(ecdf.x)/2)] > 0.05:
            ax.plot(ecdf.x, ecdf.y, label=str(junction))
            standard_splicing_junctions.add(junction)
        else:
            ax.plot(ecdf.x, ecdf.y)
        ci_lb = max(0.0, ecdf.x[int(len(ecdf.x)*0.5/100)])
        ci_ub = min(1.0, ecdf.x[int(len(ecdf.x)*99.5/100)])
        conf_intvl_99[junction] = (ci_lb, ci_ub)
        log.info(f'{junction}: {ci_lb:.2f}-{ci_ub:.2f}')
    ax.legend()
    ax.set_xlabel('Fraction of umis with junction')
    ax.set_ylabel('CDF')
    out_fpath = os.path.join(fig_dir, f'{arguments.run_name}_fraction_umis_cdf.pdf')
    fig.savefig(out_fpath)
    log.info(f'Fraction umis CDF saved to {out_fpath}')

    out_fpath = os.path.join(arguments.results_dir, f'{arguments.run_name}_standard_splicing_junctions.yml')
    with open(out_fpath, 'w') as out:
        yaml.dump(standard_splicing_junctions, out)
    log.info(f'Standard splice junctions saved to {out_fpath}')

    out_fpath = os.path.join(arguments.results_dir, f'{arguments.run_name}_splicing_junction_frac_umis_per_cell.yml')
    with open(out_fpath, 'w') as out:
        yaml.dump(junction_fracs_per_cell, out)
    log.info(f'Junction fractions per cell saved to {out_fpath}')


    if not standard_splicing_junctions:
        log.warning('No standard splicing junctions found. Gaps will be reported as deletions.')

    log.info('Processing non-standard junctions')
    # Find the non-standard splicing junctions by transitive matching with standard junctions
    splicing_junction_str = {}
    prev_len = -1
    best_standard_start, best_standard_end = {}, {}
    matching_standard_part = {start: end for start, end in standard_splicing_junctions} | {end: start for start, end in standard_splicing_junctions}
    while len(splicing_junction_str) > prev_len:
        prev_len = len(splicing_junction_str)
        for junction in splicing_junctions:
            start, end = junction
            if junction not in splicing_junction_str:
                if junction in standard_splicing_junctions:
                    splicing_junction_str[junction] = f'J{start}_{end}'
                    best_standard_start[start] = start
                    best_standard_end[end] = end
                elif start in best_standard_start or end in best_standard_end:
                    if start in best_standard_start:
                        bss = best_standard_start[start]
                        bse = matching_standard_part[bss]
                        if end not in best_standard_end:
                            best_standard_end[end] = bse
                    else:
                        bse = best_standard_end[end]
                        bss = matching_standard_part[bse]
                        best_standard_start[start] = bss
                    splicing_junction_str[junction] = f'J{bss}_{bse}>{start-bss:+d}_{end-bse:+d}'

    # Catch the rest by finding maximum overlapping standard junction
    for junction in splicing_junctions - set(splicing_junction_str.keys()):
        start, end = junction
        overlap = {}
        for standard_junction in standard_splicing_junctions:
            s_start, s_end = standard_junction
            if start <= s_start < s_end <= end:
                overlap[standard_junction] = s_end - s_start
            elif s_start <= start < end <= s_end:
                overlap[standard_junction] = end - start
            elif s_start <= start <= s_end:
                overlap[standard_junction] = s_end - start
            elif s_start <= end <= s_end:
                overlap[standard_junction] = end - s_start
            else:
                overlap[standard_junction] = 0
        best_overlap_junction = max(standard_splicing_junctions, key=lambda j: overlap[j], default=None)
        if best_overlap_junction is not None and overlap[best_overlap_junction] > 0:
            bss, bse = best_overlap_junction
            splicing_junction_str[junction] = f'J{bss}_{bse}>{start-bss:+d}_{end-bse:+d}'
            log.debug(f'{junction}, {splicing_junction_str[junction]}, {overlap}')
        else:
            log.debug(f'None {junction}, {overlap}')
    log.debug(f'{len(splicing_junction_str)}, {len(splicing_junctions)}')

    out_fpath = os.path.join(arguments.results_dir, f'{arguments.run_name}_splicing_junction_str.yml')
    with open(out_fpath, 'w') as out:
        yaml.dump(splicing_junction_str, out)
    log.info(f'Splice junction representations saved to {out_fpath}')


def haplotyped_mutation_preprocessing(arguments):
    """Haplotyped mutation statistics counting"""
    # We haplotype each read according to the segregating site present in the read that is closest
    # to the target t1 cutsite.

    log.info('Loading splicing info')
    fpath = os.path.join(arguments.results_dir, f'{arguments.run_name}_splicing_junction_str.yml')
    if not os.path.exists(fpath):
        raise ValueError(f'{fpath} does not exist. results-dir must match refsplice')
    with open(fpath) as f:
        splicing_junction_str = yaml.load(f, Loader=yaml.FullLoader)


    # Load genome
    log.info('Loading genome')
    genome = SeqIO.to_dict(SeqIO.parse(open(arguments.genome_file), 'fasta'))


    # Load annotation and build annotation-based functions
    log.info('Loading annotation...')
    target_info = TargetInfo(arguments.target_info_file, arguments.gene_name)


    def get_target_adj_mut_max_frac_and_cov(frac_cntr, cov_cntr, total_reads, max_dist=10):
        """Find the haplotype-defining mutation and report its fraction and coverage.

        The defining mutation is the highest-fraction non-splice mutation lying within
        `max_dist` of any cut site.

        Args:
            frac_cntr: Maps mutation to the fraction of its covering reads that carry it.
            cov_cntr: Maps mutation to the number of reads covering its position.
            total_reads: Total reads for this haplotype.
            max_dist: Maximum distance in bases from a cut site.

        Returns:
            `(mut, max_frac, cov)` if mutation exists
            `(None, 0, total_reads)` otherwise, for wt calling coverage
        """
        if cov_cntr:
            for mut, max_frac in frac_cntr.most_common():
                if not mut.startswith('J') and target_info.mut_min_dist_to_cutsite(mut) <= max_dist:
                    cov = cov_cntr[mut]
                    return mut, max_frac, cov
        # Reach here either because no cov_cntr (no reads have mutations) or because the loop
        # finished without finding a mutation near the cutsite
        max_frac = 0
        cov = total_reads
        return None, max_frac, cov

    def write_all_muts(bam_fpath, out_dir):
        bc = bc_from_fpath(bam_fpath)
        out_fpath = os.path.join(out_dir, f'{bc}.mutations.yml')
        output = []
        for read in pysam.AlignmentFile(bam_fpath).fetch(
                target_info.gene_chrm,
                target_info.gene_start,
                target_info.gene_end):
            if read.query is None:
                continue
            seg_site_bases, muts = get_haplotyped_mutations(
                    read,
                    genome,
                    splicing_junction_str,
                    target_info.seg_sites,
                    target_info.skip_sites
                    )
            haplotype = target_info.maternal_or_paternal(seg_site_bases)
            output.append([read.qname, (read.pos, read.aend), haplotype, seg_site_bases, muts])
        with open(out_fpath, 'w') as out:
            yaml.dump(output, out)

    mutations_dir = f'{arguments.perturbed_bam_dir}/mutations/'
    os.makedirs(mutations_dir, exist_ok=True)


    ### Look in bam files

    log.info('Processing perturbed bam files')
    bam_files = glob.glob(os.path.join(arguments.perturbed_bam_dir, '*bam'))
    bam_files.sort()
    log.info(f'Found {len(bam_files)} files')

    for i, bam_fpath in enumerate(bam_files):
        if i % 100 == 0:
            log.info(f'{i}/{len(bam_files)}')
        write_all_muts(bam_fpath, mutations_dir)
    log.info('Mutation files written')


    mut_fpaths = glob.glob(os.path.join(mutations_dir, '*.mutations.yml'))

    log.info('Loading all mutation information')
    all_muts_with_info_by_bc = {}
    for fpath in mut_fpaths:
        bc = bc_from_fpath(fpath)
        all_muts_with_info_by_bc[bc] = load_all_muts_with_info(fpath)

    all_bcs = list(all_muts_with_info_by_bc.keys())
    log.info(f'Cells with mutation information: {len(all_bcs)}')

    # We output results for all cells based on cleaned results from top mutation by target site. We
    # then compare replicates to find regions of high accuracy.

    # Output:
    # * Stats file with
    #     * Barcode
    #     * By haplotype:
    #         * Defining mutation
    #             * with coverage and frac
    #         * Before and after:
    #             * Total reads
    #             * Total umis
    #             * For each mutation
    #                 * Mutation identity
    #                 * For reads and umi:
    #                     * coverage
    #                     * count with mutation
    #                     * fraction with mutation


    log.info('Building stat counter struct')
    stat_cntrs_before_and_after = {}
    hap_reads_thresh = 1
    for i, bc in enumerate(all_bcs):
        if i % 100 == 0:
            log.info(f'{i}/{len(all_bcs)}')
        all_muts_with_info = all_muts_with_info_by_bc[bc]
        if not all_muts_with_info:
            continue

        stat_cntrs_before_and_after[bc] = {}

        for haplotype in haplotypes:
            if haplotype not in all_muts_with_info or len(all_muts_with_info[haplotype]) < hap_reads_thresh:
                continue

            # Find all muts observed in haplotype
            all_possible_muts = {}
            for _, _, muts in all_muts_with_info[haplotype]:
                for mut in muts:
                    if mut not in all_possible_muts:
                        all_possible_muts[mut] = parse_mutation(mut)

            # Count umis and build correction map and cntr
            umi_cntr = Counter(umi_tools.umi_given_read_name(read_name) for read_name, _, _ in all_muts_with_info[haplotype])
            corrected_umi_map = umi_tools.get_umi_map_from_cntr(umi_cntr)
            corrected_umi_cntr = Counter()
            for umi, count in umi_cntr.items():
                corrected_umi_cntr[corrected_umi_map[umi]] += count

            # Find all the stats before cleaning
            before_total_reads = 0
            before_total_umis = len(corrected_umi_cntr)

            stat_types = ['mut', 'cov', 'frac']
            before_read_cntr = {stat: Counter() for stat in stat_types}
            before_umi_cntr = {stat: Counter() for stat in stat_types}

            for read_name, (alignment_start, alignment_end), muts in all_muts_with_info[haplotype]:
                corrected_umi = corrected_umi_map[umi_tools.umi_given_read_name(read_name)]
                umi_inc_val = 1/corrected_umi_cntr[corrected_umi]

                before_total_reads += 1
                for mut in muts:
                    before_read_cntr['mut'][mut] += 1
                    before_umi_cntr['mut'][mut] += umi_inc_val
                for mut, (mut_type, mut_start, mut_end, bases) in all_possible_muts.items():
                    if alignment_start <= mut_start < alignment_end or alignment_start < mut_end <= alignment_end:
                        before_read_cntr['cov'][mut] += 1
                        before_umi_cntr['cov'][mut] += umi_inc_val
            before_read_cntr['frac'] = Counter({mut: before_read_cntr['mut'][mut]/before_read_cntr['cov'][mut] if before_read_cntr['cov'][mut] > 0 else 0 for mut in all_possible_muts.keys()})
            before_umi_cntr['frac'] = Counter({mut: before_umi_cntr['mut'][mut]/before_umi_cntr['cov'][mut] if before_umi_cntr['cov'][mut] > 0 else 0 for mut in all_possible_muts.keys()})

            # Find the haplotype-defining mutation: the target-adjacent mutation with max_frac of coverage
            #    def_... here means haplotype-defining, e.g. def_mut is the haplotype-defining mutation
            def_mut, def_max_frac, def_cov = get_target_adj_mut_max_frac_and_cov(before_read_cntr['frac'], before_read_cntr['cov'], before_total_reads)

            # Find all the stats after cleaning
            after_total_reads = 0
            after_observed_umis = set()

            stat_types = ['mut', 'cov', 'frac']
            after_read_cntr = {stat: Counter() for stat in stat_types}
            after_umi_cntr = {stat: Counter() for stat in stat_types}

            for read_name, (alignment_start, alignment_end), muts in all_muts_with_info[haplotype]:
                # Filter by defining mutation
                if def_mut not in muts:
                    continue

                corrected_umi = corrected_umi_map[umi_tools.umi_given_read_name(read_name)]
                umi_inc_val = 1/corrected_umi_cntr[corrected_umi]

                after_total_reads += 1
                after_observed_umis.add(corrected_umi)
                for mut in muts:
                    after_read_cntr['mut'][mut] += 1
                    after_umi_cntr['mut'][mut] += umi_inc_val
                for mut, (mut_type, mut_start, mut_end, bases) in all_possible_muts.items():
                    if alignment_start <= mut_start < alignment_end or alignment_start < mut_end <= alignment_end:
                        after_read_cntr['cov'][mut] += 1
                        after_umi_cntr['cov'][mut] += umi_inc_val
            after_total_umis = len(after_observed_umis)
            after_read_cntr['frac'] = Counter({mut: after_read_cntr['mut'][mut]/after_read_cntr['cov'][mut] if after_read_cntr['cov'][mut] > 0 else 0 for mut in all_possible_muts.keys()})
            after_umi_cntr['frac'] = Counter({mut: after_umi_cntr['mut'][mut]/after_umi_cntr['cov'][mut] if after_umi_cntr['cov'][mut] > 0 else 0 for mut in all_possible_muts.keys()})

            # Add results to output
            stat_cntrs_before_and_after[bc][haplotype] = {
                'defining mutation': {
                    'mut': def_mut,
                    'max_frac': def_max_frac,
                    'cov': def_cov,
                },
                'before': {
                    'total_reads': before_total_reads,
                    'total_umis': before_total_umis,
                    'read_cntr': before_read_cntr,
                    'umi_cntr': before_umi_cntr,
                },
                'after': {
                    'total_reads': after_total_reads,
                    'total_umis': after_total_umis,
                    'read_cntr': after_read_cntr,
                    'umi_cntr': after_umi_cntr,
                }
            }


    out_fpath = os.path.join(arguments.results_dir, f'{arguments.run_name}_stat_cntrs_before_and_after.yml')
    with open(out_fpath, 'w') as out:
        yaml.dump(stat_cntrs_before_and_after, out)
    log.info(f'Stat counter struct written to {out_fpath}')
