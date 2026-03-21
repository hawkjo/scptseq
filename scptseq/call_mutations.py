import logging
import os
import yaml
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from collections import Counter, defaultdict
from .misc import parse_mutation, TargetInfo
from .constants import haplotypes

log = logging.getLogger(__name__)
plt.set_loglevel('critical')



def call_mutations(arguments):
    """Mutation calling with QC figure genetation"""

    if not os.path.exists(arguments.output_dir):
        os.mkdir(arguments.output_dir)
    fig_dir = os.path.join(arguments.output_dir, 'figures')
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)

    # Load annotation
    log.info('Loading annotation...')

    target_info = TargetInfo(arguments.target_info_file, arguments.gene_name)

    before_after = ['before', 'after']

    fpath = os.path.join(arguments.output_dir, f'{arguments.run_name}_standard_splicing_junctions.yml')
    with open(fpath) as f:
        standard_splicing_junctions = yaml.load(f, Loader=yaml.FullLoader)
    fpath = os.path.join(arguments.output_dir, f'{arguments.run_name}_splicing_junction_str.yml')
    with open(fpath) as f:
        splicing_junction_str = yaml.load(f, Loader=yaml.FullLoader)

    # Load data

    log.info('Loading data...')
    fpath = os.path.join(arguments.output_dir, f'{arguments.run_name}_stat_cntrs_before_and_after.yml')
    with open(fpath) as f:
        stat_cntrs_before_and_after = yaml.load(f, Loader=yaml.Loader)

    # Onward

    all_bcs = set([bc for bc in stat_cntrs_before_and_after.keys()])

    ### Read count per cell figures

    xy = {boa: [[], []] for boa in before_after}
    stat_type = 'total_reads'
    for boa in before_after:
        for i, hap in enumerate(haplotypes):
            for bc in all_bcs:
                try:
                    val = stat_cntrs_before_and_after[bc][hap][boa][stat_type]
                except KeyError:
                    val = 0
                xy[boa][i].append(val)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for i, (ax, hap) in enumerate(zip(axes, haplotypes)):
        max_val = max(max(xy['before'][i]), max(xy['after'][i]))
        bins = np.array([0] + list(np.logspace(0, np.log10(max_val+1), 20))) - 0.1
        for boa in ['after']:
            ax.hist(xy[boa][i], bins, histtype='step', label=boa)
        ax.set_xlabel('Reads per Cell')
        ax.set_ylabel('Cells')
        ax.set_title(f'{arguments.run_name} {hap}')
        ax.set_xscale('log')
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_reads_per_cell_haplotype.pdf'))

    fig, ax = plt.subplots()
    data = [el for d in xy['after'] for el in d]
    max_val = max(data)
    bins = np.array([0] + list(np.logspace(0, np.log10(max_val+1), 25))) - 0.1
    ax.hist(data, bins, histtype='step')
    ax.set_xlabel('Reads per Cell')
    ax.set_ylabel('Cells')
    ax.set_title(f'{arguments.run_name}')
    ax.set_xscale('log')
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_reads_per_cell_hist.pdf'))


    # max_frac vs cov threshold fig

    def get_threshes(hap):
        xx_lower_thresh = 0.3 if hap == 'maternal' else 0.3
        xx_upper_thresh = 0.3 if hap == 'maternal' else 0.3
        min_evidence_reads=3
        return min_evidence_reads, xx_lower_thresh, xx_upper_thresh

    x, y = [], []
    for bc in all_bcs:
        for hap in haplotypes:
            try:
                xx = stat_cntrs_before_and_after[bc][hap]['defining mutation']['max_frac']
            except KeyError:
                xx = 0
            try:
                yy = stat_cntrs_before_and_after[bc][hap]['defining mutation']['cov']
            except KeyError:
                yy = 0
            x.append(xx)
            y.append(yy)
    xy_cntr = Counter([(xx, yy) for xx, yy in zip(x, y)])
    c = [xy_cntr[(xx, yy)] if yy > 0 else 0 for xx, yy in zip(x, y)]

    fig, ax = plt.subplots()
    scp = ax.scatter(x, y, c=c, norm=matplotlib.colors.LogNorm())
    ax.set_title(f'{arguments.run_name}')
    plt.colorbar(scp, label=f'Number of chromosomes')
    ax.set_xlabel('max_frac')
    ax.set_ylabel('cov')
    ax.set_yscale('log')
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    min_evidence_reads, xx_lower_thresh, xx_upper_thresh = get_threshes(hap)
    yy_thresh = min_evidence_reads/xx_lower_thresh
    ax.plot((xlim[0], xx_lower_thresh), [yy_thresh]*2, '--', color='grey')
    ax.plot([xx_lower_thresh]*2, ylim, '--', color='grey')
    thresh_x = np.linspace(xx_lower_thresh, xlim[1], 100)
    thresh_y = [min_evidence_reads/xx for xx in thresh_x]
    ax.plot(thresh_x, thresh_y, '--', color='grey')
    ax.plot([xx_upper_thresh]*2, [min_evidence_reads/xx_upper_thresh, ylim[1]], '--', color='grey')
    ax.set_xlim(xlim), ax.set_ylim(ylim)
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_max_frac_v_cov.pdf'))


    # Find and label splicing changes

    standard_splicing_junction_strs = [splicing_junction_str[junction] for junction in standard_splicing_junctions]
    standard_splicing_junction_strs

    max_dist = 10
    all_valid_muts = set()
    for bc in all_bcs:
        for hap in haplotypes:
            try:
                def_mut = stat_cntrs_before_and_after[bc][hap]['defining mutation']['mut']
            except KeyError:
                continue
            if not def_mut:
                continue
            for mut, frac in stat_cntrs_before_and_after[bc][hap]['after']['read_cntr']['frac'].items():
                if target_info.mut_min_dist_to_cutsite(mut) <= max_dist and frac > 0.5:
                    all_valid_muts.add(mut)


    # Find deletions that are actually missplicings and replace
    #   Use transitive matching with standard junctions
    splicing_junction_str_given_del = {}
    prev_len = -1
    best_standard_start = {start: parse_mutation(junc_str)[1] for (start, end), junc_str in splicing_junction_str.items()}
    best_standard_end = {end: parse_mutation(junc_str)[2] for (start, end), junc_str in splicing_junction_str.items()}
    matching_standard_part = {start: end for start, end in standard_splicing_junctions} | {end: start for start, end in standard_splicing_junctions}
    while len(splicing_junction_str_given_del) > prev_len:
        print(prev_len)
        prev_len = len(splicing_junction_str_given_del)
        for mut in all_valid_muts:
            mut_type, start, end, bases = parse_mutation(mut)
            junction = (start, end)
            if mut_type != 'del':
                continue
            if junction not in splicing_junction_str_given_del:
                if junction in standard_splicing_junctions:
                    splicing_junction_str_given_del[junction] = f'J{start}_{end}'
                    best_standard_start[start] = start
                    best_standard_end[end] = end
                elif start in best_standard_start or end in best_standard_end:
                    if start in best_standard_start:
                        bss = best_standard_start[start]
                        bse = matching_standard_part[bss]
                    else:
                        bse = best_standard_end[end]
                        bss = matching_standard_part[bse]

                    half = (bse - bss)/2
                    if start < bss + half and end > bse - half and (start == bss or end == bse or (start in best_standard_start and end in best_standard_end)): 
                        # Require ends to be on their respective sides of the splicing junction, 
                        # and require that one of the ends be a standard splicing junction or 
                        # both ends be previously accepted missplicing junction ends
                        best_standard_start[start] = bss
                        best_standard_end[end] = bse
                        splicing_junction_str_given_del[mut] = f'J{bss}_{bse}>{start-bss:+d}_{end-bse:+d}'
    print(len(splicing_junction_str_given_del))
    splicing_junction_str_given_del

    for bc in all_bcs:
        for hap in haplotypes:
            try:
                def_mut = stat_cntrs_before_and_after[bc][hap]['defining mutation']['mut']
            except KeyError:
                continue
            if not def_mut:
                continue
            for cf in ['cov', 'frac']:
                for mut, junc_str in splicing_junction_str_given_del.items():
                    if mut in stat_cntrs_before_and_after[bc][hap]['after']['read_cntr'][cf]:
                        stat_cntrs_before_and_after[bc][hap]['after']['read_cntr'][cf][junc_str] = stat_cntrs_before_and_after[bc][hap]['after']['read_cntr'][cf][mut]
                        del stat_cntrs_before_and_after[bc][hap]['after']['read_cntr'][cf][mut]


    # Call mut vs wt

    def get_mut_sig(bc, hap, max_dist=10):
        try:
            def_mut = stat_cntrs_before_and_after[bc][hap]['defining mutation']['mut']
        except KeyError:
            return set()
        if not def_mut:
            return set()
        mut_sig = set()
        for mut, frac in stat_cntrs_before_and_after[bc][hap]['after']['read_cntr']['frac'].items():
            if target_info.mut_min_dist_to_cutsite(mut) <= max_dist and frac > 0.5 and mut not in standard_splicing_junction_strs:
                mut_sig.add(mut)
        for junc_str in standard_splicing_junction_strs:
            if any(mut.startswith(junc_str) for mut in mut_sig):
                continue
            cntr = stat_cntrs_before_and_after[bc][hap]['after']['read_cntr']
            if junc_str in cntr['frac'] and cntr['frac'][junc_str] < 0.05 and cntr['cov'][junc_str] >= 3:
                mut_sig.add(f'{junc_str}>-')
        return mut_sig

    def get_cell_haplotype_status(hap, max_frac, cov):
        min_evidence_reads, xx_lower_thresh, xx_upper_thresh = get_threshes(hap)
        if max_frac >= xx_upper_thresh and max_frac*cov >= min_evidence_reads-0.5:
            return 'mut'
        elif max_frac <= xx_lower_thresh and cov >= yy_thresh:
            return 'wt'
        elif xx_lower_thresh < max_frac < xx_upper_thresh and max_frac*cov >= min_evidence_reads-0.5:
            return 'high_cov_uncalled'
        else:
            return 'low_cov_uncalled'


    status_given_bc_hap = {bc: {} for bc in all_bcs}
    mut_cntr_given_status = defaultdict(Counter)
    for bc in all_bcs:
        for hap in haplotypes:
            try:
                xx = stat_cntrs_before_and_after[bc][hap]['defining mutation']['max_frac']
            except KeyError:
                continue
            try:
                yy = stat_cntrs_before_and_after[bc][hap]['defining mutation']['cov']
            except KeyError:
                continue

            mut_sig = get_mut_sig(bc, hap)
            status = get_cell_haplotype_status(hap, xx, yy)
            mut_cntr_given_status[status][frozenset(mut_sig)] += 1
            status_given_bc_hap[bc][hap] = status if status != 'mut' else mut_sig


    # Most common mutations fig

    mut_cntr = mut_cntr_given_status['mut']
    total_cells = sum(mut_cntr.values())

    fig, ax = plt.subplots()
    nvals = 10
    vals, labels = [], []
    for mut_sig, count in mut_cntr.most_common()[:nvals]:
        vals.append(count)
        labels.append(','.join(sorted(mut_sig, key=lambda mut: parse_mutation(mut)[1])) if mut_sig else 'None')
    if len(mut_cntr) > nvals:
        vals.append(sum(count for mut_sig, count in mut_cntr.most_common()[nvals:]))
        labels.append(f'{len(mut_cntr) - nvals} others combined')
    ax.barh(range(len(vals)), width=vals)
    ax.set_title(f'Most common mutation signatures')
    ax.set_xlabel('Chromosomes')
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    ax.set_ylim(ylim[::-1])
    ax.text(total_cells, sum(ylim)/2, 'Total cells', ha='right', va='center', rotation=90)
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_most_common_muts.pdf'))


    # Call mutations and write to file

    uncalled_set = set(['n/a', 'low_cov_uncalled', 'high_cov_uncalled'])

    fpath = os.path.join(arguments.output_dir, f'{arguments.run_name}_mutation_statuses.txt')
    with open(fpath, 'w') as out:
        out.write('\t'.join(['barcode'] + haplotypes) + '\n')
        for bc in all_bcs:
            if bc in status_given_bc_hap:
                mut_sigs = []
                for hap in haplotypes:
                    if hap in status_given_bc_hap[bc]:
                        mut_sig = status_given_bc_hap[bc][hap]
                        if isinstance(mut_sig, set):
                            mut_sig = ','.join(sorted(mut_sig, key=lambda mut: parse_mutation(mut)[1]))
                    else:
                        mut_sig = 'n/a'
                    mut_sigs.append(mut_sig)
                if set(mut_sigs) - uncalled_set:
                    out.write('\t'.join([bc] + mut_sigs) + '\n')


    # Cell zygosity breakdown pie chart

    het_hom_status_given_bc = {}
    for bc in all_bcs:
        if bc not in status_given_bc_hap:
            continue
        mut_sigs = []
        for hap in haplotypes:
            if hap in status_given_bc_hap[bc]:
                mut_sig = status_given_bc_hap[bc][hap]
                if isinstance(mut_sig, set):
                    mut_sig = ','.join(sorted(mut_sig, key=lambda mut: parse_mutation(mut)[1]))
            else:
                mut_sig = 'n/a'
            mut_sigs.append(mut_sig)
        if set(mut_sigs) - uncalled_set == set():
            continue

        if set(['wt']) == set(mut_sigs):
            het_hom_status = 'wt'
        elif set(mut_sigs) & uncalled_set:
            if 'wt' in set(mut_sigs):
                het_hom_status = 'uncalled and wt'
            else:
                het_hom_status = 'uncalled and mut'
        elif mut_sigs[0] == 'wt':
            het_hom_status = 'het paternal'
        elif mut_sigs[1] == 'wt':
            het_hom_status = 'het maternal'
        else:
            het_hom_status = 'hom'
        het_hom_status_given_bc[bc] = het_hom_status

    het_hom_statuses = [het_hom_status_given_bc[bc] for bc in all_bcs if bc in het_hom_status_given_bc]
    het_hom_cntr = Counter(het_hom_statuses)
    labels = sorted(het_hom_cntr.keys())
    vals = [het_hom_cntr[label] for label in labels]
    def clean_label(label):
        for a, b in [('het', '1PC'), ('hom', '2PC'), ('mut', 'PC')]:
            label = label.replace(a, b)
        return label
    labels = [f'{clean_label(label)} ({val})' for label, val in zip(labels, vals)]
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.pie(vals, labels=labels, autopct='%1.1f%%')
    ax.set_title(f'{arguments.run_name} ({sum(vals)} cells)')
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_mut_pie_chart.pdf'))


    # Mutated targets pie chart

    tnames = ['t1', 't2']

    def target_sites_given_mut_str(mut_str, thresh=10):
        muts = mut_str.split(',')
        targets = set()
        for mut in muts:
            mut_type, start, end, bases = parse_mutation(mut)
            for tname, cutsite in zip(tnames, target_info.cutsites):
                dist = min([abs(pos - cutsite) for pos in range(start, end+1)])
                if dist <= thresh:
                    targets.add(tname)
        return targets

    def target_sites_str_given_mut_str(mut_str):
        targets = target_sites_given_mut_str(mut_str)
        return ','.join(sorted(targets))

    target_specific_statuses = []
    for bc in all_bcs:
        if bc not in status_given_bc_hap:
            continue
        mut_sigs = []
        for hap in haplotypes:
            if hap in status_given_bc_hap[bc]:
                mut_sig = status_given_bc_hap[bc][hap]
                if isinstance(mut_sig, set):
                    mut_sig = ','.join(sorted(mut_sig, key=lambda mut: parse_mutation(mut)[1]))
            else:
                mut_sig = 'n/a'
            mut_sigs.append(mut_sig)
        if set(mut_sigs) - uncalled_set == set():
            continue

        target_specific_statuses.append(' / '.join([s if s in set(['wt']) | uncalled_set else target_sites_str_given_mut_str(s) for s in mut_sigs]))

    cntr = {}
    for i, hap in enumerate(haplotypes):
        cntr[i] = Counter()
        for s in target_specific_statuses:
            mut_str = s.split(' / ')[i]
            if mut_str.startswith('t'):
                cntr[i][mut_str] += 1

    fig, ax = plt.subplots()
    labels = sorted(cntr[0].keys())
    vals = [cntr[0][label]+cntr[1][label] for label in labels]
    labels = [f'{label} ({val})' for label, val in zip(labels, vals)]
    ax.pie(vals, labels=labels, autopct='%1.1f%%')
    ax.set_title(f'{arguments.run_name} ({sum(vals)} chromosomes)')
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_target_pie.pdf'))


    # Build mutation statistics by type and location

    def split_muts_by_target_site(mut_str, thresh=10):
        muts = mut_str.split(',')
        muts_by_target = [[], []]
        for mut in muts:
            mut_type, start, end, bases = parse_mutation(mut)
            for i, cutsite in enumerate(target_info.cutsites):
                dist = min([abs(pos - cutsite) for pos in range(start, end+1)])
                if dist <= thresh:
                    muts_by_target[i].append(mut)
        return muts_by_target

    def mut_types_given_muts(muts):
        return [parse_mutation(mut)[0] for mut in muts]

    mut_count_given_hap_target = {hap: {tname: [] for tname in tnames} for hap in haplotypes}
    mut_types_given_hap_target = {hap: {tname: Counter() for tname in tnames} for hap in haplotypes}
    del_lens_given_hap_target = {hap: {tname: [] for tname in tnames} for hap in haplotypes}
    ins_lens_given_hap_target = {hap: {tname: [] for tname in tnames} for hap in haplotypes}
    combined_del_lens_given_hap = {hap: [] for hap in haplotypes}
    combined_ins_lens_given_hap = {hap: [] for hap in haplotypes}
    combined_indel_lens_given_hap = {hap: [] for hap in haplotypes}

    for bc in all_bcs:
        if bc not in status_given_bc_hap:
            continue
        mut_sigs = []
        for hap in haplotypes:
            if hap not in status_given_bc_hap[bc]:
                continue
            mut_sig = status_given_bc_hap[bc][hap]
            if mut_sig in set(['wt']) | uncalled_set:
                continue
            mut_str = ','.join(sorted(mut_sig, key=lambda mut: parse_mutation(mut)[1]))
            muts_by_target = split_muts_by_target_site(mut_str)
            del_len, ins_len, indel_len = 0, 0, 0
            for tname, mut_list in zip(tnames, muts_by_target):
                if not mut_list:
                    continue
                mut_count_given_hap_target[hap][tname].append(len(mut_list))
                for mut_type in mut_types_given_muts(mut_list):
                    mut_types_given_hap_target[hap][tname][mut_type] += 1
                for mut in mut_list:
                    mut_type, start, end, bases = parse_mutation(mut)
                    if mut_type == 'del':
                        mut_len = end - start + 1
                        del_lens_given_hap_target[hap][tname].append(mut_len)
                        del_len += mut_len
                        indel_len -= mut_len
                    if mut_type == 'ins':
                        mut_len = len(bases)
                        ins_lens_given_hap_target[hap][tname].append(mut_len)
                        ins_len += mut_len
                        indel_len += mut_len
            if del_len:
                combined_del_lens_given_hap[hap].append(del_len)
            if ins_len:
                combined_ins_lens_given_hap[hap].append(ins_len)
            if del_len or ins_len:
                combined_indel_lens_given_hap[hap].append(indel_len)


    # Mutations by type by target fig

    possible_mut_types = ['sub', 'del', 'ins', 'splice']
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, tname in zip(axes, tnames):
        vals = [sum(mut_types_given_hap_target[hap][tname][mut_type] for hap in haplotypes) for mut_type in possible_mut_types]
        ax.bar(range(4), vals)
        ax.set_xticks(range(4))
        ax.set_xticklabels(possible_mut_types)
        ax.set_xlabel('Mutations types present')
        ax.set_ylabel('Chromosomes')
        ax.set_title(f'{tname} in {arguments.run_name}')
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_muts_by_type.pdf'))


    # Deletion lens fig

    breakpoints = [20, 30, 100]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, tname in zip(axes, tnames):
        data = [el  for hap in haplotypes for el in del_lens_given_hap_target[hap][tname]]
        bins = list(np.arange(1, breakpoints[0]+2) - 0.5)
        bins += [bins[-1]+0.1, breakpoints[1], breakpoints[2], max(data) + 1]
        idx = list(range(1, breakpoints[0]+2+len(breakpoints)))
        x = list(range(1, breakpoints[0]+2))+list(np.arange(breakpoints[0]+2, breakpoints[0]+2+len(breakpoints)) + 0.5)
        idx_labels = list(range(1, breakpoints[0]+1)) + [''] + breakpoints
        h, _ = np.histogram(data, bins)
        ax.bar(x, h)
        ax.set_xlabel('Deletion lengths')
        ax.set_ylabel('Chromosomes')
        ax.set_title(f'{tname} in {arguments.run_name}')
        ax.set_xticks(idx)
        ax.set_xticklabels(idx_labels, rotation=90)
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_del_lens.pdf'))


    # Insertion lens fig

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, tname in zip(axes, tnames):
        data = [el for hap in haplotypes for el in ins_lens_given_hap_target[hap][tname]]
        bins = list(np.arange(1, breakpoints[0]+2) - 0.5)
        bins += [bins[-1]+0.1, breakpoints[1], breakpoints[2], max(max(data + [1])+1, breakpoints[2])]
        idx = list(range(1, breakpoints[0]+2+len(breakpoints)))
        x = list(range(1, breakpoints[0]+2))+list(np.arange(breakpoints[0]+2, breakpoints[0]+2+len(breakpoints)) + 0.5)
        idx_labels = list(range(1, breakpoints[0]+1)) + [''] + breakpoints
        h, _ = np.histogram(data, bins)
        ax.bar(x, h)
        ax.set_xlabel('Insertion lengths')
        ax.set_ylabel('Chromosomes')
        ax.set_title(f'{tname} in {arguments.run_name}')
        ax.set_xticks(idx)
        ax.set_xticklabels(idx_labels, rotation=90)
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_ins_lens.pdf'))


    # Combined indel lens fig

    fig, ax = plt.subplots(figsize=(9, 4))
    data = [el for hap in haplotypes for el in combined_indel_lens_given_hap[hap]]
    bins = list(np.arange(1, breakpoints[0]+2) - 0.5)
    bins += [bins[-1]+0.1, breakpoints[1], breakpoints[2], max(max(data + [1])+1, breakpoints[2])]
    bins = [-bb for bb in bins][::-1] + bins
    idx = list(range(1, breakpoints[0]+2+len(breakpoints)))
    idx = [-i for i in idx][::-1] + [0] + idx
    x = list(range(1, breakpoints[0]+2))+list(np.arange(breakpoints[0]+2, breakpoints[0]+2+len(breakpoints)) + 0.5)
    x = [-xx for xx in x][::-1] + [0] + x
    idx_labels = list(range(1, breakpoints[0]+1)) + [''] + breakpoints
    idx_labels = [-il if isinstance(il, int) else '' for il in idx_labels][::-1] + [0] + idx_labels
    h, _ = np.histogram(data, bins)
    ax.bar(x, h)
    ax.set_xlabel('Combined Indel lengths')
    ax.set_ylabel('Chromosomes')
    ax.set_title(f'{arguments.run_name}')
    ax.set_xticks(idx)
    ax.set_xticklabels(idx_labels, rotation=90)
    fig.savefig(os.path.join(fig_dir, f'{arguments.run_name}_comb_indel_lens.pdf'))
